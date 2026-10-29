"""从发布原始 CSV 与词典包复算心电归档标签一致性。仅使用 Python 标准库。"""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime
from decimal import Decimal
import argparse
import ast
import csv
import io
import json
import math
import zipfile
import sys

TARGETS = {
    "AFIB": "Atrial fibrillation",
    "CLBBB": "Complete left bundle branch block",
    "CRBBB": "Complete right bundle branch block",
}
STATES = ("positive", "negative", "equivocal", "unavailable")
Z = 1.959963984540054


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return [dict(row, _record_no=i) for i, row in enumerate(csv.DictReader(stream), 2)]


def read_dictionaries(path):
    with zipfile.ZipFile(path) as archive:
        return {
            name: list(csv.DictReader(io.StringIO(archive.read(name).decode("utf-8-sig"))))
            for name in ("scp_statements.csv", "ptbxlToSNOMED.csv", "12slv23ToSNOMED.csv", "snomed_description.csv")
        }


def key(value):
    number = Decimal(str(value).strip())
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError(f"无效匿名标识：{value!r}")
    return str(int(number))


def number(value):
    if value is None or not str(value).strip():
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def target_mappings(dictionaries):
    ptb = {r["Acronym"]: r for r in dictionaries["ptbxlToSNOMED.csv"]}
    descriptions = {int(r["snomed_id"]): r for r in dictionaries["snomed_description.csv"]}
    result = {}
    for code, label in TARGETS.items():
        row = ptb[code]
        candidates = [int(row[f"id{i}"]) for i in range(1, 5)
                      if row.get(f"id{i}") and row.get(f"name{i}", "").casefold() == label.casefold()]
        if len(candidates) != 1:
            raise ValueError(f"目标映射不唯一：{code}")
        concept = candidates[0]
        definition = descriptions[concept]
        if definition["in_both"] != "True" or definition["informative"] != "True":
            raise ValueError(f"目标不满足双方共有、信息性条件：{code}")
        automatic_codes = sorted({r["Acronym"] for r in dictionaries["12slv23ToSNOMED.csv"]
                                  if any(r.get(f"id{i}") == str(concept) for i in range(1, 6))})
        if not automatic_codes:
            raise ValueError(f"缺少自动对应条目：{code}")
        result[code] = {"concept_key": concept, "description": definition["description"],
                        "automatic_codes": automatic_codes, "in_both": True, "informative": True}
    return result


def parse_labels(raw, kind):
    if raw is None or not str(raw).strip():
        return None
    try:
        parsed = ast.literal_eval(raw)
    except (SyntaxError, ValueError):
        return None
    return parsed if isinstance(parsed, kind) else None


def reference_status(labels, target, categories):
    if labels is None:
        return "unavailable", None
    if target not in labels:
        return "negative", None
    value = number(labels[target])
    meta = categories[target]
    # 原节律/形态条目的 0 与诊断 certainty 语义不同。不得将其统一布尔化。
    nondiagnostic = number(meta.get("rhythm")) == 1 or (
        number(meta.get("form")) == 1 and number(meta.get("diagnostic")) != 1)
    if nondiagnostic and value is not None and 0 <= value <= 100:
        return "positive", value
    if value == 100:
        return "positive", value
    if value is not None and 0 < value < 100:
        return "equivocal", value
    return "unavailable", value


def automatic_status(labels, concept):
    if labels is None:
        return "unavailable", None
    values = []
    try:
        for item in labels:
            ident, certainty = item
            if int(ident) == concept:
                value = number(certainty)
                if value is None or not 0 <= value <= 100:
                    return "unavailable", value
                values.append(value)
    except (ValueError, TypeError):
        return "unavailable", None
    if not values:
        return "negative", None
    value = max(values)
    if value == 100:
        return "positive", value
    return ("equivocal", value) if value > 0 else ("unavailable", value)


def pair_state(reference, automatic, scenario="main"):
    if "unavailable" in (reference, automatic):
        return "不可评价"
    if "equivocal" in (reference, automatic):
        if scenario == "main":
            return "疑似"
        replacement = "positive" if scenario == "positive" else "negative"
        reference = replacement if reference == "equivocal" else reference
        automatic = replacement if automatic == "equivocal" else automatic
    return "a" if reference == automatic == "positive" else "b" if automatic == "positive" else "c" if reference == "positive" else "d"


def wilson(successes, denominator):
    if denominator == 0:
        return [None, None]
    proportion = successes / denominator
    divisor = 1 + Z * Z / denominator
    center = (proportion + Z * Z / (2 * denominator)) / divisor
    radius = Z * math.sqrt(proportion * (1 - proportion) / denominator + Z * Z / (4 * denominator * denominator)) / divisor
    return [max(0.0, center - radius), min(1.0, center + radius)]


def metric(rows, target, scenario="main", confidence=False):
    cells = Counter(pair_state(r[f"{target}_reference_status"], r[f"{target}_auto_status"], scenario) for r in rows)
    a, b, c, d = (cells[x] for x in "abcd")
    ref_count = Counter(r[f"{target}_reference_status"] for r in rows)
    auto_count = Counter(r[f"{target}_auto_status"] for r in rows)
    result = {
        "target": target, "n": len(rows), "evaluated": a+b+c+d,
        "missing": cells["不可评价"], "unresolved": cells["疑似"],
        "equivocal_records": sum("equivocal" in (r[f"{target}_reference_status"], r[f"{target}_auto_status"]) for r in rows),
        "reference_status_counts": {s: ref_count[s] for s in STATES},
        "automatic_status_counts": {s: auto_count[s] for s in STATES},
        "a": a, "b": b, "c": c, "d": d,
        "ref_positive": a+c, "auto_positive": a+b,
        "PPA": a/(a+c) if a+c else None, "NPA": d/(b+d) if b+d else None,
    }
    if confidence:
        result.update(PPA_CI=wilson(a, a+c), NPA_CI=wilson(d, b+d))
    return result



def compute(task, out):
    inputs = task / "input_files"
    references = read_csv(inputs / "ptbxl_database.csv")
    automatic = read_csv(inputs / "12sl_statements.csv")
    dictionaries = read_dictionaries(inputs / "statement_mappings.zip")
    mappings = target_mappings(dictionaries)
    categories = {r[""]: r for r in dictionaries["scp_statements.csv"]}
    ref_keys = [key(r["ecg_id"]) for r in references]
    auto_keys = [key(r["ecg_id"]) for r in automatic]
    if len(set(ref_keys)) != len(ref_keys) or len(set(auto_keys)) != len(auto_keys):
        raise ValueError("存在重复 ecg_id，不能静默多对多配对")
    auto_by_key = {key(r["ecg_id"]): r for r in automatic}
    fold = [r for r in references if int(r["strat_fold"]) == 10]
    adults = [r for r in fold if number(r["age"]) is not None and number(r["age"]) >= 18]
    human = [r for r in adults if r["validated_by_human"] == "True"]
    records = []
    for original in human:
        eid = key(original["ecg_id"])
        if eid not in auto_by_key:
            continue
        row = dict(original)
        auto = auto_by_key[eid]
        row.update(ecg_id=eid, patient_id=key(original["patient_id"]), age=number(original["age"]), sex=number(original["sex"]))
        row["sex_group"] = str(int(row["sex"])) if row["sex"] in (0, 1) else "缺失/其他"
        row["quality_group"] = "存在" if any(original[k].strip() for k in ("static_noise", "burst_noise", "baseline_drift", "electrodes_problems")) else "未记录"
        row["auto_record_no"] = auto["_record_no"]
        row["automatic_statements_ext"] = auto["statements_ext"]
        row["automatic_statements_ext_snomed"] = auto["statements_ext_snomed"]
        source_labels = parse_labels(original["scp_codes"], dict)
        auto_labels = parse_labels(auto["statements_ext_snomed"], list)
        for target in TARGETS:
            row[f"{target}_reference_status"], row[f"{target}_reference_value"] = reference_status(source_labels, target, categories)
            row[f"{target}_auto_status"], row[f"{target}_auto_value"] = automatic_status(auto_labels, mappings[target]["concept_key"])
            row[f"{target}_pair"] = pair_state(row[f"{target}_reference_status"], row[f"{target}_auto_status"])
        records.append(row)
    first = {}
    for row in sorted(records, key=lambda r: (datetime.fromisoformat(r["recording_date"]), int(r["ecg_id"]))):
        first.setdefault(row["patient_id"], row)
    patients = sorted(first.values(), key=lambda r: int(r["ecg_id"]))
    for row in records:
        row["is_index"] = first[row["patient_id"]]["ecg_id"] == row["ecg_id"]
    counts_by_patient = Counter(r["patient_id"] for r in records)
    flow = {
        "all_reference": len(references), "all_automatic": len(automatic),
        "reference_unique_ecg": len(set(ref_keys)), "auto_unique_ecg": len(set(auto_keys)),
        "common_ecg": len(set(ref_keys) & set(auto_keys)),
        "reference_only": len(set(ref_keys) - set(auto_keys)), "automatic_only": len(set(auto_keys) - set(ref_keys)),
        "fold10": len(fold), "age_missing": sum(number(r["age"]) is None for r in fold),
        "children": sum(number(r["age"]) is not None and number(r["age"]) < 18 for r in fold),
        "adult": len(adults), "human_adult": len(human), "paired": len(records), "index": len(patients),
        "repeat_patients": sum(n > 1 for n in counts_by_patient.values()),
        "excluded_repeat_records": len(records)-len(patients), "special_age_index": sum(r["age"] >= 300 for r in patients),
    }
    strata = {"sex": {}, "quality": {}}
    for kind, field, groups in (("sex", "sex_group", ("0", "1", "缺失/其他")), ("quality", "quality_group", ("存在", "未记录"))):
        for group in groups:
            subset = [r for r in patients if r[field] == group]
            strata[kind][group] = [metric(subset, target) for target in TARGETS]
    summary = {
        "flow": flow, "main": [metric(patients, t, confidence=True) for t in TARGETS],
        "scenarios": {s: [metric(patients, t, s) for t in TARGETS] for s in ("positive", "negative")},
        "strata": strata, "records": [metric(records, t) for t in TARGETS],
        "sex_counts": {g: sum(r["sex_group"] == g for r in patients) for g in ("0", "1", "缺失/其他")},
        "quality_counts": {g: sum(r["quality_group"] == g for r in patients) for g in ("存在", "未记录")},
        "target_values": {t: {side: dict(Counter(str(r[f'{t}_{suffix}_value']) for r in patients))
                               for side, suffix in (("ref", "reference"), ("auto", "auto"))} for t in TARGETS},
        "mappings": mappings,
        "method": {"z": Z, "interval": "two-sided 95% Wilson score", "display": "1 decimal percent; full precision retained", "comparison": "automatic rows / archived reference columns"},
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "统计结果.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["ecg_id", "patient_id", "recording_date", "age", "sex", "sex_group", "quality_group", "is_index", "_record_no", "auto_record_no", "scp_codes", "automatic_statements_ext", "automatic_statements_ext_snomed"]
    fields += [f"{t}_{suffix}" for t in TARGETS for suffix in ("reference_value", "reference_status", "auto_value", "auto_status", "pair")]
    for name, rows in (("患者配对", patients), ("记录配对", sorted(records, key=lambda r: int(r["ecg_id"])))):
        with (out / f"{name}.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    return summary, patients, sorted(records, key=lambda r: int(r["ecg_id"]))


"""从患者索引原状态计算有限患者预算的精确边界与联合见证。"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
from statistics import NormalDist

TARGET_NAMES = ('AFIB', 'CLBBB', 'CRBBB')
BITS = tuple(f'{i:03b}' for i in range(8))
NEIGHBOURS = {'FLUT', 'IRBBB', 'ILBBB', 'IVCB'}


def rv(row, side):
    suffix = 'reference_status' if side == 'reference' else 'auto_status'
    return ''.join('1' if row[f'{t}_{suffix}'] == 'positive' else '0' for t in TARGET_NAMES)


def truth(bits, target):
    return int(bits != '000') if target == '任一目标' else int(bits[TARGET_NAMES.index(target)])


def change_truth(bits, target):
    if target == '任一目标':
        return '000' if bits != '000' else '100'
    i = TARGET_NAMES.index(target)
    return bits[:i] + ('0' if bits[i] == '1' else '1') + bits[i+1:]


def witness(grouped, choices):
    used = Counter()
    answer = []
    for group, new in choices:
        row = grouped[group][used[group]]
        used[group] += 1
        answer.append({'ecg_id': row['ecg_id'], 'patient_id': row['patient_id'],
                       'reference_before': row['reference_vector'], 'reference_after': new,
                       'automatic_fixed': row['automatic_vector']})
    assert len(answer) == len({r['patient_id'] for r in answer})
    assert all(r['reference_before'] != r['reference_after'] for r in answer)
    return answer


def endpoint(value, numerator, denominator, rows, **extra):
    return {'value': float(value), 'exact': str(value), 'numerator': numerator,
            'denominator': denominator, 'changed': len(rows), 'witness': rows, **extra}


def four_cell_bounds(rows, target, selector, metric_name, budget):
    groups = defaultdict(list)
    for row in rows:
        selected = bool(selector(row))
        positive = bool(truth(row['reference_vector'], target))
        cell = 'a' if selected and positive else 'b' if selected else 'c' if positive else 'd'
        groups[cell].append(row)
    counts = {c: len(groups[c]) for c in 'abcd'}
    extremes = {}
    for amounts in product(*(range(min(budget, counts[c])+1) for c in 'abcd')):
        if sum(amounts) > budget:
            continue
        xa, xb, xc, xd = amounts
        new = {'a': counts['a']-xa+xb, 'b': counts['b']+xa-xb,
               'c': counts['c']-xc+xd, 'd': counts['d']+xc-xd}
        x, n = (new['a'], new['a']+new['c']) if metric_name == 'PPA' else (new['d'], new['b']+new['d'])
        if not n:
            raise ValueError('本任务预算出现零分母，应单列不可估计而不是丢弃。')
        p = Fraction(x, n)
        for bound in ('min', 'max'):
            old = extremes.get(bound)
            if old is None or (p < old[0] if bound == 'min' else p > old[0]):
                # 同一来源格中各参考向量可能不同；逐个匿名对象保留其实际转换。
                changed = []
                for cell, amount in zip('abcd', amounts):
                    for row in groups[cell][:amount]:
                        changed.append({'ecg_id': row['ecg_id'], 'patient_id': row['patient_id'],
                                        'reference_before': row['reference_vector'],
                                        'reference_after': change_truth(row['reference_vector'], target),
                                        'automatic_fixed': row['automatic_vector']})
                extremes[bound] = (p, endpoint(p, x, n, changed, abcd=new))
    return {'original_abcd': counts, 'min': extremes['min'][1], 'max': extremes['max'][1],
            'proof': '穷举四来源格的翻转人数；每格不超过原人数，总人数不超过K。任一多位改变对本坐标等价于该位翻转或不变；构造见证证明所有枚举点可实现。'}


def joint_dp(rows, target, budget, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row['reference_vector'], row['automatic_vector']].append(row)
    original_positive = sum(truth(r['reference_vector'], target) for r in rows)
    baseline = sum(truth(r['reference_vector'], target)*r['C0'] for r in rows)
    if mode == 'difference':
        baseline = 0
    if mode == 'workload':
        baseline = sum(r['C0'] for r in rows)
        original_positive = len(rows)
    states = {(0, 0, 0): ()}
    for group, objects in sorted(grouped.items()):
        old, auto = group
        old_truth = truth(old, target)
        old_flag = int(old != auto)
        moves = []
        for new in BITS:
            if new == old:
                continue
            new_truth = truth(new, target)
            new_flag = int(new != auto)
            if mode == 'difference':
                dn, dd = new_truth*(new_flag-old_flag), new_truth-old_truth
            elif mode == 'capture':
                dn, dd = new_truth*new_flag-old_truth*old_flag, new_truth-old_truth
            elif mode == 'workload':
                dn, dd = new_flag-old_flag, 0
            else:
                raise ValueError(mode)
            moves.append((dn, dd, new))
        local = {(0, 0, 0): ()}
        frontier = {(0, 0): ()}
        for amount in range(1, min(budget, len(objects))+1):
            nxt = {}
            for (dn, dd), chosen in frontier.items():
                for an, ad, new in moves:
                    nxt.setdefault((dn+an, dd+ad), chosen+((group, new),))
            frontier = nxt
            for (dn, dd), chosen in nxt.items():
                local[amount, dn, dd] = chosen
        nxt = {}
        for (used, dn, dd), chosen in states.items():
            for (more, an, ad), addition in local.items():
                if used+more <= budget:
                    nxt.setdefault((used+more, dn+an, dd+ad), chosen+addition)
        states = nxt
    candidates = []
    for (used, dn, dd), chosen in states.items():
        numerator, denominator = baseline+dn, original_positive+dd
        if denominator <= 0:
            raise ValueError('本任务预算出现不可估计分母。')
        candidates.append((Fraction(numerator, denominator), numerator, denominator, chosen))
    result = {'proof': '按两侧三位向量的实际配对组合计数作有界动态规划；枚举每组至多K人的七种不同参考向量替换。状态保留名额、分子改变量和分母改变量，各组容量受原人数约束。保留所有可达状态后取精确有理数极值。', 'reachable_states': len(states)}
    for bound, item in [('min', min(candidates, key=lambda x: x[0])), ('max', max(candidates, key=lambda x: x[0]))]:
        p, numerator, denominator, chosen = item
        evidence = witness(grouped, chosen)
        result[bound] = endpoint(p, numerator, denominator, evidence)
        if mode == 'workload':
            result[bound]['people'] = numerator
        if mode == 'difference':
            by_patient = {w['patient_id']: w['reference_after'] for w in evidence}
            n = c0 = c1 = 0
            for row in rows:
                new = by_patient.get(row['patient_id'], row['reference_vector'])
                positive = truth(new, target)
                n += positive
                c0 += positive*row['C0']
                c1 += positive*int(new != row['automatic_vector'])
            assert c1-c0 == numerator and n == denominator
            result[bound].update(C0_numerator=c0, C1_numerator=c1)
    return result


def simultaneous_certificate(rows, rate_bounds, budget):
    thresholds = {(t, m): Fraction(rate_bounds[t][m]['min']['exact']) for t in TARGET_NAMES for m in ('PPA', 'NPA')}
    counts = {t: rate_bounds[t]['NPA']['original_abcd'] for t in TARGET_NAMES}
    disjoint = [(a, b) for a, b in (('AFIB', 'CLBBB'), ('AFIB', 'CRBBB'), ('CLBBB', 'CRBBB'))
                if all(not (truth(r['automatic_vector'], a) and truth(r['automatic_vector'], b)) for r in rows)]
    lower_checks = []
    for total in range(2*budget):
        required = {}
        for t in TARGET_NAMES:
            c = counts[t]
            possibles = []
            for additions in range(min(total, c['a'])+1):
                removed = min(total-additions, c['d'])
                denominator = c['b']+c['d']+additions-removed
                if denominator and Fraction(c['d']-removed, denominator) <= thresholds[t, 'NPA']:
                    possibles.append(additions)
            required[t] = min(possibles) if possibles else total+1
        excluded = any(required[a]+required[b] > total for a, b in disjoint)
        lower_checks.append({'patients': total, 'necessary_new_false_positive': required, 'excluded_by_disjoint_auto_positive': excluded})
    assert all(x['excluded_by_disjoint_auto_positive'] for x in lower_checks)
    grouped = defaultdict(list)
    for row in rows:
        if row['reference_vector'] == row['automatic_vector']:
            grouped[row['reference_vector']].append(row)
    choices = None
    for bits in BITS:
        opposite = ''.join('0' if b == '1' else '1' for b in bits)
        if bits in ('000', '111') or bits >= opposite:
            continue
        if min(len(grouped[bits]), len(grouped[opposite])) >= budget:
            choices = [(r, opposite) for r in grouped[bits][:budget]] + [(r, bits) for r in grouped[opposite][:budget]]
            break
    if choices is None:
        raise ValueError('需要另行求解联合整数优化，本数据特定证书不成立。')
    changed = [{'ecg_id': r['ecg_id'], 'patient_id': r['patient_id'], 'reference_before': r['reference_vector'], 'reference_after': new, 'automatic_fixed': r['automatic_vector']} for r, new in choices]
    after = {r['patient_id']: r['reference_after'] for r in changed}
    six = []
    for t in TARGET_NAMES:
        c = Counter()
        for row in rows:
            reference = truth(after.get(row['patient_id'], row['reference_vector']), t)
            automatic = truth(row['automatic_vector'], t)
            c['a' if reference and automatic else 'b' if automatic else 'c' if reference else 'd'] += 1
        for metric, x, n in [('PPA', c['a'], c['a']+c['c']), ('NPA', c['d'], c['b']+c['d'])]:
            assert Fraction(x, n) <= thresholds[t, metric]
            six.append({'target': t, 'metric': metric, 'numerator': x, 'denominator': n, 'value': x/n, 'exact': str(Fraction(x, n)), 'threshold_exact': str(thresholds[t, metric])})
    return {'minimum_patients': len(changed), 'witness': changed, 'six': six,
            'disjoint_automatic_positive_pairs': disjoint, 'lower_bound_checks': lower_checks,
            'proof': '给更少名额M最有利的假设：每目标先最多增加x个自动阳性/新参考阴性，再把余下M−x名额全部用于去掉自动阴性/原参考阴性。该放宽仍达不到NPA阈值时，实际情景也不能达到。自动阳性集合不相交的两个目标，所需新增假阳性患者不能共享；逐M排除后以同一清单验证六指标，达到上下界相等。'}


def compute_robust(task, out, patients, primary, budget=5):
    import ast
    rows = []
    for p in patients:
        assert all(p[f'{t}_reference_status'] in ('positive', 'negative') and p[f'{t}_auto_status'] in ('positive', 'negative') for t in TARGET_NAMES)
        reference, automatic = rv(p, 'reference'), rv(p, 'auto')
        neighbours = any(code in NEIGHBOURS for code, value in ast.literal_eval(p['automatic_statements_ext']))
        rows.append({'ecg_id': p['ecg_id'], 'patient_id': p['patient_id'], 'reference_vector': reference,
                     'automatic_vector': automatic, 'A': int(automatic != '000'),
                     'B': int(automatic != '000' or neighbours or p['quality_group'] == '存在'),
                     'C0': int(reference != automatic)})
    rows.sort(key=lambda r: int(r['ecg_id']))
    n = len(rows)
    matrix = Counter((r['reference_vector'], r['automatic_vector']) for r in rows)
    targets = TARGET_NAMES+('任一目标',)
    rates = {t: {m: four_cell_bounds(rows, t, lambda r, t=t: truth(r['automatic_vector'], t), m, budget)
                  for m in ('PPA', 'NPA')} for t in targets}
    policies = {}
    fixed = {}
    for policy in ('A', 'B', 'C0'):
        selected = sum(r[policy] for r in rows)
        captures = {}
        for t in targets:
            positive = sum(truth(r['reference_vector'], t) for r in rows)
            captured = sum(truth(r['reference_vector'], t)*r[policy] for r in rows)
            captures[t] = {'positive': positive, 'captured': captured, 'missed': positive-captured, 'fraction': captured/positive}
        policies[policy] = {'n': selected, 'fraction': selected/n, 'targets': captures}
        fixed[policy] = {t: four_cell_bounds(rows, t, lambda r, p=policy: r[p], 'PPA', budget) for t in targets}
    dynamic = {t: joint_dp(rows, t, budget, 'capture') for t in targets}
    contrast = {t: joint_dp(rows, t, budget, 'difference') for t in targets}
    workload = joint_dp(rows, '任一目标', budget, 'workload')
    z = NormalDist().inv_cdf(1-0.05/(2*6))
    family = []
    for m in primary['main']:
        item = {'target': m['target']}
        for name, x, den in [('PPA', m['a'], m['a']+m['c']), ('NPA', m['d'], m['b']+m['d'])]:
            p=x/den; center=(p+z*z/(2*den))/(1+z*z/den)
            half=z*((p*(1-p)/den+z*z/(4*den*den))**0.5)/(1+z*z/den)
            item[name+'_CI']=[max(0, center-half), min(1, center+half)]
        family.append(item)
    union = rates['任一目标']['PPA']['original_abcd']
    result = {'budget': budget, 'n': n, 'matrix': [[matrix[a,b] for b in BITS] for a in BITS],
              'hidden_records': [r for r in rows if (r['reference_vector'] != '000') == (r['automatic_vector'] != '000') and r['reference_vector'] != r['automatic_vector']],
              'union': {**union, 'PPA': union['a']/(union['a']+union['c']), 'NPA': union['d']/(union['b']+union['d'])},
              'rates': rates, 'policies': policies, 'fixed_policy_bounds': fixed,
              'C1': dynamic, 'C1_workload': workload, 'paired_C1_minus_C0': contrast,
              'simultaneous': simultaneous_certificate(rows, rates, budget),
              'intersections': {a+b: sum(r[a] and r[b] for r in rows) for a,b in (('A','B'),('A','C0'),('B','C0'))},
              'ABC0': sum(r['A'] and r['B'] and r['C0'] for r in rows), 'capacity': int(n*0.15),
              'family': {'alpha': 0.05, 'm': 6, 'z': z, 'targets': family}, 'source_rows': rows}
    out.mkdir(parents=True, exist_ok=True)
    (out/'联合稳健性.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    return result


if __name__ == "__main__":
    task_default = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", type=Path, default=task_default)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    output = args.out or Path(__file__).resolve().parent / "data"
    result, patients, records = compute(args.task.resolve(), output)
    robust = compute_robust(args.task.resolve(), output, patients, result)
    print(json.dumps({"flow": result["flow"], "main": result["main"], "shared_budget": robust["budget"], "minimum_joint_patients": robust["simultaneous"]["minimum_patients"]}, ensure_ascii=False, indent=2))
