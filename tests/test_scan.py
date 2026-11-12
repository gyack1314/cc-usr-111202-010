"""条件最大扫描的候选族、精确枚举、蒙特卡洛校正与边界测试。"""
import csv
import json
import math
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from research_data.io import load
from research_data.scan import (
    ExactBudgetExceeded,
    Region,
    analyze,
    build_family,
    compare_stat,
    exact_enumeration,
    exact_likelihood_ratio,
    family_fingerprint,
    llr_float,
    monte_carlo_numpy,
    monte_carlo_python,
    parse_material,
    probability_vector,
    region_statistics,
    save_result,
    scan_max,
    _python_multinomial_draw,
)


def _material(observed, expected=None, edges=None, radii=None, **overrides):
    expected = expected or ['2', '3', '2', '4', '1']
    edges = edges or [['c0', 'c1'], ['c1', 'c2'], ['c2', 'c3'], ['c3', 'c4']]
    value = {
        'centers': [{'id': f'c{i}', 'expected': e, 'observed': observed[i]}
                    for i, e in enumerate(expected)],
        'edges': edges,
        'radii': radii if radii is not None else [0, 1, 2],
        'draws': overrides.get('draws', 999),
        'seed': overrides.get('seed', 313),
        'exact_state_budget': overrides.get('exact_state_budget', 100000),
        'comparison_tolerance': overrides.get('comparison_tolerance', 1e-12),
    }
    value.update(overrides)
    return value


class FamilyTests(unittest.TestCase):
    def test_balls_dedup_and_exclusions(self):
        value = _material([0, 0, 0, 0, 0], radii=[0, 1, 2])
        family = build_family(parse_material(value))
        members = [r.members for r in family]
        # 去重：无重复中心集合
        self.assertEqual(len(members), len(set(members)))
        # 空集与全体中心不参加检验
        self.assertNotIn(tuple(), members)
        self.assertNotIn(tuple(f'c{i}' for i in range(5)), members)
        # 半径0闭球就是单点
        singletons = {r.members for r in family if len(r.members) == 1}
        self.assertEqual(singletons, {(f'c{i}',) for i in range(5)})
        # 生成器溯源：c1 半径1 的球是 c0,c1,c2；c2 半径1 是 c1,c2,c3
        by_members = {r.members: r for r in family}
        self.assertIn(('c1', 1), by_members[('c0', 'c1', 'c2')].generators)
        self.assertIn(('c2', 1), by_members[('c1', 'c2', 'c3')].generators)
        # 同一集合只保留一次但记录全部生成器
        self.assertEqual(
            set(by_members[('c0', 'c1', 'c2')].generators),
            {('c0', 2), ('c1', 1)})

    def test_disconnected_graph_balls_stay_in_component(self):
        value = _material(
            [0, 0, 0, 0], expected=['1', '1', '1', '1'],
            edges=[['c0', 'c1']], radii=[0, 5])
        family = build_family(parse_material(value))
        members = {r.members for r in family}
        self.assertIn(('c0', 'c1'), members)
        self.assertIn(('c2',), members)
        self.assertIn(('c3',), members)
        # 大半径也不能跨越连通分量，且不出现全体集合
        self.assertNotIn(('c0', 'c1', 'c2', 'c3'), members)

    def test_family_independent_of_observed(self):
        value = _material([0, 1, 5, 4, 0])
        f1 = build_family(parse_material(value))
        value2 = _material([9, 8, 7, 6, 5])
        f2 = build_family(parse_material(value2))
        self.assertEqual([r.members for r in f1], [r.members for r in f2])
        self.assertEqual(family_fingerprint(f1), family_fingerprint(f2))

    def test_single_center_has_no_valid_region(self):
        value, _ = load('data/single-center.json')
        result = analyze(value)
        self.assertEqual(result['status'], 'no_candidate_regions')
        self.assertIsNone(result['mode'])
        self.assertEqual(result['family_size'], 0)


class StatisticTests(unittest.TestCase):
    def test_one_sided_llr_and_zero_convention(self):
        p = Fraction(1, 4)
        # 不升高（含恰好相等）一律为零
        self.assertEqual(llr_float(2, 8, p), 0.0)
        self.assertEqual(llr_float(0, 8, p), 0.0)
        # 升高时与手算一致
        c, n = 5, 8
        expected = c * math.log(c / (n * p)) + (n - c) * math.log(
            (n - c) / (n * (1 - p)))
        self.assertAlmostEqual(llr_float(c, n, p), expected)
        # 区内计数为 N、区外为 0 时 0*log0 按零，不抛异常
        p_small = Fraction(1, 10)
        value = llr_float(8, 8, p_small)
        self.assertTrue(math.isfinite(value))

    def test_tie_exact_check_not_forced_by_tolerance(self):
        k1 = (501, 1000, Fraction(1, 2))
        k2 = (500, 1000, Fraction(499, 1000))
        f1, f2 = llr_float(*k1), llr_float(*k2)
        self.assertLess(abs(f1 - f2), 1e-6)  # 浮点极近
        cmp = compare_stat(f1, k1, f2, k2, 1e-12)
        self.assertNotEqual(cmp, 0)          # 精核识别出真正不等
        self.assertEqual(
            compare_stat(f1, k1, f1, k1, 1e-12), 0)  # 真并列不被排除

    def test_ranking_tiebreak_by_id_tuple(self):
        value = _material([0] * 5)
        result = analyze(value)
        regions = [tuple(x['region']) for x in result['observed']['ranking']]
        self.assertEqual(regions, sorted(regions))  # 全零统计量按元组排序


class ExactTests(unittest.TestCase):
    def test_exact_weights_are_true_probabilities_not_uniform(self):
        value, _ = load('data/normal.json')
        material = parse_material(value)
        family = build_family(material)
        exact = exact_enumeration(material, family)
        # 概率总和为 1
        total = sum(Fraction(h['probability']) for h in exact['histogram'])
        self.assertEqual(total, Fraction(1))
        # 等权枚举会给出 1/1001；真实概率显然不是均匀的
        probs = [float(Fraction(h['probability'])) for h in exact['histogram']]
        self.assertFalse(any(abs(p - 1 / exact['states']) < 1e-12 for p in probs))

    def test_exact_matches_independent_naive_enumeration(self):
        value, _ = load('data/normal.json')
        material = parse_material(value)
        family = build_family(material)
        id_index = {c: i for i, c in enumerate(material.ids)}
        n_total = sum(material.observed)
        e_total = sum(material.expected, Fraction(0))

        def compositions(n, k):
            if k == 1:
                yield (n,)
                return
            for first in range(n + 1):
                for rest in compositions(n - first, k - 1):
                    yield (first,) + rest

        observed_max = max(
            region_statistics(material.observed, id_index, family, n_total),
            key=lambda r: r.statistic)
        tail = Fraction(0)
        for vector in compositions(n_total, len(material.ids)):
            weight = Fraction(math.factorial(n_total), 1)
            for xi, ei in zip(vector, material.expected):
                weight /= math.factorial(xi)
                weight *= (ei / e_total) ** xi
            s, key = scan_max(vector, id_index, family, n_total,
                              material.comparison_tolerance)
            if compare_stat(s, key, observed_max.statistic, observed_max.key,
                            material.comparison_tolerance) >= 0:
                tail += weight
        exact = exact_enumeration(material, family)
        self.assertEqual(exact['tail_main'], tail)

    def test_budget_enforced_no_partial_exact_report(self):
        value = _material([0, 1, 5, 4, 0], exact_state_budget=10)
        with self.assertRaises(ExactBudgetExceeded):
            analyze(value, mode='exact')
        result = analyze(value, mode='auto')
        self.assertEqual(result['mode'], 'monte_carlo')
        self.assertIsNone(result['exact'])
        self.assertTrue(any('不是精确检验' in n for n in result['notes']))


class MonteCarloTests(unittest.TestCase):
    def test_multinomial_marginals_match_expected_shares(self):
        import random
        value, _ = load('data/normal.json')
        material = parse_material(value)
        q = probability_vector(material)
        cumulative, running = [], 0.0
        for x in q:
            running += x
            cumulative.append(running)
        cumulative[-1] = 1.0
        rng = random.Random(42)
        n_total = sum(material.observed)
        totals = [0] * len(material.ids)
        reps = 30000
        for _ in range(reps):
            draw = _python_multinomial_draw(rng, cumulative, n_total,
                                            len(material.ids))
            self.assertEqual(sum(draw), n_total)
            for i, x in enumerate(draw):
                totals[i] += x
        grand = sum(totals)
        for empirical, target in zip(totals, q):
            self.assertAlmostEqual(empirical / grand, target, delta=0.01)

    def test_mc_converges_to_exact_and_uses_shared_null_maxima(self):
        value, _ = load('data/normal.json')
        exact = analyze(value)
        big = dict(value, draws=15000, seed=7)
        mc = analyze(big, mode='monte-carlo')
        self.assertAlmostEqual(mc['observed']['main_p'],
                               exact['observed']['main_p_float'], delta=0.02)
        # 区域校正 p 全部对照同一份 null 最大值：不小于主 p，且计数公式 (k+1)/(B+1)
        b_total = big['draws']
        for row in mc['observed']['ranking']:
            self.assertGreaterEqual(row['adjusted_p'], mc['observed']['main_p'])
            denominator_inverse = 1 / (b_total + 1)
            self.assertAlmostEqual(
                row['adjusted_p'] / denominator_inverse,
                round(row['adjusted_p'] / denominator_inverse))
        self.assertEqual(mc['monte_carlo']['min_resolution'],
                         1 / (b_total + 1))
        self.assertEqual(len(mc['_null_maxima']), b_total)

    def test_reproducible_seed_and_distinct_seed(self):
        value, _ = load('data/normal.json')
        material = parse_material(value)
        family = build_family(material)
        q = probability_vector(material)
        first = monte_carlo_python(material, family, q)
        second = monte_carlo_python(material, family, q)
        self.assertEqual(first, second)
        other = parse_material(dict(value, seed=999))
        self.assertNotEqual(first, monte_carlo_python(other, family, q))

    def test_numpy_engine_matches_family_and_converges(self):
        try:
            import numpy  # noqa: F401
        except ImportError:
            self.skipTest('环境未安装 numpy，加速引擎为可选依赖')
        value, _ = load('data/normal.json')
        material = parse_material(value)
        family = build_family(material)
        q = probability_vector(material)
        results = monte_carlo_numpy(material, family, q)
        self.assertEqual(len(results), material.draws)
        # 加速引擎仍输出完整固定候选族上的最大统计量
        self.assertTrue(all(math.isfinite(s) for s, _ in results))
        big = parse_material(dict(value, draws=12000, seed=11))
        exact = analyze(value)['observed']['main_p_float']
        null_big = monte_carlo_numpy(big, family, q)
        observed = max(
            region_statistics(material.observed,
                              {c: i for i, c in enumerate(material.ids)},
                              family, sum(material.observed)),
            key=lambda r: r.statistic)
        k = sum(1 for s, key in null_big
                if compare_stat(s, key, observed.statistic, observed.key,
                                material.comparison_tolerance) >= 0)
        self.assertAlmostEqual((1 + k) / (big.draws + 1), exact, delta=0.02)

    def test_n_zero_gives_all_zero_stats_and_p_one(self):
        value, _ = load('data/zero-total.json')
        for mode in ('auto', 'monte-carlo'):
            with self.subTest(mode=mode):
                result = analyze(dict(value, draws=64), mode=mode)
                self.assertEqual(float(result['observed']['main_p']), 1.0)
                self.assertTrue(all(
                    float(row['adjusted_p']) == 1.0 for row in
                    result['observed']['ranking']))
                if mode == 'monte-carlo':
                    self.assertTrue(all(
                        row['max_statistic'] == 0.0
                        for row in result['_null_maxima']))


class PersistenceAndContractTests(unittest.TestCase):
    def test_save_result_layout_and_json_safety(self):
        value, _ = load('data/normal.json')
        result = analyze(value)
        with tempfile.TemporaryDirectory() as tmp:
            written = save_result(result, tmp)
            payload = json.loads(Path(written['result']).read_text('utf8'))
            self.assertEqual(payload['candidate_family_sha256'],
                             result['candidate_family_sha256'])
            self.assertIn('seed', payload['config'])
            self.assertIn('exact_null_distribution', written)
            null = json.loads(Path(written['exact_null_distribution'])
                              .read_text('utf8'))
            self.assertTrue(all('probability' in h for h in null))

    def test_mc_save_writes_full_null_maxima_csv(self):
        value, _ = load('data/normal.json')
        result = analyze(dict(value, draws=137), mode='monte-carlo')
        with tempfile.TemporaryDirectory() as tmp:
            written = save_result(result, tmp)
            with Path(written['null_maxima']).open(encoding='utf8') as fh:
                rows = list(csv.reader(fh))
            self.assertEqual(rows[0], ['draw', 'max_statistic'])
            self.assertEqual(len(rows) - 1, 137)
            # 结果 JSON 已移除内部临时字段且可严格序列化（无 NaN/Infinity）
            json.loads(Path(written['result']).read_text('utf8'))
            self.assertNotIn('_null_maxima',
                             json.loads(Path(written['result'])
                                        .read_text('utf8')))

    def test_invalid_inputs_rejected(self):
        good = _material([0, 1, 2, 3, 4])
        bad_cases = [
            dict(good, centers=[]),
            dict(good, radii=[]),
            dict(good, radii=[-1]),
            dict(good, draws=0),
            dict(good, seed='x'),
            dict(good, exact_state_budget=0),
            dict(good, comparison_tolerance=0),
            dict(good, edges=[['c0', 'c0']]),
            dict(good, edges=[['c0', 'nope']]),
        ]
        for bad in bad_cases:
            with self.subTest(bad=list(bad.keys())):
                with self.assertRaises((ValueError, KeyError, TypeError)):
                    analyze(bad)

    def test_observed_cannot_prune_family(self):
        # 极端 observed 与零 observed 产出同一候选族，证明族不看数据裁剪
        hot = analyze(_material([10, 0, 0, 0, 0]))
        cold = analyze(_material([0, 0, 0, 0, 0]))
        self.assertEqual(hot['candidate_family'], cold['candidate_family'])
        self.assertEqual(hot['candidate_family_sha256'],
                         cold['candidate_family_sha256'])


if __name__ == '__main__':
    unittest.main()
