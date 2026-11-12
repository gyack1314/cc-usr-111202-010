"""条件最大扫描分析：候选族固定在先，校正整个搜索过程。

候选区域只由无向图与 radii 决定（各中心的固定半径闭球，去重，排除空集与
全体中心），绝不依据 observed 裁剪。条件 null 固定总异常数 N，计数向量服从
概率正比于 expected 的多项分布；每次 null 抽样都在全中心重分配 N，并扫描
完整固定候选族取最大统计量，主检验与区域校正 p 值都对照同一份 null 最大值，
不使用单区域自身尾部冒充多重比较校正。
"""
from __future__ import annotations

import csv
import json
import math
import random
import hashlib
from bisect import bisect_left
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from pathlib import Path

from .io import finite, integer, rational

# ---------------------------------------------------------------------------
# 输入解析与校验
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Material:
    ids: tuple[str, ...]
    expected: tuple[Fraction, ...]
    observed: tuple[int, ...]
    adjacency: dict[str, frozenset[str]]
    radii: tuple[int, ...]
    draws: int
    seed: int
    exact_state_budget: int
    comparison_tolerance: float


def _strict_int(x, minimum: int, name: str) -> int:
    if type(x) is not int or x < minimum:
        raise ValueError(f'{name}需要不小于{minimum}的整数')
    return x


def parse_material(value: dict) -> Material:
    """按资料约定严格解析，邻接重复方向合并，自环与未知端点拒绝。"""
    if not isinstance(value, dict):
        raise ValueError('顶层需要对象')
    rows = value.get('centers')
    if not isinstance(rows, list) or not rows:
        raise ValueError('缺少中心')
    ids: list[str] = []
    seen: set[str] = set()
    expected: list[Fraction] = []
    observed: list[int] = []
    for row in rows:
        cid = row['id']
        if not isinstance(cid, str) or not cid.strip() or cid in seen:
            raise ValueError('中心id缺失或重复')
        seen.add(cid)
        e = rational(row['expected'])
        if e <= 0:
            raise ValueError('期望必须正')
        expected.append(e)
        observed.append(integer(row['observed']))
        ids.append(cid)
    known = set(ids)
    edge_set: set[tuple[str, str]] = set()
    for edge in value.get('edges', []):
        if len(edge) != 2 or edge[0] == edge[1] or any(x not in known for x in edge):
            raise ValueError('边引用错误')
        edge_set.add(tuple(sorted(edge)))
    adjacency = {cid: set() for cid in ids}
    for a, b in edge_set:
        adjacency[a].add(b)
        adjacency[b].add(a)
    radii_raw = value.get('radii')
    if not isinstance(radii_raw, list) or not radii_raw:
        raise ValueError('radii需要非空非负整数列表')
    radii = tuple(sorted({_strict_int(r, 0, 'radius') for r in radii_raw}))
    draws = _strict_int(value.get('draws'), 1, 'draws')
    seed = value.get('seed')
    if type(seed) is not int:
        raise ValueError('seed需要整数')
    budget = _strict_int(value.get('exact_state_budget'), 1, 'exact_state_budget')
    tol = finite(value.get('comparison_tolerance'))
    if tol <= 0:
        raise ValueError('comparison_tolerance需要正数')
    return Material(
        ids=tuple(ids),
        expected=tuple(expected),
        observed=tuple(observed),
        adjacency={k: frozenset(v) for k, v in adjacency.items()},
        radii=radii,
        draws=draws,
        seed=seed,
        exact_state_budget=budget,
        comparison_tolerance=float(tol),
    )


# ---------------------------------------------------------------------------
# 候选族：只依赖图与 radii，与 observed/expected 无关
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Region:
    members: tuple[str, ...]
    generators: tuple[tuple[str, int], ...]
    p: Fraction = Fraction(0)


def _closed_balls(adjacency: dict[str, frozenset[str]], center: str,
                  radii: tuple[int, ...]) -> dict[int, frozenset[str]]:
    """BFS 求各半径闭球；图不连通时球只覆盖可达分量。"""
    radius_set = set(radii)
    reached = {center}
    frontier = {center}
    out: dict[int, frozenset[str]] = {}
    if 0 in radius_set:
        out[0] = frozenset(reached)
    for distance in range(1, max(radii) + 1):
        nxt: set[str] = set()
        for node in frontier:
            nxt.update(v for v in adjacency[node] if v not in reached)
        reached |= nxt
        frontier = nxt
        if distance in radius_set:
            out[distance] = frozenset(reached)
    return out


def build_family(material: Material) -> list[Region]:
    """构建固定候选族：相同中心集合去重，空集与全体中心不参加检验。"""
    v_total = len(material.ids)
    collected: dict[frozenset[str], list[tuple[str, int]]] = {}
    for center in material.ids:
        for radius, ball in _closed_balls(material.adjacency, center, material.radii).items():
            if not ball or len(ball) == v_total:
                continue
            collected.setdefault(ball, []).append((center, radius))
    expected_total = sum(material.expected, Fraction(0))
    id_index = {cid: i for i, cid in enumerate(material.ids)}
    family: list[Region] = []
    for ball in sorted(collected, key=lambda s: tuple(sorted(s))):
        members = tuple(sorted(ball))
        share = sum((material.expected[id_index[cid]] for cid in members),
                    Fraction(0)) / expected_total
        family.append(Region(
            members=members,
            generators=tuple(sorted(collected[ball])),
            p=share,
        ))
    return family


def family_fingerprint(family: list[Region]) -> str:
    """候选族规范指纹，用于证明族在查看 observed 之前已固定。"""
    canonical = json.dumps([list(r.members) for r in family],
                           ensure_ascii=False, sort_keys=True)
    return hashlib.sha256(canonical.encode('utf8')).hexdigest()


# ---------------------------------------------------------------------------
# 单侧对数似然比：浮点值 + 有理精确似然比（容差仅触发精核）
# ---------------------------------------------------------------------------

def llr_float(count: int, n_total: int, p: Fraction) -> float:
    """条件多项下的单侧 LLR；count <= N*p 时为 0，0*log0 按 0。"""
    if n_total == 0 or count <= n_total * p:
        return 0.0
    n_p = n_total * p
    value = count * math.log(count / float(n_p)) if count > 0 else 0.0
    rest = n_total - count
    if rest > 0:
        value += rest * math.log(rest / float(n_total * (1 - p)))
    return value


def exact_likelihood_ratio(count: int, n_total: int, p: Fraction):
    """LLR 对应的精确幂积（Fraction）；数过大时退回高精度 Decimal 对数值。

    小整数情形直接用有理概率的似然比幂积独立比较，不依赖浮点；大整数下
    幂积分母分子位数可能失控，改用 60 位 Decimal 对数，仍远严于容差。
    """
    if n_total == 0 or count <= n_total * p:
        return Fraction(1)
    n_p = n_total * p
    rest = n_total - count
    digit_estimate = 0
    if count > 0:
        ratio = Fraction(count, 1) / n_p
        digit_estimate += count * math.log10(max(ratio.numerator, ratio.denominator))
    if rest > 0:
        ratio = Fraction(rest, 1) / (n_total * (1 - p))
        digit_estimate += rest * math.log10(max(ratio.numerator, ratio.denominator))
    if digit_estimate <= 5000:
        product = Fraction(1)
        if count > 0:
            product *= (Fraction(count, 1) / n_p) ** count
        if rest > 0:
            product *= (Fraction(rest, 1) / (n_total * (1 - p))) ** rest
        return product
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    log_value = Decimal(0)
    if count > 0:
        log_value += Decimal(count) * (Decimal(count) / Decimal(n_p.numerator)
                                       * Decimal(n_p.denominator)).ln()
    if rest > 0:
        q = n_total * (1 - p)
        log_value += Decimal(rest) * (Decimal(rest) / Decimal(q.numerator)
                                      * Decimal(q.denominator)).ln()
    return log_value


def compare_stat(float_a: float, key_a, float_b: float, key_b,
                 tolerance: float) -> int:
    """统计量比较：浮点足够分明时直接用，落入容差窗口才做精确精核。

    容差只用于识别需要精核的近似相等，不把任意相近但不等的值当成并列；
    精核确认相等才返回 0，保证临界比较不会因浮点误差排除真正并列。
    """
    scale = max(1.0, abs(float_a), abs(float_b))
    if abs(float_a - float_b) <= tolerance * scale:
        exact_a = exact_likelihood_ratio(*key_a)
        exact_b = exact_likelihood_ratio(*key_b)
        if isinstance(exact_a, Fraction) and isinstance(exact_b, Fraction):
            if exact_a == exact_b:
                return 0
            return -1 if exact_a < exact_b else 1
        # Decimal 高精度对数值路径
        if exact_a == exact_b:
            return 0
        return -1 if exact_a < exact_b else 1
    return -1 if float_a < float_b else 1


# ---------------------------------------------------------------------------
# 扫描
# ---------------------------------------------------------------------------

@dataclass
class RegionStat:
    region: tuple[str, ...]
    count: int
    p: Fraction
    statistic: float
    key: tuple
    elevated: bool


def region_statistics(counts: tuple[int, ...], id_index: dict[str, int],
                      family: list[Region], n_total: int) -> list[RegionStat]:
    out: list[RegionStat] = []
    for region in family:
        c = sum(counts[id_index[cid]] for cid in region.members)
        s = llr_float(c, n_total, region.p)
        out.append(RegionStat(
            region=region.members,
            count=c,
            p=region.p,
            statistic=s,
            key=(c, n_total, region.p),
            elevated=n_total > 0 and c > n_total * region.p,
        ))
    return out


def scan_max(counts: tuple[int, ...], id_index: dict[str, int],
             family: list[Region], n_total: int,
             tolerance: float) -> tuple[float, tuple]:
    """完整固定候选族上的最大统计量（并列保留先遇到者，精确值相同）。"""
    stats = region_statistics(counts, id_index, family, n_total)
    best = stats[0]
    for item in stats[1:]:
        if compare_stat(item.statistic, item.key, best.statistic, best.key,
                        tolerance) > 0:
            best = item
    return best.statistic, best.key


# ---------------------------------------------------------------------------
# 条件 null 抽样：两种引擎共用同一概率向量与同一候选族
# ---------------------------------------------------------------------------

def probability_vector(material: Material) -> list[float]:
    total = sum(material.expected, Fraction(0))
    raw = [float(e / total) for e in material.expected]
    drift = sum(raw)
    return [x / drift for x in raw]  # 归一化修正浮点漂移，两引擎共用


def _weighted_index(rng: random.Random, cumulative) -> int:
    """按累计概率抽取一个下标；必须用 bisect_left 落在对应概率区间。"""
    return bisect_left(cumulative, rng.random())


def _python_multinomial_draw(rng: random.Random, cumulative, n_total: int,
                             v_total: int) -> list[int]:
    counts = [0] * v_total
    for _ in range(n_total):
        counts[_weighted_index(rng, cumulative)] += 1
    return counts


def monte_carlo_python(material: Material, family: list[Region],
                       probabilities: list[float]) -> list[tuple[float, tuple]]:
    """基准纯 Python 引擎：逐项二项分配，无第三方依赖。"""
    rng = random.Random(material.seed)
    cumulative = []
    running = 0.0
    for q in probabilities:
        running += q
        cumulative.append(running)
    cumulative[-1] = 1.0
    id_index = {cid: i for i, cid in enumerate(material.ids)}
    n_total = sum(material.observed)
    out: list[tuple[float, tuple]] = []
    for _ in range(material.draws):
        counts = tuple(_python_multinomial_draw(
            rng, cumulative, n_total, len(material.ids)))
        out.append(scan_max(counts, id_index, family, n_total,
                            material.comparison_tolerance))
    return out


def monte_carlo_numpy(material: Material, family: list[Region],
                      probabilities: list[float]) -> list[tuple[float, tuple]]:
    """批量矩阵加速引擎：候选族与概率向量必须与基准完全一致。

    使用独立的 numpy RNG 流（同种子可复现），多项分布与基准同参数；
    行内最大值附近的并列仍回到有理精确比较，避免 argmax 浮点偏差。
    """
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - 环境无 numpy 时
        raise RuntimeError('请求 numpy 引擎但环境未安装 numpy') from exc
    rng = np.random.default_rng(material.seed)
    v_total = len(material.ids)
    n_total = sum(material.observed)
    id_index = {cid: i for i, cid in enumerate(material.ids)}
    incidence = np.zeros((v_total, len(family)), dtype=np.int64)
    p_array = np.empty(len(family))
    for j, region in enumerate(family):
        for cid in region.members:
            incidence[id_index[cid], j] = 1
        p_array[j] = float(region.p)
    q = np.asarray(probabilities, dtype=np.float64)
    chunk = 4096
    results: list[tuple[float, tuple]] = []
    remaining = material.draws
    while remaining > 0:
        size = min(chunk, remaining)
        remaining -= size
        draws = rng.multinomial(n_total, q, size=size)          # (size, V)
        region_counts = draws @ incidence                       # (size, R)
        n_p = n_total * p_array
        elevated = region_counts > n_p
        with np.errstate(divide='ignore', invalid='ignore'):
            safe_np = np.where(n_p > 0, n_p, 1.0)
            term1 = np.where(region_counts > 0,
                             region_counts * np.log(region_counts / safe_np), 0.0)
            safe_q = np.where(n_total * (1 - p_array) > 0,
                              n_total * (1 - p_array), 1.0)
            rest = n_total - region_counts
            term2 = np.where(rest > 0, rest * np.log(rest / safe_q), 0.0)
        stats = np.where(elevated, term1 + term2, 0.0)
        for row_c, row_s in zip(region_counts, stats):
            j_max = int(np.argmax(row_s))
            # 容差窗口内可能存在真正并列，逐个精核后取精确最大
            scale = max(1.0, float(np.max(row_s)))
            near = np.where(row_s >= row_s[j_max]
                            - material.comparison_tolerance * scale)[0]
            if len(near) > 1:
                best_j = j_max
                for j in near:
                    key_j = (int(row_c[j]), n_total, family[j].p)
                    key_best = (int(row_c[best_j]), n_total, family[best_j].p)
                    if compare_stat(float(row_s[j]), key_j,
                                    float(row_s[best_j]), key_best,
                                    material.comparison_tolerance) > 0:
                        best_j = int(j)
                j_max = best_j
            c = int(row_c[j_max])
            results.append((float(row_s[j_max]), (c, n_total, family[j_max].p)))
    return results


# ---------------------------------------------------------------------------
# 精确模式：枚举全部弱组合，按多项真实概率（Fraction）加权
# ---------------------------------------------------------------------------

class ExactBudgetExceeded(RuntimeError):
    def __init__(self, states: int, budget: int):
        super().__init__(f'精确枚举需要{states}个状态，超过预算{budget}')
        self.states = states
        self.budget = budget


def _weak_compositions(n_total: int, bins: int):
    if bins == 1:
        yield (n_total,)
        return
    for first in range(n_total + 1):
        for rest in _weak_compositions(n_total - first, bins - 1):
            yield (first,) + rest


def exact_enumeration(material: Material, family: list[Region]) -> dict:
    """枚举全部多项计数向量并按真实概率加权，绝不把各向量等权。"""
    v_total = len(material.ids)
    n_total = sum(material.observed)
    states = math.comb(n_total + v_total - 1, v_total - 1)
    if states > material.exact_state_budget:
        raise ExactBudgetExceeded(states, material.exact_state_budget)
    factorial = [1]
    for k in range(1, n_total + 1):
        factorial.append(factorial[-1] * k)
    expected_total = sum(material.expected, Fraction(0))
    base = Fraction(factorial[n_total], 1) / (expected_total ** n_total)
    id_index = {cid: i for i, cid in enumerate(material.ids)}
    observed_stats = region_statistics(material.observed, id_index, family, n_total)
    obs_max = max(observed_stats, key=lambda r: r.statistic)
    tail_main = Fraction(0)
    tail_region = [Fraction(0)] * len(family)
    histogram: dict[str, dict] = {}
    weight_total = Fraction(0)
    for vector in _weak_compositions(n_total, v_total):
        # 多项真实概率：N!/Π xi! · Π (ei/E)^xi；阶乘必须按计数 xi 索引，
        # 不能按中心位置 i 取 factorial[i]
        weight = base
        for x_i, e_i in zip(vector, material.expected):
            weight *= e_i ** x_i
            weight /= factorial[x_i]
        weight_total += weight
        max_s, max_key = scan_max(tuple(vector), id_index, family, n_total,
                                  material.comparison_tolerance)
        lr_key = str(exact_likelihood_ratio(*max_key))
        bucket = histogram.setdefault(lr_key, {
            'likelihood_ratio': lr_key,
            'log_likelihood_ratio': max_s,
            'weight': Fraction(0),
        })
        bucket['weight'] += weight
        for j, obs in enumerate(observed_stats):
            if compare_stat(max_s, max_key, obs.statistic, obs.key,
                            material.comparison_tolerance) >= 0:
                tail_region[j] += weight
        if compare_stat(max_s, max_key, obs_max.statistic, obs_max.key,
                        material.comparison_tolerance) >= 0:
            tail_main += weight
    serialized_histogram = []
    for bucket in sorted(histogram.values(),
                         key=lambda b: -b['log_likelihood_ratio']):
        w: Fraction = bucket['weight']
        serialized_histogram.append({
            'likelihood_ratio': bucket['likelihood_ratio'],
            'log_likelihood_ratio': bucket['log_likelihood_ratio'],
            'weight': str(w),
            'probability': str(w / weight_total),
            'probability_float': float(w / weight_total),
        })
    return {
        'states': states,
        'weight_total': str(weight_total),
        'tail_main': tail_main / weight_total,
        'tail_region': [w / weight_total for w in tail_region],
        'histogram': serialized_histogram,
    }


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def _count_exceed(null_results, statistic, key, tolerance) -> int:
    return sum(1 for s_null, key_null in null_results
               if compare_stat(s_null, key_null, statistic, key, tolerance) >= 0)


def analyze(value: dict, mode: str = 'auto', engine: str = 'python') -> dict:
    """执行条件扫描，返回可直接 JSON 序列化的结果字典。"""
    material = parse_material(value)
    family = build_family(material)
    n_total = sum(material.observed)
    probabilities = probability_vector(material)
    id_index = {cid: i for i, cid in enumerate(material.ids)}
    config = {
        'centers': len(material.ids),
        'radii': list(material.radii),
        'draws': material.draws,
        'seed': material.seed,
        'exact_state_budget': material.exact_state_budget,
        'comparison_tolerance': material.comparison_tolerance,
        'observed_total': n_total,
        'expected_total': str(sum(material.expected, Fraction(0))),
        'null_probabilities': {cid: probabilities[i]
                               for i, cid in enumerate(material.ids)},
        'undirected_edges': sorted(
            [tuple(sorted(e)) for e in
             {tuple(sorted((a, b))) for a, ns in material.adjacency.items() for b in ns}]),
        'engine': engine,
    }
    family_payload = [{
        'region': list(region.members),
        'member_count': len(region.members),
        'expected_share': str(region.p),
        'generators': [list(g) for g in region.generators],
    } for region in family]
    result: dict = {
        'mode': None,
        'status': 'ok',
        'config': config,
        'candidate_family': family_payload,
        'candidate_family_sha256': family_fingerprint(family),
        'family_size': len(family),
        'observed': None,
        'monte_carlo': None,
        'exact': None,
        'notes': [],
    }
    if not family:
        result['status'] = 'no_candidate_regions'
        result['notes'].append('图与 radii 未产生任何有效候选区域（空集与全体中心除外），不执行检验。')
        return result

    observed_stats = region_statistics(material.observed, id_index, family, n_total)

    def cmp_observed(a: RegionStat, b: RegionStat) -> int:
        cmp = compare_stat(a.statistic, a.key, b.statistic, b.key,
                           material.comparison_tolerance)
        if cmp != 0:
            return -cmp  # 统计量降序
        if a.region == b.region:
            return 0
        return -1 if a.region < b.region else 1  # 同分按中心ID元组排序

    ordered = sorted(observed_stats, key=cmp_to_key(cmp_observed))
    max_observed = ordered[0]

    exact_states = math.comb(n_total + len(material.ids) - 1, len(material.ids) - 1)
    use_exact = mode == 'exact' or (mode == 'auto'
                                    and exact_states <= material.exact_state_budget)
    if mode == 'monte-carlo':
        use_exact = False

    if use_exact:
        if mode != 'exact' and exact_states <= material.exact_state_budget:
            result['notes'].append('auto 模式：状态数在预算内，使用精确枚举。')
        exact = exact_enumeration(material, family)
        result['mode'] = 'exact'
        result['exact'] = {
            'states': exact['states'],
            'exact_state_budget': material.exact_state_budget,
            'weight_total': exact['weight_total'],
            'main_p': str(exact['tail_main']),
            'main_p_float': float(exact['tail_main']),
            'null_distribution': exact['histogram'],
        }
        ranking = []
        tail_by_region = {stat.region: tail for stat, tail in
                          zip(observed_stats, exact['tail_region'])}
        for rank, stat in enumerate(ordered, start=1):
            tail = tail_by_region[stat.region]
            ranking.append({
                'rank': rank,
                'region': list(stat.region),
                'count': stat.count,
                'expected_share': str(stat.p),
                'elevated': stat.elevated,
                'statistic': stat.statistic,
                'adjusted_p_exact': str(tail),
                'adjusted_p': float(tail),
            })
        result['observed'] = {
            'maximum_region': list(max_observed.region),
            'maximum_statistic': max_observed.statistic,
            'main_p': str(exact['tail_main']),
            'main_p_float': float(exact['tail_main']),
            'ranking': ranking,
        }
        result['notes'].append('精确模式：枚举全部多项计数向量并按其真实概率加权，非等权计数。')
        return result

    # ---------------- 蒙特卡洛模式 ----------------
    if mode == 'exact':
        raise ExactBudgetExceeded(exact_states, material.exact_state_budget)
    if mode == 'auto':
        result['notes'].append(
            f'auto 模式：精确枚举需要{exact_states}个状态，超过预算'
            f'{material.exact_state_budget}，改用蒙特卡洛抽样；抽样结果不是精确检验。')
    result['mode'] = 'monte_carlo'
    if n_total == 0:
        null_results = [(0.0, (0, 0, Fraction(1)))] * material.draws
        result['notes'].append('N=0：条件 null 退化，所有统计量为 0，所有 p 值为 1。')
    elif engine == 'numpy':
        null_results = monte_carlo_numpy(material, family, probabilities)
        result['notes'].append('numpy 批量引擎与基准引擎共用同一候选族与同一 null 概率向量。')
    else:
        null_results = monte_carlo_python(material, family, probabilities)

    b_total = material.draws
    k_main = _count_exceed(null_results, max_observed.statistic,
                           max_observed.key, material.comparison_tolerance)
    main_p = (1 + k_main) / (b_total + 1)

    def mc_standard_error(p_hat: float) -> float:
        return math.sqrt(p_hat * (1 - p_hat) / (b_total + 1))

    ranking = []
    for rank, stat in enumerate(ordered, start=1):
        k = _count_exceed(null_results, stat.statistic, stat.key,
                          material.comparison_tolerance)
        p_adj = (1 + k) / (b_total + 1)
        ranking.append({
            'rank': rank,
            'region': list(stat.region),
            'count': stat.count,
            'expected_share': str(stat.p),
            'elevated': stat.elevated,
            'statistic': stat.statistic,
            'adjusted_p': p_adj,
            'adjusted_p_standard_error': mc_standard_error(p_adj),
        })
    result['monte_carlo'] = {
        'B': b_total,
        'seed': material.seed,
        'estimator': '(1 + 超过或等于观察值的null次数) / (B + 1)',
        'min_resolution': 1 / (b_total + 1),
        'main_p_standard_error': mc_standard_error(main_p),
        'precision_note': '有限 B 的蒙特卡洛估计，存在二项抽样误差，不是无误差的精确概率；'
                          '标准误差按 sqrt(p(1-p)/(B+1)) 估计，最小可分辨 p 值为 1/(B+1)。',
    }
    result['observed'] = {
        'maximum_region': list(max_observed.region),
        'maximum_statistic': max_observed.statistic,
        'main_p': main_p,
        'main_p_standard_error': mc_standard_error(main_p),
        'ranking': ranking,
    }
    result['_null_maxima'] = [{'draw': i, 'max_statistic': s}
                              for i, (s, _key) in enumerate(null_results)]
    return result


def save_result(result: dict, outdir: str | Path) -> dict[str, Path]:
    """保存候选族/配置摘要、种子、B、全部 null 最大统计量与观察排名。"""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    if result['mode'] == 'monte_carlo':
        null_path = out / 'null_maxima.csv'
        with null_path.open('w', newline='', encoding='utf8') as fh:
            writer = csv.writer(fh)
            writer.writerow(['draw', 'max_statistic'])
            for row in result['_null_maxima']:
                writer.writerow([row['draw'], repr(row['max_statistic'])])
        written['null_maxima'] = null_path
        result.pop('_null_maxima', None)
    elif result['mode'] == 'exact':
        null_path = out / 'exact_null_distribution.json'
        null_path.write_text(json.dumps(
            result['exact']['null_distribution'], ensure_ascii=False,
            indent=2, allow_nan=False), encoding='utf8')
        written['exact_null_distribution'] = null_path
    json_path = out / 'scan_result.json'
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                    allow_nan=False), encoding='utf8')
    written['result'] = json_path
    return written
