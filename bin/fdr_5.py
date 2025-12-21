#!/usr/bin/env python3

import math
from typing import List, Tuple


def benjamini_hochberg(p_values: List[float]) -> List[float]:
    m = len(p_values)
    if m == 0:
        return []

    indexed = sorted(enumerate(p_values), key=lambda x: x[1])

    q_temp = [0.0] * m
    for rank, (idx, p) in enumerate(indexed, start=1):
        p = min(1.0, max(0.0, p))
        q_temp[rank - 1] = min(1.0, p * m / rank)

    for i in range(m - 2, -1, -1):
        if q_temp[i] > q_temp[i + 1]:
            q_temp[i] = q_temp[i + 1]

    result = [0.0] * m
    for rank, (idx, _) in enumerate(indexed, start=1):
        result[idx] = q_temp[rank - 1]

    return result


def _median(xs: List[float]) -> float:
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return math.nan
    mid = n // 2
    if n % 2 == 1:
        return xs[mid]
    return 0.5 * (xs[mid - 1] + xs[mid])


def storey_qvalues(p_values: List[float], lambda_grid=None) -> Tuple[List[float], float]:
    m = len(p_values)
    if m == 0:
        return [], 1.0
    if lambda_grid is None:
        lambda_grid = [i / 100.0 for i in range(50, 95, 5)]

    pi0_estimates = []
    for lam in lambda_grid:
        if lam >= 1.0:
            continue
        cnt = sum(1 for x in p_values if x > lam)
        pi0_hat = cnt / ((1.0 - lam) * m)
        pi0_estimates.append(pi0_hat)
    pi0 = _median(pi0_estimates)
    if math.isnan(pi0):
        pi0 = 1.0
    pi0 = max(0.0, min(1.0, pi0))

    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_temp = [0.0] * m
    for rank, (idx, p) in enumerate(indexed, start=1):
        p = min(1.0, max(0.0, p))
        q_temp[rank - 1] = min(1.0, pi0 * m * p / rank)

    for i in range(m - 2, -1, -1):
        if q_temp[i] > q_temp[i + 1]:
            q_temp[i] = q_temp[i + 1]

    result = [0.0] * m
    for rank, (idx, _) in enumerate(indexed, start=1):
        result[idx] = q_temp[rank - 1]

    return result, pi0
