"""Exact conditional Fisher test and confidence interval for a 2x2 odds ratio.

The interval inverts one-sided noncentral hypergeometric tail tests at 0.025.
No continuity correction is applied. This is an interval for the conditional OR.
"""

import math


def _distribution(a: int, b: int, c: int, d: int, log_or: float):
    row = a + b
    col = a + c
    total = row + c + d
    lo, hi = max(0, row + col - total), min(row, col)
    xs = range(lo, hi + 1)
    logs = [math.lgamma(col + 1) - math.lgamma(x + 1) - math.lgamma(col - x + 1)
            + math.lgamma(total - col + 1) - math.lgamma(row - x + 1)
            - math.lgamma(total - col - row + x + 1) + x * log_or for x in xs]
    peak = max(logs)
    weights = [math.exp(v - peak) for v in logs]
    norm = sum(weights)
    return lo, [w / norm for w in weights]


def fisher_or_ci(a: int, b: int, c: int, d: int) -> tuple[float, tuple[float, float], float]:
    if any(type(x) is not int or x < 0 for x in (a, b, c, d)):
        raise ValueError("2x2 counts must be nonnegative integers")
    if a + b == 0 or c + d == 0:
        raise ValueError("case and control denominators must be positive")
    or_value = math.inf if b * c == 0 and a * d > 0 else (a * d / (b * c) if b * c else 0.0)
    lo, probs = _distribution(a, b, c, d, 0)
    observed = probs[a - lo]
    p_value = min(1.0, sum(p for p in probs if p <= observed * (1 + 1e-12)))
    minimum, maximum = lo, lo + len(probs) - 1

    def tail(log_or: float, upper: bool) -> float:
        offset, values = _distribution(a, b, c, d, log_or)
        return sum(values[a - offset:] if upper else values[:a - offset + 1])

    def bound(upper_tail: bool) -> float:
        left, right = -700.0, 700.0
        for _ in range(90):
            middle = (left + right) / 2
            prob = tail(middle, upper_tail)
            if (prob > 0.025) == upper_tail:
                right = middle
            else:
                left = middle
        return math.exp((left + right) / 2)

    lower = 0.0 if a == minimum else bound(True)
    upper = math.inf if a == maximum else bound(False)
    return or_value, (lower, upper), p_value
