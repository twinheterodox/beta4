"""Certify the first seed for the ninth power and twentieth derivative.

Uses the two Laurent coefficient algorithms in kernel_coefficients.py.
Run with python-flint installed for --cross-check. No floating-point step.
"""

import argparse
import json
import math
from kernel_coefficients import N, ORDERS, coefficients_by_products, coefficients_by_recurrence


def certificate(power, derivative, p, cross_check):
    assert all(p % k for k in range(2, math.isqrt(p) + 1))
    highest = max(ORDERS) * power + derivative
    assert p > max(480, highest)
    rows = [coefficients_by_recurrence(a, power, p) for a in range(N + 1)]
    if cross_check:
        for a, row in enumerate(rows):
            assert row == coefficients_by_products(a, power, p)
    binomials = [0] + [math.comb(i + derivative - 1, derivative) % p
                       for i in range(1, highest - derivative + 1)]

    # First sum the polar coefficients against all preceding harmonic terms.
    harmonic = [0] * (highest + 1)
    normalized = 0
    for a, row in enumerate(rows):
        for i, coeff in enumerate(row, 1):
            assert rows[N - a][i - 1] == (-1)**i * coeff % p
            normalized -= (-1)**a * binomials[i] * coeff * harmonic[i + derivative]
        inverse = pow(2 * a + 1, -1, p)
        for exponent in range(derivative + 1, highest + 1):
            harmonic[exponent] = (harmonic[exponent] + (-1)**a * pow(inverse, exponent, p)) % p
    normalized %= p

    # Independently group by l=a-h as in the generating function.
    by_contour = 0
    for ell in range(1, N + 1):
        for a in range(ell, N + 1):
            inv = pow(1 + 2 * (a - ell), -1, p)
            value = 0
            for i in range(len(rows[a]), 0, -1):
                value = (value + rows[a][i - 1] * binomials[i]) * inv % p
            by_contour -= (-1)**ell * pow(inv, derivative, p) * value
    assert by_contour % p == normalized and normalized
    return {"power": power, "derivative": derivative, "prime": p,
            "normalized_seed": normalized, "seed": math.factorial(derivative) * normalized % p,
            "local_coefficients": sum(map(len, rows))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cross-check", action="store_true")
    options = parser.parse_args()
    checks = [certificate(9, 20, p, options.cross_check) for p in (10007, 10009)]
    print(json.dumps({"local_products_replayed": options.cross_check, "checks": checks}, indent=2))
