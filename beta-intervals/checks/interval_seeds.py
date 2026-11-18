"""Rebuild the powered N=96 seed constants at three primes.

The first calculation uses the logarithmic derivative of the seed at a pole.
With --cross-check, python-flint also multiplies the local factors directly.
The algorithms adapt the two published fifteen-beta-values seed checkers.
"""

import argparse
import json
import math


N = 96
ROOTS = range(-239, 240, 2)
ORDERS = [min(a, N - a, 30) + 1 for a in range(N + 1)]


def unit_factors(a):
    x = N - 2 * a
    factors = [(x - v, 3) for v in ROOTS]
    factors.extend((2 * (b - a), -order)
                   for b, order in enumerate(ORDERS) if b != a)
    if x:
        factors.append((x, 1))
    return factors


def coefficients_by_recurrence(a, q, p):
    order = q * ORDERS[a]
    shift = q if 2 * a == N else 0
    length = order - shift
    constant = 1
    logarithmic = [0] * (length - 1)
    for value, exponent in unit_factors(a):
        exponent *= q
        constant = constant * pow(value, exponent, p) % p
        inverse = pow(value, -1, p)
        term = inverse
        for j in range(length - 1):
            logarithmic[j] = (logarithmic[j] + exponent * term) % p
            term = -term * inverse % p
    coefficients = [constant]
    for k in range(1, length):
        derivative = sum(logarithmic[j] * coefficients[k - 1 - j] for j in range(k))
        coefficients.append(derivative * pow(k, -1, p) % p)
    return list(reversed([0] * shift + coefficients))


def coefficients_by_products(a, q, p):
    from flint import fmpz_mod_poly_ctx

    ctx = fmpz_mod_poly_ctx(p)
    order = q * ORDERS[a]
    x = N - 2 * a
    numerator = ctx.one()
    for v in ROOTS:
        numerator = numerator.mul_low(ctx([x - v, 1]), order)
    numerator = numerator.pow_trunc(3, order).mul_low(ctx([x, 1]), order)
    denominator = ctx.one()
    for b, pole_order in enumerate(ORDERS):
        if b != a:
            factor = ctx([2 * (b - a), 1]).pow_trunc(pole_order, order)
            denominator = denominator.mul_low(factor, order)
    unit = numerator.mul_low(denominator.inverse_series_trunc(order), order)
    jet = unit.pow_trunc(q, order)
    return [int(jet[order - i]) for i in range(1, order + 1)]


def constant_by_evaluation_address(rows, q, p):
    derivative = 2 * q
    top = max(map(len, rows))
    factors = [0] + [math.factorial(derivative) * math.comb(i + derivative - 1, derivative) % p
                     for i in range(1, top + 1)]
    suffix = [0] * (top + 1)
    constant = 0
    for h in reversed(range(N)):
        for i, partial in enumerate(rows[h + 1], 1):
            suffix[i] = (suffix[i] + (-1)**(h + 1) * partial) % p
        inverse = pow(2 * h + 1, -1, p)
        value = sum(factors[i] * suffix[i] * pow(inverse, i + derivative, p)
                    for i in range(1, top + 1))
        constant = (constant - (-1)**h * value) % p
    return constant


def seed(q, p, cross_check):
    assert p > 2 * N + 240 and all(p % k for k in range(2, math.isqrt(p) + 1))
    rows = []
    for a in range(N + 1):
        row = coefficients_by_recurrence(a, q, p)
        if cross_check:
            assert row == coefficients_by_products(a, q, p)
        rows.append(row)
    largest = q * (31 + 2)
    harmonic = [0] * (largest + 1)
    beta = [0] * (largest + 1)
    constant = 0
    for a, row in enumerate(rows):
        for i, partial in enumerate(row, 1):
            assert rows[N - a][i - 1] == (-1)**i * partial % p
            s = i + 2 * q
            coefficient = (-1)**a * math.prod(range(i, s)) * partial % p
            beta[s] = (beta[s] + coefficient) % p
            constant = (constant - coefficient * harmonic[s]) % p
        inverse = pow(2 * a + 1, -1, p)
        for s in range(2 * q + 1, largest + 1):
            harmonic[s] = (harmonic[s] + (-1)**a * pow(inverse, s, p)) % p
    assert all(beta[s] == 0 for s in range(1, largest + 1, 2))
    assert constant
    if cross_check:
        assert constant == constant_by_evaluation_address(rows, q, p)
    if q == 1:
        assert constant == {10007: 3738, 10009: 5333, 10037: 938}[p]
    return {"power": q, "prime": p, "seed_residue": constant,
            "coefficients": sum(map(len, rows))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cross-check", action="store_true")
    args = parser.parse_args()
    results = [seed(q, p, args.cross_check)
               for q in (1, 2, 3) for p in (10007, 10009, 10037)]
    print(json.dumps({"direct_products_checked": args.cross_check,
                      "seeds": results}, indent=2))
