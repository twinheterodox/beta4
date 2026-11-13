#!/usr/bin/env python3
"""Fixed exact checks for the monic-substitution sections; no numerical determinant search.

Run: python3 checks/monic.py --output review/monic-checks.json
Only Python's standard library is used. All identity checks use Fraction
and integers. Printed transcendental constants are explanatory decimals;
their required sign is proved by rational inequalities in the manuscript.
"""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from functools import reduce
from itertools import combinations, permutations
from math import atan, comb, factorial, gcd, lcm, log, pi, prod
import json
from pathlib import Path


MAX_K = 8


def valuation2(x):
    x = Q(x)
    if not x:
        raise AssertionError("valuation requested for zero")
    n, d = abs(x.numerator), x.denominator
    return (n & -n).bit_length() - (d & -d).bit_length()


def det(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    result = Q(1)
    for i in range(len(a)):
        pivot = next((r for r in range(i, len(a)) if a[r][i]), None)
        if pivot is None:
            return Q(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        p = a[i][i]
        result *= p
        for r in range(i + 1, len(a)):
            t = a[r][i] / p
            for c in range(i + 1, len(a)):
                a[r][c] -= t * a[i][c]
            a[r][i] = Q(0)
    return result


def principal(a, indices):
    return [[a[i][j] for j in indices] for i in indices]


def peval(p, x):
    y = 0
    for c in reversed(p):
        y = y * x + c
    return y


def orthogonal_polynomials(n):
    ps = [[1], [0, 1]]
    for j in range(1, n + 1):
        nxt = [0] + ps[-1]
        a = j * (n - j + 1)
        for i, c in enumerate(ps[-2]):
            nxt[i] -= a * c
        ps.append(nxt)
    return ps


def polynomial_add(a, b):
    out = dict(a)
    for k, x in b.items():
        out[k] = out.get(k, Q(0)) + x
        if not out[k]:
            del out[k]
    return out


def polynomial_mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = (i + k, j + l)
            out[key] = out.get(key, Q(0)) + x * y
    return {k: v for k, v in out.items() if v}


def polynomial_det(a):
    d = len(a)
    out = {}
    for perm in permutations(range(d)):
        inversions = sum(perm[i] > perm[j] for i in range(d) for j in range(i + 1, d))
        term = {(0, 0): Q((-1) ** inversions)}
        for i, j in enumerate(perm):
            term = polynomial_mul(term, a[i][j])
        out = polynomial_add(out, term)
    return out


def matrices(k):
    nodes = [4 * m * m for m in range(1, k + 1)]
    odds = list(range(1, 2 * k, 2))
    e = [Q(prod(nodes[j] - nodes[i] for j in range(k) if j != i),
           prod(l * l - nodes[i] for l in odds)) for i in range(k)]
    rho = [Q(prod(a - l * l for a in nodes),
             prod(h * h - l * l for h in odds if h != l)) for l in odds]
    beta_weights = [Q((-1) ** j * l ** 3, 2) for j, l in enumerate(odds)]
    c = [Q(7 * (-1) ** (m - 1) * m ** 3, 32) for m in range(1, k + 1)]
    re = [-Q((-1) ** (m - 1) * m ** 3, 4)
          * sum((Q((-1) ** (v - 1), v ** 4) for v in range(1, m)), Q(0))
          - Q(1, 8 * m) for m in range(1, k + 1)]
    ro = [-beta_weights[j]
          * sum((Q((-1) ** v, (2 * v + 1) ** 4) for v in range(j)), Q(0))
          - Q(1, 4 * l) for j, l in enumerate(odds)]
    b = [[sum((rho[t] * beta_weights[t] /
               ((nodes[i] - l * l) * (nodes[j] - l * l))
               for t, l in enumerate(odds)), Q(0))
          for j in range(k)] for i in range(k)]
    r = [[sum((rho[t] * ro[t] /
               ((nodes[i] - l * l) * (nodes[j] - l * l))
               for t, l in enumerate(odds)), Q(0))
          + (e[i] * re[i] if i == j else Q(0))
          for j in range(k)] for i in range(k)]
    return nodes, odds, e, rho, b, [e[i] * c[i] for i in range(k)], r


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = Counter()
    subset_results = []
    monic_results = []

    for k in range(1, MAX_K + 1):
        nodes, odds, e, rho, b, c, r = matrices(k)
        n = 2 * k - 1
        support = [(n - 2 * a, Q(comb(n, a), 2 ** n)) for a in range(n + 1)]
        ps = orthogonal_polynomials(n)

        # Direct residue formula versus finite binomial expectation.
        for i in range(k):
            for j in range(i, k):
                avg = sum((prob * Q(ell ** 4 * prod(a - ell * ell for a in nodes),
                                       (nodes[i] - ell * ell) * (nodes[j] - ell * ell))
                           for ell, prob in support), Q(0)) / (2 * factorial(n))
                assert b[i][j] == avg
                counts["residue_binomial_entries"] += 1
                assert valuation2(2 * factorial(n) * b[i][j]) == 0
                counts["scaled_entry_units"] += 1

                for u in (Q(0), Q(1), Q(7, 3)):
                    direct = Q(prod(u + a for a in nodes), prod(u + l * l for l in odds))
                    direct /= (u + nodes[i]) * (u + nodes[j])
                    residues = (e[i] / (u + nodes[i]) if i == j else Q(0))
                    residues += sum((rho[t] / ((nodes[i] - l * l) *
                                                (nodes[j] - l * l) * (u + l * l))
                                     for t, l in enumerate(odds)), Q(0))
                    assert direct == residues
                    counts["partial_fraction_reconstructions"] += 1

        # Recurrence-generated polynomials versus the literal generating function.
        for ell, prob in support:
            plus, minus = (n + ell) // 2, (n - ell) // 2
            for j in range(n + 2):
                generated = factorial(j) * sum(
                    (-1) ** (j - h) * comb(plus, h) * comb(minus, j - h)
                    for h in range(j + 1)
                    if h <= plus and 0 <= j - h <= minus)
                assert peval(ps[j], ell) == generated
                counts["generating_function_values"] += 1
            assert peval(ps[n + 1], ell) == 0
            counts["finite_support_closures"] += 1

        for i in range(n + 1):
            for j in range(i, n + 1):
                inner = sum((prob * peval(ps[i], ell) * peval(ps[j], ell)
                             for ell, prob in support), Q(0))
                expected = factorial(i) * factorial(n) // factorial(n - i) if i == j else 0
                assert inner == expected
                counts["orthogonal_inner_products"] += 1

        for j in range(k):
            diag = (4 * j + 1) * n - 8 * j * j
            lower = (2 * j) * (2 * j - 1) * (n - 2 * j + 1) * (n - 2 * j + 2)
            assert diag % 2 == 1 and lower % 4 == 0
            for ell, prob in support:
                rhs = diag * peval(ps[2 * j], ell)
                if j + 1 < k:
                    rhs += peval(ps[2 * j + 2], ell)
                if j:
                    rhs += lower * peval(ps[2 * j - 2], ell)
                assert ell ** 2 * peval(ps[2 * j], ell) == rhs
                counts["multiplication_matrix_values"] += 1

        # All nonempty principal subsets in the fixed range.
        for d in range(1, k + 1):
            for s in combinations(range(k), d):
                actual = det(principal(b, s))
                assert actual
                delta = prod(nodes[j] - nodes[i] for i, j in combinations(s, 2))
                expected_v = 2 * valuation2(delta) - d + sum(
                    valuation2(factorial(2 * j)) - valuation2(factorial(n - 2 * j))
                    for j in range(d))
                assert valuation2(actual) == expected_v
                counts["all_subset_exact_valuations"] += 1
                subset_results.append({"K": k, "S": [i + 1 for i in s], "v2_det_B": expected_v})

        full_ratio = det(b) / prod(c)
        expected_full = Q((-1) ** k * 2 ** k * comb(2 * k, k) ** 3, 7 ** k)
        assert full_ratio == expected_full
        counts["full_coordinate_determinant_identity"] += 1
        for i, m in enumerate(range(1, k + 1)):
            closed = Q(7 * (-1) ** m * 2 ** (4 * k) * m,
                       256 * comb(2 * k - 2 * m, k - m) * comb(2 * k + 2 * m, k + m))
            assert c[i] == closed
            counts["zeta_coefficient_closed_forms"] += 1

        s = tuple(range(1, k, 2))
        if s:
            d = len(s)
            delta = prod(nodes[j] - nodes[i] for i, j in combinations(s, 2))
            delta_closed = Q(2 ** (2 * d * (d - 1))) * prod(
                Q(factorial(2 * j - 1), j) for j in range(1, d + 1))
            assert delta == delta_closed
            counts["even_subset_vandermonde"] += 1

            # A direct weighted moment Gram independently checks the unit factor.
            norms = [factorial(2 * j) * factorial(n) // factorial(n - 2 * j)
                     for j in range(d)]
            weighted = [[sum((prob * Q(ell ** 4 * prod(a - ell * ell for a in nodes),
                                       prod(nodes[m] - ell * ell for m in s) ** 2)
                              * peval(ps[2 * i], ell) * peval(ps[2 * j], ell)
                              for ell, prob in support), Q(0))
                         for j in range(d)] for i in range(d)]
            factor = [[weighted[i][j] / norms[i] for j in range(d)] for i in range(d)]
            assert all(x.denominator % 2 for row in factor for x in row)
            assert valuation2(det(factor)) == 0
            counts["selected_weighted_unit_determinants"] += 1

            if k <= 6:
                pencil = [[{(0, 0): r[i][j] / c[i],
                            (1, 0): b[i][j] / c[i],
                            **({(0, 1): Q(1)} if i == j else {})}
                           for j in s] for i in s]
                polynomial = polynomial_det(pencil)
                assert polynomial[(0, d)] == 1
                top = det(principal(b, s)) / prod(c[i] for i in s)
                assert polynomial[(d, 0)] == top
                clearer = lcm(*(a.denominator for a in polynomial.values()))
                integers = [int(clearer * a) for a in polynomial.values()]
                assert reduce(gcd, integers) == 1
                assert valuation2(clearer) >= max(0, -valuation2(top))
                counts["complete_small_monic_polynomials"] += 1
                monic_results.append({"K": k, "dimension": d,
                    "clearer": str(clearer), "top_X_coefficient": str(top),
                    "coefficients": {f"X^{i} Z^{j}": str(a)
                                     for (i, j), a in sorted(polynomial.items())}})

        if k == 2:
            assert b[1][1] == Q(-167, 140) and c[1] == Q(1, 5)
            counts["worked_example"] += 1

    # This exact rational comparison establishes the sign without decimal reliance.
    assert Q(1, 3) + Q(11, 14) == Q(47, 42)
    assert Q(2, 3) - Q(47, 84) == Q(3, 28) > 0
    counts["rational_sign_comparisons"] += 2
    dstar = -log(2) + log(2.5) / 2 + atan(.5) + atan(2) / 4 + pi / 8
    result = {
        "status": "PASS",
        "fixed_range_K": [1, MAX_K],
        "method": "Python standard-library integers and Fraction; no imported checker",
        "counts": dict(sorted(counts.items())),
        "total_exact_checks": sum(counts.values()),
        "constants_illustrative_only": {
            "D_star": dstar, "D_star_over_two": dstar / 2,
            "necessary_denominator_rate": log(2), "gamma": log(2) - dstar / 2,
            "rigorous_positive_lower_bound_for_gamma": "3/28"},
        "all_subset_results": subset_results,
        "small_monic_polynomials": monic_results,
        "scope": "Finite identity checks support, but do not prove, the general theorems."
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(json.dumps({k: result[k] for k in ("status", "counts", "total_exact_checks",
                                            "constants_illustrative_only")}, indent=2))


if __name__ == "__main__":
    main()
