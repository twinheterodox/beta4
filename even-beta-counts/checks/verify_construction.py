"""Exact checks supporting the uniform odd-index Gram construction.

The rank proof is analytic and uniform. These computations check identities
and their implementation; they do not replace the proof or constitute an audit.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, prod, isqrt
from pathlib import Path
import json
from base import Parameters, decay_certificate, modq, chi
from gram import (gram_basis, local_multiplier, seed_matrix, moment_matrix,
                  mul_linear, poly_eval)


CASES = [(1, 64, 7), (2, 160, 7), (3, 364, 7), (4, 588, 7),
         (5, 1030, 8), (6, 1296, 8), (7, 1992, 7), (8, 2502, 7),
         (9, 3260, 7), (10, 3828, 8)]


def determinant_Z(matrix):
    a = [row[:] for row in matrix]
    if not a:
        return 1
    sign, previous = 1, 1
    for k in range(len(a)-1):
        pivot = next((j for j in range(k, len(a)) if a[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = a[i][j]*a[k][k]-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
        previous = a[k][k]
        for i in range(k+1, len(a)):
            a[i][k] = 0
    return sign*a[-1][-1]


def exact_gram_check(P, modular):
    J, scalar, info = moment_matrix(P)
    p = modular['prime']
    assert [[x*modq(scalar, p) % p for x in row] for row in J] == modular['unnormalized_matrix']
    leading = [determinant_Z([row[:k] for row in J[:k]]) for k in range(1, P.B+1)]
    assert all(d > 0 for d in leading)
    c = [(-1)**(P.B-1+j)*determinant_Z(
        [[row[k] for k in range(P.B) if k != j] for row in J[:-1]])
        for j in range(P.B)]
    assert [sum(a*b for a, b in zip(row, c)) for row in J] == [0]*(P.B-1)+[leading[-1]]
    C1 = 8**P.S*P.K
    normalized_scalar = scalar*C1
    assert all((4**P.S*normalized_scalar*x).denominator == 1 for row in J for x in row)
    return {**info, 'matrix_entries_checked_exactly': P.B**2,
            'positive_leading_principal_minors': len(leading),
            'leading_principal_minor_bits': [d.bit_length() for d in leading],
            'integer_adjugate_weights_constructed': True,
            'adjugate_seed_identity_exact': True,
            'all_seed_constant_denominators_cleared_by_4_power_S': True}


def mul_truncated(a, b, degree):
    out = [0]*min(len(a)+len(b)-1, degree+1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:len(out)-i]):
            out[i+j] += x*y
    return out


def exact_base_coefficients(P, n):
    lcm = 1
    for k in range(1, n+1):
        lcm = lcm*k//gcd(lcm, k)
    delta = 4*lcm
    C = F(4**(P.S*n)*P.K**n*factorial(n+1)**P.S, factorial(n)**(P.D*P.L))
    scale = prod(theta.denominator**(P.L*n) for theta in P.theta)
    result = {}
    for h in range(n+1):
        poly = [-4*h-P.r*n, delta]
        for theta in P.theta:
            a, b = theta.numerator, theta.denominator
            for j in range(P.L*n):
                poly = mul_linear(poly, b*(-4*h-P.r*n+j)+a, b*delta)[:P.S]
        cs = [4*(k-h) for k in range(n+1) if k != h]
        initial = C/(scale*prod(cs)**P.S)
        for v in cs:
            assert delta % v == 0
            alpha = delta//v
            factor = [(-alpha)**j*comb(P.S+j-1, j) for j in range(P.S)]
            poly = mul_truncated(poly, factor, P.s)
        assert poly[0] == 0
        for i in range(1, P.S):
            cleared = initial*poly[P.S-i]
            assert cleared.denominator == 1
            result[i, h] = cleared.numerator
    return result, delta, C


def exact_modified_coefficients(P, n, poly, J, base, delta):
    result = {}
    for h in range(n+1):
        w = local_multiplier(poly, J, n, h)
        for i in range(1, P.S):
            result[i, h] = sum(c*delta**k*base.get((i+k, h), 0)
                               for k, c in enumerate(w))
    assert all(result[i, n-h] == (-1)**i*result[i, h]
               for h in range(n+1) for i in range(1, P.S))
    assert all(sum(result[i, h] for h in range(n+1)) == 0
               for i in range(1, P.S, 2))
    return result


def exact_constants(P, n, a, delta):
    q, groups = [], 0
    for b in range(1, P.B+1):
        total = 0
        for h in range(1, n+1):
            for k in range(1, 4*b*h+1, 2):
                value, power = a[P.s, h], k
                for i in range(P.s-1, -1, -1):
                    value = value*(delta*b)+(a[i, h] if i else 0)*power
                    power *= k
                assert value % k**P.s == 0
                total -= chi(k)*(value//k**P.s)
                groups += 1
        q.append(F(total, delta**P.S))
    return q, groups


def valuation(x, p):
    if not x:
        return float('inf')
    a, b = x.numerator, x.denominator
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def full_dilation_checks():
    """All coefficients and constants at two small parameter sets.

    These parameters check the identities, not the decay requirement.
    """
    records = []
    for B, s, r in ((1, 20, 3), (2, 24, 3)):
        P = Parameters(B, s, r)
        nodes, _, hs = gram_basis(B)
        J = len(nodes)-1
        tau = s-2*J-1
        base1, delta1, C1 = exact_base_coefficients(P, 1)
        seed_coeff, seed_q = [], []
        for hpoly in hs:
            a = exact_modified_coefficients(P, 1, hpoly, J, base1, delta1)
            q, _ = exact_constants(P, 1, a, delta1)
            seed_coeff.append(a)
            seed_q.append(q)
        gram_J, gram_scalar, _ = moment_matrix(P)
        assert [[q/C1 for q in column] for column in seed_q] == [
            [gram_scalar*gram_J[b][j] for b in range(B)] for j in range(B)]
        for p in (17, 19):
            assert all(p % d for d in range(2, isqrt(p)+1))
            assert p > max(4*B, P.L, 2*r+4)
            basep, deltap, _ = exact_base_coefficients(P, p)
            selected, other, constants, groups = 0, 0, 0, 0
            for j, hpoly in enumerate(hs):
                a = exact_modified_coefficients(P, p, hpoly, J, basep, deltap)
                for h in range(p+1):
                    for i in range(1, P.S):
                        value = F(a[i, h], deltap**(P.S-i))
                        if h % p == 0:
                            scaled = value*F(p)**(tau-i)
                            assert valuation(scaled, p) >= 0
                            target = F(seed_coeff[j][i, h//p], delta1**(P.S-i))
                            assert modq(scaled, p) == pow(2, -P.S, p)*modq(target, p) % p
                            selected += 1
                        else:
                            assert valuation(value, p) >= i
                            other += 1
                for i in range(1, P.S):
                    rho = F(sum(a[i, h] for h in range(p+1)), deltap**(P.S-i))
                    assert valuation(p**tau*rho, p) >= i
                q, count = exact_constants(P, p, a, deltap)
                groups += count
                for b in range(B):
                    scaled = p**tau*q[b]
                    assert valuation(scaled, p) >= 0
                    expected = chi(p)*pow(2, -P.S, p)*modq(seed_q[j][b], p) % p
                    assert modq(scaled, p) == expected
                    constants += 1
            records.append({'B': B, 's': s, 'r': r, 'prime': p, 'seed_index': 1,
                            'dilation_exponent': tau,
                            'selected_coefficient_congruences': selected,
                            'nonselected_coefficient_valuations': other,
                            'constant_congruences': constants,
                            'grouped_constant_divisibilities': groups,
                            'all_rho_valuation_bounds': True,
                            'scope': 'Exact full rational coefficient and constant checks; these small parameters are not theorem endpoint instances.'})
    return records


def uniform_endpoint_checks():
    e0 = F(68, 25)
    f66 = F(901)**4*F(351, 551)**66
    log_bound = F(901, 50)**10/F(163, 60)**29
    slope_bound = ((4*e0)**4*F(7)**7/F(11)**11)**5*e0**16
    assert f66 < 1 and log_bound < 1 and slope_bound < 1
    assert F(163, 60)+F(7, 4320) < e0
    return {'f_at_u_1_over_50_with_S_66D_lt_1': True,
            'log_18_02_lt_29_over_10': True,
            'decay_slope_gt_4_over_5': True,
            'e_upper_68_over_25': True,
            'three_rational_ratios_illustrative': [float(f66), float(log_bound), float(slope_bound)],
            'arithmetic': 'Exact rational comparisons; decimals are illustrative.'}


def main():
    if not __debug__:
        raise SystemExit('Run without -O; assertions perform the checks.')
    results = []
    for B, s, r in CASES:
        print(f'Checking the uniform seed at B={B}, s={s}...', flush=True)
        P = Parameters(B, s, r)
        modular = seed_matrix(P)
        gram = exact_gram_check(P, modular)
        results.append({'parameters': {'B': B, 's': s, 'r': r},
                        'decay': decay_certificate(P),
                        'seed': modular, 'exact_gram': gram})
    print('Checking complete prime dilations over Q at small parameters...', flush=True)
    dilation = full_dilation_checks()
    uniform = uniform_endpoint_checks()
    record = {'status': 'All exact checks passed; the uniform proof is in the manuscript. Not an independent audit or formal verification.',
              'cases': results, 'full_dilation_checks': dilation,
              'uniform_endpoint_checks': uniform}
    Path(__file__).with_name('construction-verification.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'cases': len(results),
                      'positive_leading_principal_minors': sum(z['exact_gram']['positive_leading_principal_minors'] for z in results),
                      'exact_gram_entries': sum(z['exact_gram']['matrix_entries_checked_exactly'] for z in results),
                      'full_dilation_checks': dilation}, indent=2))


if __name__ == '__main__':
    main()
