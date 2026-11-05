"""Exact finite-field partial fractions for a three-node beta construction.

No normalization is included: multiplying the entire rational function by
a nonzero rational scalar does not affect whether a seed vanishes over Q.
The prime must exceed all relevant seed factors and the series degree.
"""
from fractions import Fraction
from math import comb, prod
from itertools import combinations
import argparse
import json

THETA = (Fraction(1), Fraction(1, 2), Fraction(1, 3), Fraction(2, 3))


def constants(n, r, h=0):
    return [Fraction(-4*h-r*n)] + [Fraction(-4*h-r*n+j)+theta
            for theta in THETA for j in range((2*r+4)*n)]


def modq(x, p):
    return x.numerator * pow(x.denominator, -1, p) % p


def partials_newton(n, s, r, p):
    """Recover the regular part at each pole by its log derivative."""
    out = {}
    for h in range(n+1):
        cs = constants(n, r, h)
        assert cs.count(0) == 1
        num = [modq(c, p) for c in cs if c]
        den = [4*(k-h) % p for k in range(n+1) if k != h]
        assert all(num) and all(den)
        initial = prod(num) % p * pow(prod(den) % p, -(s+1), p) % p
        powers = [0]*s
        for factors, multiplicity in ((num, 1), (den, -(s+1))):
            for c in factors:
                inv, v = pow(c, -1, p), 1
                for j in range(1, s):
                    v = v*inv % p
                    powers[j] = (powers[j] + multiplicity*v) % p
        series = [initial]
        for j in range(1, s):
            value = sum((1 if k % 2 else -1)*powers[k]*series[j-k]
                        for k in range(1, j+1))
            series.append(value * pow(j, -1, p) % p)
        for i in range(1, s+1):
            out[i, h] = series[s-i]
    return out


def mul(a, b, degree, p):
    c = [0]*min(len(a)+len(b)-1, degree+1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:len(c)-i]):
            c[i+j] = (c[i+j]+x*y) % p
    return c


def partials_product(n, s, r, p):
    """A second algorithm: multiply truncated numerator/denominator series."""
    out = {}
    for h in range(n+1):
        cs = constants(n, r, h)
        assert cs.count(0) == 1
        v = [1]
        for c in cs:
            if c:
                cm = modq(c, p)
                assert cm
                v = mul(v, [cm, 1], s-1, p)
        for k in range(n+1):
            if k == h:
                continue
            inv = pow(4*(k-h), -1, p)
            d = [(-1)**j * comb(s+j, j) * pow(inv, s+1+j, p) % p
                 for j in range(s)]
            v = mul(v, d, s-1, p)
        for i in range(1, s+1):
            out[i, h] = v[s-i]
    return out


def chi(k):
    return 0 if k % 2 == 0 else (1 if k % 4 == 1 else -1)


def constant_term(a, n, s, b, p):
    H, answer, last = [0]*(s+1), 0, 0
    bs = [pow(b, i, p) for i in range(s+1)]
    for h in range(n+1):
        end = 4*b*h
        for k in range(last+1, end+1):
            ch = chi(k)
            if ch:
                inv, v = pow(k, -1, p), 1
                for i in range(1, s+1):
                    v = v*inv % p
                    H[i] = (H[i]+ch*v) % p
        answer = (answer-sum(a[i, h]*bs[i]*H[i]
                            for i in range(1, s+1))) % p
        last = end
    return answer


def evaluate(t, n, s, r, p):
    num = prod((t+modq(c, p)) % p for c in constants(n, r)) % p
    den = prod(pow(t+4*k, s+1, p) for k in range(n+1)) % p
    return num*pow(den, -1, p) % p


def run(n=2, s=364, r=7, p=1000003, crosscheck=False, all_pairs=False):
    assert n % 2 == 0 and s % 2 == 0 and s+1 > 4*(2*r+4)
    assert all(p % d for d in range(2, int(p**0.5)+1))
    a = partials_newton(n, s, r, p)
    if crosscheck:
        assert a == partials_product(n, s, r, p)
    for t in (7, 251, 503):
        assert evaluate(t, n, s, r, p) == sum(
            value*pow(t+4*h, -i, p) for (i, h), value in a.items()) % p
    assert all(sum(a[i, h] for h in range(n+1)) % p == 0
               for i in range(1, s+1, 2))
    q = {b: constant_term(a, n, s, b, p) for b in (1, 2, 3)}
    result = {'m': n, 's': s, 'r': r, 'prime': p,
              'q_unnormalized_mod_prime': q,
              'low_weight_elimination_residue': (15*q[1]-6*q[2]+q[3]) % p,
              'crosschecked_all_coefficients': crosscheck,
              'reconstruction_checks': 3, 'odd_rho_vanish': True}
    if all_pairs:
        zero_pairs = []
        for i, j in combinations(range(2, s+1, 2), 2):
            # det of [q, b^i, b^j], rows indexed by b=1,2,3.
            det = (q[1]*(pow(2, i, p)*pow(3, j, p)-pow(2, j, p)*pow(3, i, p))
                   - q[2]*(pow(3, j, p)-pow(3, i, p))
                   + q[3]*(pow(2, j, p)-pow(2, i, p))) % p
            if not det:
                zero_pairs.append((i, j))
        result['elimination_pairs_checked'] = comb(s//2, 2)
        result['zero_determinant_pairs'] = zero_pairs
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--m', type=int, default=2)
    parser.add_argument('--s', type=int, default=364)
    parser.add_argument('--r', type=int, default=7)
    parser.add_argument('--prime', type=int, default=1000003)
    parser.add_argument('--crosscheck', action='store_true')
    parser.add_argument('--all-pairs', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.m, args.s, args.r, args.prime,
                         args.crosscheck, args.all_pairs), indent=2))
