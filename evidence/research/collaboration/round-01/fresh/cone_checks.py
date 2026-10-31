"""Exact checks for the scalar beta(4) positive-polynomial cone.

Run from any directory with beta4/.venv/bin/python.  This file writes only
beside itself.  The mathematical proofs are in findings.txt; the checks
are finite consistency checks by the same investigator.
"""
from pathlib import Path
from math import gcd, lcm
from functools import lru_cache
import json
import random
import sys

from flint import arb, ctx, fmpq, fmpq_mat, fmpq_poly

HERE = Path(__file__).resolve().parent
t = fmpq_poly([0, 1])
one = fmpq_poly([1])


def nu(poly):
    return sum((c / (2*k+1)**4 for k, c in enumerate(poly)), fmpq(0))


def pair(poly):
    """Formal coefficients, independent of any irrationality assumption."""
    a = poly(-1)
    quotient, remainder = divmod(poly - a, one + t)
    assert not remainder
    return a, nu(quotient)


def extremizer(w, n):
    """Minimize I(w*p^2) subject to deg(p)<=n and p(-1)=1."""
    c = w(-1)
    assert c
    _, b0 = pair(w)
    if n:
        v = fmpq_mat([[nu(w*t**i)] for i in range(n)])
        gram = fmpq_mat([[nu(w*(one+t)*t**(i+j))
                          for j in range(n)] for i in range(n)])
        h = gram.solve(v)
        p = one - (one+t)*fmpq_poly([h[i, 0] for i in range(n)])
        b = b0 - sum((v[i, 0]*h[i, 0] for i in range(n)), fmpq(0))
        for j in range(n):
            assert nu(w*p*t**j) == 0
    else:
        p, b = one, b0
    assert p(-1) == 1
    assert pair(w*p*p) == (c, b)
    r = -b/c
    canonical = w*p*p/abs(c)
    assert pair(canonical) == (fmpq(1 if c > 0 else -1),
                               -r if c > 0 else r)
    return r, canonical, p


@lru_cache(maxsize=None)
def cone(degree):
    assert degree >= 1
    n = degree//2
    if degree % 2:
        low = extremizer(one-t, n)
        high = extremizer(t, n)
    else:
        low = extremizer(one, n)
        high = extremizer(t*(one-t), n-1)
    assert low[0] < high[0]
    return low, high


def chebyshev_shifted(n):
    if n == 0:
        return one
    a, b = one, 2*t-one
    for _ in range(1, n):
        a, b = b, 2*(2*t-one)*b-a
    return b


def primitive(a, b):
    d = lcm(int(a.q), int(b.q))
    aa, bb = int(a*d), int(b*d)
    g = gcd(aa, bb)
    return aa//g, bb//g


def beta_ball(bits):
    with ctx.workprec(bits):
        return (arb(4).zeta(arb(1)/4)-arb(4).zeta(arb(3)/4))/256


def check_selected_endpoints():
    data = []
    for degree in [1, 2, 3, 4, 5, 6, 8, 12, 16]:
        low, high = cone(degree)
        rec = dict(degree=degree, lower=str(low[0]), upper=str(high[0]),
                   lower_polynomial=list(map(str, low[1])),
                   upper_polynomial=list(map(str, high[1])))
        bits = max(int(low[0].q).bit_length(), int(high[0].q).bit_length())+512
        with ctx.workprec(bits):
            beta = beta_ball(bits)
            assert arb(low[0]) < beta < arb(high[0])
            for name, item, sign in [('lower', low, 1), ('upper', high, -1)]:
                a, b = primitive(*pair(item[1]))
                value = a*beta+b
                assert value > 0
                rec[name+'_primitive'] = [str(a), str(b)]
                rec[name+'_coefficient_bits'] = max(abs(a).bit_length(), abs(b).bit_length())
                rec[name+'_value'] = value.str(24)
                rec[name+'_log_value'] = value.log().str(24)
        data.append(rec)
    assert data[0]['lower'] == '1/2'
    assert data[0]['upper'] == '1'
    assert data[1]['lower'] == '81/82'
    assert data[1]['upper'] == '161/162'
    return data


def check_sos_cone():
    rng = random.Random(20261028)
    count = 0
    for degree in range(1, 13):
        (lo, pp, _), (hi, pm, _) = cone(degree)
        if degree > 1:
            previous = cone(degree-1)
            assert previous[0][0] <= lo and hi <= previous[1][0]
        n = degree//2
        for _ in range(16):
            p = fmpq_poly([rng.randrange(-8, 9) for _ in range(n+1)])
            q_degree = n if degree % 2 else n-1
            q = fmpq_poly([rng.randrange(-8, 9) for _ in range(q_degree+1)])
            if degree % 2:
                poly = t*p*p+(one-t)*q*q
            else:
                poly = p*p+t*(one-t)*q*q
            a, b = pair(poly)
            assert b+lo*a >= 0
            assert b+hi*a >= 0
            lam = (b+hi*a)/(hi-lo)
            eta = (b+lo*a)/(hi-lo)
            assert lam >= 0 and eta >= 0
            reconstructed = lam*pp+eta*pm
            assert pair(reconstructed) == (a, b)
            assert reconstructed.degree() <= degree
            count += 1
        # Rational points in the cone, including both rays and A=0.
        for lam, eta in [(0, 1), (1, 0), (1, 1), (fmpq(2, 7), fmpq(5, 9))]:
            a, b = pair(lam*pp+eta*pm)
            assert a == lam-eta and b == -lam*lo+eta*hi
            count += 1
    return count


def check_witness_lifts():
    records = []
    # All signs and A=0.  Positivity is independently enclosed by Arb.
    targets = [(0, 1), (1, 0), (-1, 1), (2, -1), (82, -81)]
    # Obtain some small, certified forms from ordinary continued fractions.
    with ctx.workprec(2048):
        beta = beta_ball(2048)
        x = beta
        p0, p1, q0, q1 = 0, 1, 1, 0
        for _ in range(10):
            ai = int(x.floor().unique_fmpz())
            p0, p1 = p1, ai*p1+p0
            q0, q1 = q1, ai*q1+q0
            err = q1*beta-p1
            if err > 0:
                targets.append((q1, -p1))
            elif err < 0:
                targets.append((-q1, p1))
            x = 1/(x-ai)
        targets = list(dict.fromkeys(targets))
        for a, b in targets:
            value = a*beta+b
            assert value > 0
            if a == 0:
                poly, n, multiplier = b*(one+t), 0, 1
            else:
                n = 0
                while True:
                    cheb = chebyshev_shifted(n)
                    norm = abs(int(cheb(-1)))
                    if arb(abs(a))/(norm*norm) < value/2:
                        break
                    n += 1
                q = cheb/cheb(-1)
                assert sum(abs(c) for c in q) == 1
                base = abs(a)*q*q*(one if a > 0 else t)
                aa, bb = pair(base)
                assert aa == a
                delta = fmpq(b)-bb
                assert delta > 0
                poly = base+delta*(one+t)
                odd_lcm = 1
                for k in range(poly.degree()):
                    odd_lcm = lcm(odd_lcm, 2*k+1)
                multiplier = norm*norm*odd_lcm**4
            assert pair(poly) == (fmpq(a), fmpq(b))
            integer_poly = multiplier*poly
            assert all(c.q == 1 for c in integer_poly)
            assert primitive(*pair(integer_poly)) == primitive(fmpq(a), fmpq(b))
            max_coefficient = max(abs(int(c)) for c in integer_poly)
            assert max_coefficient <= 4*max(abs(a), abs(b))*multiplier
            if b == 0 and a >= 0:
                minimum_degree = 0
            else:
                minimum_degree = 1
                while True:
                    low, high = cone(minimum_degree)
                    if b+low[0]*a >= 0 and b+high[0]*a >= 0:
                        break
                    minimum_degree += 1
            assert minimum_degree <= poly.degree()
            records.append(dict(a=str(a), b=str(b), chebyshev_n=n,
                                degree=poly.degree(), minimum_degree=minimum_degree,
                                error=value.str(24),
                                multiplier=str(multiplier),
                                integer_polynomial=list(map(str, integer_poly))))
    return records


def main():
    sys.set_int_max_str_digits(0)
    endpoints = check_selected_endpoints()
    checks = check_sos_cone()
    lifts = check_witness_lifts()
    last = lifts[-1]
    a, b, degree = int(last['a']), int(last['b']), last['minimum_degree']
    (lo, pp, _), (hi, pm, _) = cone(degree)
    lam, eta = (b+hi*a)/(hi-lo), (b+lo*a)/(hi-lo)
    assert lam >= 0 and eta >= 0
    witness = lam*pp+eta*pm
    assert pair(witness) == (fmpq(a), fmpq(b))
    scale = 1
    for coefficient in witness:
        scale = lcm(scale, int(coefficient.q))
    integral_witness = scale*witness
    assert all(c.q == 1 for c in integral_witness)
    assert primitive(*pair(integral_witness)) == (a, b)
    example = dict(a=str(a), b=str(b), degree=degree, lower_multiplier=str(lam),
                   upper_multiplier=str(eta), clearing_multiplier=str(scale),
                   integer_polynomial=list(map(str, integral_witness)))
    result = dict(status='Proved cone/degree statements in findings.txt; finite exact consistency checks here. No irrationality proof.',
                  endpoints=endpoints, sos_and_cone_checks=checks,
                  integer_witnesses=lifts, interior_cone_example=example)
    (HERE/'cone_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(endpoint_degrees=[r['degree'] for r in endpoints],
                          sos_and_cone_checks=checks, integer_witnesses=len(lifts)), indent=2))
    for r in endpoints:
        print(json.dumps({k: r[k] for k in ['degree', 'lower_coefficient_bits',
                                           'upper_coefficient_bits', 'lower_log_value',
                                           'upper_log_value']}))


if __name__ == '__main__':
    main()
