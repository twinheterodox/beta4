"""Exact checks for the real-rate investigation; not an independent audit.

Run with python3 -B, without -O. SymPy is used only for exact symbolic
kernel identities. All coefficient and logarithm certificates use Fraction.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, gcd, lcm
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE
PROOF_HASH = '0a4007b2a92e99ebb479c56511b836e66d2ab072d066aeac474995564f481d46'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def odd(a):
    return a >> ((a & -a).bit_length()-1)


def v2(a):
    return (a & -a).bit_length()-1


def multiplier(n, ell, N):
    J = 2*n+2
    R = odd(gcd(factorial(n), lcm(*range(1, J+1))))
    L = lcm(R, odd(lcm(*range(1, 2*N))))
    return 2**(4*n+2*ell-4+4*v2(factorial(n))) * L**4


def cubic(w, r):
    return (64*w**3+(64*r*r-16)*w*w
            -(160*r*r+20)*w+20*r*r-3)


def root_interval(r):
    # The proof establishes uniqueness above 1/4. The signs are exact.
    lo, hi = F(1, 4), F(4)
    assert cubic(lo, r) < 0 < cubic(hi, r)
    for _ in range(108):
        mid = (lo+hi)/2
        sign = cubic(mid, r)
        if sign == 0:
            return mid, mid
        if sign < 0:
            lo = mid
        else:
            hi = mid
    assert cubic(lo, r) < 0 < cubic(hi, r)
    return lo, hi


@lru_cache(None)
def elementary_log_interval(x):
    """Enclose log(x) for 1<=x<=2 by a positive atanh series."""
    assert 1 <= x <= 2
    t = (x-1)/(x+1)
    count = 48
    total = F(0)
    power = t
    for j in range(count):
        total += 2*power/(2*j+1)
        power *= t*t
    tail = 2*power/((2*count+1)*(1-t*t))
    return total, total+tail


def log_interval(x):
    assert x > 0
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1
    lo, hi = elementary_log_interval(x)
    l2, h2 = elementary_log_interval(F(2))
    return (lo+k*(l2 if k >= 0 else h2),
            hi+k*(h2 if k >= 0 else l2))


def outward(interval, places=24):
    scale = 10**places
    lo, hi = interval
    floor = (lo*scale).numerator//(lo*scale).denominator
    ceil = -((-hi*scale).numerator//(-hi*scale).denominator)
    return [str(F(floor, scale)), str(F(ceil, scale))]


def rate_interval(r, roots):
    lo = hi = F(0)
    terms = [(F(5, 2), F(1, 4)),
             (-(r+F(1, 2)), (r+F(1, 2))**2),
             (r-F(1, 2), (r-F(1, 2))**2)]
    for c, shift in terms:
        if not c:
            continue
        low = log_interval(roots[0]+shift)[0]
        high = log_interval(roots[1]+shift)[1]
        lo += c*(low if c > 0 else high)
        hi += c*(high if c > 0 else low)
    return lo, hi


def main():
    rows = json.loads((PREVIOUS/'exploration.json').read_text())['rows']
    assert len(rows) == 72

    # Check the closed adjoint kernel as an identity in TWO variables.
    import sympy as sp
    t, z = sp.symbols('t z')
    def K(j):
        return (sp.factorial(j)*(2*t-1)
                /sp.rf(z-t-sp.Rational(j, 2), j+1)
                /sp.rf(z+t-1-sp.Rational(j, 2), j+1))
    kernels = 0
    for j in range(7):
        expr = (K(j).subs(t, t+sp.Rational(1, 2))/(2*t)
                - K(j).subs(t, t-sp.Rational(1, 2))/(2*t-2)
                - K(j+1))
        assert sp.cancel(expr) == 0
        kernels += 1

    # The coefficient pairs were constructed exactly in the preserved package.
    clear_checks = divisibility_checks = 0
    for row in rows:
        n, ell, N = (row[k] for k in ('n', 'ell', 'N'))
        M = multiplier(n, ell, N)
        for name in ('constant', 'coefficient'):
            assert (M*F(row[name])).denominator == 1
            clear_checks += 1
        J = 2*n+2
        old_M = 2**(4*J+1)*lcm(*range(1, max(J, 2*N-1)+1))**4
        assert old_M % M == 0
        divisibility_checks += 1

    # Check the new integer-polynomial description on a small full set.
    sys.path.insert(0, str(PREVIOUS))
    import explore
    polynomial_checks = 0
    for n in (2, 4, 6):
        for ell in (0, 1, 2):
            P, _, _ = explore.construct(n, ell)
            assert all((factorial(n)**4*c).denominator == 1 for c in P)
            polynomial_checks += 1
    moments = explore.moments(4, 120)
    assert all((2**(j-3)*moments[j]).denominator == 1 for j in range(3, 121))

    certificates = []
    l2, h2 = log_interval(F(2))
    for r in (F(1, 2), F(1), F(2)):
        roots = root_interval(r)
        rate = rate_interval(r, roots)
        gap = (rate[0]+8*l2+8*r, rate[1]+8*h2+8*r)
        assert rate[1] < 0 < gap[0]
        assert rate[1]-rate[0] < F(1, 10**27)
        assert gap[1]-gap[0] < F(1, 10**27)
        certificates.append(dict(
            r=str(r), w_enclosure=[str(x) for x in roots],
            rate_enclosure=outward(rate), gap_enclosure=outward(gap)))
        if r == F(1, 2):
            assert gap[0] > F(933, 100)
        if r == 2:
            assert cubic(F(7, 4), r) == 0
            assert roots[0] <= F(7, 4) <= roots[1]
            # In (3.7), the three logarithm arguments are 2, 8 and 4.
            assert F(5, 2)-F(5, 2)*3+F(3, 2)*2 == -2
            assert rate[0] <= -2*l2 and -2*h2 <= rate[1]

    result = dict(
        status='passed; constructing-assistant exact checks, not an independent audit',
        counts=dict(symbolic_kernel_identities=kernels,
                    rational_coefficient_integrality=clear_checks,
                    old_multiplier_divisible_by_new=divisibility_checks,
                    integer_polynomial_identities=polynomial_checks,
                    dyadic_moment_bounds=118,
                    saddle_and_rate_certificates=len(certificates),
                    exact_r_two_rate_identity=1),
        rate_certificates=certificates,
        previous_proof_sha256=PROOF_HASH,
        previous_coefficients_sha256=sha(PREVIOUS/'exploration.json'),
        script_sha256=sha(Path(__file__)))
    (HERE/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['counts'], indent=2))
    for c in certificates:
        print('r='+c['r'], 'rate', c['rate_enclosure'], 'gap', c['gap_enclosure'])


if __name__ == '__main__':
    main()
