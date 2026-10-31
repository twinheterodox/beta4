"""Exact checks for the Euler--Boole obstruction; no numerical search.

Run with system Python from beta4. Writes only beside this script.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import comb, factorial, gcd, isqrt
import json
import sympy as sp

OUT = Path(__file__).resolve().parent


def euler_numbers(count):
    e = [1]
    for n in range(1, count):
        e.append(-sum(comb(2*n, 2*j)*e[j] for j in range(n)))
    return e


def inverse_cosh_series(count):
    c = [Q(1)]
    for n in range(1, count):
        c.append(-sum((c[n-j]/factorial(2*j) for j in range(1, n+1)), Q(0)))
    return c


def prefix(n):
    return sum((Q((-1)**k, (2*k+1)**4) for k in range(n)), Q(0))


def correction(n, m, e):
    return sum((Q(e[j]*(2*j+1)*(2*j+2)*(2*j+3),
                  12*(2*n)**(2*j+4)) for j in range(m)), Q(0))


def approx(n, m, e):
    return prefix(n)+(-1)**n*correction(n, m, e)


def vp(x, p):
    x = Q(x)
    assert x
    a, b, out = x.numerator, x.denominator, 0
    while a % p == 0:
        a //= p
        out += 1
    while b % p == 0:
        b //= p
        out -= 1
    return out


def residue(x, p):
    x = Q(x)
    assert x.denominator % p
    return x.numerator*pow(x.denominator, -1, p) % p


def primes_below(n):
    return [p for p in range(2, n)
            if all(p % d for d in range(2, isqrt(p)+1))]


def main():
    e = euler_numbers(160)
    inverse = inverse_cosh_series(160)
    for j in range(160):
        assert inverse[j]*factorial(2*j) == e[j]
        assert sp.euler(2*j) == e[j]
    x = sp.Symbol('x')
    small = sp.series(1/sp.cosh(x), x, 0, 26).removeO()
    for j in range(13):
        assert small.coeff(x, 2*j) == sp.Rational(e[j], factorial(2*j))

    # Check the finite geometric identity without trigonometric constants.
    t, a = sp.symbols('t a', positive=True)
    for m in range(1, 9):
        finite = sum((-1)**j*t**(2*j)/a**(2*j+2) for j in range(m))
        remainder = (-1)**m*t**(2*m)/(a**(2*m)*(a*a+t*t))
        assert sp.cancel(1/(a*a+t*t)-finite-remainder) == 0

    prime_checks = []
    for n in [3, 4, 7, 8, 12, 17, 24, 32, 48, 64, 80, 96]:
        s = prefix(n)
        choices = sorted(set([0, 1, 3, max(1, (11*n)//7)]))
        values = {m: approx(n, m, e) for m in choices}
        for p in primes_below(2*n):
            if p == 2 or p*p <= 2*n:
                continue
            j = (((2*n-1)//p)+1)//2
            target = (-1)**((p-1)//2)*residue(prefix(j), p) % p
            assert residue(p**4*s, p) == target
            if (6*n) % p == 0 or target == 0:
                continue
            for m, value in values.items():
                assert vp(value, p) == -4
                assert vp(Q(value.denominator), p) == 4
                assert residue(p**4*value, p) == target
                prime_checks.append(dict(N=n, m=m, p=p, J=j, valuation=-4,
                                         residue=target))

    cases = []
    for n in [1, 2, 3, 4, 8, 16, 32, 48, 64]:
        for m in sorted(set([0, 1, 2, (11*n)//7, 2*n])):
            value = approx(n, m, e)
            assert gcd(value.numerator, value.denominator) == 1
            # Rebuild the correction from formal series integration.
            laplace = sum((inverse[j]*factorial(2*j+3)/
                           (12*(2*n)**(2*j+4)) for j in range(m)), Q(0))
            assert value == prefix(n)+(-1)**n*laplace
            cases.append(dict(N=n, m=m, numerator=str(value.numerator),
                              denominator=str(value.denominator),
                              denominator_bits=value.denominator.bit_length()))

    # Demonstrate that primitive products stay primitive (Gauss's lemma).
    product_cases = []
    for parameters in [[(1,0),(2,0)], [(1,1)]*4,
                       [(3,2),(4,6),(8,12)], [(16,1),(16,2)]]:
        polynomial = sp.Poly(1, x)
        for n, m in parameters:
            z = approx(n, m, e)
            polynomial *= sp.Poly(z.denominator*x-z.numerator, x)
        assert polynomial.content() == 1
        product_cases.append(dict(parameters=parameters,
                                  degree=polynomial.degree(), content=1))

    data = dict(status='all checks passed',
                euler_coefficients_three_methods=160,
                small_symbolic_sech_coefficients=13,
                finite_geometric_identities=8,
                prime_valuation_checks=len(prime_checks),
                actual_primitive_cases=len(cases),
                prime_checks=prime_checks, cases=cases,
                product_checks=product_cases,
                limitations='Finite checks do not prove asymptotic statements; see findings.txt.')
    (OUT/'exact_checks.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in
                      ['prime_checks','cases','product_checks']}, indent=2))


if __name__ == '__main__':
    main()
