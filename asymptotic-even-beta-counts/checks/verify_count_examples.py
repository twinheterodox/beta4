#!/usr/bin/env python3
"""Exact scalar, grid and finite-kernel checks for the restricted-grid manuscript."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

import sympy as sp


HERE = Path(__file__).resolve().parent
t, z = sp.symbols("t z")


def arctan_interval(q, terms):
    value = sum((F((-1)**k, (2*k+1)*q**(2*k+1)) for k in range(terms)), F(0))
    next_term = F((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(value, value+next_term), max(value, value+next_term)


def constant_enclosure():
    a0, a1 = arctan_interval(5, 30)
    b0, b1 = arctan_interval(239, 10)
    pi0, pi1 = 16*a0-4*b1, 16*a1-4*b0
    n = 2000
    partial = sum((F(1, k**3) for k in range(1, n+1)), F(0))
    # Half the difference of the two squared reciprocals at k-1/2 and
    # k+1/2 equals k^-3 / (1-1/(4k^2))^2. These differences telescope.
    upper_tail = F(2, (2*n+1)**2)
    lower_tail = (1-F(1, 4*(n+1)**2))**2 * upper_tail
    # (2 sqrt(kappa/3))^2 = 210 zeta(3)/pi^4.
    lower_square = 210*(partial+lower_tail)/pi1**4
    upper_square = 210*(partial+upper_tail)/pi0**4
    lower, upper = F("1.6098018249961"), F("1.6098018249962")
    assert lower**2 < lower_square < upper_square < upper**2
    return [str(lower), str(upper)]


def grid_check(cutoff):
    # phi(q)^2/q is a product of factors p^(e-2)(p-1)^2. Each is
    # at least one, except for the possible factor 1/2 at p=2, e=1.
    denominators = [q for q in range(4, 2*cutoff**2+1, 4)
                    if sp.totient(q) <= cutoff]
    reduced = {sp.Rational(j, q) for q in denominators for j in range(1, q)
               if math.gcd(j, q) == 1}
    active = {sp.Rational(j, q) for q in denominators for j in range(1, q, 2)}
    assert reduced == active == {1-th for th in active}
    assert len(reduced) == sum(sp.totient(q) for q in denominators)
    return denominators, sorted(reduced)


def poly(expression):
    return sp.Poly(expression, t, domain=sp.QQ)


def taylor_quotient(numerator, denominator, terms):
    """Coefficients at zero of a polynomial quotient with unit constant over Q."""
    leading = denominator.nth(0)
    assert leading != 0
    coefficients = []
    for i in range(terms):
        product = sum(denominator.nth(j)*coefficients[i-j] for j in range(1, i+1))
        coefficients.append((numerator.nth(i)-product)/leading)
    return coefficients


def kernel_check(n, denominators, theta):
    s = 8
    dn = math.lcm(*range(1, n+1))
    A1 = math.prod(q**int(sp.totient(q)) for q in denominators)
    A2 = math.prod(p**(int(sp.totient(q))//(p-1))
                   for q in denominators for p in sp.factorint(q))
    base = poly((A1*A2)**n * math.factorial(n)**(s-len(theta))
                * math.prod(sp.rf(t+th, n) for th in theta))
    denominator = poly(sp.rf(t, n+1))**s
    R = base.as_expr()/denominator.as_expr()
    assert sp.cancel(R.subs(t, -n-t)-R) == 0
    pieces = {(k, i): denominator.exquo(poly(t+k)**i)
              for k in range(n+1) for i in range(1, s+1)}
    local = [pieces[k, s].shift(-k) for k in range(n+1)]
    basis = [sp.rf(t, i)/(math.factorial(n)**(i//n)*math.factorial(i%n))
             for i in range(3)]
    rows = []
    cancellation = None
    for i in range(3):
        for j in range(i, 3):
            P = basis[i]*basis[j]
            numerator = base*poly(P.subs(t, t+n//2)+P.subs(t, -t-n//2))
            coeff = {}
            reconstructed = poly(0)
            for k in range(n+1):
                series = taylor_quotient(numerator.shift(-k), local[k], s)
                for order in range(1, s+1):
                    value = series[s-order]
                    coeff[order, k] = value
                    assert (dn**(s-order)*value).q == 1
                    reconstructed += pieces[k, order]*value
            assert reconstructed == numerator
            assert all(coeff[order, n-k] == (-1)**order*coeff[order, k]
                       for order in range(1, s+1) for k in range(n+1))
            period = {}
            for order in range(1, s+1):
                rho = sum(coeff[order, k] for k in range(n+1))
                if order % 2:
                    assert rho == 0
                else:
                    scaled = dn**s * 4**order * rho
                    assert scaled.q == 1
                    period[str(order)] = str(scaled)
            constant = sp.Integer(0)
            for th, sign in ((sp.Rational(1, 4), 1), (sp.Rational(3, 4), -1)):
                for k in range(1, n+1):
                    for ell in range(k):
                        terms = [dn**s*coeff[order, k]/(ell+th)**order
                                 for order in range(1, s+1)]
                        assert sum(terms).q == 1
                        if cancellation is None and any(v.q != 1 for v in terms):
                            cancellation = dict(pair=[i, j], theta=str(th), pole=k, ell=ell,
                                                term_denominators=[int(v.q) for v in terms],
                                                grouped_sum=str(sum(terms)))
                        constant -= sign*sum(terms)
            assert constant.q == 1
            rows.append(dict(basis_product=[i, j], scaled_constant=str(constant),
                             scaled_period_coefficients=period))
    assert cancellation is not None
    Rh = R.subs(t, z-n//2)
    A = (z+1+n//2)**(s+1)*(z-n//2)*math.prod(z+th-n//2 for th in theta)
    B = (z-n//2)**(s+1)*(z+1+n//2)*math.prod(z+th+n//2 for th in theta)
    assert sp.cancel(A*Rh.subs(z, z+1)-B*Rh) == 0
    for th in theta:
        assert A.subs(z, n//2-1+th) == 0
        assert B.subs(z, -n//2-th) == 0
    return dict(B=2, n=n, s=s, A1=A1, A2=A2, basis_products=rows,
                grouped_tail_witness=cancellation, recurrence_and_transitions="PASS")


def main():
    enclosure = constant_enclosure()
    q2, t2 = grid_check(2)
    q4, t4 = grid_check(4)
    assert q2 == [4] and t2 == [sp.Rational(1, 4), sp.Rational(3, 4)]
    assert q4 == [4, 8, 12] and len(t4) == 10
    assert sp.Rational(3, 12) in t4 and sp.Rational(2, 12) not in t4
    kernels = [kernel_check(n, q2, t2) for n in (2, 4)]
    result = dict(status="PASS", constant_enclosure=enclosure,
                  scalar_method="Exact Machin and telescoping zeta(3) tail bounds.",
                  grid_checks={str(c): dict(denominators=q, fractions=[str(v) for v in th])
                               for c, q, th in ((2, q2, t2), (4, q4, t4))},
                  kernels=kernels, sympy_version=sp.__version__,
                  scope="Finite algebra and decimal enclosure. The large-parameter result uses the proof.",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/"count-examples.json").write_text(json.dumps(result, indent=2)+"\n")
    print("PASS. Constant lies in (1.6098018249961, 1.6098018249962). "
          "Restricted grids at B=2,4 and six basis products at each n=2,4 verified.")


if __name__ == "__main__":
    main()
