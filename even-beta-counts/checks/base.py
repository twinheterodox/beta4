"""Arithmetic and partial-fraction helpers for the uniform even-beta construction.

Python standard library only. Floating-point search proposes parameters;
decay_certificate checks the proposed parameters exactly.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, floor, gcd, isqrt, log, prod


def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


@dataclass(frozen=True)
class Parameters:
    B: int
    s: int
    r: int = 7

    def __post_init__(self):
        assert self.B >= 1 and self.r >= 3 and self.s % 2 == 0
        assert self.S > self.D*self.L and self.s > 2*self.d

    @property
    def S(self):
        return self.s+1

    @property
    def L(self):
        return 2*self.r+4

    @property
    def d(self):
        return self.B-1

    @property
    def theta(self):
        return tuple(F(a, b) for b in range(1, self.B+1)
                     for a in range(1, b+1) if gcd(a, b) == 1)

    @property
    def D(self):
        return len(self.theta)

    @property
    def K_factors(self):
        out = {}
        for theta in self.theta:
            for p, v in factor(theta.denominator).items():
                out[p] = out.get(p, 0)+self.L*v+(self.L+p-2)//(p-1)
        return out

    @property
    def K(self):
        return prod(p**v for p, v in self.K_factors.items())

    def constants(self, n, h=0):
        return [F(-4*h-self.r*n)]+[
            F(-4*h-self.r*n+j)+theta
            for theta in self.theta for j in range(self.L*n)]

    def normalization(self, n, p):
        # Reduction of C_n when p>n+1 and all block denominators.
        f = prod(range(1, n+1)) % p
        return (pow(4, self.S*n, p)*pow(self.K, n, p)
                *pow(f*(n+1), self.S, p)*pow(f, -self.D*self.L, p)) % p


def log_f(P, x):
    return 4*P.D*log((x+P.L)/x)+P.S*log((x+P.r)/(x+P.r+4))


def saddle(P):
    lo, hi = 1e-300, 1.0
    while log_f(P, hi) > 0:
        hi *= 2
    for _ in range(100):
        mid = (lo+hi)/2
        if log_f(P, mid) > 0:
            lo = mid
        else:
            hi = mid
    return hi


def rational_saddle_upper(P):
    x = saddle(P)
    denominator = 10**max(0, 1-floor(log(x)/log(10)))
    return F(int(x*denominator)+1, denominator)






def decay_certificate(P):
    u = rational_saddle_upper(P)
    f = ((u+P.L)/u)**(4*P.D)*((u+P.r)/(u+P.r+4))**P.S
    assert f < 1
    assert F(163, 60)+F(7, 4320) < F(68, 25)
    U = (F(68, 25)**(4*P.S)*P.K**4*(P.L+u)**(4*P.D*P.L)
         *4**(4*P.S)*(F(P.r)**P.r/F(P.r+4)**(P.r+4))**P.S)
    assert U < 1
    return {'saddle_upper': str(u), 'f_at_upper_lt_1': True,
            'fourth_power_rate_upper_lt_1': True,
            'fourth_power_rate_upper_illustrative': float(U),
            'rate_upper_illustrative': float(U)**0.25,
            'cleared_form_nth_root_limsup_lt': '1',
            'arithmetic': 'Exact rational; decimals are illustrative.'}


def modq(x, p):
    return x.numerator*pow(x.denominator, -1, p) % p


def chi(j):
    return 0 if j % 2 == 0 else (1 if j % 4 == 1 else -1)


def partials_newton(P, p, n=1):
    out = {}
    for h in range(n+1):
        cs = P.constants(n, h)
        assert cs.count(0) == 1
        num = [modq(c, p) for c in cs if c]
        den = [4*(k-h) % p for k in range(n+1) if k != h]
        assert all(num) and all(den)
        powers = [0]*P.s
        for factors, multiplicity in ((num, 1), (den, -P.S)):
            for c in factors:
                inv, v = pow(c, -1, p), 1
                for j in range(1, P.s):
                    v = v*inv % p
                    powers[j] = (powers[j]+multiplicity*v) % p
        series = [prod(num) % p*pow(prod(den) % p, -P.S, p) % p]
        for j in range(1, P.s):
            value = sum((1 if k % 2 else -1)*powers[k]*series[j-k]
                        for k in range(1, j+1))
            series.append(value*pow(j, -1, p) % p)
        for i in range(1, P.S):
            out[i, h] = series[P.s-i]
    return out


def mul(a, b, degree, p):
    out = [0]*min(len(a)+len(b)-1, degree+1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:len(out)-i]):
            out[i+j] = (out[i+j]+x*y) % p
    return out


def partials_product(P, p, n=1):
    """Independent direct product of local series (no Newton recurrence)."""
    out = {}
    for h in range(n+1):
        cs = P.constants(n, h)
        assert cs.count(0) == 1
        series = [1]
        for c in cs:
            if c:
                series = mul(series, [modq(c, p), 1], P.s-1, p)
        for k in range(n+1):
            if k != h:
                inv = pow(4*(k-h), -1, p)
                factor = [(-1)**j*comb(P.S+j-1, j)*pow(inv, P.S+j, p) % p
                          for j in range(P.s)]
                series = mul(series, factor, P.s-1, p)
        for i in range(1, P.S):
            out[i, h] = series[P.s-i]
    return out




def constant_term(P, a, b, p, n=1):
    harmonic, result, last = [0]*P.S, 0, 0
    bs = [pow(b, i, p) for i in range(P.S)]
    for h in range(n+1):
        for j in range(last+1, 4*b*h+1):
            if chi(j):
                inv, v = pow(j, -1, p), 1
                for i in range(1, P.S):
                    v = v*inv % p
                    harmonic[i] = (harmonic[i]+chi(j)*v) % p
        result = (result-sum(a[i, h]*bs[i]*harmonic[i]
                            for i in range(1, P.S))) % p
        last = 4*b*h
    return result


def constant_term_direct(P, a, b, p, n=1):
    """Sum each principal part by Horner evaluation at each grid point."""
    result = 0
    for h in range(1, n+1):
        for j in range(1, 4*b*h+1, 2):
            invx = b*pow(j, -1, p) % p
            value = 0
            for i in range(P.s, 0, -1):
                value = (value+a[i, h])*invx % p
            result = (result-chi(j)*value) % p
    return result


def determinant(matrix, p):
    a = [row[:] for row in matrix]
    det = 1
    for k in range(len(a)):
        pivot = next((j for j in range(k, len(a)) if a[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            det = -det
        v = a[k][k]
        det = det*v % p
        iv = pow(v, -1, p)
        for j in range(k+1, len(a)):
            w = a[j][k]*iv % p
            for l in range(k+1, len(a)):
                a[j][l] = (a[j][l]-w*a[k][l]) % p
    return det % p




def evaluate(P, t, p, n=1):
    t = modq(F(t), p)
    num = prod((t+modq(c, p)) % p for c in P.constants(n)) % p
    den = prod(pow(t+4*h, P.S, p) for h in range(n+1)) % p
    return num*pow(den, -1, p) % p


