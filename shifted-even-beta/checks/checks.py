"""Exact checks for differentiated shifted-interval constructions.

Research checks, not an independent audit. Standard library only. Polynomial
and scaled-series routines adapt even-beta-counts/checks at commit 06bd732c.
All new derivative, multiplicity and Gram identities are checked here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, factorial, gcd, log, prod, ceil, isqrt
from pathlib import Path
import json


def factors(n):
    out, p = {}, 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def chi(j):
    return 0 if j % 2 == 0 else (-1)**((j-1)//2)


def mul(a, b, limit=None):
    length = len(a)+len(b)-1
    if limit is not None:
        length = min(length, limit+1)
    out = [0]*length
    for i, x in enumerate(a[:length]):
        for j, y in enumerate(b[:length-i]):
            out[i+j] += x*y
    return out


def power(a, m):
    out = [1]
    for _ in range(m):
        out = mul(out, a)
    return out


def evalp(a, x):
    out = 0
    for c in reversed(a):
        out = out*x+c
    return out


def det(a):
    a = [list(map(F, row)) for row in a]
    result = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        pivot = a[i][i]
        result *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for h in range(i+1, len(a)):
                a[j][h] -= ratio*a[i][h]
    return result


@dataclass(frozen=True)
class Params:
    B: int
    a: int
    S: int
    r: int = 7

    def __post_init__(self):
        assert self.a >= 1 and self.B >= 1 and self.r >= 3
        assert self.S % 2 and self.S > self.m*self.D*self.L

    @property
    def k(self):
        return 2*self.a-2

    @property
    def m(self):
        return self.k+1

    @property
    def L(self):
        return 2*self.r+4

    @property
    def theta(self):
        return [F(v, d) for d in range(1, self.B+1)
                for v in range(1, d+1) if gcd(v, d) == 1]

    @property
    def D(self):
        return len(self.theta)

    @property
    def K(self):
        return prod(t.denominator**self.L * prod(
            p**((self.L+p-2)//(p-1)) for p in factors(t.denominator))
            for t in self.theta)

    def constants(self, n, h=0):
        return [F(-4*h-self.r*n)]+[
            -4*h-self.r*n+j+t for t in self.theta for j in range(self.L*n)]

    def C(self, n):
        return F(4**(self.S*n)*self.K**(self.m*n)*factorial(n+1)**self.S,
                 factorial(n)**(self.m*self.D*self.L))


def basis(B, a):
    m = 2*a-1
    nodes = sorted({F(2*j-1, b) for b in range(1, B+1)
                    for j in range(1, b+1)})
    Q = [1]
    for v in nodes:
        Q = mul(Q, [v.numerator**2, v.denominator**2])
    hs = []
    J = m*len(nodes)-a
    for b in range(1, B+1):
        h = [0]*(J+1)
        for j in range(1, b+1):
            v = F(2*j-1, b)
            num, den = v.numerator, v.denominator
            term = [-2*(-1)**(b-j+l)*comb(m, 2*l)
                    *num**(m-2*l)*den**(m+2*l) for l in range(a)]
            for other in nodes:
                if other != v:
                    term = mul(term, power([other.numerator**2,
                                            other.denominator**2], m))
            assert len(term) <= len(h)
            for i, c in enumerate(term):
                h[i] += c
        hs.append(h)
    assert all(isinstance(c, int) for h in hs for c in h)
    # Directly evaluate H_b^(k)/k! at z=iy in Q[i] with y rational.
    # The real part of (c+iy)^(-m) is computed independently by binomial expansion.
    for b, h in enumerate(hs, 1):
        for y in (F(0), F(1, 2), F(2), F(7, 3)):
            direct = F(0)
            for j in range(1, 4*b, 2):
                c = F(j, b)-2
                real = sum(comb(m, 2*l)*c**(m-2*l)*(-y*y)**l
                           for l in range(a))
                direct += chi(j)*real/(c*c+y*y)**m
            assert evalp(h, y*y)/F(evalp(Q, y*y))**m == direct
    return nodes, Q, hs, J


def local_W(poly, J, n, h):
    c = 2*n-4*h
    out = [0]*(2*J+2)
    for j, v in enumerate(poly):
        for i in range(2*j+2):
            out[i] += v*(-1)**j*n**(2*J-2*j)*comb(2*j+1, i)*c**(2*j+1-i)
    return out


def base_coefficients(P, n):
    """Return delta^(S-i) a_i,h for the undifferentiated base R_n."""
    d = 1
    for j in range(1, n+1):
        d = d*j//gcd(d, j)
    delta = 4*d
    out = {}
    for h in range(n+1):
        poly, scale = [1], 1
        cs = P.constants(n, h)
        assert cs.count(0) == 1
        for c in cs:
            num, den = c.numerator, c.denominator
            for _ in range(P.m):
                poly = mul(poly, [num, delta*den], P.S-1)
                scale *= den
        other = [4*(j-h) for j in range(n+1) if j != h]
        initial = P.C(n)/(scale*prod(other)**P.S)
        for c in other:
            assert delta % c == 0
            poly = mul(poly, [(-delta//c)**j*comb(P.S+j-1, j)
                             for j in range(P.S)], P.S-1)
        for i in range(1, P.S-P.m+1):
            cleared = initial*poly[P.S-i]
            assert cleared.denominator == 1
            out[i, h] = cleared.numerator
        assert all(poly[i] == 0 for i in range(P.m))
    return out, delta


def differentiated(P, n, poly, J, base, delta):
    """Return delta^(S+k-nu) A_nu,h for V_n=(W_n R_n)^(k)/k!."""
    out = {(nu, h): 0 for h in range(n+1) for nu in range(1, P.S)}
    for h in range(n+1):
        w = local_W(poly, J, n, h)
        for i in range(1, P.S-P.m+1):
            cleared = sum(c*delta**j*base.get((i+j, h), 0)
                          for j, c in enumerate(w))
            out[i+P.k, h] = comb(i+P.k-1, P.k)*cleared
    assert all(out[nu, n-h] == (-1)**nu*out[nu, h]
               for nu in range(1, P.S) for h in range(n+1))
    assert all(sum(out[nu, h] for h in range(n+1)) == 0
               for nu in range(1, P.S, 2))
    return out


def constants(P, n, coeff, delta):
    qs, groups = [], 0
    for b in range(1, P.B+1):
        total = 0
        for h in range(1, n+1):
            for j in range(1, 4*b*h, 2):
                # delta^(S+k) sum_nu A_nu,h (b/j)^nu, via integer Horner.
                value, jpower = coeff[P.S-1, h], j
                for nu in range(P.S-2, -1, -1):
                    value = value*(delta*b)+(coeff[nu, h] if nu else 0)*jpower
                    jpower *= j
                assert value % j**(P.S-1) == 0
                total -= chi(j)*(value//j**(P.S-1))
                groups += 1
        qs.append(F(total, delta**(P.S+P.k)))
    return qs, groups


def gram_moments(P, nodes, hs, J):
    alpha = sorted([2-c for c in P.constants(1) if c < 2])
    assert len(alpha) == P.D*P.L//2 and len(set(alpha)) == len(alpha)
    assert set(nodes) <= set(alpha)
    scale, T = 1, [1]
    for v in alpha:
        scale *= v.denominator**2
        if v not in nodes:
            T = mul(T, power([v.numerator**2, v.denominator**2], P.m))
    degree = P.a+len(T)-1+2*J
    assert degree < P.S
    moments = [comb(2*P.S-2, P.S-1)]
    for j in range(degree):
        num, den = 4*(2*j+1)*moments[-1], 2*P.S-2*j-3
        assert num % den == 0
        moments.append(num//den)
    matrix = []
    for h in hs:
        row = []
        for g in hs:
            poly = [0]*P.a+mul(T, mul(h, g))
            row.append(sum(c*d for c, d in zip(poly, moments)))
        matrix.append(row)
    scalar = F((-1)**(P.D*P.L//2+P.a), scale**P.m*2**(4*P.S-2))
    return matrix, scalar


def valuation(x, p):
    if not x:
        return float('inf')
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v


def modq(x, p):
    return x.numerator*pow(x.denominator, -1, p) % p


def kernel_and_gram_check(P):
    nodes, Q, hs, J = basis(P.B, P.a)
    E = 2*J+1
    assert E == P.m*(2*len(nodes)-1) and P.S-P.m-E > 0
    G, scalar = gram_moments(P, nodes, hs, J)
    minors = [det([r[:j] for r in G[:j]]) for j in range(1, P.B+1)]
    assert all(v > 0 for v in minors)
    weights = [(-1)**(P.B-1+j)*det([[row[l] for l in range(P.B) if l != j]
                                  for row in G[:-1]]) for j in range(P.B)]
    assert all(w.denominator == 1 for w in weights)
    assert [sum(x*y for x, y in zip(row, weights)) for row in G] == [0]*(P.B-1)+[minors[-1]]
    base, delta = base_coefficients(P, 1)
    qs, grouped = [], 0
    for col, poly in enumerate(hs):
        coeff = differentiated(P, 1, poly, J, base, delta)
        q, groups = constants(P, 1, coeff, delta)
        grouped += groups
        assert q == [scalar*G[b][col]*P.C(1) for b in range(P.B)]
        # Independent residue integration by parts: Res T H^(k)/k!.
        for b in range(1, P.B+1):
            w = local_W(poly, J, 1, 0)
            aq = {i: sum(F(c*base.get((i+j, 0), 0), delta**(P.S-i-j))
                         for j, c in enumerate(w) if i+j < P.S)
                  for i in range(1, P.S-P.m+1)}
            residue = sum(aq[i]*(-1)**(i-1)*comb(P.k+i-1, P.k)
                          *sum(chi(j)*F(b, j)**(P.k+i) for j in range(1, 4*b, 2))
                          for i in aq)
            assert residue == q[b-1]
        qs.append(q)
    return {'a': P.a, 'B': P.B, 'S': P.S, 'r': P.r, 'E': E,
            'rational_kernel_checks': 4*P.B, 'exact_Gram_entries': P.B**2,
            'positive_principal_minors': P.B, 'integer_seed_vector': True,
            'grouped_integrality_checks': grouped}


def dilation_check(P, p):
    assert all(p % d for d in range(2, isqrt(p)+1)) and p > max(4*P.B, P.L, P.k)
    nodes, _, hs, J = basis(P.B, P.a)
    E = 2*J+1
    tau = P.S-1-E
    seeds, dilated = {}, {}
    for n, target in ((1, seeds), (p, dilated)):
        base, delta = base_coefficients(P, n)
        for col, poly in enumerate(hs):
            coeff = differentiated(P, n, poly, J, base, delta)
            q, groups = constants(P, n, coeff, delta)
            target[col] = coeff, q, delta, groups
    selected, other, constants_count, groups_count = 0, 0, 0, 0
    for col in range(P.B):
        ac, qc, dc, _ = seeds[col]
        ap, qp, dp, groups = dilated[col]
        groups_count += groups
        for h in range(p+1):
            for nu in range(P.k+1, P.S):
                value = F(ap[nu, h], dp**(P.S+P.k-nu))
                if h % p == 0:
                    scaled = value*F(p)**(tau-nu)
                    target = F(ac[nu, h//p], dc**(P.S+P.k-nu))
                    assert valuation(scaled, p) >= 0
                    assert modq(scaled, p) == pow(2, -P.S, p)*modq(target, p) % p
                    selected += 1
                else:
                    assert valuation(value, p) >= nu-P.k
                    other += 1
        for nu in range(P.k+1, P.S):
            rho = F(sum(ap[nu, h] for h in range(p+1)), dp**(P.S+P.k-nu))
            assert valuation(p**tau*rho, p) >= nu
        for b in range(P.B):
            assert modq(p**tau*qp[b], p) == chi(p)*pow(2, -P.S, p)*modq(qc[b], p) % p
            constants_count += 1
    return {'a': P.a, 'B': P.B, 'S': P.S, 'r': P.r, 'prime': p,
            'selected_coefficients': selected, 'other_coefficients': other,
            'constant_congruences': constants_count, 'grouped_integrality': groups_count}


def decay(P, u):
    f = ((u+P.L)/u)**(4*P.m*P.D)*((u+P.r)/(u+P.r+4))**P.S
    U = (F(68, 25)**(4*(P.S+P.k))*P.K**(4*P.m)
         *(P.L+u)**(4*P.m*P.D*P.L)*4**(4*P.S)
         *(F(P.r)**P.r/F(P.r+4)**(P.r+4))**P.S)
    assert f < 1 and U < 1
    return {'a': P.a, 'B': P.B, 's': P.S-1, 'r': P.r, 'u': str(u),
            'f_lt_1': True, 'U_lt_1': True,
            'f_illustrative': float(f), 'U_illustrative': float(U)}


def main():
    assert __debug__, 'Run without -O'
    gram = []
    for a in range(1, 4):
        for B in range(1, 5):
            D = sum(gcd(v, d) == 1 for d in range(1, B+1) for v in range(1, d+1))
            P = Params(B, a, (2*a-1)*10*D+1, 3)
            print(f'Gram a={a}, B={B}, S={P.S}', flush=True)
            gram.append(kernel_and_gram_check(P))
    dilation = []
    for a, B in ((1, 1), (2, 1), (3, 1), (2, 2)):
        D = sum(gcd(v, d) == 1 for d in range(1, B+1) for v in range(1, d+1))
        P = Params(B, a, (2*a-1)*10*D+1, 3)
        for p in (17, 19):
            print(f'Dilation a={a}, B={B}, p={p}', flush=True)
            dilation.append(dilation_check(P, p))
    # These three exact comparisons underlie the uniform endpoint proof.
    e0 = F(68, 25)
    assert F(901)**4*F(351, 551)**66 < 1
    assert F(901, 50)**10 < F(163, 60)**29
    assert ((4*e0)**4*F(7)**7/F(11)**11)**5*e0**16 < 1
    assert F(163, 60)+F(7, 4320) < e0
    # e0 < e^(11/10), using the first six positive terms of that series.
    assert e0 < sum(F(11, 10)**j/factorial(j) for j in range(6))
    record = {'status': 'Exact supporting checks; self-review, not an independent audit.',
              'gram': gram, 'dilation': dilation, 'uniform_rational_inequalities': 5}
    Path(__file__).with_name('checks.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'Gram_entries': sum(v['exact_Gram_entries'] for v in gram),
                      'positive_minors': sum(v['positive_principal_minors'] for v in gram),
                      'dilation_selected': sum(v['selected_coefficients'] for v in dilation),
                      'dilation_other': sum(v['other_coefficients'] for v in dilation),
                      'constant_congruences': sum(v['constant_congruences'] for v in dilation),
                      'grouped_integrality': sum(v['grouped_integrality'] for v in dilation)}, indent=2))


if __name__ == '__main__':
    main()
