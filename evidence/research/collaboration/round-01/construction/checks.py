"""Exact checks for the round-01 construction findings.

Run from beta4 with:
  .venv/bin/python collaboration/round-01/construction/checks.py

The kernel is rebuilt here from its defining finite sums, without importing
the prior experiment implementations. All output stays beside this script.
The local and determinant identities use exact rationals. Finite numerical
comparisons are Arb interval checks, not an asymptotic proof.
"""
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path
import json
import time

from flint import arb, ctx, fmpq, fmpq_mat


HERE = Path(__file__).resolve().parent


def primes_through(n):
    sieve = [True] * (n + 1)
    sieve[:2] = [False, False]
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            for j in range(p*p, n+1, p):
                sieve[j] = False
    return [j for j in range(2, n+1) if sieve[j]]


def val(q, p):
    q = Q(str(q))
    if not q:
        return float('inf')
    a, b = abs(q.numerator), q.denominator
    out = 0
    while a % p == 0:
        a //= p
        out += 1
    while b % p == 0:
        b //= p
        out -= 1
    return out


def mod(q, p):
    q = Q(str(q))
    return q.numerator * pow(q.denominator, -1, p) % p


def fq(q):
    q = Q(q)
    return fmpq(q.numerator, q.denominator)


def aq(q):
    q = Q(str(q))
    return arb(q.numerator) / q.denominator


def product(xs):
    out = Q(1)
    for x in xs:
        out *= x
    return out


def cauchy_det(rows, cols):
    if not rows:
        return Q(1)
    vr = product(b-a for i, a in enumerate(rows) for b in rows[i+1:])
    vc = product(b-a for i, a in enumerate(cols) for b in cols[i+1:])
    return vr * vc / product(a+b for a in rows for b in cols)


def rational_kernel(nodes):
    bound = max(nodes)
    s3, s4 = [Q(0)], [Q(0)]
    for j in range(bound):
        s3.append(s3[-1] + Q(1, (4*j+1)**3) - Q(1, (4*j+3)**3))
        s4.append(s4[-1] + Q(1, (4*j+1)**4) - Q(1, (4*j+3)**4))
    out = []
    for a in nodes:
        line = []
        for b in nodes:
            if a == b:
                line.append(-s4[a] + s3[a]/(12*a) if a else Q(0))
            else:
                line.append((b*s3[b]-a*s3[a]) / (6*(b*b-a*a)))
        out.append(line)
    return out


def polynomial(nodes, R):
    # det(X A+R) = det(A) charpoly(-A^-1 R)(X), including A_00=2.
    A = [2 if a == 0 else 1 for a in nodes]
    neg = fmpq_mat([[fq(-q/A[i]) for q in row] for i, row in enumerate(R)])
    return [Q(str(q))*product(A) for q in neg.charpoly()]


def primitive(coeffs):
    den = lcm(*(q.denominator for q in coeffs))
    ints = [int(q*den) for q in coeffs]
    content = 0
    for q in ints:
        content = gcd(content, abs(q))
    return [q//content for q in ints], Q(den, content)


def horner(coeffs, x):
    y = arb(0)
    for q in reversed(coeffs):
        y = y*x + aq(q)
    return y


def rank_one_checks():
    frames = [([1], [3]), ([1, 2], [2, 5]), ([1, 3, 7], [2, 4, 9]),
              ([2, 5, 8, 11], [1, 3, 7, 14]),
              ([0], [0]), ([0, 1, 4], [0, 2, 7]),
              ([0, 2, 6, 9], [0, 1, 4, 13])]
    checks = []
    for rows, cols in frames:
        for tau in (Q(1, 8), Q(1, 3), Q(2)):
            R = []
            for a in rows:
                line = []
                for b in cols:
                    q = 1/tau if a == b == 0 else (a*b+tau*(a+b))/((a+b)*(a+tau)*(b+tau))
                    line.append(fq(q))
                R.append(line)
            actual = Q(str(fmpq_mat(R).det()))
            if rows[0] == cols[0] == 0:
                expect = cauchy_det(rows[1:], cols[1:]) / tau
                expect *= product(Q(a)/(a+tau) for a in rows[1:])
                expect *= product(Q(b)/(b+tau) for b in cols[1:])
            else:
                expect = cauchy_det(rows, cols)
                expect *= product(Q(a)/(a+tau) for a in rows)
                expect *= product(Q(b)/(b+tau) for b in cols)
                expect *= 1+tau*(sum(Q(1, a) for a in rows)+sum(Q(1, b) for b in cols))
            assert actual == expect
            checks.append(dict(rows=rows, cols=cols, tau=str(tau)))
    return checks


def principal_checks():
    ctx.prec = 4096
    pi = arb.pi()
    B = (arb(4).zeta(arb(1)/4)-arb(4).zeta(arb(3)/4))/256
    tau = 1/(arb(2).sqrt()*pi)
    c0 = pi**3/192
    M = 1+2/(3*arb(6).sqrt())
    frames = []
    for n in (2, 4, 6, 8, 12):
        for offset in (0, n//2, n, 2*n):
            frames.append(list(range(offset+1, offset+n+1)))
        frames.append(list(range(n)))
    examples, pchecks, entrychecks = [], 0, 0
    for nodes in frames:
        n, largest = len(nodes), max(nodes)
        R = rational_kernel(nodes)
        coeffs = polynomial(nodes, R)
        primitive_coeffs, H = primitive(coeffs)
        forced = 1
        local = []
        for p in primes_through(4*largest):
            if p <= max(3, 2*largest):
                continue
            high = sum(4*a > p for a in nodes)
            chi = 1 if p % 4 == 1 else -1
            for i, a in enumerate(nodes):
                for j, b in enumerate(nodes):
                    scaled = p**4 * R[i][j]
                    assert val(scaled, p) >= 0
                    expect = -chi if i == j and 4*a > p else 0
                    assert mod(scaled, p) == expect % p
                    entrychecks += 1
            actual = val(coeffs[n-high], p)
            assert actual == -4*high, (nodes, p, high, actual)
            assert val(H, p) >= 4*high
            forced *= p**(4*high)
            local.append(dict(p=p, rank=high, coefficient_valuation=actual))
            pchecks += 1
        value = horner(primitive_coeffs, B)
        det = horner(coeffs, B)
        assert value > 1 and det > 0
        if nodes[0] != 0:
            C = aq(cauchy_det(nodes, nodes))
            lower = c0**n * C
            for a in nodes:
                lower *= (a/(a+tau))**2
            lower *= 1+2*tau*aq(sum(Q(1, a) for a in nodes))
            upper = (c0*M)**n*C
        else:
            C = aq(cauchy_det(nodes[1:], nodes[1:]))
            lower = c0**n*C/tau
            for a in nodes[1:]:
                lower *= (a/(a+tau))**2
            tau_upper = arb(1)/8
            c_upper = c0*(M+10*pi*pi*tau_upper*tau_upper)
            upper = c_upper**n*C/tau_upper
            for a in nodes[1:]:
                upper *= (a/(a+tau_upper))**2
        assert det > lower and det < upper
        examples.append(dict(nodes=nodes, prime_layers=local,
                             primitive_log_value=value.log().str(40),
                             primitive_log_multiplier=aq(H).log().str(40),
                             forced_log_multiplier=arb(forced).log().str(40),
                             log_real_determinant=det.log().str(40),
                             lower_log_real_determinant=lower.log().str(40),
                             upper_log_real_determinant=upper.log().str(40)))
    return dict(cases=examples, coefficient_prime_checks=pchecks,
                entry_prime_checks=entrychecks,
                arb_precision_bits=ctx.prec)


def branch_location_checks():
    # At t^k= i*(sqrt(2)-1), v=i*pi/4. This diagnostic only checks
    # the listed arguments numerically; logarithmic nonremovability
    # and the exclusions for every odd k are proved in findings.txt.
    import mpmath as mp
    mp.mp.dps = 70
    q = mp.sqrt(2)-1
    rows = []
    for k in (3, 5, 7, 11):
        for j in range(k):
            t = q**(mp.mpf(1)/k)*mp.exp(1j*(mp.pi/2+2*mp.pi*j)/k)
            u, v = 2*mp.atanh(t), 2*mp.atanh(t**k)
            assert abs(v-1j*mp.pi/4) < mp.mpf('1e-65')
            other = [abs(mp.cosh(x)) for x in (u, v, u+v, u+2*v)]
            assert min(other) > mp.mpf('1e-5')
            assert abs(mp.cosh(2*v)) < mp.mpf('1e-64')
            rows.append(dict(k=k, root=j, z=mp.nstr(t*t, 30),
                             minimum_other_cosh=mp.nstr(min(other), 20)))
    return rows


def main():
    started = time.monotonic()
    output = dict(
        status='Local unpublished same-agent checks; proofs and scope are in findings.txt.',
        rank_one_identities=rank_one_checks(),
        principal=principal_checks(),
        branch_diagnostics=branch_location_checks(),
        elapsed_seconds=time.monotonic()-started,
    )
    (HERE/'checks.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(dict(rank_one_identities=len(output['rank_one_identities']),
                          principal_cases=len(output['principal']['cases']),
                          coefficient_prime_checks=output['principal']['coefficient_prime_checks'],
                          entry_prime_checks=output['principal']['entry_prime_checks'],
                          branch_diagnostics=len(output['branch_diagnostics']),
                          elapsed_seconds=output['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
