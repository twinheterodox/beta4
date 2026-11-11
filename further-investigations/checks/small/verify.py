"""Same-session checks by product division and Euler sums, without probe imports."""
from pathlib import Path
from math import comb, factorial, gcd, lcm
from functools import reduce, lru_cache
from collections import Counter
import json
from flint import fmpq as Q, arb, ctx

HERE = Path(__file__).resolve().parent
COUNTS = Counter()


def mul(a, b, size):
    c = [0]*min(size, len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:len(c)-i]):
            c[i+j] += x*y
    return c


@lru_cache(None)
def product_partials(n, A, r, v):
    """Truncated integer products followed by exact power-series division."""
    rows = []
    for h in range(n+1):
        numerator = [(2**n*factorial(n))**(A-r*v)]
        for j in range(r*n):
            z = 2*j-(r-1)*n-2*h+1
            numerator = mul(numerator, [comb(v, k)*z**(v-k)
                                       for k in range(v+1)], A)
        denominator = [1]
        for k in range(n+1):
            if k == h:
                continue
            z = 2*(k-h)
            denominator = mul(denominator, [comb(A, j)*z**(A-j)
                                           for j in range(A)], A)
        numerator += [0]*(A-len(numerator))
        denominator += [0]*(A-len(denominator))
        quotient = []
        for j in range(A):
            quotient.append((Q(numerator[j])-sum((denominator[k]*quotient[j-k]
                for k in range(1, j+1)), Q(0)))/denominator[0])
        rows.append([Q(0)]+list(reversed(quotient)))
    COUNTS['rational_functions_reconstructed_by_product_division'] += 1
    return rows


def complete_form(rows, start, derivative):
    n, A = len(rows)-1, len(rows[0])-1
    sums = {s: [Q(0)] for s in range(derivative+1, A+derivative+1)}
    for s in sums:
        for j in range(start+n):
            sums[s].append(sums[s][-1]+Q((-1)**j, (2*j+1)**s))
    q, beta = Q(0), {}
    for i in range(1, A+1):
        s = i+derivative
        c = [(-1)**h*comb(s-1, derivative)*row[i] for h, row in enumerate(rows)]
        coefficient = sum(c, Q(0))
        q -= sum((c[h]*sums[s][start+h] for h in range(n+1)), Q(0))
        if s % 2:
            assert coefficient == 0
        else:
            beta[s] = coefficient
    return q, beta


def check_coefficients(recorded, q, beta):
    assert q == Q(recorded['constant'])
    assert {str(s): str(z) for s, z in beta.items()} == recorded['beta_coefficients']
    vector = [q]+list(beta.values())
    denominator = lcm(*(int(z.q) for z in vector))
    integer_vector = [int(denominator*z) for z in vector]
    content = abs(reduce(gcd, integer_vector, 0))
    primitive = [z//content for z in integer_vector]
    assert list(map(str, primitive)) == recorded['primitive_vector']
    assert Q(denominator, content) == Q(recorded['primitive_multiplier'])
    COUNTS['complete_forms_and_primitive_vectors'] += 1


def prime_list(limit):
    return [p for p in range(2, limit+1) if all(p % d for d in range(2, int(p**0.5)+1))]


def valuation(z, p):
    assert z
    a, b = int(z.p), int(z.q)
    value = 0
    while a % p == 0:
        value += 1
        a //= p
    while b % p == 0:
        value -= 1
        b //= p
    return value


def smooth_denominator(q, n):
    d = int(q.q)
    for p in prime_list(max(n, 2)):
        while d % p == 0:
            d //= p
    assert d == 1
    COUNTS['constant_denominators_supported_at_2_and_primes_at_most_n'] += 1


def euler_beta_enclosures(max_s, terms):
    """Delta^j (2k+1)^-s = b_j h_(s-1)(1,1/3,...,1/(2j+1))."""
    homogeneous = [arb(1)]+[arb(0)]*(max_s-1)
    b, divisor = arb(1), arb(2)
    values = {s: arb(0) for s in range(2, max_s+1, 2)}
    for j in range(terms):
        odd_reciprocal = arb(1)/(2*j+1)
        for k in range(1, max_s):
            homogeneous[k] += odd_reciprocal*homogeneous[k-1]
        for s in values:
            values[s] += b*homogeneous[s-1]/divisor
        b *= arb(2*(j+1))/(2*j+3)
        divisor *= 2
    remainder = arb(2)**(-terms)
    return {s: value+arb(remainder/2, remainder/2) for s, value in values.items()}


def main():
    cases = [json.loads((HERE/name).read_text())
             for name in json.loads((HERE/'cases.json').read_text())['files']]
    elimination = json.loads((HERE/'elimination.json').read_text())['records']
    all_forms, prime_records = [], []
    for d in cases:
        n, A, r, v = d['n'], d['pole_order'], d['grid_length_ratio'], d['zero_multiplicity']
        rows = product_partials(n, A, r, v)
        assert [[str(z) for z in row] for row in rows] == d['partial_fraction_coefficients']
        COUNTS['recorded_partial_fraction_coefficients'] += (n+1)*A
        for key, start in [('full_sum', 0), ('tail_sum', (r-1)*n//2)]:
            q, beta = complete_form(rows, start, 2)
            check_coefficients(d[key], q, beta)
            all_forms.append((n, key, d[key]))
            if v >= 3:
                smooth_denominator(q, n)
        if v != 1:
            continue
        q = Q(d['tail_sum']['constant'])
        multiplier = Q(d['tail_sum']['primitive_multiplier'])
        for p in prime_list((r+1)*n-1):
            if p <= max(r*n, A+2):
                continue
            t = ((r+1)*n-1-p)//2
            total = sum(factorial(j)*factorial(r*n-1-j)*pow(
                pow(factorial(n-t+j)*factorial(t-j), A, p), -1, p)
                for j in range(t+1)) % p
            predicted = (A*pow(2, -1, p)*(-1)**((p-1)//2+n-t)
                         *pow(factorial(n), A-r, p)*total) % p
            scaled = q*p**(A+1)
            assert scaled.q % p
            actual = int(scaled.p)*pow(int(scaled.q), -1, p) % p
            assert predicted == actual
            if predicted:
                assert valuation(q, p) == -A-1
                assert valuation(multiplier, p) == A+1
            prime_records.append(dict(n=n, A=A, r=r, prime=p, residue=predicted,
                constant_valuation=valuation(q, p), primitive_multiplier_valuation=valuation(multiplier, p)))
    for d in elimination:
        n, A, r = d['n'], d['pole_order'], d['grid_length_ratio']
        rows = product_partials(n, A, r, 1)
        weighted = []
        for h, row in enumerate(rows):
            padded = row+[Q(0), Q(0)]
            x = n-2*h
            weighted.append([Q(0)]+[x*x*padded[i]+2*x*padded[i+1]+padded[i+2]
                                    for i in range(1, A+1)])
        q0, b0 = complete_form(rows, 0, 0)
        q2, b2 = complete_form(weighted, 0, 0)
        assert (q0, b0) == complete_form(rows, (r-1)*n//2, 0)
        assert (q2, b2) == complete_form(weighted, (r-1)*n//2, 0)
        COUNTS['undifferentiated_full_tail_identities'] += 2
        for key, q, beta in [('base', q0, b0), ('weighted', q2, b2)]:
            check_coefficients(d[key], q, beta)
            smooth_denominator(q, n)
            all_forms.append((n, key, d[key]))
        w0, w2 = map(int, d['primitive_integer_weights'])
        assert gcd(w0, w2) == 1
        assert w0*b0[2]+w2*b2[2] == 0
        q = w0*q0+w2*q2
        beta = {s: w0*b0[s]+w2*b2[s] for s in b0 if s != 2}
        assert all(beta.values())
        check_coefficients(d['combined'], q, beta)
        all_forms.append((n, 'combined', d['combined']))
        COUNTS['beta2_eliminations'] += 1
    max_bits = max(d['max_primitive_coefficient_bits'] for _, _, d in all_forms)
    terms = max_bits+160
    ctx.prec = max_bits+512
    beta = euler_beta_enclosures(12, terms)
    for n, key, d in all_forms:
        vector = list(map(int, d['primitive_vector']))
        support = list(map(int, d['beta_coefficients']))
        actual = arb(vector[0])+sum((z*beta[s] for z, s in zip(vector[1:], support)), arb(0))
        assert not actual.contains(0)
        assert actual.overlaps(arb(d['primitive_value']))
        assert (abs(actual) < 1) == d['primitive_abs_less_than_one']
        assert abs(actual) > 1 or abs(actual) < 1
        if n >= 3 and key in ('full_sum', 'tail_sum', 'combined'):
            assert abs(actual) > 1
            COUNTS['final_forms_at_indices_at_least_three_exceed_one'] += 1
        COUNTS['primitive_values_checked_by_independent_Euler_formula'] += 1
    COUNTS['large_prime_residue_formula'] = len(prime_records)
    COUNTS['nonzero_large_prime_residues'] = sum(bool(r['residue']) for r in prime_records)
    result = dict(all_checks_passed=True, counts=dict(COUNTS),
        Euler_terms=terms, Euler_precision_bits=ctx.prec,
        large_prime_records=prime_records,
        zero_residue_cases=[r for r in prime_records if not r['residue']],
        scope='Same-session alternative computations, not an independent audit. No asymptotic or irrationality conclusion.')
    (HERE/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'large_prime_records'}, indent=2))


if __name__ == '__main__':
    main()
