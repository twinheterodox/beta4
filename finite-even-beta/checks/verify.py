"""Exact certificates and consistency checks for the finite beta result.

Only Python's standard library is required.  The general arguments are
in the accompanying manuscript; finite computations do not replace them.
"""
from fractions import Fraction as F
from math import factorial
import json
from pathlib import Path
from seed import run, constants, evaluate, chi


def decay_certificate():
    u, S = F(3, 5000), 365
    f = ((u+18)/u)**16 * ((u+7)/(u+11))**S
    assert f < 1
    e_bound = F(68, 25)
    assert F(163, 60)+F(7, 4320) < e_bound
    upper4 = (e_bound**(4*S)*2**144*3**216*(18+u)**288
              *4**(4*S)*(F(7)**7/F(11)**11)**S)
    assert upper4 < F(3, 4)**4
    return {'saddle_upper_bound': '3/5000', 'f_at_upper_bound_lt_1': True,
            'f_at_upper_bound_illustrative_decimal': float(f),
            'e_upper_bound': '68/25', 'rate_fourth_power_upper_lt_81_over_256': True,
            'rate_fourth_power_upper_illustrative_decimal': float(upper4),
            'rate_upper_illustrative_decimal': float(upper4)**0.25,
            'limsup_integer_scaled_form_nth_root_lt': '3/4',
            'comparison_arithmetic': 'Exact rational; decimals are illustrative only.'}


def split_product(values, p):
    unit, valuation = 1, 0
    for value in values:
        assert value
        while value % p == 0:
            value //= p
            valuation += 1
        unit = unit*(value % p) % p
    return unit, valuation


def large_evaluation(n, p, t, s=364, r=7):
    """Leading p-adic digit of normalized R_n(t), directly from factors."""
    S, L = s+1, 2*r+4
    fu, fv = split_product(range(1, n+1), p)
    gu, gv = split_product(range(1, n+2), p)
    unit = (pow(4, S*n, p)*pow(2, 2*L*n, p)*pow(3, 3*L*n, p)
            *pow(gu, S, p)*pow(fu, -4*L, p)) % p
    valuation = S*gv-4*L*fv
    # The integer block includes the additional factor t-r*n.
    u, v = split_product(range(t-r*n, t+(r+4)*n+1), p)
    unit, valuation = unit*u % p, valuation+v
    for a, b in ((1, 2), (1, 3), (2, 3)):
        u, v = split_product((b*(t-r*n+j)+a for j in range(L*n)), p)
        unit = unit*u*pow(b, -L*n, p) % p
        valuation += v
    u, v = split_product((t+4*k for k in range(n+1)), p)
    return unit*pow(u, -S, p) % p, valuation-S*v


def small_evaluation(m, p, z, s=364, r=7):
    S, L = s+1, 2*r+4
    C = (pow(4, S*m, p)*pow(2, 2*L*m, p)*pow(3, 3*L*m, p)
         *pow(factorial(m+1), S, p)*pow(factorial(m), -4*L, p)) % p
    return C*evaluate(z, m, s, r, p) % p


def dilation_checks():
    rows, other_classes = [], []
    # All four possible residue classes modulo 12; in particular the
    # third-denominator numerator blocks interchange when p=2 (mod 3).
    for p in (1009, 1013, 1019, 1039):
        assert all(p % d for d in range(2, 33))
        for z in (29, 31, 37):
            unit, valuation = large_evaluation(2*p, p, p*z)
            target = pow(3, -365, p)*small_evaluation(2, p, z) % p
            assert valuation == -364 and unit == target and target
            rows.append({'p': p, 'p_mod_12': p % 12, 'z': z,
                         'valuation': valuation, 'leading_unit': unit,
                         'lambda_times_seed': target})
        for x in (1, p//2, p-1):
            _, valuation = large_evaluation(2*p, p, p*37-4*x)
            assert valuation == 0
            other_classes.append({'p': p, 'x': x, 'valuation': valuation})
    return {'zero_class': rows, 'other_classes': other_classes,
            'scope': 'Evaluations are consistency checks; the factor proof establishes the rational-function congruence.'}


def exact_arithmetic_check(n=2, s=364, r=7):
    """Check the denominator lemmas over Q at the actual finite seed."""
    S, L = s+1, 2*r+4
    C = F(4**(S*n)*2**(2*L*n)*3**(3*L*n)*factorial(n+1)**S,
          factorial(n)**(4*L))
    delta = 8  # 4*lcm(1,2), with n=2 here.
    assert n == 2 and S > 4*L
    coeff = {}
    for h in range(n+1):
        cs = constants(n, r, h)
        assert cs.count(0) == 1
        num = [c for c in cs if c]
        den = [F(4*(k-h)) for k in range(n+1) if k != h]
        initial = C
        for c in num:
            initial *= c
        for c in den:
            initial /= c**S
        powers = [F(0)]*s
        for factors, multiplicity in ((num, 1), (den, -S)):
            for c in factors:
                inv, v = 1/c, F(1)
                for j in range(1, s):
                    v *= inv
                    powers[j] += multiplicity*v
        series = [initial]
        for j in range(1, s):
            series.append(sum((1 if k % 2 else -1)*powers[k]*series[j-k]
                              for k in range(1, j+1))/j)
        for i in range(1, s+1):
            coeff[i, h] = series[s-i]
            assert (delta**(S-i)*coeff[i, h]).denominator == 1
    reductions = {}
    for b in (1, 2, 3):
        H, q, last = [F(0)]*(s+1), F(0), 0
        for h in range(n+1):
            for k in range(last+1, 4*b*h+1):
                if chi(k):
                    for i in range(1, s+1):
                        H[i] += F(chi(k), k**i)
            q -= sum(coeff[i, h]*b**i*H[i] for i in range(1, s+1))
            last = 4*b*h
        assert (delta**S*q).denominator == 1
        unnormalized = q/C
        reductions[b] = unnormalized.numerator*pow(unnormalized.denominator, -1, 1000003) % 1000003
    assert reductions == {1: 614660, 2: 346158, 3: 44336}
    return {'n': n, 's': s, 'r': r, 'coefficient_integrality_checks': len(coeff),
            'constant_integrality_checks': 3,
            'unnormalized_seed_reductions_from_exact_rationals': reductions,
            'arithmetic': 'Exact rational, normalized coefficients',
            'scope': 'Actual finite seed; not the general denominator proof.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: assertions are part of this certificate.')
    result = {'status': 'Exact finite certificates and consistency checks; not formal verification',
              'seed': run(crosscheck=True, all_pairs=True),
              'decay': decay_certificate(),
              'dilation': dilation_checks(),
              'arithmetic': exact_arithmetic_check()}
    assert result['seed']['q_unnormalized_mod_prime'] == {1: 614660, 2: 346158, 3: 44336}
    assert result['seed']['elimination_pairs_checked'] == 16471
    assert result['seed']['zero_determinant_pairs'] == []
    assert result['seed']['low_weight_elimination_residue'] == 187267
    Path(__file__).with_name('verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'seed': result['seed'], 'decay': result['decay'],
                      'dilation_evaluations': len(result['dilation']['zero_class']),
                      'other_class_evaluations': len(result['dilation']['other_classes']),
                      'arithmetic': result['arithmetic']}, indent=2))


if __name__ == '__main__':
    main()
