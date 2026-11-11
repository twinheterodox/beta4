from fractions import Fraction as Rational
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
G = Rational(221,16)

def primes_up_to(bound):
    sieve = bytearray(b'\x01')*(bound+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2,isqrt(bound)+1):
        if sieve[p]:
            sieve[p*p:bound+1:p] = b'\x00'*((bound-p*p)//p+1)
    return [p for p in range(2,bound+1) if sieve[p]]

def factorial_valuations(bound,p):
    result = [0]*(bound+1)
    for a in range(1,bound+1):
        value, exponent = a,0
        while value%p==0:
            value //= p
            exponent += 1
        result[a] = result[a-1]+exponent
    return result

def floor(q):
    return q.numerator//q.denominator

def ceil(q):
    return -floor(-q)

report = dict(small_prime_cases=0, large_prime_cases=0, residue_classes=0,
              sample_sizes=[96,192,384,768])
for n in report['sample_sizes']:
    scale = n//96
    blocks = [(j*scale,n-j*scale) for j in range(31)]
    orders = [1+min(30,a//scale,(n-a)//scale) for a in range(n+1)]
    for p in primes_up_to(n):
        if p==2 or p*p<=7*n:
            vp = factorial_valuations(7*n,p)
            bound = n if p==2 else 7*n
            lam, power = 0,p
            while power<=bound:
                lam += 1
                power *= p
            candidates = []
            for a in range(n+1):
                poles = 0
                for left,right in blocks:
                    if a<left:
                        poles += vp[right-a]-vp[left-a-1]
                    elif a>right:
                        poles += vp[a-left]-vp[a-right-1]
                    else:
                        poles += vp[a-left]+vp[right-a]
                if p==2:
                    nu = 1-(15*n//2+6)-poles
                else:
                    negative,positive = 7*n//4-a+1,a+3*n//4+1
                    nu = 3*(vp[2*negative]-vp[negative]
                            +vp[2*positive]-vp[positive])-poles
                candidates.append(-nu+(31 if p==2 else 33)*lam)
            exponent = max(0,max(candidates))
            remainder = exponent-G*vp[n]
            upper = 15*n+36+107*(lam+1) if p==2 else 128*(lam+1)
            assert remainder<=upper
            report['small_prime_cases'] += 1
            continue
        address_counts = [0]*p
        pole_counts = [0]*p
        zero_counts = [0]*p
        omitted_counts = [0]*p
        for a,order in enumerate(orders):
            address_counts[a%p] += 1
            pole_counts[a%p] += order
        for ell in range(-7*n//4-1,3*n//4+1):
            residue = (-(2*ell+1)*pow(2,-1,p))%p
            zero_counts[residue] += 1
            if -7*n//4<=ell<3*n//4:
                omitted_counts[residue] += 1
        candidates = []
        envelope_candidates = []
        x = Rational(n,p)
        for a in range(p):
            if address_counts[a]<2:
                continue
            y = Rational(a,p)
            assert x>=2 or y<=x-1
            predicted_poles = sum(floor(y-Rational(j,96)*x)
                                  -ceil(y-Rational(96-j,96)*x)+1
                                  for j in range(31))
            predicted_zeros = floor(y+Rational(3,4)*x+Rational(1,2)) \
                              -floor(y-Rational(7,4)*x+Rational(1,2))
            assert predicted_poles==pole_counts[a]
            assert predicted_zeros==omitted_counts[a]
            prefactor = int((n-2*a)%p==0)
            F = pole_counts[a]-3*zero_counts[a]-prefactor
            candidates.append(F+2)
            envelope_candidates.append(2+predicted_poles-3*predicted_zeros)
            report['residue_classes'] += 1
        assert max(0,max(candidates))<=max(0,max(envelope_candidates))
        report['large_prime_cases'] += 1
(ROOT/'independent-local-bounds.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
