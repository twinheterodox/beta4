"""Independent global-polynomial reconstruction of the powered seed.

No manuscript program is imported. Construct full polynomials first, then
Taylor-shift them globally at every pole, divide out the exact pole power,
and invert by formal power-series division. No reflection shortcut,
logarithmic derivative recurrence, or factorial/prefix local product is used.
"""
from pathlib import Path
import hashlib
import json
import math
import time
from flint import nmod_poly

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'results'

def product(polys, prime):
    while len(polys) > 1:
        polys = [polys[i] * polys[i+1] if i+1 < len(polys) else polys[i]
                 for i in range(0, len(polys), 2)]
    return polys[0] if polys else nmod_poly([1], prime)

def audit_seed(prime):
    start = time.monotonic()
    assert all(prime % d for d in range(2, math.isqrt(prime)+1))
    n = 480
    blocks = [(0, n, 8)] + [(3*k, n-3*k, 1)
        for k in [*range(1,9), *range(10,47,2)]]
    orders = [sum(mult for lo,hi,mult in blocks if lo <= j <= hi)
              for j in range(n+1)]
    poles = [n-2*j for j in range(n+1)]
    zeros = list(range(-1599,1600,2))
    numerator = nmod_poly([0,1],prime) * product(
        [nmod_poly([-z,1],prime) for z in zeros], prime)**3
    denominator = product([nmod_poly([-a,1],prime)**r
                           for a,r in zip(poles,orders)], prime)
    assert numerator.degree() == 4801
    assert denominator.degree() == 13427
    assert sum(orders) == 13392+35
    probes = [5000,5001,7777]
    assert all(int(denominator(x)) for x in probes)
    q_values = [1,2,3]
    states = {}
    for q in q_values:
        cap = 37*q
        states[q] = dict(beta=[0]*(cap+1),tail=[0]*(cap+1),constant=0,
                         probe_sums=[0]*len(probes),partials=[])
    for j,(pole,r) in enumerate(zip(poles,orders)):
        cap = 3*r
        linear = nmod_poly([pole,1],prime)
        full_den = denominator.compose(linear)
        assert all(int(full_den[k]) == 0 for k in range(r))
        assert int(full_den[r]) != 0
        den = full_den.right_shift(r).truncate(cap)
        num = numerator.compose(linear).truncate(cap)
        base = (num * den.inverse_series_trunc(cap)).truncate(cap)
        for q,state in states.items():
            size = q*r
            jet = base.pow_trunc(q,size)
            coeffs = [0]+[int(jet[size-i]) for i in range(1,size+1)]
            state['partials'].append(coeffs)
            for i in range(1,size+1):
                s = i+2*q
                weight = ((-1)**j * math.prod(range(i,s)) * coeffs[i]) % prime
                state['beta'][s] = (state['beta'][s]+weight) % prime
                state['constant'] = (state['constant']-weight*state['tail'][s]) % prime
            for at,x in enumerate(probes):
                inverse = pow((x-pole)%prime,-1,prime)
                state['probe_sums'][at] = (state['probe_sums'][at] + sum(
                    coeffs[i]*pow(inverse,i,prime) for i in range(1,size+1))) % prime
            inverse = pow(2*j+1,-1,prime)
            for s in range(2*q+1,37*q+1):
                state['tail'][s] = (state['tail'][s]+(-1)**j*pow(inverse,s,prime)) % prime
        if j in (0,120,240,360,480):
            print('seed prime',prime,'completed pole',j,flush=True)
    reports = []
    for q,state in states.items():
        assert all(state['beta'][s] == 0 for s in range(2*q+1,37*q+1,2))
        assert state['constant'] != 0
        for j,coeffs in enumerate(state['partials']):
            opposite = state['partials'][n-j]
            assert all(coeffs[i] == ((-1)**i*opposite[i]) % prime
                       for i in range(1,len(coeffs)))
        for x,total in zip(probes,state['probe_sums']):
            exact = pow(int(numerator(x))*pow(int(denominator(x)),-1,prime)%prime,q,prime)
            assert total == exact
        even = {str(s):state['beta'][s] for s in range(2*q+2,37*q+1,2)}
        record = SOURCE/f'power_seed_k{q}_p{prime}.json'
        compared = False
        if record.exists():
            published = json.loads(record.read_text())
            assert state['constant'] == published['constant_mod_prime']
            assert even == published['even_beta_coefficients']
            compared = True
        reports.append(dict(power=q,prime=prime,constant=state['constant'],
            even_beta_coefficients=even,poles=481,global_reconstruction_points=probes,
            reflection_coefficients_checked=sum(len(a)-1 for a in state['partials']),
            matched_published_record=compared))
    return dict(prime=prime,seconds=time.monotonic()-start,results=reports)

if __name__ == '__main__':
    output=[]
    for prime in (10007,10009):
        output.append(audit_seed(prime))
        (HERE/'fresh-seeds.json').write_text(json.dumps(output,indent=2)+'\n')
        print('Finished prime',prime,'seconds',output[-1]['seconds'],
              'constants',[(r['power'],r['constant']) for r in output[-1]['results']],flush=True)
