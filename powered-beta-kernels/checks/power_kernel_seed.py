#!/usr/bin/env python3
"""Exact finite-field seed certificates for the proposed powered-kernel route.

This script does not itself establish the all-prime dilation lemma.
It computes the rational constant coordinate of the precisely specified
trimmed seed, reduced modulo an auxiliary prime exceeding every denominator.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import nmod_poly, fmpz
from beta_determinants import profile


def product(factors, cap, prime):
    factors=list(factors)
    if not factors:return nmod_poly([1],prime)
    while len(factors)>1:
        rows=[factors[i].mul_low(factors[i+1],cap) for i in range(0,len(factors)-1,2)]
        if len(factors)%2:rows.append(factors[-1])
        factors=rows
    return factors[0]


def seed(n,rn,orders,power,prime,trim=True):
    assert n%2==0 and fmpz(prime).is_prime() and prime>2*n
    assert orders==orders[::-1] and sum(orders)%2==1
    bound=n+2*rn+1-(2 if trim else 0)
    zeros=range(-bound,bound+1,2)
    poles=[n-2*k for k in range(n+1)]
    maximum=max(orders)*power
    derivative=2*power
    harmonic=[[0]*(n+1) for _ in range(maximum+derivative+1)]
    for s in range(derivative+1,maximum+derivative+1):
        for k in range(n):
            harmonic[s][k+1]=(harmonic[s][k]+(-1)**k*pow(pow(2*k+1,s,prime),-1,prime))%prime
    constants=0
    beta=[0]*(maximum+derivative+1)
    principal=[]
    jet_hash=hashlib.sha256()
    for k,a in enumerate(poles):
        cap=orders[k]*power
        if k>n//2:
            old=principal[n-k]
            aa=[0]+[(-1)**i*old[i]%prime for i in range(1,cap+1)]
        else:
            numerator=product((nmod_poly([a-z,1],prime) for z in zeros),cap,prime)
            numerator=numerator.pow_trunc(3*power,cap).mul_low(nmod_poly([a,1],prime).pow_trunc(power,cap),cap)
            denominator=product((nmod_poly([a-b,1],prime).pow_trunc(orders[kk]*power,cap)
                                 for kk,b in enumerate(poles) if kk!=k),cap,prime)
            jet=numerator.mul_low(denominator.inverse_series_trunc(cap),cap)
            aa=[0]+[int(jet[cap-i]) for i in range(1,cap+1)]
        principal.append(aa)
        jet_hash.update(json.dumps(aa,separators=(',',':')).encode())
        for i in range(1,cap+1):
            rising=1
            for j in range(derivative):rising=rising*(i+j)%prime
            weight=(-1)**k*rising*aa[i]%prime
            beta[i+derivative]=(beta[i+derivative]+weight)%prime
            constants=(constants-weight*harmonic[i+derivative][k])%prime
    assert all(beta[s]==0 for s in range(derivative+1,maximum+derivative+1,2))
    return {'n':n,'rn':rn,'power':power,'prime':prime,'trim_outer_zero_pair':trim,
            'zero_bound':bound,'pole_orders':orders,'derivative_order':derivative,
            'constant_mod_prime':constants,'constant_nonzero_certified':constants!=0,
            'even_beta_coefficients':{str(s):beta[s] for s in range(derivative+2,maximum+derivative+1,2)},
            'principal_parts_sha256':jet_hash.hexdigest()}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--power',type=int,default=2)
    p.add_argument('--prime',type=int,default=10007)
    p.add_argument('--n',type=int,default=480)
    p.add_argument('--rn',type=int,default=560)
    p.add_argument('--order',type=int,default=35)
    p.add_argument('--profile',default='krisman35')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.monotonic()
    result=seed(a.n,a.rn,profile(a.n,a.order,a.profile),a.power,a.prime)
    result['seconds']=time.monotonic()-start
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pole_orders','even_beta_coefficients']},indent=2))
