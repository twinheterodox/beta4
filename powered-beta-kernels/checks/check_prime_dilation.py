#!/usr/bin/env python3
"""Exact rational tests of the proposed prime-dilation lemma.

These finite examples test formulas; they do not replace the general proof.
Every pole is computed independently, without reflection.
"""
import hashlib
import json
import math
from pathlib import Path
import time

from flint import fmpq, fmpq_series, fmpz_poly, ctx
from beta_determinants import truncated_product


def valuation(q, p):
    if not q:
        return math.inf
    a,b=abs(int(q.numerator)),int(q.denominator)
    v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v


def residue(q,p):
    assert int(q.denominator)%p
    return int(q.numerator)%p*pow(int(q.denominator)%p,-1,p)%p


def decompose(n, rn, blocks, power, trim):
    orders=[sum(m for a,b,m in blocks if a<=j<=b) for j in range(n+1)]
    bound=n+2*rn+1-2*trim
    poles=[n-2*j for j in range(n+1)]
    zeros=list(range(-bound,bound+1,2))
    assert 1+3*len(zeros)<sum(orders)
    cap=power*max(orders);ctx.cap=cap
    beta=[fmpq(0)]*(cap+2*power+1)
    const=fmpq(0);partials=[]
    tails=[[fmpq(0)]*(n+1) for _ in beta]
    for s in range(2*power+1,len(beta)):
        for j in range(n):
            tails[s][j+1]=tails[s][j]+fmpq((-1)**j,(2*j+1)**s)
    for j,a in enumerate(poles):
        size=power*orders[j]
        num=truncated_product((fmpz_poly([a-z,1]) for z in zeros),size).pow_trunc(3*power,size)
        num=num.mul_low(fmpz_poly([a,1]).pow_trunc(power,size),size)
        den=truncated_product((fmpz_poly([a-b,1]).pow_trunc(power*orders[h],size)
                               for h,b in enumerate(poles) if h!=j),size)
        cn,cd=num.content(),den.content();num,den=num/cn,den/cd
        jet=(fmpq_series(num.coeffs(),prec=size)/fmpq_series(den.coeffs(),prec=size))*fmpq(cn,cd)
        aa=[fmpq(0)]+[jet[size-i] for i in range(1,size+1)]
        partials.append(aa)
        for i in range(1,size+1):
            s=i+2*power
            weight=(-1)**j*math.prod(range(i,s))*aa[i]
            beta[s]+=weight;const-=weight*tails[s][j]
    assert all(not beta[s] for s in range(2*power+1,len(beta),2))
    # Check the full rational-function identity, independently of local jets.
    for x in [fmpq(1,3),fmpq(5,7)]:
        actual=x**power
        for z in zeros:actual*=(x-z)**(3*power)
        for a,b in zip(poles,orders):actual/=(x-a)**(power*b)
        reconstructed=sum((aa[i]/(x-a)**i for a,aa in zip(poles,partials)
                           for i in range(1,len(aa))),fmpq(0))
        assert actual==reconstructed
    return {'partials':partials,'beta':beta,'constant':const,'orders':orders,
            'zeros':zeros,'poles':poles}


def check(N, rN, blocks, power, p):
    start=time.monotonic();L=N+2*rN
    assert p>N+L+1 and p>2*N
    S=sum(m for a,b,m in blocks);W=sum((b-a)*m for a,b,m in blocks)
    GN=W-3*L;A=GN+S-1;tau=power*(A+2)
    seed=decompose(N,rN,blocks,power,True)
    big=decompose(N*p,rN*p,[(a*p,b*p,m) for a,b,m in blocks],power,False)
    unit=1
    for z in big['zeros']:
        if z%p:unit=unit*pow(-z,3,p)%p
    for a,order in zip(big['poles'],big['orders']):
        if a%p:unit=unit*pow(pow(-a,order,p),-1,p)%p
    selected=nonselected=0
    for j,aa in enumerate(big['partials']):
        for i in range(1,len(aa)):
            if j%p==0:
                scaled=aa[i]*fmpq(p)**(power*A-i)
                assert valuation(scaled,p)>=0
                assert residue(scaled,p)==pow(unit,power,p)*residue(seed['partials'][j//p][i],p)%p
                selected+=1
            else:
                assert valuation(aa[i],p)>=-power*GN+i
                nonselected+=1
    for q in big['beta']:
        assert valuation(q*fmpq(p)**tau,p)>=1
    lhs=residue(big['constant']*fmpq(p)**tau,p)
    rhs=pow(unit,power,p)*(-1)**((p-1)//2)*residue(seed['constant'],p)%p
    assert lhs==rhs
    return {'status':'passed','seed_n':N,'seed_rn':rN,'blocks':blocks,'power':power,
            'prime':p,'A':A,'tau':tau,'unit_mod_prime':unit,
            'seed_constant_mod_prime':residue(seed['constant'],p),
            'scaled_dilated_constant_mod_prime':lhs,'selected_coefficients_checked':selected,
            'nonselected_coefficient_bounds_checked':nonselected,
            'all_scaled_beta_coefficients_vanish':True,'full_partial_fraction_spot_checks':4,
            'seconds':time.monotonic()-start}


if __name__=='__main__':
    results=[]
    cases=[(2,2,[(0,2,17)],1,11),(2,2,[(0,2,17)],2,11),
           (2,2,[(0,2,17)],2,13),(4,4,[(0,4,9),(1,3,8)],2,19)]
    for case in cases:
        result=check(*case);results.append(result);print(json.dumps(result),flush=True)
    output={'checks':results,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Finite exact rational checks, not a proof of the all-prime assertion.'}
    Path(__file__).with_name('PRIME_DILATION_VALIDATION.json').write_text(json.dumps(output,indent=2)+'\n')
