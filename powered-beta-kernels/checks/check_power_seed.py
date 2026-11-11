#!/usr/bin/env python3
"""Independent standard-library-only checker for powered-kernel seed residues.

Uses logarithmic derivatives and harmonic/factorial prefix sums, not FLINT
polynomial products or reflected jets. Every pole is computed separately.
"""
import argparse
import json
import math
from pathlib import Path
import time


def is_prime(p):
    return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))


def check(record):
    started=time.monotonic()
    n,q,p=record['n'],record['power'],record['prime']
    bound=record['zero_bound'];orders=record['pole_orders'];cap=max(orders)*q
    # Reconstruct the published profile from its defining parameters, rather
    # than trusting the input record's pole orders or endpoint convention.
    assert n==480 and record['rn']==560 and bound==1599
    assert record['trim_outer_zero_pair'] is True and record['derivative_order']==2*q
    kappas=list(range(1,9))+list(range(10,47,2))
    expected_orders=[8+sum(3*k<=j<=480-3*k for k in kappas) for j in range(481)]
    assert orders==expected_orders
    assert is_prime(p) and p>max(n+bound,2*n,cap+2*q)
    assert n%2==0 and bound%2==1 and orders==orders[::-1]
    blocks=[];previous=0
    for a in range(n//2+1):
        weight=orders[a]-previous
        assert weight>=0
        if weight:blocks.append((a,n-a,weight))
        previous=orders[a]

    def prefixes(values,maximum):
        table=[[0]*(len(values)+1) for _ in range(maximum+1)]
        for j,value in enumerate(values,1):
            inv=pow(value,-1,p);power=1
            for m in range(1,maximum+1):
                power=power*inv%p
                table[m][j]=(table[m][j-1]+power)%p
        return table

    positive=prefixes(list(range(1,n+1)),cap)
    odd_count=(n+bound+1)//2
    odd=prefixes(list(range(1,2*odd_count,2)),cap)
    factorial=[1]
    for j in range(1,n+1):factorial.append(factorial[-1]*j%p)
    odd_factorial=[1]
    for j in range(odd_count):odd_factorial.append(odd_factorial[-1]*(2*j+1)%p)
    half_powers=[1]
    for j in range(cap):half_powers.append(half_powers[-1]*pow(2,-1,p)%p)

    maximum_s=cap+2*q
    assert set(record['even_beta_coefficients'])==set(map(str,range(2*q+2,maximum_s+1,2)))
    tails=[[0]*(n+1) for _ in range(maximum_s+1)]
    for j in range(n):
        inv=pow(2*j+1,-1,p);power=1
        for s in range(1,maximum_s+1):
            power=power*inv%p
            tails[s][j+1]=(tails[s][j]+(-1)**j*power)%p
    constant=0;beta=[0]*(maximum_s+1)
    for k in range(n+1):
        a=n-2*k;size=q*orders[k];shift=q if a==0 else 0
        pos=(a+bound+1)//2;neg=(bound-a+1)//2
        numerator=pow((-1)**neg*odd_factorial[pos]*odd_factorial[neg]%p,3*q,p)
        if a:numerator=numerator*pow(a,q,p)%p
        denominator=1
        for left,right,multiplicity in blocks:
            block_value=1
            if right>k:
                lo=max(left-k,1);hi=right-k
                block_value=block_value*factorial[hi]*pow(factorial[lo-1],-1,p)%p
            if left<k:
                lo=max(k-right,1);hi=k-left
                block_value=block_value*((-1)**(hi-lo+1))*factorial[hi]*pow(factorial[lo-1],-1,p)%p
            count=right-left+1-(left<=k<=right)
            block_value=block_value*pow(2,count,p)%p
            denominator=denominator*pow(block_value,multiplicity*q,p)%p
        g0=numerator*pow(denominator,-1,p)%p
        log_derivative=[0]*size
        for m in range(1,size-shift):
            num_sum=(odd[m][pos]+(-1)**m*odd[m][neg])%p
            den_sum=0
            for left,right,multiplicity in blocks:
                if left<=k<=right:
                    subtotal=positive[m][right-k]+(-1)**m*positive[m][k-left]
                elif k<left:
                    subtotal=positive[m][right-k]-positive[m][left-k-1]
                else:
                    subtotal=(-1)**m*(positive[m][k-left]-positive[m][k-right-1])
                den_sum=(den_sum+multiplicity*subtotal)%p
            ee=q*pow(pow(a,m,p),-1,p) if a else 0
            log_derivative[m]=(-1)**(m-1)*(ee+3*q*num_sum-q*half_powers[m]*den_sum)%p
        coefficients=[g0]
        for j in range(1,size-shift):
            coefficients.append(sum(log_derivative[m]*coefficients[j-m] for m in range(1,j+1))*pow(j,-1,p)%p)
        jet=[0]*shift+coefficients
        for i in range(1,size+1):
            rising=math.prod(range(i,i+2*q))%p
            weight=(-1)**k*rising*jet[size-i]%p
            beta[i+2*q]=(beta[i+2*q]+weight)%p
            constant=(constant-weight*tails[i+2*q][k])%p
    assert all(beta[s]==0 for s in range(2*q+1,maximum_s+1,2))
    assert constant==record['constant_mod_prime'],(constant,record['constant_mod_prime'])
    for s,v in record['even_beta_coefficients'].items():assert beta[int(s)]==v,(s,beta[int(s)],v)
    assert constant!=0
    return {'status':'passed','n':n,'power':q,'prime':p,'constant_mod_prime':constant,
            'checked_all_even_coefficients':True,'computed_each_pole_without_reflection':True,
            'method':'pure Python logarithmic-derivative recurrence and prefix sums',
            'seconds':time.monotonic()-started}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    results=[]
    for path in a.records:
        result=check(json.loads(path.read_text()));result['source']=str(path)
        results.append(result);print(json.dumps(result),flush=True)
    a.output.write_text(json.dumps(results,indent=2)+'\n')
