#!/usr/bin/env python3
"""Exact polynomial and logarithm certificate for the manuscript's phase.

The polynomial test was developed in the second research session and is
specialized here to r=3/4. This program imports no other checker.
"""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def evaluate(coeff,x):
    answer=Q(0)
    for c in reversed(coeff):answer=answer*x+c
    return answer

def logs(x):
    k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    def unit(y):
        w=(y-1)/(y+1);power=w;s=Q(0)
        for j in range(48):s+=2*power/(2*j+1);power*=w*w
        return s,s+2*power/(97*(1-w*w))
    a,b=unit(x);c,d=unit(Q(2))
    return (a+k*c,b+k*d) if k>=0 else (a+k*d,b+k*c)

def h_bounds(a,b):
    lo=logs(a)[0];hi=logs(b)[1]
    products=[a*lo,a*hi,b*lo,b*hi]
    return min(products),max(products)

def run():
    a=[1];b=[8]
    for _ in range(3):a=mul(a,[5,2]);b=mul(b,[0,1])
    for j in range(31):a=mul(a,[72+j,96]);b=mul(b,[168-j,96])
    coeff=[x-y for x,y in zip(a,b)]
    while coeff[-1]==0:coeff.pop()
    assert len(coeff)==34 and all(c>0 for c in coeff[:3]) and all(c<0 for c in coeff[3:])
    lo=Q(0);hi=Q(1)
    assert evaluate(coeff,lo)>0 and evaluate(coeff,hi)<0
    for _ in range(40):
        mid=(lo+hi)/2
        if evaluate(coeff,mid)>0:lo=mid
        else:hi=mid
    scale=10**9;r=Q(3,4)
    left=Q(((lo+r)*scale).__floor__(),scale)
    right=Q(((hi+r)*scale).__ceil__(),scale)
    assert evaluate(coeff,left-r)>0 and evaluate(coeff,right-r)<0
    upper=3*(h_bounds(left+1+r,right+1+r)[1]-h_bounds(left-r,right-r)[0])
    for j in range(31):
        d=Q(j,96)
        upper-=h_bounds(left+1-d,right+1-d)[0]-h_bounds(left+d,right+d)[1]
    assert upper<Q(-1883,100)
    rounded=Q((upper*10**12).__ceil__(),10**12)
    tail=Q(3)*Q(5,2)*3/(3-r)-Q(341,16)*Q(3,4)
    assert tail==Q(-383,64) and tail<-3
    out={'status':'passed','degree':len(coeff)-1,'coefficients_ascending':[str(c) for c in coeff],
         'sign_pattern':'+++-------------------------------','sign_changes':1,
         'critical_point_bracket':[str(left),str(right)],
         'maximum_upper':str(rounded),'maximum_upper_decimal':float(rounded),
         'simplified_phase_upper':'-1883/100','tail_x_phi_prime_upper':str(tail),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('phase-polynomial-certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='coefficients_ascending'},indent=2))

if __name__=='__main__':run()
