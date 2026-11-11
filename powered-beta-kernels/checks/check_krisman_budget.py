#!/usr/bin/env python3
"""Independent exact checks of the arithmetic and real bounds used here.

Standard library only. Reconstructs every rational arithmetic event, rounds
each integral term upward exactly, and checks the elementary phase proof
with rational bounds for logarithms. No downloaded code is imported.
"""
from fractions import Fraction as Q
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
import time

BLOCKS=[(Q(0),8)]+[(Q(a,160),1) for a in list(range(1,9))+list(range(10,47,2))]
R=Q(7,6);G=Q(179,10);V=Q(-744,25)
A=Q(12507029707866288048422972493354331,790300921079786767495600904601600)
U=1243345293752876677540040193981


def events():
    offsets=set([(d,Q(0)) for d,m in BLOCKS]+[(1-d,Q(0)) for d,m in BLOCKS]
                +[(-R,Q(-1,2)),(1+R,Q(-1,2)),(Q(0),Q(1,2))])
    out={Q(k) for k in range(1,401)}
    out.update(Q(1)+Q(j,20) for j in range(40))
    out.update(Q(3)+Q(j,2) for j in range(795))
    for (s,t),(ss,tt) in combinations(sorted(offsets),2):
        ds,dt=ss-s,t-tt
        if not ds:continue
        for k in range(math.floor(ds-dt)-1,math.ceil(400*ds-dt)+2):
            x=(dt+k)/ds
            if 1<x<400:out.add(x)
    return sorted(out)


def e_star(x):
    a,b=x.numerator,x.denominator;den=480*b
    base=2;offsets={}
    def add_offset(num,jump):
        residue=num%den
        if residue:offsets[residue]=offsets.get(residue,0)+jump
    for d,m in BLOCKS:
        low=int(480*d)*a;high=(480-int(480*d))*a
        base+=m*((-low)//den-(-high)//den)
        add_offset(low,m);add_offset(high,-m)
    base-=3*((560*a+240*b)//den-(-1040*a+240*b)//den)
    add_offset(-560*a-240*b,-3)
    add_offset(1040*a-240*b,3)
    current=gain=0
    for pos,jump in sorted(offsets.items()):
        current+=jump;gain=max(gain,current)
    return Q(base+gain)-G*(a//b)


def log_bounds(x,terms=30):
    x=Q(x);assert x>0
    k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    def basic(y):
        z=(y-1)/(y+1)
        lo=2*sum((z**(2*j+1)/Q(2*j+1) for j in range(terms)),Q(0))
        return lo,lo+2*z**(2*terms+1)/(Q(2*terms+1)*(1-z*z))
    low,high=basic(x);ll,hh=basic(Q(2))
    if k>=0:return low+k*ll,high+k*hh
    return low+k*hh,high+k*ll


def linear_logs(terms):
    low=high=Q(0)
    for c,x in terms:
        ll,hh=log_bounds(x)
        if c>=0:low+=c*ll;high+=c*hh
        else:low+=c*hh;high+=c*ll
    return low,high


def compact_interval(low,high):
    """Outward dyadic rounding for a small, still rigorous saved record."""
    scale=2**96
    a=low*scale;b=high*scale
    return [str(Q(a.numerator//a.denominator,scale)),
            str(Q(-((-b.numerator)//b.denominator),scale))]


def phase_certificate():
    intervals=[(R,Q(6,5)),(Q(6,5),Q(5,4)),(Q(5,4),Q(4,3))]
    curvature=[]
    for a,b in intervals:
        upper=3/(a+1+R)-3/(b-R)-sum(m*(1/(b+1-d)-1/(a+d)) for d,m in BLOCKS)
        assert upper<0;curvature.append(str(upper))
    x=Q(58911,50000)
    phi=[(3*(x+1+R),x+1+R),(-3*(x-R),x-R)]
    slope=[(Q(3),x+1+R),(Q(-3),x-R)]
    for d,m in BLOCKS:
        phi.extend([(-m*(x+1-d),x+1-d),(m*(x+d),x+d)])
        slope.extend([(-m,x+1-d),(m,x+d)])
    pl,ph=linear_logs(phi);dl,dh=linear_logs(slope)
    assert ph<Q(-297613,10000) and Q(-1,500)<dl and dh<0
    tangent=Q(-297613,10000)+Q(1,500)*(x-R)
    assert tangent<V
    sb2=sum(m*(1-2*d)/(2+1-d) for d,m in BLOCKS)
    sb3=sum(m*(1-2*d)/(3+1-d) for d,m in BLOCKS)
    assert sb2>Q(191,20) and 3*log_bounds(24)[1]<Q(477,50)
    assert sb3>7 and 3*log_bounds(6)[1]<Q(27,5)
    assert Q(180,11)-Q(837,40)<-4
    qlo,qhi=linear_logs([(Q(20),Q(2))]);qlo+=V+A;qhi+=V+A
    assert qhi<Q(-714,10000)
    return {'status':'passed','curvature_upper_bounds':curvature,
            'phi_x0_interval':compact_interval(pl,ph),'slope_x0_interval':compact_interval(dl,dh),
            'tangent_max_upper':str(tangent),'certified_V':str(V),
            'Q_interval':compact_interval(qlo,qhi),'Q_strict_upper':str(Q(-714,10000)),
            'method':'Rational atanh series with explicit positive remainder, 30 terms'}


if __name__=='__main__':
    start=time.monotonic();ee=events();assert len(ee)==62132
    print('Reconstructed 62,132 exact event points.',flush=True)
    total=0;digest=hashlib.sha256()
    for a,b in zip(ee,ee[1:]):
        v=e_star((a+b)/2);term=(1/a-1/b)*v*2**96
        ceil=-((-term.numerator)//term.denominator)
        total+=ceil;digest.update(f'{a}|{b}|{v}|{ceil}\n'.encode())
    assert total==U,(total,U)
    assert Q(U,2**96)+Q(44,400)+G/Q(798)==A
    phase=phase_certificate()
    result={'status':'passed','event_points':len(ee),'intervals':len(ee)-1,
            'rounded_arithmetic_sum':total,'rounding_denominator':str(2**96),
            'A_cert':str(A),'event_value_sha256':digest.hexdigest(),'phase':phase,
            'seconds':time.monotonic()-start,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact finite certificates; the prime-number-theorem transfer and elementary inequalities are mathematical arguments in the cited paper and the draft.'}
    Path(__file__).with_name('KRISMAN_BUDGET_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='phase'},indent=2),flush=True)
