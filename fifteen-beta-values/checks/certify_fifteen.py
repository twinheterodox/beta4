#!/usr/bin/env python3
"""Exact arithmetic and real-phase certificates for beta(4) through beta(32).

The mathematical denominator improvement and asymptotic passage are
proved in the accompanying manuscript. This script certifies the finite inputs.
"""
from fractions import Fraction as Q
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
DS=[Q(j,96) for j in range(31)]
R=Q(3,4);G=Q(221,16);X=400;COMMON=96;BITS=64

def events():
    offsets=set([(d,Q(0)) for d in DS]+[(1-d,Q(0)) for d in DS]
                +[(-R,Q(-1,2)),(1+R,Q(-1,2))])
    out={Q(k) for k in range(1,X+1)}
    for (s,t),(ss,tt) in combinations(sorted(offsets),2):
        ds,dt=ss-s,t-tt
        if not ds:continue
        for k in range(math.floor(ds-dt)-1,math.ceil(X*ds-dt)+2):
            z=(dt+k)/ds
            if 1<z<X:out.add(z)
    return sorted(out)

def e_star(x):
    a,b=x.numerator,x.denominator;den=COMMON*b
    base=0;jumps={}
    def add(num,jump):
        residue=num%den
        if residue:jumps[residue]=jumps.get(residue,0)+jump
    for j in range(31):
        lo=j*a;hi=(COMMON-j)*a
        base+=(-lo)//den-(-hi)//den
        add(lo,1);add(hi,-1)
    base-=3*((72*a+48*b)//den-(-168*a+48*b)//den)
    add(-72*a-48*b,-3);add(168*a-48*b,3)
    cap=COMMON*(a-b) if a<2*b else den
    current=maximum=base
    for pos,jump in sorted(jumps.items()):
        if pos>cap:break
        current+=jump;maximum=max(maximum,current)
    return Q(max(0,maximum+2))-G*(a//b)

def log_bounds(x,terms=40):
    x=Q(x);assert x>0;k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    def basic(y):
        z=(y-1)/(y+1)
        lo=2*sum((z**(2*j+1)/Q(2*j+1) for j in range(terms)),Q(0))
        return lo,lo+2*z**(2*terms+1)/(Q(2*terms+1)*(1-z*z))
    lo,hi=basic(x);ll,hh=basic(Q(2))
    return (lo+k*ll,hi+k*hh) if k>=0 else (lo+k*hh,hi+k*ll)

def linear_logs(terms):
    low=high=Q(0)
    for c,x in terms:
        a,b=log_bounds(x)
        low+=c*(a if c>=0 else b);high+=c*(b if c>=0 else a)
    return low,high

def outward(a,b):
    scale=10**20
    return [str(Q((a*scale).numerator//(a*scale).denominator,scale)),
            str(Q(-((-b*scale).numerator//(-b*scale).denominator),scale))]

def phase_certificate():
    cut=[R,Q(31,40),Q(4,5),Q(17,20)];curvature=[]
    for a,b in zip(cut,cut[1:]):
        bound=3/(a+1+R)-3/(b-R)-sum(1/(b+1-d)-1/(a+d) for d in DS)
        assert bound<0;curvature.append(str(bound))
    x=Q(7573,10000)
    value=[(3*(x+1+R),x+1+R),(-3*(x-R),x-R)]
    slope=[(Q(3),x+1+R),(Q(-3),x-R)]
    for d in DS:
        value.extend([(-(x+1-d),x+1-d),(x+d,x+d)])
        slope.extend([(Q(-1),x+1-d),(Q(1),x+d)])
    vl,vh=linear_logs(value);sl,sh=linear_logs(slope)
    assert vh<Q(-3767,200) and sl>Q(-1,50) and sh<Q(1,50)
    V=Q(-1883,100)
    tangent=Q(-3767,200)+Q(1,50)*max(x-R,Q(17,20)-x)
    assert tangent<V
    gamma=sum(1-2*d for d in DS)
    # g(x)=gamma/(x+1)-3 log((x+1+r)/(x-r)). Its derivative has
    # numerator -G*x^2+(6*(1+2r)-gamma)*x+3*(1+2r)+gamma*r*(1+r).
    # That numerator is strictly decreasing for x>=0. Thus g has no
    # interior minimum there; g(17/20)>0 and g(infinity)=0 imply g>0.
    assert 6*(1+2*R)-gamma<0
    a=Q(17,20)
    gl,gh=linear_logs([(Q(-3),(a+1+R)/(a-R))]);gl+=gamma/(a+1);gh+=gamma/(a+1)
    assert gl>0
    tail_slope=3*(1+2*R)*3/(3-R)-gamma*3/4
    assert tail_slope<-3
    ql,qh=linear_logs([(Q(15),Q(2))]);ql+=V+Q(41,5);qh+=V+Q(41,5)
    assert qh<Q(-23,100)
    return dict(V=str(V),curvature_upper_bounds=curvature,x0=str(x),
                phase_at_x0=outward(vl,vh),slope_at_x0=outward(sl,sh),
                tangent_upper=str(tangent),g_at_seventeen_twentieths=outward(gl,gh),
                tail_slope_upper=str(tail_slope),cleared_exponent=outward(ql,qh),
                cleared_exponent_strict_upper='-23/100')

def run():
    start=time.monotonic();ee=events();total=0;digest=hashlib.sha256()
    trace=HERE/'fifteen-arithmetic-midpoints.tsv'
    with trace.open('w') as f:
        for a,b in zip(ee,ee[1:]):
            mid=(a+b)/2;v=e_star(mid);term=v*(1/a-1/b)*2**BITS
            rounded=-((-term.numerator)//term.denominator);total+=rounded
            line=f'{mid.numerator}\t{mid.denominator}\t{v.numerator}\t{v.denominator}\n'
            f.write(line);digest.update(line.encode())
    body=Q(total,2**BITS);tail=Q(36,X)+G/Q(2*(X-1));A=body+tail
    assert A<Q(41,5),(float(A),A)
    phase=phase_certificate()
    result=dict(status='passed',event_points=len(ee),intervals=len(ee)-1,
                arithmetic_body_upper=str(body),arithmetic_tail_upper=str(tail),
                arithmetic_total_upper=str(A),arithmetic_total_decimal=float(A),
                arithmetic_simplified_upper='41/5',rounding_bits=BITS,
                midpoint_trace_sha256=digest.hexdigest(),phase=phase,
                seconds=time.monotonic()-start,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'fifteen-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':run()
