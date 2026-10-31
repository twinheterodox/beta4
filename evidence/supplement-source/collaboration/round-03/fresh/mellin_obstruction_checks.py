#!/usr/bin/env python3
"""Standard-library exact certificate for the Mellin positive-source obstruction.
Sturm remainders and polylogarithm tails use Fraction, not floating point.
"""
from fractions import Fraction as F
from math import gcd,lcm
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def add(a,b):
    return trim([(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))])

def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):c[j+k]+=x*y
    return trim(c)

def divmodpoly(a,b):
    a=list(map(F,trim(a)));b=list(map(F,trim(b)))
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while a!=[0] and len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]+=v
        for k,z in enumerate(b):a[j+k]-=v*z
        a=trim(a)
    return trim(q),a

def evalp(a,x):
    v=F(0)
    for c in a[::-1]:v=v*x+c
    return v

def primitive_positive_scale(a):
    den=lcm(*(x.denominator for x in a));a=[int(x*den) for x in a]
    g=gcd(*a)
    return [F(x,g) for x in a]

def variations(signs):
    signs=[x for x in signs if x]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

f=[256,-800,1088,-1160,880,-464,260,-232,220,-145,68,-25,4]
D=mul(mul([2,-2,1],[2,2,1]),mul([4,-4,2,-2,1],[16,-16,8,0,-4,0,2,-2,1]))
Ds=[[4,0,0,0,1],[8,0,0,-4,0,0,1],[32,0,0,0,0,-8,0,0,0,0,1]]
N17=[F(0)]
for den,degree,scale in zip(Ds,[2,3,5],[900,-2025,-5625]):
    q,r=divmodpoly(D,den);assert r==[0]
    N17=add(N17,[0]*degree+[scale*x for x in q])
assert N17==[0,0]+[225*x for x in f]
seq=[list(map(F,f)),[F(j*f[j]) for j in range(1,len(f))]]
while len(seq[-1])>1:
    _,r=divmodpoly(seq[-2],seq[-1]);r=[-x for x in r]
    assert r!=[0]
    seq.append(primitive_positive_scale(r))
records=[]
for x in [F(0),F(3,5),F(5,8),F(1)]:
    vals=[evalp(a,x) for a in seq]
    signs=[(v>0)-(v<0) for v in vals]
    records.append(dict(x=str(x),f=str(vals[0]),sturm_signs=signs,variations=variations(signs)))
assert records[0]['variations']-records[-1]['variations']==1
assert records[1]['variations']-records[2]['variations']==1
assert F(records[0]['f'])>0 and F(records[-1]['f'])<0
assert F(records[1]['f'])>0 and F(records[2]['f'])<0

params=[(2,F(225,136),(F(0),F(1,2)),F(1,2)),
        (3,F(-25,34),(F(1,4),F(1,4)),F(3,8)),
        (5,F(9,68),(F(1,8),F(-1,8)),F(3,16))]
N=32
moments={}
for s in [2,3]:
    center=F(0);error=F(0)
    for m,c,z,rho in params:
        assert z[0]*z[0]+z[1]*z[1]<=rho*rho<1
        a,b=F(1),F(0);partial=F(0)
        for k in range(1,N+1):
            a,b=a*z[0]-b*z[1],a*z[1]+b*z[0]
            partial+=b/F(k**s)
        weight=c*m*m if s==2 else 2*c*m
        center+=weight*partial
        error+=abs(weight)*rho**(N+1)/(F((N+1)**s)*(1-rho))
    lo,hi=center-error,center+error
    moments['M1' if s==2 else 'M2']=dict(lower=str(lo),upper=str(hi),
        lower_float=float(lo),upper_float=float(hi),tail_bound=str(error))
assert F(moments['M1']['lower'])>F(9,10)
assert F(moments['M2']['lower'])>F(193,100)
lower=(F(2,3)/4+F(6,13)/6)*F(193,100)+F(2,3)**2/12*F(9,10)
assert lower==F(1309,2600)>F(1,2)
out=dict(status='Exact certificates passed; this concerns only nonnegative constrained polynomial sources.',
   rational_kernel_numerator='225*x^2*f(x)/17',f_coefficients_ascending=f,
   denominator_coefficients_ascending=[str(x) for x in D],
   sturm=records,sturm_degrees=[len(a)-1 for a in seq],
   unique_root_bracket=['3/5','5/8'],moment_series_terms=N,moments=moments,
   elementary_log_bounds=['log(2)>2/3','-log(x0)>log(8/5)>6/13'],
   proved_positive_constant_lower=str(lower),primitive_value_lower='strictly greater than 1/2')
(HERE/'mellin_obstruction_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k in ['status','sturm','sturm_degrees','unique_root_bracket','moment_series_terms','proved_positive_constant_lower','primitive_value_lower']},indent=2))
print('Moment enclosures:',{k:(v['lower_float'],v['upper_float']) for k,v in moments.items()})
