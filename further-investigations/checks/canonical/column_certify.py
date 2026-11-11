"""Rational trial coefficients, exact Sturm counts, and Arb maximum bounds."""
from pathlib import Path
from fractions import Fraction
import argparse,json,time
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
HERE=Path(__file__).resolve().parent

def primitive(P):
    z=P.numer();c=z.content()
    return Poly([Q(v,c) for v in z])

def sturm(P):
    ss=[primitive(P),primitive(P.derivative())]
    while not ss[-1].is_zero():
        rem=-(ss[-2]%ss[-1])
        if rem.is_zero():break
        ss.append(primitive(rem))
    assert ss[-1].degree()==0,'Repeated roots require separate treatment.'
    return ss

def isolate(P,domains,bits=80):
    ss=sturm(P);cache={}
    def variation(x):
        if x not in cache:
            signs=[1 if z>0 else -1 for f in ss if (z:=f(x))!=0]
            cache[x]=sum(a!=b for a,b in zip(signs,signs[1:]))
        return cache[x]
    roots=[];counts=[]
    for left,right in domains:
        assert P(left)!=0 and P(right)!=0
        count=variation(left)-variation(right);counts.append(count)
        stack=[(left,right,count)]
        while stack:
            lo,hi,c=stack.pop()
            if c==0:continue
            if c==1 and hi-lo<Q(1,2**bits):roots.append((lo,hi));continue
            mid=(lo+hi)/2
            for j in range(3,P.degree()+5):
                if P(mid)!=0:break
                mid=(lo*(j-1)+hi)/j
            assert P(mid)!=0
            lc=variation(lo)-variation(mid);rc=variation(mid)-variation(hi)
            assert lc+rc==c
            stack.extend([(lo,mid,lc),(mid,hi,rc)])
    return sorted(roots),counts,len(cache)

def cheb_sum(co,xx):
    # co[k-1] multiplies T_k / k.
    angle=xx.acos();val=arb(0)
    for k,c in enumerate(co,1):
        val+=arb(c)*(k*angle).cos()/k
    return val

def run(degree):
    start=time.monotonic();ctx.prec=512
    data=json.loads((HERE/f'column_opt_d{degree}.json').read_text())
    def rat(v):
        f=Fraction(f'{v:.12f}');return Q(f.numerator,f.denominator)
    p=[rat(z) for z in data['p']];v=[rat(z) for z in data['v']]
    B,H,g=Q(7,48),Q(2,48),Q(4,48)
    t=Poly([0,1]);U=[Poly([1]),2*t]
    for j in range(2,degree):U.append(2*t*U[-1]-U[-2])
    dp=sum((2*p[k-1]*U[k-1] for k in range(1,degree+1)),Poly([]))
    dv=Poly(v)
    NX=B*(1-t)*(1+t*t)-2*H*t*(1+t*t)-2*(B+H)*t*t*(1-t)-(dp+dv)*t*(1-t)*(1+t*t)
    NY=B*(1-t)-g*t+2*sum((v[k-1]*U[k-1] for k in range(1,degree+1)),Poly([]))*t*(1-t)
    xr,xc,xev=isolate(NX,[(Q(-1),Q(0)),(Q(0),Q(1))]);yr,yc,yev=isolate(NY,[(Q(0),Q(1))])
    (HERE/f'column_roots_d{degree}.json').write_text(json.dumps(dict(
        p=[str(z) for z in p],v=[str(z) for z in v],x_root_counts=xc,s_root_counts=yc,
        x_root_brackets=[[str(a),str(b)] for a,b in xr],s_root_brackets=[[str(a),str(b)] for a,b in yr]),indent=2)+'\n')
    print(json.dumps(dict(stage='Sturm isolation',x_counts=xc,s_counts=yc,seconds=time.monotonic()-start)),flush=True)
    def fx(x):
        pp=sum((arb(vv)*x**(i+1)/(i+1) for i,vv in enumerate(v)),arb(0))
        return arb(B)*abs(x).log()+2*arb(H)*(1-x).log()-arb(B+H)*(1+x*x).log()-2*cheb_sum(p,x)-pp
    def fy(s):return arb(B)*s.log()+arb(g)*(1-s).log()+2*cheb_sum(v,s)
    xv=[fx(arb(-1))]+[fx(arb(lo).union(arb(hi))) for lo,hi in xr]
    yv=[fy(arb(lo).union(arb(hi))) for lo,hi in yr]
    # Global maxima occur at these critical points or the regular x=-1 endpoint.
    # The other endpoint limits are -infinity because B,H,g>0.
    mx=max(float(z.upper()) for z in xv);my=max(float(z.upper()) for z in yv)
    # Round upward on a rational 10^-10 grid and verify every interval against it.
    import math
    X=Q(math.ceil(mx*10**10)+1,10**10);Y=Q(math.ceil(my*10**10)+1,10**10)
    assert all(z<arb(X) for z in xv) and all(z<arb(Y) for z in yv)
    norm=sum((p[k-1]**2/Q(k)+v[k-1]**2/Q(2*k) for k in range(1,degree+1)),Q(0))
    value=arb(B-1)*arb(2).log()+arb(norm+X+Y)
    out=dict(degree=degree,p=[str(z) for z in p],v=[str(z) for z in v],
        derivative_x=[str(z) for z in primitive(NX)],derivative_s=[str(z) for z in primitive(NY)],
        x_root_counts=xc,s_root_counts=yc,
        x_root_brackets=[[str(a),str(b)] for a,b in xr],s_root_brackets=[[str(a),str(b)] for a,b in yr],
        x_value_enclosures=[z.str(30) for z in xv],s_value_enclosures=[z.str(30) for z in yv],
        rational_x_upper=str(X),rational_s_upper=str(Y),rational_norm=str(norm),
        half_column_logdet_upper=value.str(40),full_column_logdet_upper=(2*value).str(40),
        precision_bits=ctx.prec,sturm_evaluation_points=xev+yev,seconds=time.monotonic()-start,
        scope='Exact rational dual functions; exhaustive critical-point counts; rigorous Arb evaluations. Uniform energy theorem is justified separately.')
    (HERE/f'column_certificate_d{degree}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['degree','rational_x_upper','rational_s_upper','half_column_logdet_upper','full_column_logdet_upper','seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degree',type=int,default=100);a=p.parse_args();run(a.degree)
