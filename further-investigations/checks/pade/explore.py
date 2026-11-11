"""Exact small-weight Kawashima--Poels polynomials and real beta forms.

Standard-library arithmetic constructs the polynomials. Optional mpmath is
used only for exploratory decimal diagnostics, never for accepting identities.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, gcd, lcm
from pathlib import Path
import hashlib
import json


def trim(a):
    a = list(a)
    while len(a)>1 and not a[-1]:
        a.pop()
    return a


def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return trim(c)


def scale(a,c): return trim([x*c for x in a])


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return trim(c)


def power(a,k):
    out=[F(1)]
    for _ in range(k): out=mul(out,a)
    return out


def ev(a,x):
    out=0
    for c in reversed(a): out=out*x+c
    return out


def shift(a,x):
    out=[F(0)]*len(a)
    for d,c in enumerate(a):
        for j in range(d+1): out[j]+=c*comb(d,j)*x**(d-j)
    return trim(out)


@lru_cache(None)
def T(k):
    out=[F(1)]
    for j in range(k): out=mul(out,[F(j)-F(k,2),F(1)])
    return tuple(out)


def Rplus(a):
    b=add(shift(a,F(1,2)),shift(a,F(-1,2)))
    q=[F(0)]*(len(b)-1)
    for j in range(len(b)-1,0,-1):
        q[j-1]=b[j]/2
        b[j]=0
        b[j-1]+=q[j-1]
    assert b[0]==0
    return trim(q)


def construct(n,ell=0,m=5):
    assert n>=0 and m%2==1 and m>=3
    H=scale(mul(power(T(n),m),T(2*ell)),F(1,factorial(n)**(m-1)))
    rem=H[:]
    a={}
    for k in range(len(H)-1,n-1,-2):
        c=rem[k] if len(rem)>k else F(0)
        if c:
            a[k]=c
            rem=add(rem,scale(T(k),-c))
    assert rem==[0], (n,ell,m,rem)
    P=[F(0)]
    for k,c in a.items(): P=add(P,scale(T(k-n),c))
    return P,H,a


@lru_cache(None)
def moments(s,N):
    # Unique formal solution f(z)+f(z+1)=(s-1)!/z^s for even s.
    assert s%2==0
    c=[]
    for j in range(N+1):
        rhs=F(factorial(s-1) if j==s-1 else 0)
        c.append((rhs-sum((-1)**(j-r)*comb(j,r)*c[r] for r in range(j)))/2)
    return tuple(c)


def Qpoly(P,s):
    c=moments(s,len(P))
    q=[F(0)]*max(1,len(P)-1)
    for d,p in enumerate(P):
        for j in range(d): q[d-j-1]+=p*c[j]
    return trim(q)


def phi(P,s):
    return sum(x*y for x,y in zip(P,moments(s,len(P))))


def primitive(v):
    den=lcm(*(x.denominator for x in v))
    ints=[int(x*den) for x in v]
    g=gcd(*ints)
    if not g: return ints
    if ints[-1]<0: g=-g
    return [x//g for x in ints]


def beta_pair(P,s,N):
    """Return constant, beta coefficient for P(z)Theta_s^+(z)-Q_s(z).

    z=N+1/2, Theta_s^+(z)=(-1)^N (s-1)! 2^s (beta(s)-S_N).
    """
    z=F(2*N+1,2)
    q=Qpoly(P,s)
    coeff=(-1)**N*factorial(s-1)*2**s*ev(P,z)
    tail=sum((F((-1)**j,(2*j+1)**s) for j in range(N)),F(0))
    const=-coeff*tail-ev(q,z)
    return const,coeff


def dump_fraction(x): return str(x)


def main():
    import mpmath as mp
    mp.mp.dps=220
    beta={s:mp.dirichlet(s,[0,1,0,-1]) for s in [2,4]}
    rows=[]
    for n in [2,4,6,8,10,12]:
        for ell in [0,1,2]:
            P,H,a=construct(n,ell)
            for s in [2,4]:
                c=moments(s,len(P)+2*n)
                assert all(sum(P[d]*c[j+d] for d in range(len(P)))==0 for j in range(2*n))
            for N in sorted({n//2,n,2*n,3*n}):
                s=4
                const,coeff=beta_pair(P,s,N)
                p,q=primitive([const,coeff])
                raw=mp.mpf(const.numerator)/const.denominator+mp.mpf(coeff.numerator)/coeff.denominator*beta[s]
                val=mp.mpf(p)+mp.mpf(q)*beta[s]
                rows.append(dict(n=n,ell=ell,N=N,constant=str(const),coefficient=str(coeff),
                    primitive_constant=str(p),primitive_coefficient=str(q),
                    log_abs_raw_approx=str(mp.log(abs(raw))),
                    log_abs_primitive_approx=str(mp.log(abs(val))),
                    log_abs_rational_error_approx=str(mp.log(abs(val/q))) if q else None))
        print('finished n=',n,flush=True)
    out=dict(scope='Exploratory decimals only; exact coefficient identities asserted.',
             normalization='Definition 4.4; Theta_s^+=(s-1)! sum (-1)^j/(z+j)^s.',
             rows=rows,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_name('exploration.json').write_text(json.dumps(out,indent=2)+'\n')
    for n in [2,4,6,8,10,12]:
        choices=[r for r in rows if r['n']==n]
        best=min(choices,key=lambda r:float(r['log_abs_primitive_approx']))
        print(n,'best:',{k:best[k] for k in ['ell','N','log_abs_raw_approx','log_abs_primitive_approx']})


if __name__=='__main__': main()
