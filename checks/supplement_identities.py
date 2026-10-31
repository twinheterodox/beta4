"""Finite exact checks of rewritten formulas; no asymptotic extrapolation."""
from fractions import Fraction as F
from math import comb, factorial, gcd, isqrt
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def primes(n):
    return [p for p in range(2,n) if all(p % q for q in range(2,isqrt(p)+1))]

def vp(x,p):
    x=F(x)
    if not x:
        return None
    a,b=x.numerator,x.denominator
    out=0
    while a%p==0: a//=p; out+=1
    while b%p==0: b//=p; out-=1
    return out

def pair(poly):
    """Formal pair of the scalar Mellin polynomial, coefficients ascending."""
    A=sum(c*(-1)**j for j,c in enumerate(poly))
    q=[]
    for j in range(len(poly)-1):
        q.append(poly[j]-(A if j==0 else q[-1]))
    assert len(poly)==1 or q[-1]==poly[-1]
    return F(A), sum((q[j]*F(1,(2*j+1)**4) for j in range(len(q))),F(0))

def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def solve(a,b):
    n=len(b)
    a=[list(map(F,row))+[F(y)] for row,y in zip(a,b)]
    for k in range(n):
        i=next(i for i in range(k,n) if a[i][k])
        a[k],a[i]=a[i],a[k]
        z=a[k][k]; a[k]=[v/z for v in a[k]]
        for i in range(n):
            if i!=k:
                z=a[i][k]; a[i]=[v-z*w for v,w in zip(a[i],a[k])]
    return [row[-1] for row in a]

def continuant_jets(n,a):
    """Integer Taylor coefficients through order three, not numerical derivatives."""
    q0,q1=[1,0,0,0],[a,1,0,0]
    r0,r1=[0,0,0,0],[1,0,0,0]
    for k in range(1,n):
        def step(y,z):
            return [a*y[j]+(y[j-1] if j else 0)+k*k*z[j] for j in range(4)]
        q0,q1=q1,step(q1,q0)
        r0,r1=r1,step(r1,r0)
    jet=[]
    for j in range(4):
        jet.append((F(r1[j])-sum(F(q1[k])*jet[j-k] for k in range(1,j+1)))/q1[0])
    return q1,r1,jet

def falling(x,j):
    z=1
    for k in range(j): z*=x-k
    return z

def AB(N,x,odd=False):
    ans=0
    for j in range(N+(not odd)):
        z=comb(2*N,2*j+int(odd))*2**j*falling(x,j)
        for k in range(N-j): z*=2*x+2*k+1
        ans+=z
    return ans

def run():
    out={}
    reductions=prime_count=0
    for n in (2,4,6,8,12,16,24):
        for u,v in ((1,1),(2,3),(5,6),(7,2)):
            poly=[0]*n+[comb(n,j)*u**(n-j)*(-v)**j for j in range(n+1)]
            A,B=pair(poly)
            Q=[comb(n,j)*v**(n-j)*u**j for j in range(n+1)]
            coeff=[sum(F(Q[j]*(-1)**(k-j),(2*(k-j)+1)**4)
                       for j in range(min(k,n)+1)) for k in range(2*n)]
            assert A==(u+v)**n
            assert B/A==-sum(coeff,F(0))/sum(Q)
            reductions+=1
            for p in primes(2*n):
                if 3*p>4*n and (u+v)%p:
                    assert vp(B,p)==-4
                    prime_count+=1
    out['quadratic']={'reversal_identities':reductions,'exact_surviving_prime_valuations':prime_count}
    # Cone witnesses at both signs have the exact prescribed leading coefficient.
    cone=0
    for n in range(1,10):
        for A in (-7,-1,1,5):
            poly=[F(abs(A)*comb(2*n,j)*(-1)**j,4**n) for j in range(2*n+1)]
            if A<0: poly=[F(0)]+poly
            assert pair(poly)[0]==A
            cone+=1
    out['cone_formal_pairs']=cone
    cf=pc=0
    for N in range(1,5):
        h=1
        for k in range(1,2*N,2): h*=k
        prefix=sum((F((-1)**k,(2*k+1)**4) for k in range(N)),F(0))
        for m in range(1,21):
            n=2*m; q,r,jet=continuant_jets(n,2*N)
            A=AB(N,m); C=AB(N,F(2*m-1,2),True)
            assert C.denominator==1
            assert F(q[0],factorial(n))==F(comb(2*m,m),4**m)*A/h
            alpha=prefix+F((-1)**(N+1),2)*jet[3]
            cf+=1
            for p in primes(2*m):
                if p>max(m,2*N) and C%p:
                    assert vp(alpha,p)==-(8+4*vp(A,p))
                    assert q[1]%p and r[0]%p
                    pc+=1
    out['continuants']={'closed_coefficient_and_quotient_cases':cf,'exact_reduced_prime_exponents':pc}
    # Independently reconstruct the 24 already archived forms with rational Schur solves.
    archived=json.loads((ROOT/'evidence/research/collaboration/round-18/primitive-check.json').read_text())
    expected={(r['n'],r['model']):(int(r['A']),int(r['B'])) for r in archived['records']}
    counts=0
    for n in range(1,13):
        d=n-1
        c=[F(1,(2*j+1)**2) for j in range(d)]
        cp=[x*x for x in c]
        B=[[F(1,(2*(i+j)+1)**2)+F(1,(2*(i+j+1)+1)**2)
            for j in range(d)] for i in range(d)]
        Bp=[[F(1,(2*(i+j)+1)**4)+F(1,(2*(i+j+1)+1)**4)
             for j in range(d)] for i in range(d)]
        vec=solve(B,c)
        derivative=2*sum((x*y for x,y in zip(cp,vec)),F(0))-sum(
            (vec[i]*Bp[i][j]*vec[j] for i in range(d) for j in range(d)),F(0))
        direct_vec=solve(Bp,cp)
        direct=sum((x*y for x,y in zip(cp,direct_vec)),F(0))
        for model,value in (('first_variation',derivative),('direct',direct)):
            # Archived B is the constant coefficient, so the form is A theta + B.
            assert (value.denominator,-value.numerator)==expected[n,model],(n,model)
            counts+=1
        P=[F(1)]+[F(0)]*d
        for j,z in enumerate(vec): P[j]-=z; P[j+1]-=z
        coeffA,coeffB=pair(mul(P,P))
        assert coeffA==1 and coeffB==-derivative
    out['ordinary_moments']={'exact_primitive_pairs_reconstructed':counts,
                             'positive_integral_formal_identities':12,
                             'new_indices':0,'new_interval_checks':0}
    # Residue in the monodromy proof, divided by pi^3.
    # -( (5i/2)^3-(i/2)^3 )/(12i) = (125-1)/(8*12).
    assert F(125-1,8*12)==F(31,24)
    out['monodromy_residue_over_pi_cubed']='31/24'
    out['status']='passed'
    out['scope']='Finite exact identities and already selected forms; not an asymptotic proof or independent model audit.'
    return out

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
