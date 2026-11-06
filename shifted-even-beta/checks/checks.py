"""Exact self-checks for multiple odd grids. Standard library only; no -O."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, gcd, prod
from pathlib import Path
import json


def chi(j):
    return 0 if j % 2 == 0 else (-1)**((j-1)//2)


def lcm(x, y):
    return x*y//gcd(x,y)


def lcm_to(n):
    out=1
    for j in range(1,n+1): out=lcm(out,j)
    return out


def mul(p,q,limit=None):
    length=len(p)+len(q)-1
    if limit is not None: length=min(length,limit+1)
    out=[0]*length
    for i,x in enumerate(p[:length]):
        for j,y in enumerate(q[:length-i]): out[i+j]+=x*y
    return out


def add(p,q):
    out=[0]*max(len(p),len(q))
    for i,x in enumerate(p):out[i]+=x
    for i,x in enumerate(q):out[i]+=x
    return out


def power(p,n,limit=None):
    out=[1]
    for _ in range(n):out=mul(out,p,limit)
    return out


def evaluate(p,x):
    out=0
    for c in reversed(p):out=out*x+c
    return out


def primes(n):
    out=[]
    for p in range(2,n+1):
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
    return out


@lru_cache(None)
def parameters(B,r):
    L=2*r+1
    theta=tuple(F(v,d) for d in range(1,2*B,2)
                for v in range(1,2*d,2) if gcd(v,d)==1)
    K=prod(((2*t.denominator)**L
            *prod(p**((L+p-2)//(p-1)) for p in primes(t.denominator)))
           for t in theta)
    integer_K=K//prod(t.denominator**L for t in theta)
    assert K==integer_K*prod(t.denominator**L for t in theta)
    return theta,len(theta),L,K,integer_K


@lru_cache(None)
def kernels(B,m):
    """Return integral h_b and the positive kernel distances."""
    c=(m+1)//2
    distances=tuple(sorted({F(2*j,b) for b in range(1,2*B,2)
                            for j in range(1,(b-1)//2+1)}))
    blocks={t:power([t.numerator**2,t.denominator**2],m) for t in distances}
    full=[1]
    for block in blocks.values():full=mul(full,block)
    out=[]
    for b in range(1,2*B,2):
        h=[chi(b)*(-1)**c*x for x in full]
        for j in range(1,(b-1)//2+1):
            t=F(2*j,b);v,d=t.numerator,t.denominator
            numerator=[-2*chi(b-2*j)*(-1)**ell*comb(m,2*ell+1)
                       *v**(m-2*ell-1)*d**(m+2*ell+1) for ell in range(c)]
            term=[0]*c+numerator
            for alpha,block in blocks.items():
                if alpha!=t:term=mul(term,block)
            h=add(h,term)
        assert len(h)<=m*len(distances)+1
        out.append(tuple(h))
    return tuple(out),distances


def cmul(z,w):
    return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]


def cpow(z,n):
    out=(F(1),F(0))
    for _ in range(n):out=cmul(out,z)
    return out


def cinv(z):
    norm=z[0]**2+z[1]**2
    return z[0]/norm,-z[1]/norm


@lru_cache(None)
def base_partials(B,m,r,a,n):
    assert m%2==n%2==1 and a%2==0
    theta,D,L,K,integer_K=parameters(B,r)
    assert a>m*D*L
    delta=2*lcm_to(n);rows=[]
    divisor=factorial(n)**(m*D*L)
    scale=integer_K**(m*n)
    for h in range(n+1):
        local=[1]
        for t in theta:
            v,d=t.numerator,t.denominator
            for j in range(L*n):
                constant=d*(-2*h-2*r*n+2*j)+v
                block=[comb(m,ell)*constant**(m-ell)*(d*delta)**ell
                       for ell in range(m+1)]
                local=mul(local,block,a-1)
        local += [0]*(a-len(local))
        normalized=[]
        for x in local:
            quotient,remainder=divmod(scale*x,divisor)
            assert remainder==0
            normalized.append(quotient)
        local=normalized
        for ell in range(n+1):
            if ell==h:continue
            ratio=delta//(2*(ell-h))
            assert ratio*2*(ell-h)==delta
            local=mul(local,[comb(a+j-1,j)*(-ratio)**j for j in range(a)],a-1)
        rows.append(tuple(x*comb(n,h)**a for x in local))
    return tuple(rows)


def multiplier_local(B,m,n,h,hpoly,delta,limit):
    J=m*(parameters(B,1)[1]-1)//2
    out=[0]*(min(2*J,limit)+1)
    for j,c in enumerate(hpoly):
        for ell in range(min(2*j,limit)+1):
            out[ell]+=c*(-1)**j*n**(2*J-2*j)*comb(2*j,ell)*(n-2*h)**(2*j-ell)*delta**ell
    return out


def partials(B,m,r,a,n,hpoly):
    delta=2*lcm_to(n);rows=[];integer_rows=[]
    for h,base in enumerate(base_partials(B,m,r,a,n)):
        local=mul(list(base),multiplier_local(B,m,n,h,hpoly,delta,a-1),a-1)
        local += [0]*(a-len(local))
        integer_rows.append(local)
        rows.append(tuple(F(local[a-i],delta**(a-i)) for i in range(1,a+1)))
    return tuple(rows),integer_rows


def forms(B,m,r,a,n,hpoly):
    rows,integer_rows=partials(B,m,r,a,n,hpoly)
    k=m-1;delta=2*lcm_to(n);qexp=a+k
    A=tuple(tuple(comb(i+k-1,k)*row[i-1] for i in range(1,a+1)) for row in rows)
    rho=tuple(sum((-1)**h*A[h][i-1] for h in range(n+1)) for i in range(1,a+1))
    for i,value in enumerate(rho,1):
        assert (delta**(a-i)*value).denominator==1
        if (i+k)%2:assert value==0
    constants=[];groups=0
    for b in range(1,2*B,2):
        q=F(0)
        for h in range(n+1):
            for j in range(1,2*b*h,2):
                group=sum(A[h][i-1]*F(b,j)**(i+k) for i in range(1,a+1))
                assert (delta**qexp*group).denominator==1
                q-=(-1)**h*chi(j)*group
                groups+=1
        assert (delta**qexp*q).denominator==1
        constants.append(q)
    return tuple(constants),rho,A,rows,integer_rows,groups


def moment(a,j):
    assert 0<=j<a
    return F(factorial(2*j)*factorial(2*a-2*j-2),
             2**(2*a-1)*factorial(j)*factorial(a-j-1)*factorial(a-1))


def integral_matrix(B,m,r,a):
    theta,D,L,K,_=parameters(B,r)
    hs,distances=kernels(B,m)
    roots={F(2*r+1)-2*j-t for t in theta for j in range(L)}
    positive=sorted(t for t in roots if t>0)
    assert len(roots)==D*L and 0 in roots and len(positive)==(D*L-1)//2
    assert set(distances)<=set(positive)
    weight=[F(1,prod(t.denominator**(2*m) for t in distances))]
    for t in positive:
        if t not in distances:weight=mul(weight,power([t*t,1],m))
    sigma=(-1)**((m+1)//2+(D*L-1)//2+1)
    C=2**a*K**m
    return [[sigma*C*sum(c*moment(a,j) for j,c in enumerate(mul(weight,mul(hb,hj))))
             for hj in hs] for hb in hs],sigma


def determinant(M):
    A=[list(map(F,row)) for row in M];out=F(1)
    for i in range(len(A)):
        pivot=next(j for j in range(i,len(A)) if A[j][i])
        if pivot!=i:A[i],A[pivot]=A[pivot],A[i];out=-out
        value=A[i][i];out*=value
        for j in range(i+1,len(A)):
            ratio=A[j][i]/value
            for k in range(i+1,len(A)):A[j][k]-=ratio*A[i][k]
    return out


@lru_cache(None)
def seed(B,m,r,a):
    hs,_=kernels(B,m)
    columns=[forms(B,m,r,a,1,h)[0] for h in hs]
    M=[list(row) for row in zip(*columns)]
    direct,sigma=integral_matrix(B,m,r,a)
    assert M==direct
    for j in range(1,B+1):
        assert determinant([[sigma*M[u][v] for v in range(j)] for u in range(j)])>0
    system=[list(map(F,row))+[F(int(i==B-1))] for i,row in enumerate(M)]
    for i in range(B):
        value=system[i][i];assert value
        system[i]=[x/value for x in system[i]]
        for j in range(B):
            if j==i:continue
            value=system[j][i]
            system[j]=[x-value*y for x,y in zip(system[j],system[i])]
    coefficients=[row[-1] for row in system];scale=1
    for c in coefficients:scale=lcm(scale,c.denominator)
    integers=[int(c*scale) for c in coefficients]
    hpoly=[0]
    for c,h in zip(integers,hs):hpoly=add(hpoly,[c*x for x in h])
    q=forms(B,m,r,a,1,tuple(hpoly))[0]
    assert q==tuple([F(0)]*(B-1)+[F(scale)])
    return tuple(hpoly),M,sigma


def raw_value(B,m,r,a,n,hpoly,t):
    theta,D,L,K,_=parameters(B,r)
    C=2**(a*n)*K**(m*n)*factorial(n)**(a-m*D*L)
    P=prod(t-2*r*n+2*j+theta0 for theta0 in theta for j in range(L*n))
    J=m*(D-1)//2
    W=n**(2*J)*evaluate(hpoly,-(t+n)**2/n**2)
    return C*P**m*W/prod(t+2*h for h in range(n+1))**a


def modp(x,p):
    x=F(x);assert x.denominator%p
    return x.numerator%p*pow(x.denominator,-1,p)%p


def valuation(x,p):
    x=F(x)
    if not x:return 10**9
    num,den=abs(x.numerator),x.denominator;value=0
    while num%p==0:num//=p;value+=1
    while den%p==0:den//=p;value-=1
    return value


def main():
    assert __debug__
    counts=dict(kernel_evaluations=0,gram_entries=0,positive_principal_minors=0,
                seed_targets=0,reconstructions=0,reflection_coefficients=0,
                differentiated_grid_zeros=0,grouped_integrality=0,
                selected_pole_congruences=0,other_pole_bounds=0,
                constant_congruences=0,beta_coefficient_bounds=0)
    cases=[(1,1,2),(1,3,2),(1,5,1),(2,1,1),(2,3,1),
           (2,1,2),(3,1,1),(3,3,1),(4,1,1)]
    for B,m,r in cases:
        theta,D,L,K,_=parameters(B,r);a=2*((m*D*L+2)//2)
        hs,distances=kernels(B,m);c=(m+1)//2
        for b,h in zip(range(1,2*B,2),hs):
            for y in [F(1,2),F(1),F(3,2)]:
                direct=(F(0),F(0))
                for j in range(1,2*b,2):
                    term=cpow(cinv((F(j,b)-1,y)),m)
                    direct=direct[0]+chi(j)*term[0],direct[1]+chi(j)*term[1]
                x=y*y
                Q=prod(t.denominator**2*x+t.numerator**2 for t in distances)
                predicted=(F(0),y*evaluate(h,x)/(x**c*Q**m))
                assert direct==predicted
                counts['kernel_evaluations']+=1
        hp,M,sigma=seed(B,m,r,a)
        counts['gram_entries']+=B*B
        counts['positive_principal_minors']+=B
        counts['seed_targets']+=1
        for n in [1,3]:
            q,rho,A,rows,integer_rows,groups=forms(B,m,r,a,n,hp)
            counts['grouped_integrality']+=groups
            for h in range(n+1):
                for i in range(1,a+1):
                    assert rows[n-h][i-1]==(-1)**(i+1)*rows[h][i-1]
                    counts['reflection_coefficients']+=1
            for t in [F(1,2),F(2*r*n+1),F(6*r*n+1,3)]:
                val=sum(rows[h][i-1]/(t+2*h)**i for h in range(n+1) for i in range(1,a+1))
                assert val==raw_value(B,m,r,a,n,hp,t)
                counts['reconstructions']+=1
            grid={F(j,b) for b in range(1,2*B,2) for j in range(-2*b*n+1,2*r*b*n,2)}
            for t in grid:
                assert sum(A[h][i-1]/(t+2*h)**(i+m-1) for h in range(n+1) for i in range(1,a+1))==0
                counts['differentiated_grid_zeros']+=1
        print('Gram and structure passed',B,m,r,a,flush=True)
    prime_cases=[(1,1,2,3),(1,1,2,7),(1,3,2,5),(2,1,1,7),
                 (2,1,1,11),(2,3,1,7),(3,1,1,11),(3,3,1,11),(4,1,1,17)]
    for B,m,r,p in prime_cases:
        theta,D,L,K,_=parameters(B,r);a=2*((m*D*L+2)//2)
        hp,_,_=seed(B,m,r,a)
        E=m*(D-1);tau=a-E+m-1
        q1,_,A1,_,_,g1=forms(B,m,r,a,1,hp)
        q,rho,A,rows,integer_rows,groups=forms(B,m,r,a,p,hp)
        counts['grouped_integrality']+=groups
        for y in [0,1]:
            for i in range(1,a+1):
                nu=i+m-1
                assert modp(F(p)**(tau-nu)*A[p*y][i-1],p)==modp(A1[y][i-1],p)
                counts['selected_pole_congruences']+=1
        for h in range(1,p):
            for coefficient in integer_rows[h]:
                assert coefficient%(p**a)==0
                counts['other_pole_bounds']+=1
        for qv,qseed in zip(q,q1):
            assert modp(p**tau*qv,p)==chi(p)*modp(qseed,p)%p
            counts['constant_congruences']+=1
        for i,value in enumerate(rho,1):
            assert valuation(p**tau*value,p)>=i+m-1
            counts['beta_coefficient_bounds']+=1
        print('Prime dilation passed',B,m,r,p,flush=True)
    record={'status':'All exact self-checks passed; not an independent audit.',
            'cases':cases,'prime_cases':prime_cases,'counts':counts,
            'scope':'Finite identity checks supplement the general proof; they are not a proof by extrapolation.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(counts,indent=2),flush=True)


if __name__=='__main__':main()
