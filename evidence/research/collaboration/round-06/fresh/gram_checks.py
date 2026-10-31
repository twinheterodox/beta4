"""Exact checks for difference-variable Gram determinants. Writes only this folder."""
from fractions import Fraction as F
from math import comb, gcd, lcm, factorial
from pathlib import Path
import json
from flint import fmpq, fmpq_mat, arb, ctx

HERE=Path(__file__).resolve().parent

def prime_list(M):
    sieve=[True]*(M+1)
    if M>=0:sieve[0]=False
    if M>=1:sieve[1]=False
    for p in range(2,int(M**.5)+1):
        if sieve[p]:
            for k in range(p*p,M+1,p):sieve[k]=False
    return [p for p in range(2,M+1) if sieve[p]]

def val_int(a,p):
    a=abs(int(a))
    if a==0:return 10**9
    v=0
    while a%p==0:a//=p;v+=1
    return v

def valuation(a,p):
    a=F(a)
    return val_int(a.numerator,p)-val_int(a.denominator,p)

def modq(x,p):
    x=F(x)
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p

def tails(s,M):
    out=[F(0)]
    for k in range(M):out.append(out[-1]+F((-1)**k,(2*k+1)**s))
    return out

def moments(M):
    t3=tails(3,2*M);t4=tails(4,M)
    A=[];B=[]
    for m in range(M+1):
        c=comb(2*m,m)
        b=-c*t4[m]+sum((F(comb(2*m,j),6*(j-m))*t3[j] for j in range(2*m+1) if j!=m),F(0))
        A.append(c);B.append(b)
    return A,B

def mono_pair(i,j,t3,t4):
    if i==j:return F((-1)**i),-F((-1)**i)*t4[i]
    assert (i+j)%2==0
    return F(0),-F((-1)**i,6*(j-i))*(t3[i]-t3[j])

def moments_direct(M):
    t3=tails(3,2*M);t4=tails(4,2*M)
    out=[]
    for m in range(M+1):
        a=b=F(0)
        for i in range(2*m+1):
            j=2*m-i;c=(-1)**j*comb(2*m,i)
            aa,bb=mono_pair(i,j,t3,t4)
            a+=c*aa;b+=c*bb
        out.append((a,b))
    return out

def q(x):
    x=F(x);return fmpq(x.numerator,x.denominator)

def primitive(coeff):
    coeff=[F(str(c)) for c in coeff]
    d=lcm(*(c.denominator for c in coeff))
    ints=[int(c*d) for c in coeff]
    g=gcd(*ints)
    if ints[-1]<0:g=-g
    return [x//g for x in ints],F(d,g)

def determinant_fraction(mat):
    a=[[F(x) for x in row] for row in mat];out=F(1);N=len(a)
    for i in range(N):
        j=next((j for j in range(i,N) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];out=-out
        d=a[i][i];out*=d
        for j in range(i+1,N):
            c=a[j][i]/d
            for k in range(i+1,N):a[j][k]-=c*a[i][k]
    return out

def rank_mod(mat,p):
    a=[[int(x)%p for x in row] for row in mat];r=0
    for j in range(len(a[0])):
        z=next((i for i in range(r,len(a)) if a[i][j]),None)
        if z is None:continue
        a[r],a[z]=a[z],a[r];w=pow(a[r][j],-1,p)
        a[r]=[(v*w)%p for v in a[r]]
        for i in range(len(a)):
            if i!=r:
                c=a[i][j]
                a[i]=[(x-c*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def theta():return (arb(4).zeta(arb(1)/4)-arb(4).zeta(arb(3)/4))/256

def horner(coeff,x):
    out=arb(0)
    for c in reversed(coeff):out=out*x+int(c)
    return out

def main():
    max_n=12;max_m=2*max_n-2
    A,B=moments(max_m)
    direct=moments_direct(max_m)
    assert direct==list(zip(A,B))
    prime_checks=[];rows=[];residue_checks=0;rank_checks=0
    for n in range(1,max_n+1):
        AA=fmpq_mat([[A[i+j] for j in range(n)] for i in range(n)])
        BB=fmpq_mat([[q(B[i+j]) for j in range(n)] for i in range(n)])
        assert AA.det()==2**(n-1)
        p=(-AA.inv()*BB).charpoly()*AA.det()
        coeff=list(p);ints,kappa=primitive(coeff)
        assert len(ints)==n+1 and ints[-1]>0
        assert gcd(*ints)==1
        for x in [0,1,2]:
            check=determinant_fraction([[A[i+j]*x+B[i+j] for j in range(n)] for i in range(n)])
            assert check==F(str(p(x)))
        # Actual primitive exponents, and independently reconstructed residue matrices.
        p_records=[]
        for pp in prime_list(8*n):
            if pp<=4*n:continue
            m0=pp//4+1
            r=max(0,2*n-1-m0)
            residue=[]
            for m in range(2*n-1):
                v=modq(pp**3*B[m],pp)
                h=(pp+1)//2
                expected=modq(F((-1)**((pp-1)//2),6)*sum((F(comb(2*m,j),j-m) for j in range(h,2*m+1)),F(0)),pp)
                assert v==expected
                if m<m0:assert v==0
                residue.append(v);residue_checks+=1
            if m0<=2*n-2:
                mm=m0
                first=F(1,mm) if pp==4*mm-1 else F((2*mm-1)*(mm+1),mm*(mm-1))
                expected=modq(F((-1)**((pp-1)//2),6)*first,pp)
                assert residue[m0]==expected and expected!=0
            mat=[[residue[i+j] for j in range(n)] for i in range(n)]
            assert rank_mod(mat,pp)==r;rank_checks+=1
            actual=valuation(kappa,pp)
            assert actual==3*r,(n,pp,r,actual)
            p_records.append(dict(p=pp,m0=m0,rank=r,actual_primitive_exponent=actual))
        ctx.prec=max(max(abs(z).bit_length() for z in ints)+512,1024)
        tt=theta();vv=horner(ints,tt)
        assert vv>0
        raw=vv/(arb(kappa.numerator)/kappa.denominator)
        # Elementary rigorous determinant bounds; theta<1 in upper bound.
        lower=arb(1)/(arb(6144)**n*factorial(n)*arb(2)**(2*n*n+n))
        upper=arb(4)**(-(n-1)**2)
        assert raw>lower and raw<upper if n>1 else raw>lower and raw<=upper
        # Store a short logarithmic enclosure as well as exact coefficients.
        rows.append(dict(n=n,primitive_coefficients=ints,kappa=str(kappa),primitive_log=str(vv.log()),raw_log=str(raw.log()),prime_checks=p_records))
        print('n',n,'degree',len(ints)-1,'log primitive',float(vv.log()),'large-prime certificates',len(p_records),flush=True)
    result=dict(moment_formula_checks=len(direct),determinant_evaluation_checks=3*max_n,residue_checks=residue_checks,rank_checks=rank_checks,rows=rows)
    (HERE/'gram_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print('All exact and Arb checks passed.',result['moment_formula_checks'],residue_checks,rank_checks,flush=True)
if __name__=='__main__':main()
