"""Explicit shifted hypergeometric certificate problems, from the archived construction."""
from math import comb
import algebra as c

def factstep(a,alpha=0):
    out=c.ONE
    for j in range(1,a+1):out*=a*c.N+alpha+j
    return out

def recurrence(raw,alpha,forward_kind=False):
    g=[c.parse(raw[i][2]) for i in range(3)]
    if forward_kind:
        rho=c.Rat(factstep(15)*factstep(12,alpha)*factstep(13,alpha)*factstep(14,alpha),
             factstep(10)**3*factstep(1)*factstep(14)*factstep(2,alpha)*factstep(3,alpha)*factstep(4,alpha))
        for j in range(3):
            for h in range(j,3):g[j]*=rho.shift(h)
    B=c.primitive([c.Rat(1),-g[2],-g[1],-g[0]])
    for j in range(3):assert not (c.Rat(B[3-j],B[0])+g[j]).p
    assert c.nonzero(B[0])
    return B

def fproblem(B,alpha):
    A=c.linear_product([(12,1,1+alpha),(13,1,1+alpha),(14,1,1+alpha),(14,-1,42),(15,-1,45)])
    Bm=c.linear_product([(0,1,0),(1,1,3),(2,1,6+alpha),(3,1,9+alpha),(4,1,12+alpha)])
    D=c.linear_product([(a,1,r+(alpha if a>=2 else 0)) for a in range(1,5) for r in range(1,3*a+1)]
                     +[(14,-1,r) for r in range(1,43)]+[(15,-1,r) for r in range(1,46)])
    K=c.ONE
    for r in range(1,31):K*=(10*c.N+r)**3
    F=[]
    for j in range(4):
        normalization=c.ONE
        for r in range(1,15*j+1):normalization*=(15*c.N+r)**2
        for r in range(10*j+1,31):normalization*=(10*c.N+r)**3
        factors=([(a,1,r+alpha) for a in range(12,15) for r in range(1,a*j+1)]
                +[(a,1,r+(alpha if a>=2 else 0)) for a in range(1,5) for r in range(a*j+1,3*a+1)]
                +[(14,-1,r) for r in range(14*j+1,43)]+[(15,-1,r) for r in range(15*j+1,46)])
        F=c.add(F,c.scale(c.linear_product(factors),B[3-j]*normalization))
    return A,Bm,F,D,K

def dproblem(h,alpha):
    D=c.linear_product([(a,1,r+alpha) for a in range(1,3) for r in range(1,h*a+1)]
                     +[(b,-1,r) for b in range(11,13) for r in range(1,h*b+1)]
                     +[(10,-1,1+r) for r in range(1,10*h+1)])
    A=c.scale(c.linear_product([(16,1,1+alpha),(17,1,1+alpha),(11,-1,11*h),
                              (12,-1,12*h),(10,-1,1+10*h)]),-1)
    Bm=c.linear_product([(0,1,0),(-1,1,0),(-2,1,0),(1,1,h+alpha),(2,1,2*h+alpha)])
    return A,Bm,D

def dshift(j,h,alpha):
    return c.linear_product([(a,1,r+alpha) for a in range(16,18) for r in range(1,a*j+1)]
           +[(a,1,r+alpha) for a in range(1,3) for r in range(a*j+1,h*a+1)]
           +[(b,-1,r) for b in range(11,13) for r in range(b*j+1,h*b+1)]
           +[(10,-1,1+r) for r in range(10*j+1,10*h+1)]
           +[(-a,1,-r) for a in range(1,3) for r in range(a*j)])

def power(j,alpha):
    return [c.fmpq_poly([comb(j,i)*alpha**(j-i)]) for i in range(j+1)]

def coeffs(P):
    return [[str(x) for x in p.coeffs()] for p in P]
