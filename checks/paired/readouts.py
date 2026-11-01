"""Rational reconstruction of all forward and dual coefficients; no numerical special values."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial, prod
from collections import Counter
import json
import sys
sys.set_int_max_str_digits(0)

HERE=Path(__file__).resolve().parent

def choose(x,m):
    return prod((F(x-i) for i in range(m)),start=F(1))/factorial(m)

def h(m,r,a=F(1)):
    return sum((1/(a+j)**r for j in range(m)),F(0))

def bs(n,k,alpha):
    upper=prod((choose(k+j*n+alpha,10*n) for j in (12,13,14)),start=F(1))
    def cb(a,b): return comb(a,b) if 0<=b<=a else 0
    b=upper*cb(15*n,n+k)*cb(15*n,k)
    if 0<=k<=14*n:
        logarithm=sum((sum((1/(k+j*n+alpha-i) for i in range(10*n)),F(0)) for j in (12,13,14)),F(0))
        logarithm+=h(14*n-k,1)-h(n+k,1)+h(15*n-k,1)-h(k,1)
        derivative=b*logarithm
    elif 14*n<k<=15*n:
        derivative=upper*cb(15*n,k)*F((-1)**k,(n+k)*comb(n+k-1,15*n))
    elif -n<=k<0:
        i=-k
        derivative=upper*cb(15*n,n-i)*F((-1)**(i-1),(15*n+i)*comb(15*n+i-1,i-1))
    else:
        derivative=F(0)
    return b,derivative

def forward(n,alpha):
    a=1+alpha
    q=sum((bs(n,k,alpha)[0] for k in range(14*n+1)),F(0))
    p=[]
    for s in (3,4):
        p.append(sum((b*h(4*n+k,s,a)-derivative*h(4*n+k,s-1,a)/(s-1)
                      for k in range(-n,15*n+1) for b,derivative in [bs(n,k,alpha)]),F(0)))
    return [q,*p]

def dual(n,alpha):
    roots=[F(1+2*n+i) for i in range(14*n)]+[F(1+n+i) for i in range(16*n)]
    poles=Counter([-alpha-i for i in range(11*n+1)]+[-alpha-n-i for i in range(11*n+1)]+[-alpha-2*n-i for i in range(8*n+2)])
    rows=[[F(0) for k in range(3)] for r in range(3)]
    for pole,m in poles.items():
        j=int(-pole-alpha)
        c=F((-1)**n);l1=F(0);l2=F(0)
        for r in roots:
            d=r-pole;c*=d;l1-=1/d;l2-=1/(2*d*d)
        for other,e in poles.items():
            if other==pole:continue
            d=pole-other;c/=d**e;l1-=e/d;l2+=e/(2*d*d)
        la={m:c,m-1:c*l1,m-2:c*(l1*l1/2+l2)}
        for k in range(3):
            cf={i:sum((comb(k,l)*(-pole)**(k-l)*(-1)**l*la.get(i+l,0) for l in range(k+1)),F(0)) for i in range(1,4)}
            rows[0][k]+=3*cf[3];rows[1][k]+=2*cf[2]
            rows[2][k]-=sum((i*cf[i]*h(2*n+j,i+1,1+alpha) for i in range(1,4)),F(0))
    return rows

def determinant(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
