from fractions import Fraction as F
from math import prod
from mpmath import iv
iv.dps=65

def poly(roots):
    a=[1]
    for root in roots:
        b=[0]*(len(a)+1)
        for i,v in enumerate(a):
            b[i]+=root*v
            b[i+1]+=v
        a=b
    return a

A=poly([240]*3+list(range(72,103)))
B=poly([0]*3+[168-j for j in range(31)])
coef=[a-b for a,b in zip(A,B)]
assert coef[-1]==0 and all(c>0 for c in coef[:3]) and all(c<0 for c in coef[3:-1])
def sign(x):
    z=(x-F(3,4))*96
    n,d=z.as_integer_ratio()
    left=(n+240*d)**3*prod(n+(72+j)*d for j in range(31))
    right=n**3*prod(n+(168-j)*d for j in range(31))
    return (left>right)-(left<right)
lo,hi= F(378667493,500000000),F(757334987,1000000000)
assert sign(lo)==1 and sign(hi)==-1
for k in range(150):
    m=(lo+hi)/2
    if sign(m)>0: lo=m
    else: hi=m
def V(z):
    return iv.mpf(z.numerator)/z.denominator
def H(z):
    return z*iv.log(z)
x=iv.mpf([V(lo).a,V(hi).b])
r=V(F(3,4))
phi=3*(H(x+1+r)-H(x-r))-sum(H(x+1-V(F(j,96)))-H(x+V(F(j,96))) for j in range(31))
assert phi < V(F(-18835112315041,10**12))
tail=3*(1+2*F(3,4))*3/(3-F(3,4))-F(341,16)*F(3,4)
printed=V(F(41,5))+15*iv.log(2)-V(F(1883,100))
sharper=V(F(745969417189431202995707,92003136067626388684800))+15*iv.log(2)+phi
print('root',x,'\nphase',phi,'\ntail derivative',tail,'\nprinted rate',printed,'\nunrounded upper rate',sharper,sep=' ')
