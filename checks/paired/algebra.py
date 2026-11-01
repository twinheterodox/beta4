"""Exact univariate rational functions and polynomial telescoping.

The descending solve is the standard coefficient method for
A(k)P(k+1)-B(k)P(k)=F(k). Every result is checked by full substitution.
This compact implementation is specific to the supplement certificate data.
"""
import ast
from math import comb, lcm
from flint import fmpq, fmpq_poly
N=fmpq_poly([0,1]); ZERO=fmpq_poly([]); ONE=fmpq_poly([1])
class Rat:
    def __init__(self,p=0,q=1):
        if isinstance(p,Rat): self.p,self.q=p.p,p.q; return
        p=p if isinstance(p,fmpq_poly) else fmpq_poly([p])
        q=q if isinstance(q,fmpq_poly) else fmpq_poly([q])
        if not q: raise ZeroDivisionError()
        g=p.gcd(q);p,q=p//g,q//g
        self.p,self.q=p/q.leading_coefficient(),q/q.leading_coefficient()
    def __add__(self,b):
        b=Rat(b)
        if self.q==b.q:return Rat(self.p+b.p,self.q)
        return Rat(self.p*b.q+b.p*self.q,self.q*b.q)
    __radd__=__add__
    def __neg__(self):return Rat(-self.p,self.q)
    def __sub__(self,b):return self+-Rat(b)
    def __rsub__(self,b):return Rat(b)+-self
    def __mul__(self,b):
        b=Rat(b);return Rat(self.p*b.p,self.q*b.q)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=Rat(b);return Rat(self.p*b.q,self.q*b.p)
    def __rtruediv__(self,b):return Rat(b)/self
    def __pow__(self,e):return Rat(self.p**e,self.q**e) if e>=0 else Rat(self.q**(-e),self.p**(-e))
    def __eq__(self,b):
        b=Rat(b);return self.p==b.p and self.q==b.q
    def __bool__(self):return bool(self.p)
    def shift(self,h):return Rat(self.p(N+h),self.q(N+h))
    def at(self,n):return self.p(n)/self.q(n)
    def pack(self):return [[str(x) for x in z.coeffs()] for z in (self.p,self.q)]

def parse(s,variables=None):
    env={'n':Rat(N)};env.update(variables or {})
    def go(t):
        if isinstance(t,ast.Constant) and type(t.value) is int:return Rat(t.value)
        if isinstance(t,ast.Name) and t.id in env:return env[t.id]
        if isinstance(t,ast.UnaryOp):
            if isinstance(t.op,ast.USub):return -go(t.operand)
            if isinstance(t.op,ast.UAdd):return go(t.operand)
        if isinstance(t,ast.BinOp):
            a=go(t.left);b=go(t.right)
            if isinstance(t.op,ast.Add):return a+b
            if isinstance(t.op,ast.Sub):return a-b
            if isinstance(t.op,ast.Mult):return a*b
            if isinstance(t.op,ast.Div):return a/b
            if isinstance(t.op,ast.Pow):
                assert b.q==ONE and b.p.degree()<=0 and b.p[0].denom()==1
                return a**int(b.p[0])
        raise ValueError(ast.dump(t))
    return go(ast.parse(s,mode='eval').body)

def unpack(x):return Rat(*(fmpq_poly([fmpq(z) for z in row]) for row in x))
def transpose(a):return [list(x) for x in zip(*a)]
def mm(a,b):return [[sum((x*y for x,y in zip(r,c)),Rat()) for c in zip(*b)] for r in a]
def shifted(a,h):return [[x.shift(h) for x in r] for r in a]
def det(a):
    a=[r[:] for r in a];out=Rat(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return Rat()
        if j!=i:a[j],a[i]=a[i],a[j];out=-out
        v=a[i][i];out*=v
        for j in range(i+1,len(a)):
            t=a[j][i]/v
            for k in range(i+1,len(a)):a[j][k]-=t*a[i][k]
    return out

def nonzero(poly,start=0):
    if not poly:return False
    if poly.leading_coefficient()<0:poly=-poly
    return bool(poly(start)) and bool(poly(start+1)) and all(x>=0 for x in poly(N+start+1).coeffs())
def poly_lcm(a,b):return a*b//a.gcd(b)
def primitive(values):
    d=ONE
    for v in values:d=poly_lcm(d,v.q)
    p=[v.p*(d//v.q) for v in values];g=ZERO
    for v in p:g=g.gcd(v)
    return [v//g for v in p]
def trim(a):
    while a and not a[-1]:a.pop()
    return a

def add(a,b):
    z=[ZERO]*max(len(a),len(b))
    for j,v in enumerate(a):z[j]+=v
    for j,v in enumerate(b):z[j]+=v
    return trim(z)
def scale(a,c):return trim([v*c for v in a])
def multiply(a,b):
    if not a or not b:return []
    z=[ZERO]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return trim(z)
def linear_product(factors):
    z=[ONE]
    for a,b,c in factors:z=multiply(z,[a*N+c,fmpq_poly([b])])
    return z

def shift_k(p):
    out=[ZERO]*len(p)
    for i,v in enumerate(p):
        for j in range(i+1):out[j]+=comb(i,j)*v
    return trim(out)
def operator(a,b,p):return add(multiply(a,shift_k(p)),scale(multiply(b,p),-1))
def evaluate_k(p,k):
    out=ZERO
    for v in reversed(p):out=out*k+v
    return out

def solve(a,b,f,d):
    assert len(a)==len(b)==6 and a[5]==b[5]==ONE
    assert a[4]-b[4]==fmpq_poly([-d]) and len(f)<=d+4
    residual=f[:]+[ZERO]*(d+5-len(f))
    h=[ZERO]*d+[ONE];part=[ZERO]*(d+1)
    hr=operator(a,b,h)+[ZERO]*(d+5-len(operator(a,b,h)))
    for i in range(d-1,-1,-1):
        part[i]=residual[i+4]/(i-d);h[i]=-hr[i+4]/(i-d)
        image=operator(a,b,[ZERO]*i+[ONE])
        for j,v in enumerate(image):
            residual[j]-=part[i]*v;hr[j]+=h[i]*v
        assert not residual[i+4] and not hr[i+4]
    idx=next(j for j in range(4) if hr[j])
    top,rem=divmod(residual[idx],hr[idx]);assert not rem
    assert all(residual[j]==top*hr[j] for j in range(4))
    p=trim([part[j]+top*h[j] for j in range(d+1)])
    assert operator(a,b,p)==trim(f[:])
    return p
