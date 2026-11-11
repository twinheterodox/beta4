"""Exact finite checks of the displayed determinant and normalization identities."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial, gcd, prod
from pathlib import Path
from collections import Counter
import json
from flint import fmpq, fmpq_mat

counts=Counter()
def det(a):
    return Q(str(fmpq_mat([[fmpq(str(x)) for x in row] for row in a]).det()))
def bezout(a,b):
    if b==0:
        return abs(a),1 if a>=0 else -1,0
    g,x,y=bezout(b,a%b)
    return g,y,x-(a//b)*y
def bezout_all(values):
    g=0; weights=[]
    for v in values:
        h,x,y=bezout(g,v)
        weights=[x*w for w in weights]+[y]
        g=h
    assert sum(x*y for x,y in zip(values,weights))==g
    return g,weights

for b in range(1,5):
    for ell in range(1,9):
        for seed in (1,2,3):
            for zeros in range(b+1):
                d=[seed*2**i if i<b-zeros else 0 for i in range(b)]
                base=[[0]*b for _ in range(b)]
                for i in range(b):
                    base[i][i]=d[i]
                    if i<b-1:base[i][i+1]=ell
                cof=[(-1)**(b+j+1)*ell**(b-j-1)*prod(d[:j]) for j in range(b)]
                g,w=bezout_all(cof)
                assert g==prod(gcd(ell,x) for x in d[:-1])
                for k in (-3,-1,0,1,4):
                    mat=[row[:] for row in base]
                    for j in range(b):mat[-1][j]+=ell*k*w[j]
                    assert det(mat)==prod(d)+ell*g*k
                    counts['Smith_progression_attainments']+=1

for e in range(1,7):
    for x in map(Q,[-3,-1,0,1,2]):
        D=[[Q(int(i==j))-x*int(i==j+1) for j in range(e)] for i in range(e+1)]
        E=[[x**(i-j) if j<=i else Q(0) for j in range(e+1)] for i in range(e)]
        assert all(sum(E[i][r]*D[r][j] for r in range(e+1))==int(i==j) for i in range(e) for j in range(e))
        null=[x**(e-j) for j in range(e+1)]
        assert all(sum(null[r]*D[r][j] for r in range(e+1))==0 for j in range(e))
        counts['chain_left_inverse_and_null_row']+=1

for h in range(2,6):
    for k in range(1,min(h,3)+1):
        for jj in combinations(range(h),k):
            for rr in combinations(range(6),k):
                lhs=det([[prod(range(r+1,r+j+1)) for j in jj] for r in rr])
                rhs=det([[(r+h)**j for j in jj] for r in rr])
                assert 0<=lhs<=rhs
                counts['factorial_alternant_comparisons']+=1

for b in (1,2,3):
    d=[Q((i+1)**2,4) for i in range(b)]
    H=sum(factorial(k)*sum((prod(d[i] for i in I) for I in combinations(range(b),k)),Q(0)) for k in range(b+1))
    values=[]
    for signs in product((-1,1),repeat=b*b):
        mat=[[Q(int(i==j))+Q(i+1,2)*signs[i*b+j] for j in range(b)] for i in range(b)]
        values.append(det(mat))
    assert sum(values,Q(0))/len(values)==1
    assert sum((v*v for v in values),Q(0))/len(values)==H
    counts['rounding_determinant_outcomes']+=len(values)

segments=[(0,11,96,0),(11,Q(65,2),86,0),(Q(65,2),Q(145,4),82,0),
          (Q(145,4),52,227,-4),(52,58,182,-3),(58,59,124,-2),(59,65,65,-1)]
area=sum((Q(a)*(Q(hi)-Q(lo))+Q(b,2)*(Q(hi)**2-Q(lo)**2) for lo,hi,a,b in segments),Q(0))
assert area==Q(33079,8)
assert area/48**2==Q(33079,18432)
counts['prime_cost_integral']=1
out=dict(status='passed',counts=dict(counts),scope='Finite algebraic checks of the displayed mechanisms; general proofs are in the manuscript.')
Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
