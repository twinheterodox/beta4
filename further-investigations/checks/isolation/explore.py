"""Exact exploration of fixed and index-dependent polynomial multipliers.

No irrationality claim is made. Fractions and exact row reduction are used
for ranks, coefficient cancellation and primitive rational coefficient pairs.
"""
from fractions import Fraction as Q
from math import comb, factorial, gcd, lcm
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def mul(p, q, count):
    out = [Q(0)]*count
    for i, x in enumerate(p[:count]):
        for j, y in enumerate(q[:count-i]):
            out[i+j] += x*y
    return out


def rref(mat):
    mat = [list(map(Q, row)) for row in mat]
    if not mat: return [], []
    pivots = []
    row = 0
    for col in range(len(mat[0])):
        pivot = next((i for i in range(row, len(mat)) if mat[i][col]), None)
        if pivot is None: continue
        mat[row], mat[pivot] = mat[pivot], mat[row]
        div = mat[row][col]
        mat[row] = [x/div for x in mat[row]]
        for i in range(len(mat)):
            if i != row and mat[i][col]:
                c = mat[i][col]
                mat[i] = [x-c*y for x, y in zip(mat[i], mat[row])]
        pivots.append(col)
        row += 1
        if row == len(mat): break
    return mat, pivots


def nullspace(mat):
    reduced, pivots = rref(mat)
    dimension = len(mat[0])
    free = [j for j in range(dimension) if j not in pivots]
    out = []
    for j in free:
        v = [Q(0)]*dimension
        v[j] = 1
        for i, pivot in enumerate(pivots):
            v[pivot] = -reduced[i][j]
        out.append(v)
    return out


def primitive(values):
    den = lcm(*(Q(x).denominator for x in values))
    integers = [int(x*den) for x in values]
    div = gcd(*integers)
    if not div: return integers
    integers = [x//div for x in integers]
    if next(x for x in reversed(integers) if x) < 0:
        integers = [-x for x in integers]
    return integers


def data(m, r, a, n, E=None):
    assert n % 2 == m % 2 == 1 and a % 2 == 0
    L = 2*r+1
    assert a > m*L
    if E is None: E = a-2
    assert E % 2 == 0
    basis = E//2+1
    rho = [[Q(0)]*basis for _ in range(a)]
    q = [Q(0)]*basis
    pole_coefficients=[]
    for h in range(n+1):
        local = [Q(1)]
        for j in range(L*n):
            z = -2*h-2*r*n+2*j+1
            local = mul(local, [Q(comb(m, k)*z**(m-k)) for k in range(m+1)], a)
        for j in range(n+1):
            if j == h: continue
            d = 2*(j-h)
            local = mul(local, [Q((-1)**k*comb(a+k-1,k), d**(a+k)) for k in range(a)], a)
        at_pole=[]
        for j in range(basis):
            w = [Q(comb(2*j,k)*n**(E-2*j)*(n-2*h)**(2*j-k))
                 for k in range(min(2*j,a-1)+1)]
            product = mul(local,w,a)
            coefficients=[]
            for i in range(1,a+1):
                nu = i+m-1
                c = comb(i+m-2,m-1)*product[a-i]
                coefficients.append(c)
                rho[i-1][j] += (-1)**h*c
                tail = sum((Q((-1)**((v-1)//2),v**nu) for v in range(1,2*h,2)), Q(0))
                q[j] -= (-1)**h*c*tail
            at_pole.append(coefficients)
        pole_coefficients.append(at_pole)
    assert all(not x for i in range(1,a+1,2) for x in rho[i-1])
    unwanted = [rho[i-1] for i in range(2,a+1,2) if i+m-1 != 4]
    target = rho[4-m]
    return {'unwanted':unwanted, 'target':target, 'constant':q,
            'highest':rho[-1], 'basis':basis,'pole_coefficients':pole_coefficients}


def main():
    records = []
    for m,r,a in [(1,1,4),(1,1,6),(1,2,6),(1,2,10),(1,2,18),
                  (1,3,8),(1,3,20),(3,1,10),(3,2,16)]:
        stack=[];per_index=[]
        for n in (1,3,5):
            entry=data(m,r,a,n)
            stack.extend(entry['unwanted'])
            ns=nullspace(entry['unwanted'])
            assert len(ns)==1
            w=ns[0]
            assert all(sum(x*y for x,y in zip(row,w))==0 for row in entry['unwanted'])
            q=sum(x*y for x,y in zip(entry['constant'],w))
            target=sum(x*y for x,y in zip(entry['target'],w))
            assert target
            pair=primitive([q,target])
            per_index.append({'n':n,'rank':len(rref(entry['unwanted'])[1]),
                              'cumulative_rank':len(rref(stack)[1]),
                              'multiplier_primitive':primitive(w),
                              'primitive_pair':pair,
                              'target_only':True})
        records.append({'m':m,'r':r,'a':a,'E':a-2,'basis_dimension':entry['basis'],
                        'indices':per_index,'fixed_multiplier_space_dimension':entry['basis']-len(rref(stack)[1])})
        print(m,r,a,'joint dimension',records[-1]['fixed_multiplier_space_dimension'],
              'pair digit lengths',[(x['n'],len(str(abs(x['primitive_pair'][1])))) for x in per_index],flush=True)
    result={'scope':'Finite exact exploration, not an all-index height or irrationality theorem.',
            'normalization':'Base rational function omits its scalar kappa_n; final pairs have full rational denominators and common integer gcd removed.',
            'records':records}
    (HERE/'explore.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
