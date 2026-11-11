"""Exact coefficient congruences for the lowest homogeneous part.

This checks the newly derived general-seed extension in a small single-grid
case. The proof for arbitrary fixed seeds is in proof.txt, not inferred
from these samples.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
from hashlib import sha256
import json
from explore import data

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def base(n,E):return data(1,1,4,n,E)


def coefficients(n,poly):
    E=max(x for u,x in poly)
    entry=base(n,E)
    w=[Q() for j in range(E//2+1)]
    for (u,x),c in poly.items():
        w[x//2]+=Q(c)*Q(n)**(u-E+x)
    norm=2**(7*n)*factorial(n)
    return [norm*sum(x*y for x,y in zip(w,row))
            for row in (entry['unwanted'][0],entry['target'])]


def reduce_mod_p(value,p):
    value=Q(value)
    assert value.denominator%p
    return value.numerator*pow(value.denominator,-1,p)%p


def main():
    assert __debug__
    polynomials=[
        {(0,0):1,(1,0):1,(0,2):1},
        {(2,0):1,(0,2):2,(3,0):1,(0,4):1},
        {(1,0):1,(2,2):3},
        {(3,0):1,(1,2):1,(0,4):2},
    ]
    rows=[]
    for index,poly in enumerate(polynomials):
        e0=min(u+x for u,x in poly)
        low={key:c for key,c in poly.items() if sum(key)==e0}
        for seed in (1,3):
            wanted=coefficients(seed,low)
            for p in (7,11,17):
                actual=coefficients(seed*p,poly)
                for nu,got,want in zip((2,4),actual,wanted):
                    left=reduce_mod_p(Q(p)**(4-e0-nu)*got,p)
                    right=reduce_mod_p(want,p)
                    assert left==right
                    rows.append({'polynomial':index,'minimum_degree':e0,'seed':seed,
                                 'prime':p,'beta_index':nu,'residue':left})
    result={'status':'All exact general-seed coefficient congruences passed.',
            'scope':'48 finite checks supporting the all-parameter proof; not an independent audit.',
            'count':len(rows),'checks':rows,
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'dilation-checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(rows))


if __name__=='__main__':main()
