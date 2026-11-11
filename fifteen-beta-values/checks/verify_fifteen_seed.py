#!/usr/bin/env python3
"""Separate seed calculation by logarithmic-derivative Taylor recurrence.

Does not import the truncated-product checker. All poles are computed
individually. The saved product rows are read only after computation.
"""
import hashlib
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent

def calculate(p):
    start=time.monotonic()
    assert all(p%d for d in range(2,math.isqrt(p)+1))
    n=96
    orders=[sum(j<=a<=96-j for j in range(31)) for a in range(n+1)]
    roots=range(-239,240,2)
    rows=[]
    for a,order in enumerate(orders):
        u=n-2*a
        shift=int(u==0)
        length=order-shift
        const=1
        sums=[0]*max(0,length-1)
        factors=[(u-z,3) for z in roots]
        factors.extend((2*(b-a),-orders[b]) for b in range(n+1) if b!=a)
        if u:factors.append((u,1))
        for value,multiplicity in factors:
            assert value%p
            const=const*pow(value,multiplicity,p)%p
            inv=pow(value,-1,p);power=inv
            for j in range(len(sums)):
                sums[j]=(sums[j]+multiplicity*power)%p
                power=power*inv%p
        coeff=[const]
        for k in range(1,length):
            coeff.append(pow(k,-1,p)*sum((-1)**j*sums[j]*coeff[k-1-j] for j in range(k))%p)
        coeff=[0]*shift+coeff
        rows.append([0]+[coeff[order-i] for i in range(1,order+1)])
    beta=[0]*34;harm=[0]*34;constant=0
    for a,row in enumerate(rows):
        for i,c in enumerate(row[1:],1):
            w=(-1)**a*i*(i+1)*c%p
            beta[i+2]=(beta[i+2]+w)%p
            constant=(constant-w*harm[i+2])%p
        inv=pow(2*a+1,-1,p)
        for s in range(3,34):harm[s]=(harm[s]+(-1)**a*pow(inv,s,p))%p
    assert not any(beta[s] for s in range(3,34,2))
    assert constant
    saved=HERE/f'fifteen-seed-rows-p{p}.json'
    assert rows==json.loads(saved.read_text())
    return dict(prime=p,seed_constant_residue=constant,coefficient_comparisons=sum(orders),
                beta_residues={str(s):beta[s] for s in range(4,33,2)},
                rows_sha256=hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
                seconds=time.monotonic()-start)

if __name__=='__main__':
    out={'status':'passed','algorithm':'logarithmic-derivative recurrence',
         'results':[calculate(p) for p in (10007,10009)]}
    (HERE/'fifteen-seed-second-check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
