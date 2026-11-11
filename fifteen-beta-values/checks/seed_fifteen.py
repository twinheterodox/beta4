#!/usr/bin/env python3
"""Exact finite-field nonvanishing certificate for the trimmed seed."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import time
from flint import fmpz_mod_poly_ctx

N=96;L=240;S=31;DS=[Q(j,96) for j in range(S)]

def prefixes(ctx,factors,length):
    out=[ctx.one()]
    for f in factors:out.append(out[-1].mul_low(f,length))
    return out,[f.inverse_series_trunc(length) for f in out]

def product_rows(p):
    ctx=fmpz_mod_poly_ctx(p)
    oddlow=-N-L+1
    roots=list(range(-L+1,L,2))
    odds=range(oddlow,N+L,2)
    num,ni=prefixes(ctx,[ctx([a,1]) for a in odds],S)
    den,di=prefixes(ctx,[ctx([2*h,1]) if h else ctx.one() for h in range(-N,N+1)],S)
    blocks=[(int(d*N),int((1-d)*N)) for d in DS]
    B=[sum(a<=j<=b for a,b in blocks) for j in range(N+1)]
    rows=[]
    for j in range(N+1):
        x=N-2*j;lo=(x-L+1-oddlow)//2;hi=(x+L-1-oddlow)//2+1
        a=num[hi].mul_low(ni[lo],S).pow_trunc(3,S).mul_low(ctx([x,1]),S)
        b=ctx.one()
        for left,right in blocks:
            b=b.mul_low(den[right-j+N+1].mul_low(di[left-j+N],S),S)
        series=a.mul_low(b.inverse_series_trunc(S),S)
        rows.append([0]+[int(series[B[j]-i]) for i in range(1,B[j]+1)])
    return B,rows,roots

def collect(p,B,rows):
    beta=[0]*(S+3);harm=[0]*(S+3);rho=0
    for j in range(N+1):
        for i in range(1,B[j]+1):
            c=(-1)**j*i*(i+1)*rows[j][i]%p;s=i+2
            beta[s]=(beta[s]+c)%p;rho=(rho-c*harm[s])%p
        inv=pow(2*j+1,-1,p);power=1
        for s in range(1,S+3):
            power=power*inv%p;harm[s]=(harm[s]+(-1)**j*power)%p
    assert all(beta[s]==0 for s in range(3,S+3,2))
    return rho,beta

def check_rows(p,B,rows,roots):
    for j in range(N+1):
        for i in range(1,B[j]+1):assert rows[N-j][i]==(-1)**i*rows[j][i]%p
    checks=[]
    for y in (N+L+2,N+L+4,N+L+8):
        direct=y%p
        for z in roots:direct=direct*pow(y-z,3,p)%p
        for j in range(N+1):direct=direct*pow(pow(y-N+2*j,-1,p),B[j],p)%p
        recon=0
        for j in range(N+1):
            inv=pow(y-N+2*j,-1,p);power=1
            for i in range(1,B[j]+1):power=power*inv%p;recon=(recon+rows[j][i]*power)%p
        assert direct==recon;checks.append(dict(point=y,value=direct))
    return checks

def run(p):
    start=time.monotonic()
    assert p>2*N+L+1 and all(p%d for d in range(2,math.isqrt(p)+1))
    B,rows,roots=product_rows(p);rho,beta=collect(p,B,rows)
    checks=check_rows(p,B,rows,roots);assert rho
    canonical=json.dumps(rows,separators=(',',':')).encode()
    out=dict(status='passed',prime=p,base_index=N,trimmed_zero_endpoint=L-1,
             seed_constant_residue=rho,beta_residues={str(s):beta[s] for s in range(4,33,2)},
             rational_reconstruction_checks=checks,coefficient_count=sum(B),
             rows_sha256=hashlib.sha256(canonical).hexdigest(),seconds=time.monotonic()-start)
    here=Path(__file__).resolve().parent
    (here/f'fifteen-seed-p{p}.json').write_text(json.dumps(out,indent=2)+'\n')
    (here/f'fifteen-seed-rows-p{p}.json').write_bytes(canonical+b'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prime',type=int,default=10007)
    args=ap.parse_args();run(args.prime)
