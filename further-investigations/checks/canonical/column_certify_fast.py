"""Reuse the column certificate with exact Bernstein--Descartes isolation."""
from pathlib import Path
from math import comb,gcd,lcm
from functools import reduce
import sys,time,json,argparse
from flint import fmpq as Q,fmpq_poly as Poly
HERE=Path(__file__).resolve().parent
# Companion module is in the same replay directory.
import column_certify as inherited

def primitive_list(values):
    divisor=reduce(gcd,values)
    assert divisor
    return [x//divisor for x in values]

def bernstein(P,lo,hi):
    d=P.degree();t=Poly([lo,hi-lo]);C=P(t)
    values=[sum((C[j]*Q(comb(k,j),comb(d,j)) for j in range(k+1)),Q(0)) for k in range(d+1)]
    denominator=lcm(*(int(x.q) for x in values))
    return primitive_list([int(x.p)*(denominator//int(x.q)) for x in values])

def variation(bb):
    signs=[1 if v>0 else -1 for v in bb if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def split(bb):
    d=len(bb)-1;layer=bb;left=[bb[0]<<d];right=[bb[-1]<<d]
    for j in range(1,d+1):
        layer=[a+b for a,b in zip(layer,layer[1:])]
        left.append(layer[0]<<(d-j));right.append(layer[-1]<<(d-j))
    return primitive_list(left),primitive_list(right[::-1])

def isolate(P,domains,bits=80):
    roots=[];counts=[];evaluations=0
    for lo,hi in domains:
        assert P(lo)!=0 and P(hi)!=0
        stack=[(lo,hi,bernstein(P,lo,hi),0)];count=0
        while stack:
            a,b,bb,depth=stack.pop();v=variation(bb);evaluations+=1
            if v==0:continue
            assert depth<80
            if v==1:
                fa=P(a);fb=P(b);assert fa*fb<0
                while b-a>Q(1,2**bits):
                    mid=(a+b)/2;fm=P(mid);evaluations+=1
                    assert fm!=0
                    if fa*fm<0:b=mid;fb=fm
                    else:a=mid;fa=fm
                roots.append((a,b));count+=1
            else:
                mid=(a+b)/2;assert P(mid)!=0
                left,right=split(bb)
                stack.extend([(a,mid,left,depth+1),(mid,b,right,depth+1)])
        counts.append(count)
    return sorted(roots),counts,evaluations

def main(degree):
    inherited.HERE=HERE;inherited.isolate=isolate
    inherited.run(degree)
    path=HERE/f'column_certificate_d{degree}.json'
    data=json.loads(path.read_text())
    data['root_isolation_method']='Exact Bernstein coefficients and Descartes variations; sign bisection after a unique root is proved.'
    data['root_isolation_evaluation_points']=data.pop('sturm_evaluation_points')
    data['scope']='Exact rational dual functions; complete Bernstein--Descartes root isolation; Arb evaluation. The inherited uniform energy theorem is unchanged.'
    path.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degree',type=int,default=200);a=p.parse_args();main(a.degree)
