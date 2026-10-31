"""Coordinator checks of the construction worker's determinant identity.

Direct rational determinants versus the closed product formula, using
positive and common-zero nodes. This checks an algebraic identity, not
irrationality or the worker's entire argument.
"""
from pathlib import Path
import json,random
from fractions import Fraction as Q
from flint import fmpq,fmpq_mat

HERE=Path(__file__).resolve().parent


def rational_det(matrix):
    return fmpq_mat([[fmpq(x.numerator,x.denominator) for x in row] for row in matrix]).det()


def cauchy_product(a,b):
    num=1;den=1
    for i in range(len(a)):
        for j in range(i+1,len(a)):num*=(a[j]-a[i])*(b[j]-b[i])
        for y in b:den*=a[i]+y
    return Q(num,den)


def run():
    rng=random.Random(20261028);cases=[]
    for n in range(1,10):
        for zero in [False,True]:
            for trial in range(8):
                a=sorted(rng.sample(range(1,50),n-int(zero)))
                b=sorted(rng.sample(range(1,50),n-int(zero)))
                if zero:a=[0]+a;b=[0]+b
                tau=Q(rng.randrange(1,8),rng.randrange(2,13))
                matrix=[]
                for x in a:
                    row=[]
                    for y in b:
                        row.append(1/tau if x==y==0 else
                                   (x*y+tau*(x+y))/((x+y)*(x+tau)*(y+tau)))
                    matrix.append(row)
                det=rational_det(matrix)
                aa,bb=(a[1:],b[1:]) if zero else (a,b)
                expected=cauchy_product(aa,bb)
                for x in aa:expected*=x/(x+tau)
                for y in bb:expected*=y/(y+tau)
                if zero:expected/=tau
                else:expected*=1+tau*sum((Q(1,x) for x in a+b),Q(0))
                assert det==fmpq(expected.numerator,expected.denominator) and det>0
                cases.append(dict(n=n,common_zero=zero,a=a,b=b,tau=str(tau),determinant=str(det)))
    out=dict(status='Coordinator exact checks of the algebraic identity, not an independent mathematical audit.',cases=len(cases),results=cases)
    (HERE/'kernel_identity_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Exact determinant identities checked:',len(cases))


if __name__=='__main__':run()
