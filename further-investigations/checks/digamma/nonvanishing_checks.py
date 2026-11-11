"""Finite exact checks of the common-node deletion lemma.

Each determinant polynomial is reconstructed from direct determinants,
then checked at two additional rational points. A gcd of 1 verifies that
the main determinant and its common-node cofactors have no common root.
The all-size proof is in the supplement; finite checks do not prove that lemma.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from flint import fmpq, fmpq_mat, fmpq_poly
from loewner import family,nodes,rational_matrix

def determinant(terms,R,C,x):
    if not R:return fmpq(1)
    data=rational_matrix(terms,R,C,Q(x))
    return fmpq_mat([[fmpq(q.numerator,q.denominator) for q in row] for row in data]).det()

def polynomial(terms,R,C):
    d=len(set(R)&set(C))
    values=[determinant(terms,R,C,x) for x in range(d+1)]
    coeff=[]
    while values:
        coeff.append(values[0])
        values=[values[j+1]-values[j] for j in range(len(values)-1)]
    P=fmpq_poly();basis=fmpq_poly([1])
    for k,a in enumerate(coeff):
        if k:basis=basis*fmpq_poly([1-k,1])/k
        P+=a*basis
    for x in (d+1,d+2):assert P(x)==determinant(terms,R,C,x)
    return P

if __name__=='__main__':
    records=[];polynomials=0
    slopes=['-1','0','1','1/7','5/48','1/16','-1/16']
    for slope in slopes:
        terms=family('three',Q(slope))
        for r in (1,2,4,6):
            for m in (0,1,2):
                for prefix in range(m+1):
                    R,C=nodes(r,m,prefix)
                    P=polynomial(terms,R,C);g=P
                    polynomials+=1
                    for t in sorted(set(R)&set(C)):
                        minor=polynomial(terms,[a for a in R if a!=t],[a for a in C if a!=t])
                        g=g.gcd(minor)
                        polynomials+=1
                    assert g.degree()==0,(slope,r,m,prefix,g)
                    records.append(dict(slope=slope,r=r,m=m,prefix=prefix,gcd_degree=0))
    result=dict(status='Finite exact checks; all-size proof is separate.',
        cases=len(records),determinant_polynomials=polynomials,
        extra_evaluation_checks=2*polynomials,records=records)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'}))
