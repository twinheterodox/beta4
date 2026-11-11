"""Exact cross-checks of the multiplier study; not an independent audit."""
from fractions import Fraction as Q
from math import comb, gcd
from pathlib import Path
from hashlib import sha256
import json
from explore import data, mul, primitive, rref

HERE=Path(__file__).resolve().parent


def direct_derivative_value(m,r,a,n,w,t):
    """Taylor coefficient at an ordinary point, not at a pole."""
    count=m
    jet=[Q(1)]
    for j in range((2*r+1)*n):
        z=t-2*r*n+2*j+1
        jet=mul(jet,[comb(m,k)*z**(m-k) for k in range(count)],count)
    E=2*(len(w)-1)
    wjet=[sum((w[j]*n**(E-2*j)*comb(2*j,k)*(t+n)**(2*j-k)
               for j in range(len(w)) if 2*j>=k),Q()) for k in range(count)]
    jet=mul(jet,wjet,count)
    for h in range(n+1):
        z=t+2*h
        jet=mul(jet,[(-1)**k*comb(a+k-1,k)/z**(a+k) for k in range(count)],count)
    return jet[m-1]


def main():
    assert __debug__
    records=json.loads((HERE/'explore.json').read_text())['records']
    reconstruction_count=moment_count=0
    for record in records:
        m,r,a,E=record['m'],record['r'],record['a'],record['E']
        stack=[]
        for row in record['indices']:
            n=row['n'];d=data(m,r,a,n,E)
            stack.extend(d['unwanted'])
            w=list(map(Q,row['multiplier_primitive']))
            assert gcd(*(int(v) for v in w))==1
            assert all(sum(x*y for x,y in zip(w,constraint))==0 for constraint in d['unwanted'])
            pair=primitive([sum(x*y for x,y in zip(w,d['constant'])),
                            sum(x*y for x,y in zip(w,d['target']))])
            assert pair==row['primitive_pair']
            for t in (Q(2*r*n)+Q(1,2),Q(2*r*n+1),Q(2*r*n+2)):
                partial=sum((sum(w[j]*d['pole_coefficients'][h][j][i-1]
                                 for j in range(len(w)))/(t+2*h)**(i+m-1)
                             for h in range(n+1) for i in range(1,a+1)),Q())
                assert partial==direct_derivative_value(m,r,a,n,w,t)
                reconstruction_count+=1
            # Independent positive-weight recurrence for the highest coefficient.
            weights=[Q(1)]
            for h in range(n):
                ratio=(Q(r*n+h)+Q(1,2))/(Q((r+1)*n-h)-Q(1,2))
                weights.append(weights[-1]*ratio**m*Q(n-h,h+1)**a)
            assert weights==weights[::-1]
            mass=sum(weights)
            for j in range(d['basis']):
                moment=sum((weights[h]*Q(n-2*h,n)**(2*j) for h in range(n+1)),Q())/mass
                assert moment>0
                assert d['highest'][j]/d['highest'][0]==moment
                moment_count+=1
        assert len(rref(stack)[1])==record['basis_dimension']
        assert record['fixed_multiplier_space_dimension']==0

    harmonic_checks=0
    for n in (1,3,5):
        H=[Q()];H2=[Q()];O=[Q()];O2=[Q()]
        for j in range(1,2*n+1):
            H.append(H[-1]+Q(1,j));H2.append(H2[-1]+Q(1,j*j))
            O.append(O[-1]+Q(1,2*j-1));O2.append(O2[-1]+Q(1,(2*j-1)**2))
        weights=[Q(1)]
        for h in range(n):
            weights.append(weights[-1]*(Q(n+h)+Q(1,2))/(Q(2*n-h)-Q(1,2))*Q(n-h,h+1)**4)
        mass=sum(weights);moments=[Q(),Q()]
        for h,g in enumerate(weights):
            x=n-2*h
            ell1=O[2*n-h]-O[n+h]-2*(H[n-h]-H[h])
            ell2=-O2[2*n-h]-O2[n+h]+H2[n-h]+H2[h]
            A=(ell1*ell1+ell2)/2
            moments[0]+=g*A
            moments[1]+=g*(x*x*A+2*x*ell1+1)
        d=data(1,1,4,n)
        assert d['unwanted'][0][0]/d['highest'][0]==moments[0]/mass
        assert d['unwanted'][0][1]/d['highest'][0]==moments[1]/(n*n*mass)
        harmonic_checks+=2
    rows=[primitive(data(1,1,4,n)['unwanted'][0]) for n in (1,3)]
    assert rows==[[-5,9],[4713857,89181]]
    det=rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]
    assert det==-42870618
    result={'status':'All exact cross-checks passed.',
            'scope':'Constructing-assistant checks; not an independent audit or formal verification.',
            'parameter_triples':len(records),'exact_target_isolations':27,
            'rational_function_reconstructions':reconstruction_count,
            'positive_weight_moment_identities':moment_count,
            'exceptional_case_harmonic_identities':harmonic_checks,
            'exceptional_case_matrix':rows,'exceptional_case_determinant':det,
            'exploration_sha256':sha256((HERE/'explore.json').read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
