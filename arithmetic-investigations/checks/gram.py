#!/usr/bin/env python3
"""Rebuild the recorded Gram cases and their complete primitive polynomials.
Portable adaptation of the original exact assessor; no new indices.
"""
import argparse, json, math
from pathlib import Path
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat, fmpq_mpoly_ctx
ROOT=Path(__file__).resolve().parents[1]

def blocks(h, r):
    labels = [(a, j) for a in range(h) for j in range(r - 1)]
    def moment(a, b, d):
        s = 4 + 2 * a + 2 * b
        return fmpq(math.factorial(s - 1), d ** s)
    C = fmpq_mat([[moment(a, b, 2*j+1) for b,j in labels] for a in range(h)])
    B = fmpq_mat([[moment(a,b,2*(i+j)+1) + moment(a,b,2*(i+j)+3)
                   for b,j in labels] for a,i in labels])
    if labels:
        Y = B.solve(C.transpose())
        assert B * Y == C.transpose()
        R = C * Y
        assert R == R.transpose()
        bd = B.det()
        assert bd > 0
    else:
        R = fmpq_mat(h,h)
        bd = fmpq(1)
    return B,C,R,bd


def det_polynomial(R, h):
    ring = fmpq_mpoly_ctx.get(tuple(f'x{k}' for k in range(2*h-1)), ordering='lex')
    xs = ring.gens()
    M = [[math.factorial(3+2*a+2*b)*xs[a+b]-R[a,b]
          for b in range(h)] for a in range(h)]
    # Subset Laplace expansion, adding the last row at each stage.
    dp = {0:ring.constant(1)}
    for mask in range(1,1<<h):
        k = mask.bit_count()
        cols = [j for j in range(h) if mask >> j & 1]
        value = ring.constant(0)
        for pos,j in enumerate(cols):
            value += (-1)**(k-1+pos)*dp[mask^(1<<j)]*M[k-1][j]
        dp[mask] = value
    return dp[(1<<h)-1], ring


def primitive(P):
    coefficients = list(P.coeffs())
    denominator = math.lcm(*(int(c.denominator) for c in coefficients))
    numerators = [int(c*denominator) for c in coefficients]
    content = math.gcd(*numerators)
    q = fmpq(denominator, content)
    assert all((q*c).denominator == 1 for c in coefficients)
    assert math.gcd(*(abs(int(q*c)) for c in coefficients)) == 1
    return q


def evaluate(R,q,h,bits):
    ctx.prec = bits
    values = []
    for k in range(2*h-1):
        s = 4+2*k
        values.append((arb(s).zeta(arb(1)/4)-arb(s).zeta(arb(3)/4))/arb(4)**s)
    S = arb_mat([[math.factorial(3+2*a+2*b)*values[a+b]-arb(R[a,b])
                  for b in range(h)] for a in range(h)])
    sd = S.det()
    if not sd > 0:
        raise ArithmeticError('Insufficient precision to enclose positive determinant')
    pv = arb(q)*sd
    return {'precision_bits':bits, 'det_S':sd.str(40),
            'log_det_S':sd.log().str(40), 'log_q':arb(q).log().str(40),
            'log_primitive_value':pv.log().str(40),
            'primitive_value':pv.str(40), 'primitive_above_one':bool(pv>1),
            'primitive_below_one':bool(pv<1)}


def main(output):
    records=[]
    for file in sorted((ROOT/'evidence/gram').glob('h*-r*.json')):
        data=json.loads(file.read_text());h=data['h'];r=data['r']
        B,C,R,bd=blocks(h,r)
        assert bd==fmpq(data['bulk_determinant'])
        assert [[str(R[a,b]) for b in range(h)] for a in range(h)]==data['R']
        P,ring=det_polynomial(R,h);q=primitive(P)
        assert q==fmpq(data['q'])
        assert len(P.coeffs())==data['polynomial_terms']
        target=[0]*(2*h-1);target[h-1]=h
        assert P.to_dict()[tuple(target)]==fmpq(data['fixed_coefficient'])
        low=[0]*(2*h-1)
        for a in range(1,h):low[2*a]+=1
        fac=math.prod(math.factorial(3+4*a) for a in range(1,h))
        assert P.to_dict().get(tuple(low),0)==-fac*R[0,0]
        if r>1:
            labels=[(a,j) for a in range(h) for j in range(r-1)]
            ell=lambda n:math.lcm(*range(1,n+1,2))
            hd=[ell(2*r-3)**(2+2*a) for a in range(h)]
            ed=[ell(2*i+2*r-1)**(2+2*a) for a,i in labels]
            M=fmpq_mat([[B[i,j]*ed[i]*ed[j] for j in range(len(labels))] for i in range(len(labels))])
            assert all(x.denominator==1 for x in M.entries())
            assert all((C[a,j]*hd[a]*ed[j]).denominator==1 for a in range(h) for j in range(len(labels)))
            multiplier=math.prod(hd)**2*M.det()
            assert multiplier>0 and (multiplier/q).denominator==1
        for offset in (0,1,3):
            xs=[fmpq(offset+k+1,k+2) for k in range(2*h-1)]
            direct=fmpq_mat([[math.factorial(3+2*a+2*b)*xs[a+b]-R[a,b] for b in range(h)] for a in range(h)]).det()
            assert P(*xs)==direct
        if (h,r)==(1,2):
            assert (q*P).to_dict()=={(1,):fmpq(82),(0,):fmpq(-81)}
        finite=evaluate(R,q,h,4096)
        for key in ('log_q','log_det_S','log_primitive_value'):
            assert arb(finite[key]).overlaps(arb(data[key]))
        if h==7:assert finite['primitive_above_one']
        record={'h':h,'r':r,'coefficients':len(P.coeffs()),'exact_primitive_normalization':True,
                'specializations_passed':3,**finite}
        records.append(record)
        print(json.dumps({'h':h,'r':r,'passed':True,'coefficients':len(P.coeffs())}),flush=True)
    Path(output).write_text(json.dumps({'all_passed':True,'scope':'Recomputed recorded finite cases only; no asymptotic inference.','records':records},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='review/gram-checks.json')
    main(p.parse_args().output)
