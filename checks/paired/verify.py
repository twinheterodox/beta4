"""Replay the all-index algebra, telescoping, boundaries, and initial values.

Run from any working directory with python-flint installed. This verifies
rational/polynomial identities, not a proof-assistant formalization. The
manuscript explains how the identities transfer to sums and derivative tails.
"""
from pathlib import Path
from fractions import Fraction
from math import prod
import json,hashlib,time,datetime,sys
import algebra as c
from problems import recurrence,fproblem,dproblem,dshift,power
from readouts import forward,dual,determinant,bs
from flint import fmpq,fmpq_poly
HERE=Path(__file__).resolve().parent
start=time.monotonic();records=[]

def krylov(m,comp):
    e=[[c.Rat(int(i==0))] for i in range(len(m))]
    v=c.mm(m,e);w=c.mm(m,c.mm(c.shifted(m,1),e))
    u=[[e[i][0],v[i][0],w[i][0]] for i in range(len(m))]
    assert c.mm(m,c.shifted(u,1))==c.mm(u,comp)
    return u

def matrix(raw):return [[c.unpack(x) for x in r] for r in raw]
def formula_lambda(n,a):
    from math import factorial
    poch=lambda a,k:prod((a+j for j in range(k)),start=Fraction(1))
    return (Fraction(factorial(8*n)*factorial(9*n)**2*factorial(11*n)**2*factorial(12*n),6*factorial(14*n)*factorial(16*n))
       *poch(a,11*n)*poch(a,12*n)**2*poch(a,13*n)**2*poch(a,14*n)
       / (poch(a,16*n)*poch(a,17*n)**2*poch(a,18*n)**2*poch(a,19*n)))

def scalar_recurrence(comp,al,forward_kind=False):
    g=[comp[i][2] for i in range(3)]
    if forward_kind:
        from problems import factstep
        rho=c.Rat(factstep(15)*factstep(12,al)*factstep(13,al)*factstep(14,al),
             factstep(10)**3*factstep(1)*factstep(14)*factstep(2,al)*factstep(3,al)*factstep(4,al))
        for j in range(3):
            for h in range(j,3):g[j]*=rho.shift(h)
    b=c.primitive([c.Rat(1),-g[2],-g[1],-g[0]])
    assert c.nonzero(b[0])
    return b

for label in ('quarter','threequarter'):
    raw=json.loads((HERE/(label+'.json')).read_text());al=fmpq(raw['alpha']);a=1+al
    M,Q,K,C,CQ=[matrix(raw[key]) for key in ('M','Q','K','C','CQ')]
    scalar=c.unpack(raw['c']);dq=c.det(Q);assert dq==c.unpack(raw['detQ'])
    assert c.mm(c.mm(M,c.shifted(K,1)),c.transpose(Q))==[[scalar*x for x in r] for r in K]
    E=c.Rat((8*c.N+1)*(9*c.N+1)*(10*c.N+1)*(11*c.N+a))
    u=[[c.Rat(1)],[c.Rat(10*c.N+a)],[c.Rat((10*c.N+a)**2)]]
    assert c.mm(K,u)==[[-E],[c.Rat()],[c.Rat()],[c.Rat()]]
    U=krylov(M,C);UQ=krylov(Q,CQ)
    for m in (M,Q,K,C,CQ,U,UQ):
        assert all(c.nonzero(x.q) for row in m for x in row)
    for m in (Q,C,CQ,UQ):
        d=c.det(m);assert c.nonzero(d.p) and c.nonzero(d.q)
    # Verify the actual initial readouts, including the degenerate n=0 case.
    af=Fraction(str(al));rho=[[1,0,0,0],[0,0,Fraction(-1,2),(2*af+1)/4],[0,0,0,Fraction(-1,6)]]
    initial=dual(0,af);k0=[[Fraction(str(x.at(0))) for x in r] for r in K]
    def cross(v,w):return [v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]]
    products=[cross(initial[1],initial[0]),cross(initial[2],initial[0]),cross(initial[2],initial[1])]
    for j in range(3):
        assert [sum(rho[j][i]*k0[i][h] for i in range(4)) for h in range(3)]==[(-1 if j!=1 else 1)*(1+af)*x/6 for x in products[j]]
    # Closed scalar has the same one-step ratio as the matrix pairing.
    from problems import factstep
    rC=c.Rat(factstep(15)*factstep(12,al)*factstep(13,al)*factstep(14,al),
         factstep(10)**3*factstep(1)*factstep(14)*factstep(2,al)*factstep(3,al)*factstep(4,al))
    ratio=rC*scalar*E/E.shift(1)/dq
    candidate=c.Rat(1)
    for j,e in {8:1,9:2,11:2,12:1,14:-1,16:-1}.items():candidate*=c.Rat(factstep(j))**e
    for j,e in {11:1,12:2,13:2,14:1,16:-1,17:-2,18:-2,19:-1}.items():candidate*=c.Rat(factstep(j,al))**e
    assert candidate==ratio==c.unpack(raw['lambda_step'])
    print(label,'matrix transport, regularity and normalization PASS',flush=True)
    degrees=[]
    B=scalar_recurrence(C,al,True);A,Bm,F,D,norm=fproblem(B,al)
    P=c.solve(A,Bm,F,114);degrees.append([max(x.degree() for x in P),len(P)-1])
    assert c.nonzero(c.evaluate_k(D,1-12*c.N-al),2)
    assert 5+len(P)-1-(len(D)-1)-2<=0
    DB=scalar_recurrence(CQ,al);A,Bm,D=dproblem(3,al);F=[]
    for j in range(4):F=c.add(F,c.scale(dshift(j,3,al),DB[3-j]))
    P=c.solve(A,Bm,F,107);degrees.append([max(x.degree() for x in P),len(P)-1])
    A,Bm,D=dproblem(1,al)
    for j in range(3):
        den=c.ONE
        for i in range(3):den=c.poly_lcm(den,Q[i][j].q)
        F=c.scale(c.multiply(dshift(1,1,al),power(j,al)),den)
        for i in range(3):F=c.add(F,c.scale(c.multiply(D,power(i,al)),-Q[i][j].p*(den//Q[i][j].q)))
        P=c.solve(A,Bm,F,35);degrees.append([max(x.degree() for x in P),len(P)-1])
        assert c.nonzero(den) and c.nonzero(c.evaluate_k(D,1-16*c.N-al),1)
        assert 5+len(P)-1-(len(D)-1)-4<=0
    # Finite exceptional indices are part of the proof, not extrapolation.
    vals={n:forward(n,af) for n in range(5)}
    for n in (0,1):
        for j in (1,2):assert sum(Fraction(str(B[i](n)))*vals[n+3-i][j] for i in range(4))==0
    d0=dual(0,af);d1=dual(1,af)
    assert [[sum(d0[r][i]*Fraction(str(Q[i][j].at(0))) for i in range(3)) for j in range(3)] for r in range(3)]==d1
    walk=[[c.Rat(int(i==j)) for j in range(4)] for i in range(4)]
    for n in range(3):
        dn=dual(n,af);un=[Fraction(1),10*n+1+af,(10*n+1+af)**2];lam=formula_lambda(n,1+af)
        assert vals[n]==[lam*determinant(dn[1],dn[0],un),-lam*determinant(dn[2],dn[0],un),lam*determinant(dn[2],dn[1],un)]
        state=[sum(rho[j][i]*Fraction(str(walk[i][0].at(0))) for i in range(4))*bs(n,0,af)[0] for j in range(3)]
        assert state==vals[n]
        walk=c.mm(walk,[[c.Rat(x.at(n)) for x in r] for r in M])
    records.append({'alpha':str(al),'polynomial_identity_degrees':degrees,'matrix_transport':True,'regularity':True,'closed_scalar_ratio':True,'forward_exceptions':[0,1],'dual_exception':0,'initial_readout_indices':[0,1,2]})
    print(label,'five polynomial certificates and all boundary/initial checks PASS',flush=True)
input_paths = sorted(HERE.glob('*.py')) + [HERE / name for name in
    ('quarter.json', 'threequarter.json', 'provenance.json')]
report={'status':'passed','run_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seconds':round(time.monotonic()-start,3),'quarters':records,'inputs':[{'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in input_paths],'scope':'Exact identities and stated boundary conditions; analytic transfer is explained in the manuscript. No independent or formal verification claim.'}
(HERE/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',report['seconds'],'seconds',flush=True)
