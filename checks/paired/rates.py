"""Exact limiting-system comparisons and rigorous numerical rate enclosures."""
from pathlib import Path
import json,datetime
from flint import fmpq,fmpq_mat,fmpq_poly,arb,ctx
from algebra import unpack,mm,det,Rat
HERE=Path(__file__).resolve().parent
ctx.dps=70
limits=[]
for label in ('quarter','threequarter'):
    data=json.loads((HERE/(label+'.json')).read_text())
    Q=[[unpack(x) for x in row] for row in data['Q']]
    out=[]
    for i in range(3):
        row=[]
        for j in range(3):
            x=Q[i][j];degree=x.p.degree()-x.q.degree()+i-j
            assert degree<=0
            row.append(x.p.leading_coefficient()/x.q.leading_coefficient() if degree==0 else fmpq(0))
        out.append(row)
    limits.append(out)
assert limits[0]==limits[1]
a=limits[0];mat=fmpq_mat(a);cp=mat.charpoly();p=Rat
u=[[p(1)],[p(10)],[p(100)]];ar=[[p(x) for x in row] for row in a]
v=mm(ar,u);w=mm(ar,v);cyclic=det([[u[i][0],v[i][0],w[i][0]] for i in range(3)])
assert cyclic and cyclic.p.degree()==0
assert cp[3]==1
# Sign brackets and exact Sturm counts isolate the unique small real root.
# The cubic discriminant is negative, so the other roots are a conjugate pair.
A,B,C=cp[2],cp[1],cp[0]
disc=A*A*B*B-4*B**3-4*A**3*C-27*C*C+18*A*B*C
assert disc<0
def enclose(poly,lo,hi):
    assert poly(lo)<0<poly(hi)
    for _ in range(250):
        mid=(lo+hi)/2
        if poly(mid)<0:lo=mid
        else:hi=mid
    return arb((lo+hi)/2, (hi-lo)/2)
r=enclose(cp,fmpq(7,10**24),fmpq(8,10**24))
lg=(arb(mat.det()).log()-r.log())/2;g=lg-r.log()
x=fmpq_poly([0,1]);poly=fmpq_poly([1]);rhs=fmpq_poly([1])
for j in (0,1,2,3,4):poly*=x+j
for j in (12,13,14):rhs*=x+j
poly-=rhs*(14-x)*(15-x)
xs=enclose(poly,fmpq(0),fmpq(14));L=lambda x:x*x.log()
sigma=sum(L(xs+j)-L(arb(10))-L(xs+j-10) for j in (12,13,14))
sigma+=2*L(arb(15))-L(xs+1)-L(14-xs)-L(xs)-L(15-xs)
c0=14*arb(15).log()+16*arb(16).log()-11*arb(2).log()-11*arb(3).log()-8*arb(4).log()
assert c0<lg and g>arb('112.80753') and g<arb('112.80754')
assert sigma>arb('65.45153') and sigma<arb('65.45154')
ceiling=json.loads((HERE/'arithmetic-result.json').read_text())
upper=arb(fmpq(ceiling['combined_upper']));needed=2*sigma+180*arb(2).log()-g
assert needed-upper>arb('54.65')
out={'status':'passed','run_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'limiting_matrix':[[str(x) for x in row] for row in a],'characteristic_polynomial':[str(x) for x in cp.coeffs()],'cyclic_vector_determinant':str(cyclic.p[0]/cyclic.q[0]),'g':g.str(30),'sigma':sigma.str(30),'crude_state_exponent':c0.str(30),'growing_mode_log_modulus':lg.str(30),'required_odd_saving':needed.str(30),'certificate_ceiling':upper.str(30),'shortfall':(needed-upper).str(30),'scope':'Exact algebra and Arb enclosures for the rate comparison. The analytic recurrence theorem and its transfer hypotheses are stated in the supplement.'}
(HERE/'rates-result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('limiting_matrix','characteristic_polynomial','cyclic_vector_determinant')},indent=2))
