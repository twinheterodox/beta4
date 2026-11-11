"""Global integral certificate for the repair energy dual."""
from pathlib import Path
from fractions import Fraction
import argparse,json,time
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
HERE=Path(__file__).resolve().parent

def run(degree,grid):
    start=time.monotonic();ctx.prec=512
    data=json.loads((HERE/f'energy_opt_d{degree}.json').read_text())
    vv=[]
    for v in data['coefficients']:
        f=Fraction(f'{v:.12f}');vv.append(Q(f.numerator,f.denominator))
    t=Poly([0,1]);T=[Poly([1]),t]
    for k in range(2,2*degree):T.append(2*t*T[-1]-T[-2])
    P=sum((16*v*T[2*i+1]/(2*i+1) for i,v in enumerate(vv)),Poly([]));D=P.derivative()
    PA=arb_poly(P);DA=arb_poly(D)
    norm=8*sum((v*v/Q(2*i+1) for i,v in enumerate(vv)),Q(0))
    curvature=sum((Q(16,3)*abs(v)*(2*i+1)*((2*i+1)**2-1) for i,v in enumerate(vv)),Q(0))
    width=arb(Q(1,grid));bump=arb(curvature)*width**2/8;integral=arb(0)
    cases=dict(positive=0,negative=0,crossing=0,uncertain=0)
    for j in range(grid):
        x=arb(Q(2*j+1,2*grid));log=(54/x).log()
        value=2*x*(log+1)-PA(x);derivative=2*log-DA(x)
        center=value+bump;left=center-derivative*width/2;right=center+derivative*width/2
        if left>=0 and right>=0:bound=width*center;cases['positive']+=1
        elif left<=0 and right<=0:bound=arb(0);cases['negative']+=1
        elif (left>0 and right<0) or (left<0 and right>0):
            positive=left if left>0 else right
            bound=width*positive**2/(2*abs(right-left));cases['crossing']+=1
        else:
            bound=width*arb(0).max(left).max(right);cases['uncertain']+=1
        integral+=bound
    total=(arb(norm)+integral)/576
    target=Q(249,50000)  # 0.00498
    assert total<arb(target),total
    out=dict(degree=degree,grid=grid,rational_coefficients=[str(v) for v in vv],
             rational_norm=str(norm),global_second_derivative_upper=str(curvature),
             positive_part_integral_upper=integral.str(40),repair_rate_upper=total.str(40),
             rounded_repair_rate_upper=str(target),precision_bits=ctx.prec,cases=cases,
             seconds=time.monotonic()-start,
             scope='Global tangent-line integral upper bound. The general determinant energy inequality is proved separately.')
    (HERE/f'energy_certificate_d{degree}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['rational_coefficients','rational_norm']},indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degree',type=int,default=40);p.add_argument('--grid',type=int,default=16384)
    a=p.parse_args();run(a.degree,a.grid)
