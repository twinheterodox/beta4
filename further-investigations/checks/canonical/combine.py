"""Recombine the certified bounds, including the explicit Jacobi integrals."""
from pathlib import Path
import json
from flint import arb, fmpq as Q, ctx

HERE=Path(__file__).resolve().parent
ctx.prec=256
def read(name): return json.loads((HERE/name).read_text())
def antiderivative(x):
    x=arb(x)
    return arb(0) if x==0 else x*x*x.log()/2-x*x/4
def integral(u,v,t):
    u,v,t=map(Q,(u,v,t))
    result=sum((antiderivative(t+a)-antiderivative(a) for a in [Q(0),u,v,u+v]),arb(0))
    return result-antiderivative(2*t+u+v)+antiderivative(u+v)
row=read('potential_certificate.json')
repair=read('energy_certificate_d40.json')
column=read('column_certificate_d200.json')
base=arb(Q(61,72))*arb(2).log()+2*integral('1/4','0','11/24')+integral('13/6','0','1/12')
main=base+arb(Q(row['extra_logdet_upper_rational']))
assert Q(row['extra_logdet_upper_rational'])==Q(112701993,2500000000)
col=2*(arb(Q(-41,48))*arb(2).log()+arb(Q(column['rational_norm'])+Q(column['rational_x_upper'])+Q(column['rational_s_upper'])))
rounded=arb(Q(249,50000))
cost=arb(Q(33079,18432))+arb(Q(505,4608))*arb(2).log()
decay=-(main+rounded+col)/2
margin=decay-cost
assert margin>arb(Q(1888937116500327,10**20))
assert decay>arb('1.87063289244902359295')
out={k:v.str(45) for k,v in dict(main_row_upper=main,rounded_repair=rounded,column_upper=col,
     canonical_decay_lower=decay,arithmetic_cost_upper=cost,margin_lower=margin).items()}
out['scope']='Recomputed from rational trials and replayed certificates; asymptotic passage is a separate proof.'
(HERE/'combined-replay.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
