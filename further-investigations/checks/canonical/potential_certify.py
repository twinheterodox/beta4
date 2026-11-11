"""Global interval certificates for the rational interpolation potentials."""
from pathlib import Path
import json,time
from flint import fmpq as Q,arb,ctx
HERE=Path(__file__).resolve().parent
def F(x):return arb(0) if x==0 else x*x.log()
def field(x):
    t=arb(Q(11,24));a=arb(Q(1,4));v=1-x*x
    if v==0:return arb(0)
    u=2*v*t*(t+a)/(v*(2*t+a)+(v*v*a*a+4*v*t*(t+a)).sqrt())
    return F(t)+F(t+a)-F(t-u)-F(t+a-u)-2*F(u)+u*v.log()
def field_derivative(x):
    if x==0:return arb(0)
    c=arb(Q(7,6));a=arb(Q(1,4))
    return (c-((c*c-a*a*x*x)/(1-x*x)).sqrt())/x
def certify(rec,kind):
    p=[(Q(x),Q(w)) for x,w in rec['p']];q=[(Q(x),Q(w)) for x,w in rec['q']]
    pa=[(arb(x),arb(w)) for x,w in p];qa=[(arb(x),arb(w)) for x,w in q]
    target=Q(rec['proposed_a' if kind=='output' else 'proposed_b'])
    edges=sorted(set(([Q(0),Q(1)]+[x for x,w in p]) if kind=='output' else ([Q(-1),Q(1)]+[x for x,w in q])))
    stack=[(lo,hi,0) for lo,hi in zip(edges,edges[1:])];leaves=evals=depth=0
    while stack:
        lo,hi,d=stack.pop();mid=(lo+hi)/2;x=arb(mid);hh=arb((hi-lo)/2)
        if kind=='output':
            val=arb(Q(13,12))*x.log();der=arb(Q(13,12))/x;curv=arb(0)
            for z,w in pa:val+=w*abs(x-z).log();der+=w/(x-z)
            for (y,w),(yr,wr) in zip(qa,q):
                val-=w*(1-y*x).log();der+=w*y/(1-y*x)
                endpoint=hi if yr>=0 else lo
                curv+=w*y*y/(1-y*arb(endpoint))**2
        else:
            val=field(x);der=field_derivative(x);curv=arb(0)
            for y,w in qa:val+=w*abs(x-y).log();der+=w/(x-y)
            for z,w in pa:
                val-=w*(1-z*x).log();der+=w*z/(1-z*x)
                curv+=w*z*z/(1-z*arb(hi))**2
        upper=val+abs(der)*hh+curv*hh*hh/2;evals+=1
        if upper<arb(target):leaves+=1;depth=max(depth,d);continue
        assert not val>arb(target),(kind,rec['r'],str(mid),str(val),str(target))
        assert d<60,(kind,rec['r'],'subdivision unresolved',str(lo),str(hi))
        stack.extend([(lo,mid,d+1),(mid,hi,d+1)])
    return dict(kind=kind,upper=str(target),leaves=leaves,interval_evaluations=evals,max_depth=depth)
def main():
    start=time.monotonic();ctx.prec=192
    data=json.loads((HERE/'potential_trials.json').read_text());records=[]
    for rec in data['records']:
        r=Q(rec['r']);assert sum(Q(w) for x,w in rec['p'])==r+Q(1,12)
        assert sum(Q(w) for x,w in rec['q'])==r
        assert all(Q(w)>0 and 0<Q(x)<1 for x,w in rec['p'])
        assert all(Q(w)>0 and -1<Q(x)<1 for x,w in rec['q'])
        a=certify(rec,'output');b=certify(rec,'input');e=Q(a['upper'])+Q(b['upper'])
        assert e>=0
        out=dict(r=str(r),output=a,input=b,envelope_upper=str(e));records.append(out)
        print(json.dumps(dict(r=str(r),envelope=str(e),checks=a['interval_evaluations']+b['interval_evaluations'],seconds=time.monotonic()-start)),flush=True)
    area=Q(0)
    for a,b in zip(records,records[1:]):area+=(Q(b['r'])-Q(a['r']))*(Q(a['envelope_upper'])+Q(b['envelope_upper']))
    area+=2*(Q(11,12)-Q(records[-1]['r']))*Q(records[-1]['envelope_upper'])
    out=dict(records=records,extra_logdet_upper_rational=str(area),extra_logdet_upper=arb(area).str(35),
        precision_bits=ctx.prec,seconds=time.monotonic()-start,
        scope='Global upper bounds using concavity and Arb subdivision. Operator and asymptotic passage are proved separately.')
    (HERE/'potential_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(area=str(area),decimal=arb(area).str(25),seconds=time.monotonic()-start)))
if __name__=='__main__':main()
