"""Separate Sturm root counts for the degree-200 column certificate."""
from pathlib import Path
import sys,json,time
from flint import fmpq as Q,fmpq_poly as Poly
HERE=Path(__file__).resolve().parent
# Companion module is in the same replay directory.
from column_certify import sturm

def main():
    start=time.monotonic();data=json.loads((HERE/'column_certificate_d200.json').read_text());out={}
    for key,domains,expected in [
        ('x',[(Q(-1),Q(0)),(Q(0),Q(1))],data['x_root_counts']),
        ('s',[(Q(0),Q(1))],data['s_root_counts'])]:
        P=Poly([Q(x) for x in data['derivative_'+key]]);seq=sturm(P)
        def variation(x):
            signs=[1 if z>0 else -1 for q in seq if (z:=q(x))!=0]
            return sum(a!=b for a,b in zip(signs,signs[1:]))
        counts=[variation(a)-variation(b) for a,b in domains]
        assert counts==expected
        brackets=data['x_root_brackets' if key=='x' else 's_root_brackets']
        previous=None
        for lo,hi in brackets:
            lo,hi=Q(lo),Q(hi)
            assert lo<hi and P(lo)*P(hi)<0
            if previous is not None:assert previous<=lo
            previous=hi
        out[key]=dict(root_counts=counts,sign_change_brackets=len(brackets))
        print(key,counts,'seconds',time.monotonic()-start,flush=True)
    out['seconds']=time.monotonic()-start
    out['scope']='Independent exact Sturm counts match the complete Bernstein--Descartes isolation; every recorded bracket is disjoint and has a sign change.'
    (HERE/'column_root_crosscheck.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
