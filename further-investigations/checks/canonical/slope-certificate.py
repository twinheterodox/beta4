"""Exclude rational roots of all three seed determinant pencils."""
from pathlib import Path
import json,time
from flint import fmpq as Q,fmpq_mat as Mat,nmod_mat as Mod,nmod_poly as Poly
HERE=Path(__file__).resolve().parent

start=time.monotonic();data=json.loads((HERE/'nonvanishing-seeds.json').read_text())
B=Mat([[Q(x) for x in row] for row in data['rational_base_matrix']])
J=Mat([[Q(x) for x in row] for row in data['rational_correction_matrix']])
records=[]
for sigma,p in ((0,103),(1,107),(-1,103)):
    A=Mat([[B[r,k]+sigma*B[r+1,k] for k in range(48)] for r in range(48)])
    C=Mat([[J[r,k]+sigma*J[r+1,k] for k in range(48)] for r in range(48)])
    relations=[]
    for r in range(4):
        partner=8-r if sigma==0 else 7-r
        sign=1 if sigma in (0,1) else -1
        assert all(C[r,k]==sign*C[partner,k] for k in range(48))
        relations.append(dict(row=r,partner=partner,sign=sign))
    def reduce(M):return Mod([[int(x.p)*pow(int(x.q),-1,p)%p for x in row] for row in M.tolist()],p)
    am,cm=reduce(A),reduce(C);det0=am.det();assert det0
    cp=(-am.inv()*cm).charpoly()
    polynomial=Poly([int(det0*cp[48-j]) for j in range(49)],p)
    assert polynomial.degree()==44
    values=[]
    for x in range(p):
        value=int((am+x*cm).det());assert value==int(polynomial(x)) and value!=0
        values.append(value)
    rec=dict(sigma=sigma,prime=p,rational_degree_bound=44,
        exact_row_relations_for_coefficient_matrix=relations,
        coefficients_mod_prime=list(map(int,polynomial)),degree_mod_prime=polynomial.degree(),
        leading_coefficient_mod_prime=int(polynomial[44]),finite_field_values=values,
        roots_mod_prime=[],no_rational_root=True)
    records.append(rec)
    print(json.dumps({k:rec[k] for k in ['sigma','prime','degree_mod_prime','leading_coefficient_mod_prime','no_rational_root']}),flush=True)
out=dict(scope='Exact certificate: every fixed rational slope has nonzero seed determinants. It supplies prime-scale nonvanishing, not a real decay or irrationality theorem.',
    certificates=records,elapsed_seconds=time.monotonic()-start)
(HERE/'slope-certificate.json').write_text(json.dumps(out,indent=2)+'\n')
