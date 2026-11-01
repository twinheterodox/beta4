"""One-time import of archived rational certificate data (not needed to replay).

Pass the original coordinator directory explicitly. Canonicalizes, records
source hashes, and retains the algebraic data rather than archived verdicts.
"""
from pathlib import Path
import json,sys,hashlib
from algebra import parse,Rat
from flint import fmpq
src=Path(sys.argv[1]);dest=Path(__file__).resolve().parent
provenance=[]
def read(name):
    p=src/name;provenance.append({'file':name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return json.loads(p.read_text())
k=read('adjoint-pairing.json')['K']
for label in ('quarter','threequarter'):
    a=read('deformed-matrices-'+label+'.json');b=read('verified-matrices-'+label+'.json')
    al=Rat(fmpq(a['alpha']));env={'alpha':al}
    out={'alpha':a['alpha']}
    for key,arr in {'M':a['M'],'Q':[r[:3] for r in a['G'][:3]],'C':b['C'],'CQ':b['CQ'],'K':k}.items():
        out[key]=[[parse(x,env).pack() for x in r] for r in arr]
    for key in ('c','detQ','lambda_step'):out[key]=parse(b[key],env).pack()
    (dest/(label+'.json')).write_text(json.dumps(out,separators=(',',':'))+'\n')
(dest/'provenance.json').write_text(json.dumps({'operation':'Canonical exact rational functions; no archived status flags imported.','sources':provenance},indent=2)+'\n')
