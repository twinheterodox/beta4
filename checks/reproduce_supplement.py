"""Replay selected supplement certificates in scratch copies; preserve evidence."""
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import mpmath
import sympy

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'checks/supplement-results'
DEST.mkdir(parents=True,exist_ok=True)
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def verify_frozen():
    count=0
    for name in ('manifest.json','supplement-manifest.json'):
        data=json.loads((ROOT/'evidence'/name).read_text())
        for item in data['files']:
            assert digest(ROOT/item['path'])==item['sha256'],item['path']
            count+=1
    return count

start=time.monotonic()
before=verify_frozen()
jobs=[
 ('mellin','collaboration/round-03/fresh/mellin_obstruction_checks.py','mellin_obstruction_checks.json'),
 ('boole','collaboration/round-05/fresh/exact_checks.py','exact_checks.json'),
 ('cyclotomic','collaboration/round-11/norms/check_sine_norms.py','checks.json'),
]
records=[]
with tempfile.TemporaryDirectory(prefix='beta4-supplement-') as temporary:
    for label,rel,outname in jobs:
        source=ROOT/'evidence/supplement-source'/rel
        folder=Path(temporary)/label; folder.mkdir()
        script=folder/source.name; shutil.copyfile(source,script)
        run=subprocess.run([sys.executable,str(script)],cwd=folder,
                           capture_output=True,text=True,timeout=180)
        if run.returncode: raise RuntimeError(run.stdout+run.stderr)
        generated=folder/outname
        actual=json.loads(generated.read_text())
        expected=json.loads(source.with_name(outname).read_text())
        assert actual==expected,label+' differs from archived output'
        target=DEST/(label+'.json'); shutil.copyfile(generated,target)
        records.append({'label':label,'source':str(source.relative_to(ROOT)),
            'source_sha256':digest(source),'result':str(target.relative_to(ROOT)),
            'sha256':digest(target),'archived_output_matches':True})
spec=importlib.util.spec_from_file_location('supplement_identities',ROOT/'checks/supplement_identities.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
identities=module.run()
target=DEST/'rewritten-identities.json'
target.write_text(json.dumps(identities,indent=2)+'\n')
records.append({'label':'rewritten-identities','source':'checks/supplement_identities.py',
    'source_sha256':digest(ROOT/'checks/supplement_identities.py'),
    'result':str(target.relative_to(ROOT)),'sha256':digest(target)})
assert verify_frozen()==before
report={'status':'passed','run_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'python_version':sys.version.split()[0],'sympy_version':sympy.__version__,
    'mpmath_version':mpmath.__version__,'seconds':round(time.monotonic()-start,3),
    'frozen_files_unchanged':before,'original_manifest_files_unchanged':120,
    'runs':records,'new_exact_checks':identities,
    'scope':'Three archived certificate replays and finite exact checks of rewritten formulas; no independent model audit or formal verification.',
    'numerical_diagnostics':'Cyclotomic mpmath values are 100-digit diagnostics, not interval certificates. No new ordinary-moment intervals were computed.'}
(ROOT/'checks/supplement-reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='runs'},indent=2))
