"""Replay selected certificates in fresh directories, preserving evidence files."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
SUITES = {
    'isolation': ['verify.py', 'check_dilation.py', 'polynomial_gate.py'],
    'pade': ['verify.py'],
    'digamma': ['nonvanishing_checks.py'],
    'canonical': ['potential_certify.py', 'energy_certify.py',
                  'column_certify_fast.py', 'column_crosscheck.py',
                  'slope-certificate.py', 'combine.py'],
    'small': ['verify.py'],
    'algebra': ['verify.py'],
}

def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are part of the checks.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--suite', choices=['all', *SUITES], default='all')
    args = parser.parse_args()
    manifest = json.loads((ROOT/'evidence/manifest.json').read_text())
    for row in manifest['files']:
        actual = hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()
        assert actual == row['sha256'], row['path']
    parent = ROOT/'build/checks'
    parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='replay-', dir=parent))
    selected = SUITES if args.suite == 'all' else {args.suite: SUITES[args.suite]}
    report = dict(status='running', scope='Finite computational replay; not an independent proof audit.',
                  suites={}, work_directory=str(work.relative_to(ROOT)))
    for suite, scripts in selected.items():
        dest = work/suite
        dest.mkdir()
        for source in [ROOT/'evidence'/suite, ROOT/'checks'/suite]:
            if source.exists():
                for path in source.iterdir():
                    if path.suffix in ('.py', '.json'):
                        shutil.copyfile(path,dest/path.name)
        report['suites'][suite] = []
        for script in scripts:
            print(f'Checking {suite}/{script}',flush=True)
            start=time.monotonic()
            run=subprocess.run([sys.executable,'-B',script],cwd=dest,
                               capture_output=True,text=True)
            (dest/(script+'.log')).write_text(run.stdout+run.stderr)
            result=dict(script=script,passed=run.returncode==0,
                        seconds=round(time.monotonic()-start,3))
            report['suites'][suite].append(result)
            if run.returncode:
                report['status']='failed'
                (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
                raise SystemExit(run.stdout[-5000:]+run.stderr[-5000:])
            print(f'  passed ({result["seconds"]} seconds)',flush=True)
    report['status']='passed'
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'review'/('replay-'+args.suite+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
