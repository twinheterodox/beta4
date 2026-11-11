"""Build this supplement without changing any earlier paper."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'paper/further-investigations.tex'
BUILD = ROOT / 'build/latex'
OUTPUT = ROOT / 'output/pdf/further-investigations.pdf'
BUILD.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
for i in range(3):
    run = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
        '-file-line-error', '-output-directory='+str(BUILD), str(SOURCE)],
        cwd=ROOT, capture_output=True, text=True, errors='replace')
    (BUILD/f'pass-{i+1}.txt').write_text(run.stdout+run.stderr)
    if run.returncode:
        raise SystemExit(run.stdout[-6000:]+run.stderr)
log = (BUILD/'further-investigations.log').read_text(errors='replace')
failures = [s for s in log.splitlines() if any(t in s.lower() for t in
            ('undefined','multiply defined','overfull','underfull'))]
if failures:
    raise SystemExit('LaTeX checks failed:\n'+'\n'.join(failures))
shutil.copyfile(BUILD/OUTPUT.name, OUTPUT)
info = subprocess.check_output(['pdfinfo', str(OUTPUT)], text=True)
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = dict(status='compiled', pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1)),
    source_files=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p))
                  for p in sorted((ROOT/'paper').rglob('*.tex'))],
    pdf=str(OUTPUT.relative_to(ROOT)), pdf_sha256=digest(OUTPUT),
    unresolved_references=False, overfull_or_underfull_boxes=False,
    scope='Build and log checks; mathematical and visual reviews recorded separately.')
(ROOT/'review/build.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='source_files'},indent=2))
