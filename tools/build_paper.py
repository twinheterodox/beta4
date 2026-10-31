"""Build the main LaTeX paper and reject unresolved references or box overflow."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'paper/beta4-determinant-obstructions.tex'
BUILD = ROOT / 'build/latex'
OUTPUT = ROOT / 'output/pdf/beta4-determinant-obstructions.pdf'
BUILD.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
for i in range(3):
    run = subprocess.run([
        'pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-file-line-error',
        '-output-directory=' + str(BUILD), str(SOURCE)], cwd=ROOT,
        capture_output=True, text=True)
    (BUILD / f'build-pass-{i+1}.txt').write_text(run.stdout + run.stderr)
    if run.returncode:
        raise SystemExit(run.stdout[-6000:] + run.stderr)
log = SOURCE.with_suffix('.log').name
text = (BUILD / log).read_text()
failures = [line for line in text.splitlines() if
            'undefined' in line.lower() or 'multiply defined' in line.lower()
            or 'Overfull' in line or 'Underfull' in line]
if failures:
    raise SystemExit('LaTeX checks failed:\n' + '\n'.join(failures))
shutil.copyfile(BUILD / OUTPUT.name, OUTPUT)
info = subprocess.check_output(['pdfinfo', str(OUTPUT)], text=True)
pages = int(re.search(r'^Pages:\s+(\d+)', info, re.M).group(1))
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = {
    'status': 'compiled', 'pages': pages,
    'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': digest(SOURCE),
    'pdf': str(OUTPUT.relative_to(ROOT)), 'pdf_sha256': digest(OUTPUT),
    'tex_engine': subprocess.check_output(['pdflatex', '--version'], text=True).splitlines()[0],
    'unresolved_references': False, 'overfull_or_underfull_boxes': False,
    'scope': 'Compilation and log checks only; visual inspection and mathematics are recorded separately.'}
(ROOT / 'paper/build-record.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
