#!/usr/bin/env python3
"""Build the canonical manuscript; keep intermediate products in build/."""
from pathlib import Path
import hashlib, json, re, shutil, subprocess

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'build'
BUILD.mkdir(exist_ok=True)
for _ in range(3):
    run = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                          '-file-line-error', '-output-directory=build',
                          'paper/arithmetic-investigations.tex'],
                         cwd=ROOT, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    (BUILD / 'terminal.txt').write_text(run.stdout)
    if run.returncode:
        print(run.stdout[-7000:])
        raise SystemExit(run.returncode)
log = (BUILD / 'arithmetic-investigations.log').read_text()
warnings = [line for line in log.splitlines()
            if any(s in line for s in ('Overfull', 'Underfull', 'undefined',
                                        'multiply defined', 'Rerun to get'))]
output = ROOT / 'output/pdf/arithmetic-investigations.pdf'
shutil.copy2(BUILD / output.name, output)
info = subprocess.check_output(['pdfinfo', str(output)], text=True)
record = {'pdf_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
          'pages': int(re.search(r'^Pages:\s+(\d+)', info, re.M).group(1)),
          'warnings': warnings,
          'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted((ROOT / 'paper').rglob('*.tex'))}}
(ROOT / 'review/build.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'}, indent=2))
