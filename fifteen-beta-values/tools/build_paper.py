#!/usr/bin/env python3
"""Build the canonical PDF and reject unresolved references or overfull boxes."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "paper/fifteen-beta-values.tex"
BUILD = ROOT / "build/latex"
OUTPUT = ROOT / "output/pdf/fifteen-beta-values.pdf"
BUILD.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
for i in range(3):
    run = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
         "-file-line-error", "-output-directory="+str(BUILD), str(SOURCE)],
        cwd=ROOT, capture_output=True, text=True)
    (BUILD / f"pass-{i+1}.txt").write_text(run.stdout+run.stderr)
    if run.returncode:
        raise SystemExit(run.stdout[-6500:]+run.stderr)
log = (BUILD / "fifteen-beta-values.log").read_text()
issues = [line for line in log.splitlines()
          if "undefined" in line.lower() or "multiply defined" in line.lower()
          or "Overfull" in line]
if issues:
    raise SystemExit("LaTeX review required:\n"+"\n".join(issues))
shutil.copyfile(BUILD / OUTPUT.name, OUTPUT)
info = subprocess.check_output(["pdfinfo", str(OUTPUT)], text=True)
pages = int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1))
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
record = {
    "status": "compiled; correctness and visual reviews recorded separately",
    "pages": pages,
    "source_files": {str(p.relative_to(ROOT)): digest(p)
                     for p in sorted((ROOT / "paper").rglob("*.tex"))},
    "pdf": str(OUTPUT.relative_to(ROOT)), "pdf_sha256": digest(OUTPUT),
    "unresolved_references": False, "overfull_boxes": False,
    "underfull_warnings": [line for line in log.splitlines() if "Underfull" in line],
    "independent_audit_record": "review/independent-audit-2026-11-11.json",
    "review_scope": "Compilation only. The independent audit identifies the preceding 14-page PDF; subsequent changes are recorded in review/audit-response-2026-11-11.json, review/exposition-update-2026-11-11.json and review/publication-preparation-2026-11-11.json.",
}
(ROOT / "paper/build-record.json").write_text(json.dumps(record, indent=2)+"\n")
print(json.dumps({"pages": pages, "pdf": record["pdf"],
                  "pdf_sha256": record["pdf_sha256"],
                  "unresolved_references": False, "overfull_boxes": False}, indent=2))
