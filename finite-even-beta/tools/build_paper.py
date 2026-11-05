"""Build the canonical PDF and reject unresolved references or overfull boxes."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "paper/finite-even-beta.tex"
BUILD = ROOT / "build/latex"
OUTPUT = ROOT / "output/pdf/finite-even-beta.pdf"
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
log = (BUILD / "finite-even-beta.log").read_text()
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
    "status": "compiled; visual and mathematical reviews recorded separately",
    "pages": pages,
    "source_files": {str(p.relative_to(ROOT)): digest(p)
                     for p in sorted((ROOT / "paper").rglob("*.tex"))},
    "pdf": str(OUTPUT.relative_to(ROOT)), "pdf_sha256": digest(OUTPUT),
    "unresolved_references": False, "overfull_boxes": False,
    "underfull_warnings": [line for line in log.splitlines() if "Underfull" in line],
    "review_scope": "Typeset exposition of the previously audited plain-text proof. No new independent audit or formal verification is claimed.",
}
(ROOT / "paper/build-record.json").write_text(json.dumps(record, indent=2)+"\n")
print(json.dumps({"pages": pages, "pdf": record["pdf"],
                  "pdf_sha256": record["pdf_sha256"],
                  "unresolved_references": False, "overfull_boxes": False}, indent=2))
