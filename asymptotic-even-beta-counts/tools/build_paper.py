"""Build the paper and check the final LaTeX log for reference and layout errors."""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STEM = "asymptotic-even-beta-counts"
source = Path("paper") / (STEM + ".tex")
build = Path("build/latex")
output = Path("output/pdf") / (STEM + ".pdf")
(ROOT / build).mkdir(parents=True, exist_ok=True)
(ROOT / output).parent.mkdir(parents=True, exist_ok=True)

# Fix the manuscript date and omit the path-dependent PDF trailer identifier.
env = dict(os.environ, SOURCE_DATE_EPOCH="1794916800", FORCE_SOURCE_DATE="1", TZ="UTC")
for run in (1, 2):
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
               "-file-line-error", f"-jobname={STEM}", f"-output-directory={build}",
               r"\pdftrailerid{}\input{" + source.as_posix() + "}"]
    result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
    (ROOT / build / f"pass-{run}.txt").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise SystemExit(result.stdout[-6500:] + result.stderr)

log = (ROOT / build / (STEM + ".log")).read_text()
issues = [line for line in log.splitlines() if "undefined" in line.lower()
          or "multiply defined" in line.lower() or "Overfull" in line
          or "Rerun to get" in line or "has changed. Rerun" in line]
if issues:
    raise SystemExit("\n".join(issues))
shutil.copyfile(ROOT / build / output.name, ROOT / output)
info = subprocess.check_output(["pdfinfo", output.as_posix()], cwd=ROOT, text=True)
record = {
    "source": source.as_posix(),
    "source_sha256": hashlib.sha256((ROOT / source).read_bytes()).hexdigest(),
    "pdf": output.as_posix(),
    "pdf_sha256": hashlib.sha256((ROOT / output).read_bytes()).hexdigest(),
    "pages": int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1)),
    "latex_engine": subprocess.check_output(["pdflatex", "--version"], text=True).splitlines()[0],
    "unresolved_references": False,
    "overfull_boxes": False,
}
(ROOT / "paper/build-record.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
