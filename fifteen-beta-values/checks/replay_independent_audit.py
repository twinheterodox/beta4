#!/usr/bin/env python3
"""Run the preserved auditor programs in a scratch directory and compare their results."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "independent-audit"
BUILD = ROOT / "build/independent-audit-replay"
BUILD.mkdir(parents=True, exist_ok=True)

for p in SOURCE.iterdir():
    if p.suffix in (".py", ".cpp"):
        shutil.copyfile(p, BUILD / p.name)
references = BUILD / "replay/checks"
references.mkdir(parents=True, exist_ok=True)
for name in ("fifteen-arithmetic-midpoints.tsv", "fifteen-seed-rows-p10007.json",
             "fifteen-seed-rows-p10009.json"):
    shutil.copyfile(HERE / name, references / name)

commands = [
    ("compile", ["g++", "-O2", "-std=c++17", str(BUILD / "independent_floor.cpp"),
                 "-o", str(BUILD / "independent_floor")]),
]
commands += [(name, [sys.executable, "-B", str(BUILD / (name + ".py"))])
             for name in ("independent_arithmetic", "independent_phase",
                          "independent_seed", "independent_local_bounds")]
for label, command in commands:
    run = subprocess.run(command, cwd=BUILD, capture_output=True, text=True)
    (BUILD / (label + ".log")).write_text(run.stdout + run.stderr)
    if run.returncode:
        raise SystemExit(label + " failed:\n" + run.stdout + run.stderr)
    print(label + ": passed", flush=True)

results = ["independent-arithmetic.json", "independent-phase.json",
           "independent-seed.json", "independent-local-bounds.json"]
for name in results:
    assert json.loads((BUILD / name).read_text()) == json.loads((SOURCE / name).read_text()), name
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = {
    "status": "passed",
    "scope": "Packaging replay of the auditor's preserved programs; not a new mathematical audit.",
    "all_four_results_equal_the_original_audit_results": True,
    "auditor_programs_sha256": {p.name: sha(p) for p in sorted(SOURCE.iterdir())
                               if p.suffix in (".py", ".cpp")},
    "result_files_sha256": {name: sha(BUILD / name) for name in results},
    "runner_sha256": sha(Path(__file__)),
}
(HERE / "independent-audit-replay.json").write_text(json.dumps(report, indent=2) + "\n")
print("All four results agree with the original audit.")
