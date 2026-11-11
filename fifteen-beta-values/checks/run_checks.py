#!/usr/bin/env python3
"""Replay the finite certificates; use --with-product for the FLINT calculation."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-product", action="store_true",
                        help="also regenerate the seed by truncated products (requires python-flint)")
    args = parser.parse_args()
    build = ROOT / "build/checks"
    build.mkdir(parents=True, exist_ok=True)
    steps = []

    def run(label, command):
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        (build / (label + ".txt")).write_text(result.stdout + result.stderr)
        if result.returncode:
            raise SystemExit(label + " failed:\n" + result.stdout + result.stderr)
        steps.append(label)
        print(label + ": passed", flush=True)
        return result.stdout

    py = [sys.executable, "-B"]
    run("arithmetic_and_concavity", py + [str(HERE / "certify_fifteen.py")])
    run("phase_polynomial", py + [str(HERE / "certify_phase_polynomial.py")])
    binary = build / "verify_fifteen_arithmetic"
    run("compile_direct_floor_checker", ["g++", "-O2", "-std=c++17",
         str(HERE / "verify_fifteen_arithmetic.cpp"), "-o", str(binary)])
    floor_result = json.loads(run("direct_floor_cells",
                             [str(binary), str(HERE / "fifteen-arithmetic-midpoints.tsv")]))
    assert floor_result["direct_floor_intervals"] == 38371
    (HERE / "fifteen-arithmetic-second-check.json").write_text(json.dumps(floor_result, indent=2)+"\n")
    if args.with_product:
        for p in (10007, 10009):
            run(f"seed_product_{p}", py + [str(HERE / "seed_fifteen.py"), "--prime", str(p)])
    run("seed_recurrence", py + [str(HERE / "verify_fifteen_seed.py")])
    seed = json.loads((HERE / "fifteen-seed-second-check.json").read_text())
    assert [(v["prime"], v["seed_constant_residue"]) for v in seed["results"]] == [(10007,3738),(10009,5333)]
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    files = sorted(p for p in HERE.iterdir() if p.suffix in (".py", ".cpp", ".json", ".tsv")
                   and p.name != "check-report.json")
    report = {
        "status": "passed", "steps": steps,
        "product_calculation_replayed": args.with_product,
        "arithmetic_intervals": floor_result["direct_floor_intervals"],
        "seed_coefficients_compared": sum(v["coefficient_comparisons"] for v in seed["results"]),
        "seed_residues": {str(v["prime"]): v["seed_constant_residue"] for v in seed["results"]},
        "scope": "Finite certificate replay. General arguments remain mathematical proofs in the manuscript; no external audit or formal verification is asserted.",
        "file_sha256": {p.name: digest(p) for p in files},
    }
    (HERE / "check-report.json").write_text(json.dumps(report, indent=2)+"\n")
    print("All finite certificate checks passed.")

if __name__ == "__main__":
    main()
