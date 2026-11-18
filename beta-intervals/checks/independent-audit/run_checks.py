"""Run the standalone finite checks in a temporary directory."""
from pathlib import Path
from shutil import copy2
from subprocess import run
from sys import executable
from tempfile import TemporaryDirectory

source = Path(__file__).resolve().parent
scripts = [
    "seed_recheck.py",
    "density_recheck.py",
    "phase_recheck.py",
    "dilation_recheck.py",
]
with TemporaryDirectory(prefix="beta-intervals-check-") as tmp:
    work = Path(tmp)
    for name in ["density_max.cpp", *scripts]:
        copy2(source / name, work / name)
    run(["g++", "-O3", "-std=c++17", "-fPIC", "-shared",
         "density_max.cpp", "-o", "density_max.so"], cwd=work, check=True)
    for name in scripts:
        print(name, flush=True)
        run([executable, name], cwd=work, check=True)
