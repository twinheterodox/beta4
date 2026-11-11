PRESERVED INDEPENDENT AUDIT CALCULATIONS
11 November 2026

The five programs and four JSON result files in this directory are
unchanged copies supplied with the independent mathematical audit.
Their hashes are recorded in review/independent-audit-2026-11-11.json
relative to the paper directory. The programs do not import the
manuscript's checkers. Comparisons with its saved data occur after
the independent calculations.

From the repository root, using an environment with the optional
python-flint dependency and a C++17 compiler:

    python3 -B fifteen-beta-values/checks/replay_independent_audit.py

See ../requirements-optional.txt for the pinned python-flint version.
The wrapper copies the auditor programs into the ignored
build/independent-audit-replay/ directory, compiles the floor checker,
supplies the three reference data files, and runs all four calculations.
It compares the generated JSON results with the auditor's preserved
results and writes checks/independent-audit-replay.json. No original
audit file is overwritten.

independent_arithmetic.py and independent_floor.cpp independently
enumerate the event partition and check both closed boundary points
and open cells at all 38,371 midpoints. They use 96-bit upward rounding.

independent_phase.py reconstructs the critical polynomial and uses
512-bit Arb arithmetic for the maximum and decay-rate enclosures.

independent_seed.py uses polynomial Bezout inversion at each pole,
followed by binomial translation, to recover every seed coefficient
at both primes. It compares all 4,154 coefficients with the supplied
rows after computing them.

independent_local_bounds.py supplements the general proof with 56
small-prime cases and 222 large-prime cases, checking 32,722 repeated
residue classes.

The audit's general mathematical conclusions and scope are in
review/independent-audit-2026-11-11.txt. Replaying these calculations
checks the computational inputs; it is not a new audit of the proof.
