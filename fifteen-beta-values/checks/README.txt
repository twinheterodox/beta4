FINITE CERTIFICATES FOR THE FIFTEEN-VALUE PAPER

From the repository root:

    python3 -B fifteen-beta-values/checks/run_checks.py

This requires Python 3.10 or later and a C++17 compiler named g++.
The Python calculations use the standard library. The command regenerates
the arithmetic partition, checks its 38,371 midpoint maxima by two
algorithms, checks the phase by polynomial and concavity arguments, and
recomputes all 2,077 seed coefficients at each of 10007 and 10009 by a
logarithmic-derivative recurrence. The saved truncated-product rows are
read only after that calculation, for comparison.

To regenerate those rows as well, install the optional pinned dependency
in a virtual environment and run:

    python3 -m pip install -r fifteen-beta-values/checks/requirements-optional.txt
    python3 -B fifteen-beta-values/checks/run_checks.py --with-product

The complete command above was run for the manuscript. The report is
check-report.json; it records whether the product calculation was replayed
and hashes the scripts and certificate data. Scratch executables and logs
go under ../build/checks/ and are not publication files.

WHAT EACH CHECK ESTABLISHES

certify_fifteen.py
  Enumerates affine-offset intersections and integer events on [1,400].
  Sweeps the jumps in the residue variable to find the maxima. Integrates
  the resulting step function using rational arithmetic with upward
  64-bit rounding, adds an explicit infinite-tail bound, and proves the
  arithmetic cost is below 41/5. Also certifies max Phi < -18.83 using
  concavity on an initial interval and a tail derivative bound.

verify_fifteen_arithmetic.cpp
  Recomputes the maximum at every supplied midpoint by evaluating the
  floor counts directly on each residue cell. It checks the maxima, not
  a separately enumerated event partition. The general completeness
  argument for that partition is in Section 4.2 of the manuscript.

certify_phase_polynomial.py
  Reconstructs the degree-33 critical polynomial, checks its coefficient
  signs, and brackets the unique positive root. Rational logarithm
  enclosures prove max Phi <= -18.835112315041. This script imports no
  other checker.

seed_fifteen.py
  Uses python-flint to multiply truncated finite-field series at all
  poles. Checks reflection and three rational-function evaluations per
  prime. Produces the coefficient rows and nonzero constant residues
  3738 modulo 10007 and 5333 modulo 10009.

verify_fifteen_seed.py
  Reconstructs the seed coefficients by logarithmic differentiation and
  a Taylor recurrence, using only Python's standard library. Checks both
  prime moduli by trial division, odd-index cancellation, and agreement
  with all 4,154 saved coefficients. It does not import seed_fifteen.py.

The finite calculations establish the numerical inputs specified above.
The denominator lemma, uniform analytic estimates, prime-dilation
argument and final irrationality contradiction are general proofs in
the manuscript. These replays are not a separate mathematical audit or
formal verification of those arguments.

INDEPENDENT AUDIT

The separate audit reported on 11 November 2026 checked the complete
argument and independently reproduced the arithmetic, phase and seed
inputs. Its original programs and results are in independent-audit/;
the report and audited-file hashes are in ../review/. To replay those
programs in a scratch directory, with the optional dependency installed:

    python3 -B fifteen-beta-values/checks/replay_independent_audit.py

The independent-audit/README.txt explains the algorithms and coverage.
