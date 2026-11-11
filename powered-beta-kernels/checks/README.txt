EXACT CHECKS FOR POWERED HYPERGEOMETRIC BETA KERNELS

Run the following commands from the powered-beta-kernels directory.
Do not use Python's -O option: assertions perform the checks.

Standard-library checks

  python3 -B checks/check_power_seed.py checks/results/power_seed_k1_p10007.json checks/results/power_seed_k2_p10007.json checks/results/power_seed_k2_p10009.json checks/results/power_seed_k3_p10007.json checks/results/power_seed_k3_p10009.json --output checks/POWER_SEED_VALIDATION.json
  python3 -B checks/check_krisman_budget.py

check_power_seed.py reconstructs every pole independently using a
logarithmic-derivative recurrence and factorial/harmonic prefix sums.
It checks the five reported seed constants and all 191 even-beta
coefficients, including the control case q=1. It reconstructs the pole
profile and verifies primality and the required denominator units.

check_krisman_budget.py reconstructs all 62,132 event points and 62,131
interval maxima, the exact upward-rounded arithmetic sum, and rational
enclosures proving Q < -0.0714. The all-index asymptotic passage is a
mathematical argument in the cited paper, not a finite computational test.

Checks requiring python-flint

Install python-flint==0.9.0 in an isolated environment if needed. The
complete dependency list is in requirements.txt. Then run:

  python3 -B checks/check_prime_dilation.py
  python3 -B checks/power_kernel_seed.py --power 2 --prime 10007 --output checks/results/power_seed_k2_p10007.json

Repeat the seed-generation command for (power, prime) equal to
(1,10007), (2,10009), (3,10007) and (3,10009). The primary generator
uses finite-field polynomial products and series inversion. The two
seed algorithms were developed in the same research session.

check_prime_dilation.py checks four exact rational examples, including
uniform and nonuniform blocks and both odd prime classes modulo 4:
393 selected-pole congruences, 3,708 other-pole bounds, every beta
coefficient and constant congruence, and rational-function spot checks.
These small examples check formulas, not the general irrationality proof.

beta_determinants.py contains only the two helper functions used by these
programs. Their definitions are unchanged from the original research
module; the separate determinant experiments are not part of this paper.

Independent audit calculations

  python3 -B checks/audit/fresh_seeds.py
  g++ -O2 -std=c++17 checks/audit/direct_estars.cpp -o checks/audit/direct_estars
  python3 -B checks/audit/fresh_budget.py
  python3 -B checks/audit/fresh_rational.py

The first and last commands require python-flint. The middle two use a
C++17 compiler and the Python standard library. The compiled executable
is a local build artifact and is Git-ignored.

These programs reproduce the separate auditor's calculations. They do
not import the manuscript's checkers. fresh_seeds.py compares its output
with the supplied JSON records after reconstructing the global
polynomials independently. The audit reports six seed constants (one
additional control), 1,644 selected-pole congruences, 24,552 other-pole
bounds and 12 constant congruences, as well as the arithmetic certificate.
The corresponding JSON files preserve the original audit outputs.

Reruns may change timings and output paths. Compare the mathematical
fields, exact bounds, check counts and residues. The assembly review
records the portable replay; rerunning an auditor's code is a self-check,
not an additional independent audit.
