EXACT CHECKS FOR SHIFTED EVEN-BETA INTERVALS, VERSION 2

Run from the shifted-even-beta directory:

  python3 -B checks/checks.py
  python3 -B checks/endpoints.py
  python3 -B checks/worked_example.py

These programs use only the Python standard library. Do not use -O:
their assertions perform the checks. Each writes the corresponding
JSON result beside its source.

The parameters now agree directly with the manuscript: m is positive
odd, B is the count, r is the integer zero-range parameter, and a is
the even pole order. The endpoint is a+m-1. The v1 checker's notation
and numerical records belong to the earlier version.

checks.py
  Reproduces 8,909 structural checks with exact rational arithmetic:
    57 differentiated-kernel evaluations;
    49 Gram entries, comparing residues with moments;
    19 positive leading principal minors;
     9 integer seed-vector identities;
    54 rational-function reconstructions;
  1308 reflection identities;
   332 differentiated grid zeros;
  4516 grouped-integrality checks;
   404 selected-pole congruences;
  1940 other-pole valuation bounds;
    19 constant congruences;
   202 beta-coefficient valuation bounds.

  The nine seed cases are (B,m,r):
    (1,1,2), (1,3,2), (1,5,1), (2,1,1), (2,3,1),
    (2,1,2), (3,1,1), (3,3,1), (4,1,1).
  All use a=mDL+1. Complete prime cases are recorded in checks.json.
  These small parameters check identities, not endpoint decay.

endpoints.py
  Reproduces all 21 rational endpoint certificates, including the 12
  displayed in Section 8.2; checks the uniform inequalities for every
  m and B; certifies the six listed numerical values of a_B using
  rational exponential bounds; and proves the seven integer-r
  comparisons in Appendix B.

  Floating point locates trial u and candidate ceilings. Every claimed
  inequality and ceiling is then certified with rational arithmetic.
  U_approx and f_approx are illustrative floating-point displays.
  The general proof supplies the all-parameter arguments.

worked_example.py
  Checks B=2, m=1, r=2, a=70 from Section 3.1. It compares the displayed
  moment matrix with the seed partial fractions, verifies the
  determinant formula, and checks the integer multiplier's seed vector.
  This example was added after the independent audit and is a self-check.

optimization.py
  Optional numerical exploration:
    python3 -B checks/optimization.py

  Reproduces the limiting profile transition for B=1 and scans integer
  r=1,...,15 for B=1,2,3,4,5,10. Its JSON is numerical evidence, not
  an exact certificate. The manuscript's 18.58 bound and integer-r
  comparisons are certified by endpoints.py.

The checkers are supplied self-checks, not an independent audit.
The independent auditor's separately written calculations and general
proof review are recorded in review/independent-audit-2026-11-06.json.
No proof-assistant verification is claimed.
