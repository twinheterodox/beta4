EXACT CHECKS FOR SHIFTED EVEN BETA INTERVALS

Run from the paper directory:

  python3 -B checks/checks.py
  python3 -B checks/endpoints.py

Python's standard library is sufficient. Do not use -O; assertions
perform the checks. The outputs checks.json and endpoints.json are
written beside the programs.

The code parameter a is related to the manuscript's m by m=2a-1;
the first even beta index is 2a=m+1. Parameters k=m-1, S and r
otherwise have the same meaning as in the paper.

checks.py verifies differentiated kernel evaluations over Q, Gram
entries through partial fractions and rational moments, positive
leading principal minors, and integer adjugate seed identities.
Full rational computations test the selected-pole coefficient
congruences, other-pole bounds, rational constants and grouped
integrality at prime indices 17 and 19. The small parameters used
for these identities are not claimed to satisfy analytic decay.

endpoints.py verifies the exact rational f and U comparisons for
the five displayed intervals and a direct all-m certificate at
B=1, S=67m. All printed decimal values are illustrative.

checks.py and its results match the supporting files accompanying the
audited plain-text proof byte for byte. On 6 November, endpoints.py
was extended by the case beta(4),...,beta(198); its four original cases
and the all-m certificate are unchanged. All five cases pass the exact
rational comparisons. The new case was checked by the constructing
assistant and was not part of the 5 November independent audit.
The independent auditor reports using independently written code,
without running or importing these programs. That code was not
provided and is not represented as part of this package.
