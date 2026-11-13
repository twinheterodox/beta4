Reproducing the supporting calculations
=======================================
Run from the manuscript directory. Python 3.11 or later is recommended.
monic.py uses only the standard library. The other scripts use
python-flint (FLINT exact arithmetic and Arb balls); the version used
for assembly is recorded in requirements.txt.

    python3 checks/monic.py --output review/monic-checks.json
    python3 checks/linear.py --output review/linear-checks.json
    python3 checks/gram.py --output review/gram-checks.json
    python3 checks/loewner.py

monic.py
--------
Recomputes the finite binomial reformulation, the generating-function
orthogonality, the finite multiplication matrix and its parity, the
all-subset valuation formula at the recorded small sizes, complete small
monic polynomials, the full-coordinate determinant identity, and the
rational positivity bound. The 2,980 counted exact checks are inherited
from the audited checker, with its calculation unchanged.

The all-K valuation and uniform analytic estimates are proved in the
manuscript. A passing finite computation is not a proof of either limit.

linear.py
---------
Reconstructs the entire endpoint partition on [1,118] (13,091 cells) and
the local Taylor-functional partition on [1,6] (582 cells). Integer
patterns and Taylor recursion independently recover the 464 template
identities and the integral 62795/199056. The checker also encloses the
real phase maximum and recomputes the constants H_* and C_*.

The input JSON files retain exact rational endpoints and certificate
data. Their original generators are not imported. The checker does not
replay the large n = 116, 232, 464 cancellation bases, and it does not
prove nonvanishing at the actual beta values.

gram.py
-------
Reconstructs the rational blocks and complete determinant polynomial for
the recorded cases r = 1,...,6. It removes the gcd of ALL coefficients
and verifies the saved primitive multiplier, both isolated coefficients,
the rational block-clearing assertion and three exact specializations.
Arb then checks the positive primitive values used in the table.

No saved large polynomial is needed: it is regenerated from the seven
by seven Schur matrix. The scalar control gives 82 beta(4) - 81. A
separate earlier Fraction/dictionary implementation checked every
coefficient for r = 2; the present replay is not a new external audit.
These cases make no assertion about asymptotic primitive growth.

loewner.py
----------
Replays the declared geometry r = 200k, m = 13k, its 200 monotonic
potential cells, tail and energy bounds. It writes
review/loewner-checks.json. The incompatible coarse rates remain
INCONCLUSIVE, not a failure theorem for the construction. This
calculation does not independently prove the moving-contour and
all-deletion estimates; those belong to the reviewed general argument.

Provenance
----------
The scripts are portable adaptations of the recorded calculations:
paths, output destinations and explanatory headers were changed;
private session records and environment details were not copied.
The evidence guide gives the distinction between general proofs,
recomputed certificates, retained observations and earlier reviews.
