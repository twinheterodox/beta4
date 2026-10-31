# Main-paper calculations

The paper's asymptotic results are proved in the text. These scripts check finite instances of its identities, coefficient valuations, primitive normalization and real bounds.

Install the single pinned dependency into a Python environment, then run:

```sh
python -m pip install -r checks/requirements.txt
python checks/reproduce.py
```

The recorded run used Python 3.14.4 and python-flint 0.9.0. The runner copies each archived script into a temporary directory before executing it. It writes fresh results to `checks/results/` and records their hashes and the environment in [reproduction.json](reproduction.json). It verifies that all 120 archived snapshots remain unchanged.

The rerun covered:

- 21 regularized determinant identities, 25 principal frames, 130 coefficient-prime checks and 9,572 entry-prime checks; positivity and real bounds used Arb at 4,096-bit precision.
- 144 further exact Cauchy determinant identities, including common-zero frames.
- Five Arb checks of the appendix's simplified common-zero upper and lower bounds, using the same regularization parameter in both.
- 23 difference-variable moment identities, 36 direct rational determinant evaluations, 1,100 residue checks and 74 exact rank/primitive-exponent checks, for matrix sizes 1 through 12. Primitive evaluations and real bounds were also checked with Arb.
- The manuscript's Rodrigues calculation: 13 monic polynomials, 78 orthogonality equalities, 13 exact squared norms and 9 direct Gram determinant products, computed with Python rational arithmetic.

The exact structural records match the archived results. Numerical interval strings are not required to match byte-for-byte; their inequalities are checked afresh. The unrelated slow-filter branch diagnostics are not part of this paper and are not rerun.

These are model-authored calculations and a replay of existing checks, with an additional rational-polynomial implementation for the norm formulas. They are not a fresh independent model audit, a human audit or formal verification.
