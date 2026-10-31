# Supplement calculations

Install the pinned dependencies into a Python environment, then run from the repository root:

```sh
python -m pip install -r checks/supplement-requirements.txt
python checks/reproduce_supplement.py
```

The recorded run used Python 3.14.4, SymPy 1.14.0 and mpmath 1.3.0. It passed three replays, matching the complete archived JSON outputs:

- The Mellin certificate: exact rational kernel division, Sturm sequence, unique-root bracket, 32-term rational polylogarithmic tail enclosures and the lower bound 1309/2600. This finite certificate is an explicit input to the supplement's kernel inequality.
- Euler–Boole: 160 Euler coefficients by three implementations, 13 symbolic sech coefficients, eight finite geometric identities, 652 exact prime valuations, 43 reduced fractions and four primitive-product checks.
- Cyclotomic normalization: exact algebraic integers and primitive minimal polynomials at primes 3, 5, 7 and 11. The accompanying 100-digit all-conjugate diagnostics at primes 11, 31 and 101 are **not interval certificates**.

The new rational-arithmetic checks reproduce 28 quadratic-source reversal identities, 57 surviving prime valuations, 36 cone coefficient pairs, 80 continuant cases and 203 exact reduced prime exponents. They also independently implement rational Schur solves to reconstruct the **same 24** archived ordinary-moment primitive pairs and verify twelve positive-integral formal identities. Here “independently implement” means a different calculation in this model's writing pass, not another model or human auditor. No new ordinary-moment index or numerical interval was evaluated.

All scripts run away from the frozen evidence. The runner verifies the original 120 snapshots and seven additional source files before and after execution. Results and hashes are in [supplement-reproduction.json](supplement-reproduction.json). These checks do not prove asymptotic statements by sampling or establish priority.
