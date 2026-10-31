# Companion supplement

[Read the supplement](../output/pdf/beta4-investigations.pdf) · [LaTeX entry point](beta4-investigations.tex) · [Audit record](audit.txt)

This 31-page draft expands the mathematical record behind the [focused determinant paper](../paper/README.md). It contains complete proofs of selected results and a guide to all twelve method groups and 57 grouped claim records.

The revised abstract and introduction explain the irrationality criterion, the denominator problem and the distinction between proved failures and unresolved constructions. A table guides readers to the proved results, and each construction now begins with its motivation and mechanism. Specialized terms are introduced where they are used. The [accessibility revision record](../planning/accessibility-review.txt) documents these changes and the final layout review; existing statements and proofs are unchanged. The assembly audit above records the earlier draft.

The rewritten material covers the positive-polynomial cone and degree law; repeated-quadratic sources; Euler–Boole truncations; fixed-endpoint differentiated continued fractions; a sine-normalized cyclotomic family; a constrained rational Mellin kernel; a scoped monodromy obstruction; and the positivity and sharp scalar rate of the ordinary-moment first variation. Two source-level reductions are also derived. The completed Mellin kernel's identification with β(4) uses the explicitly declared Broadhurst identity; the kernel inequality itself is proved independently of that identification.

The other constructions are retained with their precise archived status. A source-register entry does not certify that its full proof has been re-audited in this writing pass. In particular, finite failures, an unsuccessful sufficient multiplier, and a theorem about actual primitive values are different conclusions.

From the repository root:

```sh
python3 tools/render_supplement_register.py
python3 tools/build_supplement.py
python3 checks/reproduce_supplement.py
```

The build needs LaTeX and Poppler. The check runner needs the versions in [supplement-requirements.txt](../checks/supplement-requirements.txt). It executes preserved scripts in temporary copies and checks all 127 frozen evidence hashes. Its [report](../checks/supplement-reproduction.json) separates exact calculations from uncertified numerical diagnostics. The [coverage record](coverage.json) states the scope of each claim's treatment; the [source ledger](../sources/supplement-reading.json) records actual reading and external dependencies.

This is a local, unpublished draft prepared with Codex. Fresh independent mathematical review and further priority scrutiny remain outstanding. No new irrationality research campaign or parameter search was started for the writeup.
