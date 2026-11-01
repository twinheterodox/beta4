# Companion supplement

[Read the supplement](../output/pdf/beta4-investigations.pdf) · [LaTeX entry point](beta4-investigations.tex) · [Current revision review](../planning/supplement-self-contained-review.txt)

This 67-page draft is a mathematical starting point for continuing the investigation behind the [focused determinant paper](../paper/README.md). Its substantial constructions can now be reconstructed from the text and accompanying certificate files, without recovering definitions or missing proof steps from research logs.

The introduction explains the irrationality criterion and the denominator problem, then gives separate guides to the obstruction proofs and unfinished constructions. The case studies cover interleaved determinants, squared and unsquared quarter poles, the common-base polylogarithm ladder, endpoint interpolation, prime deletion, modular transformations, coupled period reductions, rationality-dependent sources, paired hypergeometric identities, compatible pole multipliers and slow logarithmic dilations. Each explains its motivation, supplies explicit formulas and established arguments, and identifies the remaining estimate or nonvanishing question.

The original complete proofs remain: the positive-polynomial cone and degree law; repeated quadratics; Euler–Boole corrections; fixed-endpoint differentiated continued fractions; the sine-normalized cyclotomic family; the constrained Mellin kernel; a scoped monodromy obstruction; and ordinary-moment first variation. Shorter exploratory variants and bounded experiments remain in a separate catalog, with a linked register covering all 57 grouped records.

The scope of “self-contained” is mathematical exposition with identified inputs. Classical results are cited where used, the interleaved case states the kernel comparison imported from the companion paper, and the ladder uses the declared Broadhurst identity. The paired case states its asymptotic recurrence input and includes a [standalone exact certificate replay](../checks/paired/README.txt) for its large algebraic identities. None of these featured cases requires the original research environment. The [coverage record](coverage.json) distinguishes the expanded arguments from variants retained only with their archived scope.

From the repository root, with the [check dependencies](../checks/supplement-requirements.txt) installed:

```sh
python3 tools/render_supplement_register.py
python3 tools/build_supplement.py
python3 checks/reproduce_supplement.py
python3 checks/case_studies.py
python3 checks/paired/verify.py
python3 checks/paired/arithmetic.py
python3 checks/paired/rates.py
```

The build needs LaTeX and Poppler. The existing [reproduction report](../checks/supplement-reproduction.json) distinguishes exact calculations from numerical diagnostics and checks all 127 frozen evidence hashes. The new [case-study checks](../checks/case-studies-result.json) reconstruct formulas independently of the archived calculation scripts. The paired replay checks symbolic all-index identities; its finite initial cases and numerical rate enclosures have separate roles, explained in the manuscript and replay guide.

The [source ledger](../sources/supplement-reading.json) records versions, local PDFs, checksums, actual reading scope and external dependencies. The [original assembly audit](audit.txt) and earlier editorial records describe their respective historical drafts.

This is a local, unpublished draft prepared with Codex. The revision includes source review and exact computational checks by its assembling model, not a fresh independent audit or formal verification. It does not settle the irrationality of β(4).
