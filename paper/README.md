# Main manuscript

[Read the paper](../output/pdf/beta4-determinant-obstructions.pdf) · [LaTeX source](beta4-determinant-obstructions.tex)

This 18-page draft of *Arithmetic obstructions in positive determinant approximations to β(4)* gives the two selected obstruction proofs in one notation, including their analytic representations, real lower bounds, exact local arithmetic and scope. The common-zero extension is an appendix.

The introduction explains the irrationality criterion before introducing the matrices, gives a shared three-step proof outline, and compares the real decay and arithmetic cost in a table. It derives the first matrix from the shifted cubic tail and introduces the second through its positive double integral. A three-node example at the prime seven illustrates the arithmetic argument.

The [equation-by-equation expository review](../planning/expository-review.txt) adds intermediate calculations for the integral kernel, Cauchy comparison, moment coefficients and prime residues. It also explains the offset regimes and replaces a compressed determinant argument with an explicit block matrix obtained by scaling rows. The theorem statements are unchanged; their proofs have been expanded. The [earlier accessibility review](../planning/accessibility-review.txt) records the preceding revision.

The authorship line is “Research note prepared with Codex.” The original OpenAI determinant work is credited as motivation; other literature is cited where it provides context or a standard framework. The specific arguments used in the two proofs are derived in the paper. Background inputs are stated explicitly, including the prime number theorem.

Build from the repository root:

```sh
python3 tools/build_paper.py
```

This requires `pdflatex` with the listed standard LaTeX packages, and `pdfinfo` for the build record. The script places the PDF under `output/pdf/` and rejects unresolved references and box warnings. The [audit record](audit.txt) distinguishes the mathematical review, computational checks and PDF inspection.

The [broader technical supplement](../supplement/README.md) also has a complete draft, organized from the existing [claim inventory](../planning/theorem-inventory.txt) and [method catalog](../planning/method-catalog.txt). The main-paper assembly audit records the earlier writing pass; its reference to the then-pending supplement is historical. The original research snapshots remain unchanged. Neither assembled manuscript has been published or submitted for a new independent cross-check.
