# Main manuscript

[Read the paper](../output/pdf/beta4-determinant-obstructions.pdf) · [LaTeX source](beta4-determinant-obstructions.tex)

This 16-page draft of *Arithmetic obstructions in positive determinant approximations to β(4)* gives the two selected obstruction proofs in one notation, including their analytic representations, real lower bounds, exact local arithmetic and scope. The common-zero extension is an appendix.

The accessibility revision explains the irrationality criterion before introducing the matrices, gives a shared three-step proof outline, compares the real decay and arithmetic cost in a table, and works through a three-node example at the prime seven. Short introductions explain why each construction and estimate is needed. The existing theorem statements and proofs are unchanged; see the [editorial review record](../planning/accessibility-review.txt).

The authorship line is “Research note prepared with Codex.” The original OpenAI determinant work is credited as motivation; other literature is cited where it provides context or a standard framework. The specific arguments used in the two proofs are derived in the paper. Background inputs are stated explicitly, including the prime number theorem.

Build from the repository root:

```sh
python3 tools/build_paper.py
```

This requires `pdflatex` with the listed standard LaTeX packages, and `pdfinfo` for the build record. The script places the PDF under `output/pdf/` and rejects unresolved references and box warnings. The [audit record](audit.txt) distinguishes the mathematical review, computational checks and PDF inspection.

The [broader technical supplement](../supplement/README.md) also has a complete draft, organized from the existing [claim inventory](../planning/theorem-inventory.txt) and [method catalog](../planning/method-catalog.txt). The main-paper assembly audit records the earlier writing pass; its reference to the then-pending supplement is historical. The original research snapshots remain unchanged. Neither assembled manuscript has been published or submitted for a new independent cross-check.
