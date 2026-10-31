# Main manuscript

[Read the paper](../output/pdf/beta4-determinant-obstructions.pdf) · [LaTeX source](beta4-determinant-obstructions.tex)

This 19-page draft of *Arithmetic obstructions in positive determinant approximations to β(4)* gives the two selected obstruction proofs in one notation, including their analytic representations, real lower bounds, exact local arithmetic and scope. The common-zero extension is an appendix.

The introduction starts with the integer-arithmetic contradiction behind an irrationality proof. A two-by-two example explains how determinants produce polynomials, and a simple denominator example explains primitive normalization. The two main theorems then appear with plain-language interpretations, before the detailed matrix formulas. A shared three-step proof outline and an exponent table explain why the arithmetic cost defeats the real decay.

Sections 3–4 treat the first construction and its arithmetic obstruction; Sections 5–6 do the same for the second. Section 3.1 derives the first matrix from the shifted cubic tail, while Section 5.1 introduces the second through its positive double integral. A three-node example at the prime seven illustrates the first arithmetic argument. The exact large-prime formula formerly included in Theorem 1.2 is now Proposition 6.1, immediately before its proof.

The [current editorial review](../planning/narrative-review.txt) records this restructuring and the preservation of all nine existing proof blocks. The [earlier equation-by-equation review](../planning/expository-review.txt) records the preceding expansion of the kernel, Cauchy, moment and prime-residue calculations; the [accessibility review](../planning/accessibility-review.txt) describes an earlier revision. Those records are historical, not new mathematical audits.

The authorship line is “Research note prepared with Codex.” The original OpenAI determinant work is credited as motivation; other literature is cited where it provides context or a standard framework. The specific arguments used in the two proofs are derived in the paper. Background inputs are stated explicitly, including the prime number theorem.

Build from the repository root:

```sh
python3 tools/build_paper.py
```

This requires `pdflatex` with the listed standard LaTeX packages, and `pdfinfo` for the build record. The script places the PDF under `output/pdf/` and rejects unresolved references and box warnings. The [audit record](audit.txt) distinguishes the mathematical review, computational checks and PDF inspection.

The [broader technical supplement](../supplement/README.md) also has a complete draft, organized from the existing [claim inventory](../planning/theorem-inventory.txt) and [method catalog](../planning/method-catalog.txt). The main-paper assembly audit records the earlier writing pass; its reference to the then-pending supplement is historical. The original research snapshots remain unchanged. Neither assembled manuscript has been published or submitted for a new independent cross-check.
