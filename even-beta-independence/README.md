# Linear independence of even Dirichlet beta values

**[Read the paper (PDF, 19 pages)](output/pdf/even-beta-independence.pdf)** · [LaTeX source](paper/even-beta-independence.tex) · [Exact supporting checks](checks/README.txt)

For every sufficiently large even integer $s$, the paper proves

$$
\dim_{\mathbb Q}\operatorname{span}_{\mathbb Q}\{\beta(2),\beta(4),\ldots,\beta(s)\}
\geq 0.21\sqrt{\frac{s}{\log s}}.
$$

Here $\beta(s)=1-3^{-s}+5^{-s}-\cdots$. The bound measures how many values in a growing even-beta family are linearly independent over the rationals. It does not identify any particular value, including $\beta(4)$, as irrational. No explicit threshold for “sufficiently large” is asserted.

The proof adapts Fischler's rational-function construction to sums over half-integers. Reflection removes the odd beta values, while the resulting differential system is regular at the evaluation point $-1$. A local denominator estimate controls the integer coefficients. A quadratic change of variable permits the multiplicity argument, counting vanishing at both preimages of the evaluation point. The paper states its three imported results precisely and gives the adaptation, estimates and rank proof in full.

## Review and supporting records

An [independent model proof audit](review/independent-model-review-2026-11-03.json), reported on 3 November 2026, found no substantive mathematical gap in the current 19-page revision. The [full report](review/independent-model-review-2026-11-03.txt) covers the denominator and rank arguments, the application of Fischler's multiplicity estimate, the parameter constraints and the final bound. The reviewer reported that all supplied exact checks passed, together with separate checks through Taylor degree 60 and depth 12 and numerical checks of the tail identities. The audit uses the cited results and is not formal verification. Priority has not been established.

The exposition now includes a two-integration denominator example, a worked common-row example with h=1, an explanation of all six points in the rank count, and the square-root balance behind the parameter choice. The [final exposition review](review/final-expository-review-2026-11-03.txt) recommends no further substantive changes. The exposition revisions used OpenAI Codex with the author-supplied model label Astra-6-Max; the independent reviewing model was not named.

- [Developing assistant's audit](review/audit.txt) and [original review questions](review/review-request.txt).
- [Earlier independent model cross-check](review/independent-model-review-2026-11-02.json), applying to its recorded revision.
- [Final exposition revision](review/human-reader-exposition-revision-2026-11-03-round2.json), [review-status update](review/verification-status-update-2026-11-03.json) and [repository consolidation](review/consolidation-2026-11-03.json).
- [Exact-check coverage and limitations](checks/README.txt), [earlier publication replay](checks/reproduction.json), [local-system results](checks/local-system-results.json) and [linear-form results](checks/half-shift-forms-checks.json).
- [Source versions and reading scopes](sources/reading-record.json). Archive-relative paths identify the original local copies; the papers themselves are not distributed. Source URLs and checksums are supplied.
- [Build record](paper/build-record.json) and [visual review](review/visual-review.json).

The scripts check finite identities, denominator calculations and the numerical constant using exact arithmetic. They support the written proof; they do not establish the asymptotic theorem by checking finitely many cases.

## Reproduce the checks and build

From this directory, install the pinned Python dependency and run both scripts:

```sh
python3 -m pip install -r checks/requirements.txt
python3 checks/check_local_system.py
python3 checks/check_linear_forms.py
```

The scripts write their JSON reports in `checks/`. They use Python and SymPy; no numerical approximation to a beta value is needed.

With LaTeX and Poppler (`pdflatex` and `pdfinfo`) installed:

```sh
python3 tools/build_paper.py
```

The build checks references and overfull boxes and writes the canonical PDF to `output/pdf/even-beta-independence.pdf`. This directory is the single working manuscript: edit the source and rebuild here, with earlier revisions preserved in Git. Files in `build/latex/` are generated intermediates. The package needs no files from the original research workspace.

Prepared with Codex. Original material is covered by the repository's [Apache 2.0 license](../LICENSE).
