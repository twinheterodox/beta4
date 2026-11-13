Linear forms, moment determinants, and arithmetic constraints
for even beta values
============================================================
twinheterodox, 13 November 2026

The manuscript is output/pdf/arithmetic-investigations.pdf.
Its source is paper/arithmetic-investigations.tex.

Purpose and reading order
-------------------------
The supplement records later attempts to reduce the known fifteen-value
beta set or to prove beta(4) irrational. Its main exact result determines
the power of two in every principal minor of a beta-coefficient matrix.
That necessary denominator cost excludes one precisely specified
rational-approximation certificate.

Sections 2-4 give that result and its application as a self-contained
route through the paper. Sections 5-6 develop controlled elimination of
zeta(4) and uniform determinant-content savings. Section 7 removes a
rational Gram block and proves the linear scale of the residual decay.
Section 8 proves formal rank, classifies zero multipliers, and constructs
eleven independent small integer forms in fourteen candidate beta values.
Section 9 records two further tests and distinguishes the remaining gaps.

The supplement proves no new irrationality theorem. Its obstruction is
to a globally monic-cleared certificate using the generic q^-2 error
envelope; it is not an obstruction to every approximant or moment method.
The other unsuccessful comparisons retain their stated narrower scopes.

Build
-----
From this directory:

    python3 tools/build.py

Requirements: Python 3, pdflatex with the packages named in the source,
and pdfinfo (Poppler). Intermediate TeX files stay in build/. The build
record in review/build.json identifies the PDF, source hashes and warnings.

Checks and evidence
-------------------
See checks/README.txt for commands, dependencies and precise coverage.
The scripts recompute the exact valuation identities, the fixed finite
arithmetic partitions for the eleven-form theorem, the six recorded Gram
cases, and the one recorded Loewner screen. They do not replace the
general proofs.

evidence/coverage.json maps the claims to supporting calculations and
states which numerical observations are retained historical results.
sources/ledger.json identifies cited versions, checksums and reading
scope. Source papers are not redistributed.

Review
------
The mathematical investigations were developed and cross-checked in
cooperating Codex sessions. The manuscript passed a separate
cooperating-session correctness and exposition review; its completed
scope is recorded in review/cooperating-review.txt. Finite replays, independent
implementations and mathematical reviews are distinguished there.
None is described as formal verification.
