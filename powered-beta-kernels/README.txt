IRRATIONAL EVEN BETA VALUES FROM POWERED HYPERGEOMETRIC KERNELS
11 November 2026

Each of the following sets contains an irrational number:

  beta(6), beta(8), ..., beta(74);
  beta(8), beta(10), ..., beta(110).

The earlier shifted-interval paper gives endpoints 92 and 128 for these
starting indices. The new intervals contain 35 and 52 candidate values,
respectively; each removes nine values from the earlier interval.
The two assertions need not identify distinct irrational numbers.

The proof raises Krisman's hypergeometric kernel to a fixed power and
differentiates to remove the early beta values. The denominator cost and
real decay scale together. A prime-dilation argument and a finite nonzero
seed then give infinitely many nonzero integers tending to zero under
the assumption that all values in the chosen interval are rational.

Canonical files

  paper/powered-beta-kernels.tex
  output/pdf/powered-beta-kernels.pdf

This directory is the canonical working version. The original research
workspace is retained privately as a historical record.

Build from this directory

  python3 -B tools/build_paper.py

The build requires pdfLaTeX with the usual AMS and Latin Modern packages,
and Poppler's pdfinfo. Intermediate files go in the Git-ignored build/
directory. The source has no external figures or bibliography files.

Computational checks

See checks/README.txt. The seed checker and arithmetic-budget checker
use only the Python standard library. The primary seed generator and
rational dilation checker also need python-flint 0.9.0. Two original
helper functions are retained verbatim; unrelated determinant pilots
are not included. The complete finite certificate data are supplied.

The auditor's separate implementations and numerical records are also
available in checks/audit/. They are distinct from the manuscript's
supplied checkers. The auxiliary arithmetic-budget program needs a
C++17 compiler. No executable binary is distributed.

Attribution and inputs

Krisman's kernel and quantitative arithmetic/decay estimates are imported
inputs. The prime-dilation strategy follows the earlier shifted-interval
paper. Rivoal, Lai--Sprang and Zudilin are cited directly for the relevant
methodological precedents; Sorokin gives earlier shifted-interval context.
sources/reading-record.json identifies all six sources, their hashes and
the passages actually read. Third-party PDFs are not redistributed.

Review scope

review/independent-audit-2026-11-11.txt preserves the separate AI audit,
with local paths normalized for public use. Its companion JSON identifies
the audited source and PDF hashes. The audit found no substantive error
or gap and independently reconstructed the finite certificates. It used
fresh implementations after preliminary inspection and replay of the
supplied checkers; it was not a blind review. No Lean project was rebuilt,
and this is not formal verification.

The current manuscript is a subsequent exposition and attribution
revision. review/assembly-review.json records the transfer into this
repository, preservation checks, portable replays and PDF comparison.
The general power criterion remains conditional on a nonzero seed;
the manuscript verifies that condition for powers 1, 2 and 3.

The original material follows the repository's Apache-2.0 license.
Prepared with Codex. Author: twinheterodox.
