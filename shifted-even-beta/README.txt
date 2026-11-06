IRRATIONAL EVEN BETA VALUES IN SHIFTED INTERVALS

This separate 17-page paper proves that, for every positive odd m,
at least one of β(m+1),β(m+3),...,β(67m-1) is irrational.

More generally, for every B≥1 it gives at least B irrational values
from β(m+1) through β(m(s_B+1)-1), where

  s_B=2 ceil[D_B(33+(45/2)log B)],  D_B=Σ_{d≤B}φ(d).

For each fixed count B, the upper endpoint grows linearly with the
starting index. Discarding the early values from the previous uniform
count would give an endpoint of order m² log m. Section 1.1 explains
this quadratic-versus-linear comparison and gives examples starting
at β(4), β(6), β(20) and β(100). It distinguishes Sorokin and Huttner's
direct shifted-block bounds of order m log m, and compares Krisman's
specific interval β(4) through β(36).

The parameter criterion gives these sharper examples:

  β(4),...,β(198): at least one irrational value;
  β(6),...,β(330): at least one irrational value;
  β(8),...,β(462): at least one;
  β(20),...,β(3080): at least two;
  β(20),...,β(6950): at least three.

Canonical files

  paper/shifted-even-beta.tex
  output/pdf/shifted-even-beta.pdf

The proof uses even derivatives to remove low-index values and repeated
numerator zeros to retain the grouped denominator cancellation. The
differentiated seed kernels form a positive Gram matrix. Prime dilation
propagates the prescribed seed constant and makes the final rational
combination nonzero infinitely often. After clearing denominators it
would give nonzero integers tending to zero, contradicting the assumption
that the interval contains fewer than B irrational values.

The argument is self-contained apart from standard classical inputs:
Stirling's formula, the prime number theorem, Wilson's and Fermat's
theorems, the residue theorem, Cauchy's estimate and elementary linear
algebra and real calculus. The cited papers provide context and method
attribution; no proposition from them is imported into the proof.

Reproduce from this directory

  python3 -B checks/checks.py
  python3 -B checks/endpoints.py
  python3 -B tools/build_paper.py

The checkers use only the Python standard library. Run without -O:
assertions perform the checks. The build uses pdflatex and Poppler's
pdfinfo. Intermediate files are in the Git-ignored build/ directory.

Review scope

An independent audit supplied by the user on 5 November 2026 found no
substantive mathematical error or gap in the complete underlying
plain-text proof. The auditor checked arbitrary allowed parameters and
independently reproduced the finite calculations without running or
importing the supplied checkers. An additional dilation at seed index
l=3, prime 17 and n=51 also passed.

The audit identifies its input by SHA-256 prefix c7197fc5cd341c3c.
review/independent-audit-2026-11-05.json records the full hash and scope;
the accompanying text preserves the user-supplied report with local
filesystem references normalized. Auditor identity and model were not
specified. This is not formal verification or certification of priority.

The PDF is a subsequent exposition of that audited argument. The sole
requested wording correction, naming q_w(1) in the prime exclusion,
is included. review/assembly-review.json records the source and PDF
hashes, correspondence with the audited proof, reproduced checks and
visual inspection. The mathematical and exposition self-reviews are
recorded separately. The abstract contains results only, the title has
no manual line break, and subsections are numbered.

The beta(4) through beta(198) example was added on 6 November and
checked by exact rational arithmetic using the audited general criterion.
It was not among the four examples in the 5 November audit.

sources/reading-record.json records versions, hashes and reading scopes
for the cited sources. Third-party PDFs are not redistributed.

Prepared with Codex under the research pseudonym twinheterodox.
