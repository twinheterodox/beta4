IRRATIONAL EVEN BETA VALUES IN SHIFTED INTERVALS
Version 2, 6 November 2026

For every positive odd m, this paper proves that at least one of

  beta(m+1), beta(m+3), ..., beta(19m-1)

is irrational, where beta is the Dirichlet beta function.
This improves the v1 endpoint 67m-1.

More generally, for every B>=1 it gives at least B irrational values
from beta(m+1) through beta(C_B m-1), where

  D_B = sum of phi(d) over odd d with 1<=d<=2B-1,
  C_B = 1+2 ceil[D_B(9+(40/9)log(2B-1))].

The introduction reproduces the v1 general theorem immediately before
v2, using matching endpoint notation. For B=1,2,3 the multipliers
improve from 67,197,463 to 19,85,229. In the proof, a_B=C_B-1
sets the pole order a=m a_B.

For each count B, the upper endpoint grows linearly with the starting
index. Discarding the early values from the previous uniform count
would give an endpoint of order m^2 log m. Section 1.1 gives concrete
comparisons with that count, Sorokin and Huttner's shifted bounds, and
Krisman's shorter interval beta(4) through beta(36).

Individual parameter choices give, for example:

  beta(4),...,beta(54): at least one irrational value;
  beta(6),...,beta(92): at least one;
  beta(8),...,beta(128): at least one;
  beta(20),...,beta(352): at least one;
  beta(20),...,beta(1344): at least two;
  beta(20),...,beta(3538): at least three;
  beta(100),...,beta(1838): at least one.

Starting at beta(2), the construction gives counts 2, 3, 4, 5 and 10
through indices 70, 186, 360, 554 and 2794, respectively.

Canonical files

  paper/shifted-even-beta.tex
  output/pdf/shifted-even-beta.pdf

The proof uses odd numerator grids and poles spaced by 2. Even
derivatives remove the low-index beta values; repeated numerator zeros
retain the grouped denominator cancellation. A positive Gram matrix
supplies an integer polynomial multiplier with a prescribed seed vector.
Prime dilation propagates its nonzero last constant. If the interval
contained fewer than B irrational values, a fixed combination of the
form sequences would give infinitely many nonzero integers tending
to zero.

Section 3.1 works through the two-grid example giving two irrational
values through beta(70). Section 6.2 spells out the local coefficient
scaling. Appendix B retains the sharper uniform endpoint
s(m)=m-1+2 ceil[(879m-50)/100]<18.58m and explains its multiplier
and the choice r=2.

The argument is self-contained apart from standard classical inputs:
Stirling's formula, the prime number theorem, Wilson's and Fermat's
theorems, the residue theorem, Cauchy's estimate and elementary linear
algebra and real calculus. The cited papers provide context and method
attribution; no proposition from them is imported into the proof.

Reproduce from this directory

  python3 -B checks/checks.py
  python3 -B checks/endpoints.py
  python3 -B checks/worked_example.py
  python3 -B tools/build_paper.py

For the numerical transition in Appendix B:

  python3 -B checks/optimization.py

The checkers use only the Python standard library. Run without -O:
assertions perform the checks. The build uses pdflatex and Poppler's
pdfinfo. Intermediate files are in the Git-ignored build/ directory.

Review scope

An independent audit supplied on 6 November 2026 found no substantive
gap in the v2 argument. It covered the general proof and independently
reproduced all 8,909 recorded structural checks, all 21 endpoint
certificates, six uniform parameter records and seven integer-r
comparisons. The auditor inspected the supplied checkers for coverage
but did not run or import them.

The audit identifies its input by SHA-256 prefix 603211b50609dc4a.
review/independent-audit-2026-11-06.json records the full hash and scope;
the accompanying text preserves the supplied report with local-file
references normalized. Auditor identity/model and a separate audit
date were not specified. The audit is not formal verification or a
certification of priority.

The PDF is a subsequent exposition of the audited argument.
review/assembly-review.json records the source correspondence, exact
replays and visual inspection. The worked example has a separate
post-audit exact check. Historical v1 reviews retain their original
scope; see review/README.txt.
