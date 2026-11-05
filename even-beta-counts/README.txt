EXPLICIT LOWER BOUNDS FOR THE NUMBER OF IRRATIONAL EVEN BETA VALUES

This separate 20-page paper proves, for every even integer s≥160,

  I(s) ≥ ⌊π√(2/135) √(s/log s)⌋,

where I(s) counts irrational values among β(2), β(4), …, β(s).
The constant is π√(2/135)=0.3823824806363651… . The bound gives
two at s=160 and three at s=364. The floor can be removed for
every even s≥42,170.

The underlying finite theorem gives at least B irrational values
through the endpoint

  s_B = 2⌈D_B(33+(45/2)log B)⌉,   D_B = Σ_{d=1}^B φ(d),

for every positive integer B. Exact choices of the construction
parameters also give

  B:  2    3    4     5     6     7     8     9     10
  s: 160  364  588  1030  1296  1992  2502  3260   3828.

Canonical manuscript and PDF

  paper/even-beta-counts.tex
  output/pdf/even-beta-counts.pdf

The paper is self-contained apart from standard classical inputs:
the residue theorem, Rolle's theorem, Stirling's formula, Wilson's
and Fermat's congruences, the prime number theorem, and elementary
totient identities. Its literature references supply context and
attribution; no proposition from those papers is a proof input.

Proof and exposition

The introduction explains how fixed integer weights on B sequences
cancel B−1 possible irrational values. A nonzero constant at one
seed index is propagated to infinitely many prime indices, yielding
the nonzero integers tending to zero that give the contradiction.

The new rank argument represents the seed constants as a Gram matrix
of independent polynomials against a positive weight. A two-grid
example explains the choice of multiplier polynomials. Subsequent
sections give the grouped denominator argument, the exponential and
infinite-tail estimates, and the local power-series justification
of coefficient extraction modulo a prime. The endpoint section
explains both the square-root scale and the factor 2 in its constant.

The abstract contains the result and its finite scope, with no proof
sketch. All subsections are numbered. The title has no manual line
break. The original three-through-364 paper in ../finite-even-beta/
has been preserved unchanged; this is a separate manuscript.

Reproduce from this directory

  python3 -B checks/verify_construction.py
  python3 -B checks/verify_endpoints.py
  python3 -B tools/build_paper.py

The checkers require Python 3.8 or later and only its standard library.
Run without -O, since assertions perform the checks. The build requires
pdflatex with the packages listed in the source, and Poppler's pdfinfo.
Intermediate build and rendering files are in the Git-ignored build/
directory; output/pdf/ contains the single canonical current PDF.

checks/construction-verification.json records exact decay comparisons,
385 rational Gram entries, 55 positive leading principal minors,
integer seed identities, and full rational tests of coefficient and
constant dilation at small parameters.

checks/endpoint-verification.json records 230 unrounded and 23 rounded
endpoint comparisons with rigorous rational logarithm intervals. It
also checks all constants in the infinite-range argument, the floors
at s=160 and 364, and the positive explicit error coefficient. It also
checks the footnote's alternative coefficient 0.358 at 64, 158 and 364.

Review scope

An independent mathematical audit, reported on 5 November 2026,
found no substantive mathematical error or gap. It checked the
general proof and independently reproduced the finite calculations
without running or importing the supplied checkers.

review/independent-audit-2026-11-05.txt preserves the supplied report,
with its local filesystem path normalized. The accompanying JSON
record identifies the 19-page PDF snapshot and source hashes captured
before this revision. The report did not supply a reviewer name,
model or checksum; those details are not inferred.

The present revision adds the requested explanation of the uniform
Stirling error in Section 5.2 and the audit note in Appendix A.3.
review/audit-response-2026-11-05.json records these subsequent edits
and their validation separately. The audit was not formal verification
or a fresh audit of all cited literature.

review/assembly-review.json identifies the exact source and PDF version,
the mathematical self-review, the checks, visual inspection and the
preservation of the earlier paper. review/exposition-review.txt records
the review against mathematical-exposition.txt.

review/comparison-appendix-2026-11-04.json records the earlier revision in which
the introduction states the earlier counting bound, and Appendix B
gives its deduction from the augmented-span estimate. The preceding
version is recorded in review/comparison-justification-2026-11-04.json.
The earlier manuscript and PDF retain their original versions.

The developing assistant's self-review and supplied computational
checks remain separate from the independent audit. Audit records for
the earlier papers retain their original scopes.

sources/reading-record.json records the versions, hashes and actual
reading scopes of the cited sources. Kumar and Tosi's stated asymptotic
count has a larger coefficient; the introduction distinguishes that
statement from the explicit finite range proved here. No new priority
search was performed during this assembly.

Prepared with Codex under the research pseudonym twinheterodox.
