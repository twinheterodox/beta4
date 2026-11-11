Further investigations of irrationality constructions for Dirichlet beta values
10 November 2026

Manuscript
  output/pdf/further-investigations.pdf
  paper/further-investigations.tex

This is a separate technical supplement to the earlier investigation
supplement. It covers the 29 later records from 6--8 November 2026
concerning individual beta values, small beta sets, and rational relations
among 1, pi-squared and Catalan's constant. The existing count, shifted-
interval, linear-independence and obstruction papers remain separate.

Contents
  Section 1: The integer contradiction and complete primitive normalization.
  Section 2: Fixed polynomial isolation, prime dilation and telescopers.
  Section 3: Real Pade remainders, exact decay rates and an arithmetic gap.
  Section 4: Mixed moments, uniform denominator savings and nonvanishing.
  Section 5: The audited canonical decay estimate and its positive margin.
  Section 6: Attainable determinant values, row lattices, output rounding,
             index losses and the limitations of polynomial kernel inverses.
  Section 7: Small beta sets, rational boundary costs and finite comparisons.
  Section 8: Digamma nonvanishing and other Pade variants.
  Appendix A: Complete catalog of the 29 research records.
  Appendix B: Reproduction and independent-audit scope.

The canonical decay exceeds the available arithmetic cost before the
required integer corrections. Their uniform real effect remains open.
The supplement establishes no new individual irrationality or independence
theorem. Finite unfavorable values are kept separate from general
obstruction theorems; a proved clearing-factor cost is kept separate from
a lower bound on primitive coefficient growth.

Build
  From this directory:

    python3 tools/build.py

  Requires pdflatex and pdfinfo (Poppler). The script makes three LaTeX
  passes, rejects unresolved references and overfull/underfull boxes,
  writes the canonical PDF above, and records hashes in review/build.json.
  Intermediate files go to build/, which is ignored by Git.

Checks
  Use Python with the packages in checks/requirements.txt installed.
  From this directory, without Python's -O option:

    python3 -B checks/replay.py

  Or select one suite, for example:

    python3 -B checks/replay.py --suite canonical

  The runner verifies the evidence manifest, copies selected inputs and
  scripts into a fresh build/checks/ directory, then runs there. It never
  overwrites the retained evidence. No network or private archive access
  is needed. See checks/README.txt for coverage and the distinction between
  reconstruction from definitions and checking stored rational input data.

Evidence and review
  review/coverage.json                All 29 records, outcomes and proof hashes.
  evidence/manifest.json              Selected input/script provenance and hashes.
  sources/ledger.json                Cited versions, PDF hashes and reading scope.
  review/independent-audit-summary.txt Scope of the 8 November independent audit.
  review/independent-audit-2026-11-10.txt  Audit of the 41-page assembly and response.
  review/verification.json            Current computational replay summary.
  review/assembly-review.txt          Original assembly and mathematical-scope review.
  review/exposition-revision-2026-11-10.txt  Subsequent exposition and layout review.

The 8 November audit covers the underlying canonical estimate and its
specified dependencies. A second audit of the 41-page assembly, reported
on 10 November, confirmed that estimate and identified a missing uniform
integrality hypothesis in Section 6.5. The revision states it and cites
the arithmetic results that supply it for both admissible correction
classes. The two exposition corrections were also applied. The audit
report distinguishes the independently reviewed PDF from these subsequent
edits and records the checks replayed independently. The original private
research records were preserved during preparation.
