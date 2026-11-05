AT LEAST THREE OF BETA(2), BETA(4), ..., BETA(364) ARE IRRATIONAL.

The paper proves that at least three of beta(2), beta(4), ..., beta(364)
are irrational.

Canonical manuscript:
  paper/finite-even-beta.tex
Canonical PDF:
  output/pdf/finite-even-beta.pdf

There is one editable source and one current PDF in this directory.
The historical plain-text proof, private build history and research
provenance are preserved separately from this checkout.

The current PDF has 16 pages. The introduction explains the
contradiction assumption and why three forms cancel two possible
irrational values. Sections 2 and 3 motivate the construction and grouped
denominator bound. Section 4 includes the explicit infinite-tail estimate
and its exact decay checks. Section 5 motivates prime dilation before its
worked product modulo 5. Section 6 explains why the determinant tests
nonvanishing, gives the test, and then derives the seed residues.
Appendix A compares the finite estimates behind selected asymptotic
results at endpoint 364. It includes Rivoal--Zudilin's uniform family,
Nishimoto's explicit estimate, Fischler's 2018 and 2020 size estimates,
and the full finite-parameter bound in the companion square-root paper.

Proof outline

Three character sums have coefficient columns proportional to (1,2^i,3^i).
An integer combination cancels any two chosen beta values for every n.
A grouped denominator argument and an exact decay bound make the cleared
forms tend to zero. Prime dilation transfers the nonzero constant term
of a seed at n=2 to infinitely many indices n=2p. An exact finite-field
certificate checks every one of the 16,471 possible pairs.

Reproduce from this directory

  python3 checks/verify.py
  python3 checks/check_asymptotic_comparisons.py
  python3 tools/build_paper.py

Both checkers use only Python's standard library. Use Python
3.8 or later, without -O. The build requires pdflatex with the packages
listed in the source, and Poppler's pdfinfo. Build products other than
the canonical PDF are in build/ and are ignored by Git.

checks/verification.json records the finite seed, exact denominator and
decay checks, and sampled dilation consistency checks. The general proof
is in the manuscript. The optional check_comparison.py and comparison.json
preserve a separate earlier literature calculation; it is not part of
the paper or the three-value proof.
checks/asymptotic-comparisons.json records the new Appendix A checks:
91 Rivoal--Zudilin cases, 429 Fischler 2018 cases, 111 Fischler 2020
cases, and exact constants for the Nishimoto and square-root comparisons.
These compare the specified estimates; they are separate from the proof
of the three-value theorem and from complete audits of the source papers.

Review and source scope

An independent audit of the preceding plain-text proof was supplied on
3 November 2026. It reported no substantive mathematical gap and rebuilt
the finite certificate independently. review/independent-audit.json gives
the exact audited hash, reported checks and limits. The reviewer name
and model were not specified for that report. The typeset manuscript is
a subsequent exposition revision and has not received a new independent
audit. No formal verification is claimed.

review/assembly-review.json preserves the initial 12-page assembly review
and its hashes. review/exposition-review.json records the 13-page revision
with its explanatory examples. review/clarity-review.json records the
14-page revision with the grouped-denominator and decay explanations.
review/front-matter-review.json records the earlier title and abstract edit.
review/focus-review.json records the focus revision: the abstract states
only the main theorem, the corollary and separate comparison appendix are
removed, and the alternating weights are introduced when the construction
needs them. review/title-review.json records the selected title and opening.
review/literature-comparison-review.json records the preceding version,
with the explicit Sorokin/Huttner statement and the OpenAI Catalan reference.
review/strategy-review.json records the structural revision: the proof
strategy precedes the machinery, the seed-to-nonvanishing connection is
explicit, and the infinite-tail calculation is expanded.
review/numbering-review.json records the numbering update and the explicit
contradiction in Subsection 1.2. review/wording-review.json records the
prose revision throughout all seven sections. review/explicit-title-review.json
records the explicit list in the title. review/asymptotic-context-review.json
records the brief comparison to asymptotic linear-independence results.
review/asymptotic-wording-review.json records the explicit count of linearly
independent values. review/asymptotic-attribution-review.json records the
direct attribution of the companion result. review/asymptotic-endpoints-review.json
records the preceding 13-page version, with the distinction framed in terms
of sufficiently large endpoints. review/asymptotic-comparison-appendix-review.json
records the addition of the introduction paragraph, Appendix A, primary-source
records and exact comparison checker. review/sequences-review.json records
the explicit introduction of sequences of forms and the integer sequence.
review/combination-review.json identifies the current 16-page version,
clarifying that the three forms at each index are combined to give one
sequence of rational numbers.
Each record identifies its exact PDF and source version; earlier reviews
retain their scope.

sources/reading-record.json records primary-source versions, hashes and
actual reading scopes. No third-party PDFs are included. The novelty
claim is qualified: the explicit count of three through 364 appears new among
checked sources. The introduction gives the block conditions in Sorokin
and Huttner and cites OpenAI's irrationality theorem for beta(2).

Prepared with Codex under the research pseudonym twinheterodox.
