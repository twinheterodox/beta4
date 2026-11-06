# Changelog

## 6 November 2026 — Shifted-interval comparisons and exposition

- Added the interval β(4), …, β(198), its explicit parameter choice and exact rational decay check. It improves the general theorem's endpoint 200; the earlier count's endpoint 160 and Krisman's endpoint 36 remain shorter.
- Compared the five individual endpoints with the general theorem in the opening table and kept the later comparison table consistent. This additional example is recorded separately from the four examples in the 5 November independent audit.
- Added proof-method attribution and explained the construction's design, the seed choice through a two-form example, the role of the Gram matrix, the need for grouped denominator cancellation, and the endpoint coefficients. The current manuscript is 17 pages; its theorem and exact-check results are unchanged.

## 5 November 2026 — Irrational values in shifted intervals

- Added the separate 16-page paper *Irrational even beta values in shifted intervals*, with complete LaTeX sources, standalone exact checkers and a reproducible build.
- For every positive odd m, the paper gives an irrational member of β(m+1), …, β(67m−1). Its general family guarantees any fixed number of irrational values with an upper endpoint linear in m, and includes sharper individual intervals.
- Recorded the independent audit of the underlying plain-text proof, which found no substantive error or gap and independently reproduced the finite calculations. Applied the requested wording clarification in the prime exclusion.
- Kept the audit's scope separate from the subsequent LaTeX assembly, mathematical self-review and visual inspection. The earlier manuscripts and PDFs are unchanged.
- Expanded Section 1.1 to explain linear versus quadratic endpoint growth, distinguish the earlier direct shifted bounds, and compare concrete intervals starting at β(6), β(20) and β(100).

## 5 November 2026 — Independent audit of the uniform count paper

- Recorded the supplied independent mathematical audit: no substantive error or gap found in the 19-page count-paper revision, with the finite calculations reproduced using independently written arithmetic.
- Preserved the reported scope, check counts and the PDF and source hashes captured before revision. The reviewer identity and model were not supplied.
- Explained why the Stirling estimate remains uniform down to the fixed lower endpoint t = 1/B, as recommended by the audit, and added an audit note to Appendix A.
- Updated the review status and recorded the subsequent exposition edits separately from the audited snapshot. Earlier audit records retain their original scopes.
- Rebuilt the count paper as a 20-page PDF, including the uniformity explanation and audit note.

## 4 November 2026 — Uniform irrationality counts from s = 160

- Added the separate 18-page paper *Explicit lower bounds for the number of irrational even beta values* in `even-beta-counts/`.
- The main theorem gives $I(s)\ge\lfloor\pi\sqrt{2s/(135\log s)}\rfloor$ for every even s ≥ 160. It includes sharper finite endpoints and the bound without a floor for every even s ≥ 42,170.
- Proved uniform rank through a positive Gram matrix, with complete denominator, decay and prime-dilation arguments. Added a two-grid example, the local power-series justification of coefficient extraction, and an explicit derivation of the square-root constant.
- Added standalone exact construction and endpoint checkers, their results, a reproducible build, and source and review records. The abstract states results only; all subsections are numbered.
- Explained the full 0.21 irrationality-count consequence of the earlier proof in the new paper, using its estimate for the span including 1. The earlier manuscript and PDF retain their original versions.
- Kept the earlier counting result in the introduction and moved its justification to Appendix B. The current PDF has 19 pages.
- Preserved the existing three-through-364 paper and the linear-independence paper unchanged. The new manuscript has been self-reviewed and computationally checked; earlier independent audits retain their original scopes.

## 4 November 2026 — Finite consequences of asymptotic estimates

- Added a brief introduction paragraph and Appendix A checking whether the finite estimates behind selected asymptotic results yield three irrational values through 364.
- Compared Rivoal–Zudilin's uniform family, Nishimoto's explicit estimate, Fischler's 2018 and 2020 size estimates, and all admissible parameters in the companion square-root bound. The specified estimates do not give the three-value count at this endpoint.
- Added a self-contained exact comparison checker, its results, and primary-source records. The theorem proof and its certificate are unchanged; the existing independent audit retains its original plain-text scope.
- Rebuilt and visually reviewed the 16-page PDF.
- Clarified the proof overview: three sequences of forms, a fixed integer combination at every index, and a sequence of integers tending to zero with infinitely many nonzero terms.

## 4 November 2026 — Finite-set manuscript wording

- Added a brief comparison to asymptotic linear-independence results, describing the square-root order as the strongest known to us and citing the pinned companion paper.
- Changed the title to *At least three of β(2), β(4), …, β(364) are irrational.* and updated the PDF title metadata.
- Revised prose throughout all seven sections, tightening the introduction and clarifying the denominator, dilation and seed arguments.
- Kept the theorem statements, displayed equations and certificate unchanged. The current PDF has 13 pages, with numbered subsections and a clean build.

## 3 November 2026 — Finite count of irrational even beta values

- Prepared the standalone paper *At least three even beta values through β(364) are irrational*, now 14 pages, proving that at least three of β(2), β(4), …, β(364) are irrational.
- Included exact seed and decay certificates and the prime-dilation proof.
- Recorded the independent audit of the preceding plain-text proof separately from the typeset exposition review, with exact source hashes and reading scopes.
- Revised the construction to motivate the pole spacing, reflection, zero blocks, endpoint and normalization before their formulas. Added a worked prime-dilation example modulo 5, including the exchange of the fractional shifts, and kept the certificate code with its introduction.
- Defined the decay profile before its formula and supplied the gamma-product conversion. Reorganized grouped integrality around the unique possible exceptional pole group, and clarified the comparison with the list omitting β(2).
- Removed the manual title line break and limited the abstract to the main theorem. Removed the redundant corollary and separate comparison appendix, and introduced the alternating weights where the construction needs them.
- Stated the Sorokin and Huttner block theorem with its parameter conditions, cited OpenAI's irrationality theorem for β(2), and removed unsolicited optimality disclaimers.
- Reorganized the exposition around the contradiction assumption, explained the determinant and seed propagation before the computation, and expanded the infinite-tail estimate. Moved the exact decay checks into the decay section.
- Numbered all subsections and made explicit why integers tending to zero but nonzero infinitely often yield the contradiction.

## 3 November 2026 — Exposition, independent audit and consolidation

- Updated the even-beta paper to the 19-page revision, with worked denominator and common-row examples, the six-point rank explanation and the parameter balance behind the square-root bound.
- Recorded the independent model proof audit of that revision: no substantive mathematical gap identified, with the supplied and additional computational checks reported to pass.
- Added the exposition model attribution and preserved the audit report and reviewed revision hashes.
- Consolidated the manuscript into this repository as its single working copy, with the current PDF, build instructions and review records together.

## 3 November 2026 — Author attribution

- Added twinheterodox as the author of *Linear independence of even Dirichlet beta values*, retaining the Codex preparation credit, and updated the PDF author metadata.

## 2 November 2026 — Even beta values

- Added *Linear independence of even Dirichlet beta values*, with its PDF, complete LaTeX source, exact supporting checks and source-reading records in `even-beta-independence/`.
- The paper proves a lower bound of $0.21\sqrt{s/\log s}$ for the rational span of β(2), β(4), …, β(s), for every sufficiently large even s. It does not establish the irrationality of an individual value such as β(4).
- Recorded the reported independent model cross-check separately from computational checks and the developing assistant's audit.
- Updated the main README to present the even-beta result alongside the existing determinant-obstruction paper and its supplement.

Earlier determinant papers, supporting records and their history remain available at their existing paths.
