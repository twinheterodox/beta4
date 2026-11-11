Portable checks for the further-investigations supplement

Run from the supplement directory:
  python3 -B checks/replay.py

Dependencies: Python 3.11 or newer, python-flint 0.9.0, SymPy 1.14.0.
The tested environment used Python 3.14.4. Assertions are required; do
not use -O. No external services or private records are accessed.

Suites (select with --suite NAME)

isolation
  Reconstructs the recorded coefficient cancellations, checks 48 exact
  prime-dilation congruences, and recomputes the finite polynomial-kernel
  ranks and telescoping generators through total degree 16. Some routines
  share the exact polynomial constructor; this is not an independent audit.

pade
  Checks seven symbolic two-variable adjoint-kernel identities, integrality
  of 144 recorded rational coefficients under the improved multiplier,
  72 divisibilities between clearing factors, nine integer-polynomial
  identities, 118 dyadic moment bounds and three rational rate enclosures.
  It reconstructs the nine polynomials; the 144 coefficients are stored
  rational inputs from the earlier exact construction. Historical archive
  hash assertions were replaced by the portable evidence manifest.

digamma
  Reconstructs 714 determinant polynomials by exact interpolation, checks
  1,428 extra evaluations, and verifies the common-cofactor gcd in 168
  parameter/geometry cases. The all-size proof is separate.

canonical
  Replays 26 global row-potential bounds and the repair positive-part
  integral. Reconstructs the degree-200 derivative polynomials, isolates
  every critical point, verifies its maximum enclosure, and independently
  recounts the roots with Sturm sequences. Recomputes the three seed
  determinant polynomials and all 313 field evaluations from the supplied
  rational seed matrices. Those matrices are inputs here; reconstruction
  of their moments belongs to the earlier independent audit. Finally
  recomputes the Jacobi integral constant and combines the certified bounds.

small
  Reconstructs partial fractions by truncated product division, independently
  of the original producer's logarithmic-derivative recurrence. Checks all
  222 complete forms and primitive vectors, 2,772 partial-fraction entries,
  the prime-residue formula and zero cancellation. Evaluates beta values
  with an Euler-transformed series, distinct from the producer's Hurwitz
  zeta evaluation. It imports neither original producer.

algebra
  Tests explicit Smith-progression attainment, chain left inverses,
  factorial-alternant comparisons and determinant second moments, and
  evaluates the piecewise prime-cost integral exactly. These are finite
  checks of mechanisms whose general proofs are given in the manuscript.

Each run uses a fresh ignored build/checks/ directory and saves logs and
new result files there. review/replay-NAME.json records completion. The
curated evidence files retain their original hashes unless the provenance
manifest explicitly identifies a path or metadata adaptation. Newly
written assembly checkers are identified separately in that manifest.

The complete catalog includes shorter historical tests that are summarized
but not replayed by this package. No blanket claim is made that every
calculation in every private record has been rerun. The independently
reported canonical audit and this assembling-assistant replay are distinct.
