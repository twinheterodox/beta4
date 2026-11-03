EXACT SUPPORTING CHECKS

Install the pinned dependency and run from the even-beta-independence directory:

  python3 -m pip install -r checks/requirements.txt
  python3 checks/check_local_system.py
  python3 checks/check_linear_forms.py

check_local_system.py writes local-system-results.json. It checks:

  - 52 comparisons of the distinct-product least common multiple against
    direct enumeration of subsets.
  - 2,304 exact Taylor-coefficient comparisons, including 1,204 nonzero
    coefficients, for a system of depth 7 through degree 36. One route
    solves the original differential system recursively; the other uses
    the gauge and successive integrations. The proposed denominator
    clears every tested coefficient.
  - 100 symbolic identities for the enlarged system used in Section 6.
  - Three checks of the quadratic pullback and residues.
  - 12 descending-residue identities, for top-polynomial degrees 0 to 5
    and both types of link in the triangular system.
  - A rational certificate that the limiting constant exceeds 0.21 and
    that the chosen rank-inequality slope exceeds 1.

check_linear_forms.py writes half-shift-forms-checks.json. It compares
six coefficient vectors obtained from direct partial fractions and from
the differential-system recurrence. Its fixture is

  F(t) = 1/[t^4 (t+1)^4], n=1, a=4, M=2,
  p=0,1 and derivative orders 4,5,6.

One of those vectors is zero. This does not contradict the theorem, which
uses an asymptotic weak rank condition on a family of rows, rather than
nonvanishing of each row. The fixture tests the algebraic identities; it
is not an instance satisfying all asymptotic parameter conditions.

The earlier publication replay is recorded in reproduction.json; both
check reports matched the original outputs exactly, and an isolated build
reproduced the PDF text for the revision recorded there. The current
relocation build is recorded in ../review/consolidation-2026-11-03.json.

All calculations use exact rational or symbolic arithmetic. The decimal
constant in the JSON is a diagnostic; the inequality certificate is
rational. No numerical approximation to any beta value is used.

The computational implementations were prepared and checked by the developing
assistant. Separate computational routes are not independent model review.
The current 19-page paper passed an independent model proof audit, reported
on 3 November 2026; its scope is recorded in
../review/independent-model-review-2026-11-03.json. That reviewer also
reported separate checks through Taylor degree 60 and depth 12 and
numerical tail-identity checks. Their additional code and outputs were
not supplied and are not part of these scripts.
The tests do not prove the general denominator lemma, the multiplicity
argument, or the asymptotic dimension estimate. Those require the written
proof and its imported mathematical inputs.
