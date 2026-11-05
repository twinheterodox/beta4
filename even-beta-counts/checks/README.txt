EXACT CHECKS FOR THE UNIFORM EVEN-BETA COUNT PAPER

From the paper directory:

  python3 -B checks/verify_construction.py
  python3 -B checks/verify_endpoints.py

Both scripts use only Python's standard library. Run with Python 3.8+
without -O. They are self-contained within this paper directory and
have no dependency on private research records or other checkouts.

Construction checks

base.py contains arithmetic, rational decay and partial-fraction helpers.
gram.py constructs the integer kernel polynomials and the Gram moments.
verify_construction.py checks all ten parameter rows in Section 8:

  - exact rational f(u)<1 and U<1 comparisons;
  - two finite-field partial-fraction calculations;
  - two constant-term calculations and the residue identity;
  - comparison with 385 independently constructed rational Gram entries;
  - 55 positive leading principal minors and ten integer adjugate identities;
  - all selected and other-pole coefficients at the small dilation tests;
  - complete rational constant-term congruences and grouped divisibility.

The dilation tests use (B,s,r)=(1,20,3),(2,24,3), seed n=1 and
prime indices n=17,19. These small parameters check the identities,
not additional theorem endpoints. They cover 272 selected coefficient
congruences, 2312 other-pole valuation bounds, 10 constant congruences
and 4802 grouped constant divisibilities.

Floating-point calculations only propose a rational upper bound for
the saddle point. The subsequent f and U comparisons use Fraction
arithmetic exactly; decimals in the output are illustrative.

Endpoint checks

verify_endpoints.py uses the positive series for log(y)=2 atanh(z),
with an explicit rational remainder after 18 terms and reduction to
1<=y<2. It computes totients by direct coprime counting and isolates
each necessary even endpoint between equal integer bounds.

It checks all 230 intervals needed before the analytic infinite range
in Proposition 9.2, all 23 remaining rounded intervals down to s=160,
and the rational constants in the analytic argument for every k>=256.
The floor values at 160 and 364, the positive coefficient after the
explicit rounding error, and four direct inversions of the finite
family are also certified.

The output files are construction-verification.json and
endpoint-verification.json. The general all-index proof and the
uniform positive-definiteness argument are in the manuscript.
These computations are not an independent proof audit or a formal
verification of those arguments.
