EXACT CERTIFICATES

Run python3 verify.py from this directory (Python 3.8+, standard library
only; do not use -O). It writes verification.json and fails if a required
seed residue, denominator check, decay comparison or pair determinant
check fails.

seed.py computes the unnormalized n=2 partial fractions modulo 1,000,003
by two methods: a logarithmic-derivative recurrence and multiplication of
truncated series. It checks 1,092 matching coefficients, reconstruction
at three points, and zero odd coefficient sums. The all-pairs test checks
16,471 determinants. verify.py explicitly requires that the zero list is
empty, rather than merely recording that list in its output.

verify.py separately computes the normalized seed coefficients over the
rationals, checks 1,092 coefficient divisibilities and the three constant
divisibilities, and recovers the same three modular seed residues. It
also checks the two decay comparisons by exact rational arithmetic.

The 12 dilation evaluations and 12 evaluations in other pole classes are
consistency checks. They do not replace the factor and local-power-series
arguments for all sufficiently large primes in Section 5 of the paper.

For Appendix A, run python3 check_asymptotic_comparisons.py. It writes
asymptotic-comparisons.json and checks the specified finite estimates at
endpoint 364: all 91 cases in Rivoal--Zudilin's uniform family, all 429
Fischler 2018 pairs and all 111 Fischler 2020 pairs. It also checks the
rational constants in the appendix's general inequalities for Nishimoto
and the companion square-root bound. It is self-contained and uses only
the standard library. The comparison does not audit complete source
proofs or exclude refinements of their estimates.

Optional historical calculation: python3 check_comparison.py writes
comparison.json and checks an earlier specialization of Sorokin's estimate.
The former appendix about that specialization was removed. This preserved
calculation is separate from the current Appendix A and the theorem proof.

No floating-point value is used to decide a certificate inequality or
determinant test. Decimals in the output are illustrative. These supplied
scripts are distinct from the independent auditor's reported calculation;
the independent checker source was not supplied.
