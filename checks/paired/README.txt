PAIRED QUARTER-SHIFT IDENTITIES — STANDALONE REPLAY

Purpose
The supplement's paired-identity case study gives two exact presentations
of the same rational coefficients, at a=1/4 and 3/4. This directory contains
the algebraic data and a replay that does not require the original research
workspace, the upstream notebooks, ramanujantools, or a network connection.

From the repository root, using Python with python-flint installed:

  python3 checks/paired/verify.py
  python3 checks/paired/arithmetic.py
  python3 checks/paired/rates.py

The last command reads arithmetic-result.json, so run them in this order.
Requirements are included in ../supplement-requirements.txt. The replay is
small enough to run on an ordinary workstation (about 30 seconds for the
identity checks on the writing machine). It regenerates the large polynomial
certificates in memory and checks full substitution; it does not serialize
roughly 186 MB of expanded decimal coefficients.

Files and meanings
- quarter.json, threequarter.json: rational matrices M,Q,K,C,CQ and scalars.
  Each rational function is [[numerator coefficients],[denominator
  coefficients]], with coefficient order 1,n,n^2,... . Every coefficient is
  an exact integer or rational string. These are input data, not verdicts.
- algebra.py: exact Q(n) arithmetic, polynomial products, shifted sign tests,
  and a descending-coefficient solution of A(k)P(k+1)-B(k)P(k)=F(k).
- problems.py: all linear factors, shifts and normalizations for the five
  polynomial problems at each quarter.
- readouts.py: direct finite reconstruction of Q,P3,P4 and dual residue rows.
  Simple-zero Gamma-derivative tails are included explicitly.
- verify.py: matrix transport, Krylov companions, regularity, initial cross
  products, normalizing scalar, ten full polynomial substitutions, boundary
  degree/sign checks, exceptional recurrence indices, and initial readouts.
- arithmetic.py: exact floor partitions and rational integration for the
  ceiling on the specified valuation estimates. It is not a bound on actual
  gcds or primitive height.
- rates.py: exact equality of the two limiting matrices, the cyclic-vector
  check, and Arb enclosures of the rate comparison.
- result.json, arithmetic-result.json, rates-result.json: current outputs.
- provenance.json: hashes of the original data before canonicalization.
- import_data.py: historical import utility; unnecessary for a replay.

How this certifies an all-index statement
The matrix and polynomial substitutions are identities in a symbolic n,
not checks at sampled values of n. The shifted polynomial sign tests prove
nonzero denominators at all required integer indices. The supplement
explains how the telescoping identities sum, how the derivative-tail
boundaries vanish, and why the finite exceptional steps complete induction.
Those analytic arguments are part of the proof and are not replaced by a
computer PASS flag. The numerical rate comparison additionally uses the
explicitly cited asymptotic recurrence theorem and the identifying bounds
in the manuscript.

Status
This replay was assembled by the same writing model. It is not a new
independent-model audit, a human audit or proof-assistant formalization.
The identities do not establish useful primitive size or nonvanishing of
the final character-projected form for all indices.
