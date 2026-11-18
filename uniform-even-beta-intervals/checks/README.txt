Exact checks for Uniform intervals for irrational even beta values

Run from the directory containing paper/, checks/, output/, and tools/.
Use Python 3 with python-flint installed (tested with python-flint 0.9.0).
The complete replay compares its five recomputed JSON objects with the
five reference files distributed in this directory.

  python3 -B checks/verify.py

Individual calculations
-----------------------
  python3 -B checks/derivative_seed.py --cross-check
    Rebuild the ninth-power, twentieth-derivative seed modulo 10007
    and 10009. Compare the two Laurent coefficient algorithms and the
    two orders of harmonic summation. Output: derivative-seed.json.
    Omitting --cross-check requires only the Python standard library
    and omits the direct-product comparison.

  python3 -B checks/explicit_threshold.py
    Compute the exact height bound for the fixed rational seed and both
    rational prime-band inequalities. Output: explicit-threshold.json.
    This program uses only the Python standard library.

  python3 -B checks/prime_bridge.py
    Compute the 135 prime exponents, the ninth-power seed modulo each
    prime, and the two end points of the covered range of q. Compute
    the residues in one product ring. Output: prime-bridge.json.

  python3 -B checks/initial_intervals.py
    Cover q=1,...,1890 using 89 seed evaluations modulo 65537. Check
    the defining differential identity for all the local Taylor
    coefficients used before reflection. Output: initial-intervals.json.

  python3 -B checks/interval_replay.py
    Compute the 89 seeds at a second prime, 100003. Rebuild each
    contributing address directly, sum the harmonic correction in
    reverse order, and check the interval inclusions. Replay the
    bridge in each of its 135 prime fields and the two control fields.
    This program reads the two bundled coverage tables.
    Output: interval-replay.json.

All individual commands print JSON to standard output. Kernel geometry
and cover definitions are shared between implementations. The second
calculation checks the first using another modulus or summation order.
The two Laurent routines in kernel_coefficients.py are taken from
the companion beta-intervals seed checker. No program uses an
approximation to a Dirichlet beta value.
