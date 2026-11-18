Exact checks for the interval paper

Run the commands below from the beta-intervals directory with Python
assertions enabled.

  python3 -B checks/interval_seeds.py

The standard-library recurrence prints the nine residues in
checks/recurrence-results.json. The powers are 1, 2 and 3 and the primes
are 10007, 10009 and 10037. Power 1 gives a control. A nonzero residue for
each of powers 2 and 3 supplies the finite input to the two interval proofs.

  python3 -B checks/interval_seeds.py --cross-check

This optional check uses python-flint. It compares all 37,386 principal-part
coefficients with direct truncated products and checks the harmonic sums in
two orders. The recorded output is checks/cross-check-results.json.

  python3 -B checks/independent-audit/run_checks.py

An independent AI auditor wrote these programs separately. They check kernel
arithmetic, prime dilation, seed residues, all 38,371 intervals in the
prime-density integral, and enclosures for the phase. The run requires
Python 3.9 or later, NumPy, mpmath, and a C++17 compiler named g++. Its
output is recorded in checks/independent-audit/results.txt. The audit
scope and the subsequent revisions are described in review/.
