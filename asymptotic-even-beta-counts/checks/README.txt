Exact checks for the asymptotic counting paper

From the asymptotic-even-beta-counts directory, run

  python3 -B checks/verify_count_examples.py

Use Python 3 with assertions enabled and SymPy (the recorded run used
version 1.14.0). The script writes checks/count-examples.json. That file
includes the input parameters, an exact enclosure for the leading constant,
the active grids for B=2 and 4, and six polynomial basis products at each
of n=2 and 4 for B=2 and s=8. Each of the latter two cases includes a
finite-tail sum whose individual terms have denominators but whose grouped
sum is integral. The general cancellation and asymptotic result are proved
in the paper.
