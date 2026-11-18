An asymptotic lower bound for the number of irrational even beta values
twinheterodox, 17 November 2026

Read output/pdf/asymptotic-even-beta-counts.pdf. The complete LaTeX source
is paper/asymptotic-even-beta-counts.tex.

For even s, let I(s) count the irrational values among beta(2), beta(4),
..., beta(s). The paper proves

  liminf I(s)/sqrt(s log log s/log s) >= 2 sqrt(kappa/3),
  kappa = zeta(2) zeta(3)/zeta(6).

The constant rounds to 1.6098018250. The bound also holds after omitting any
fixed finite prefix. The proof gives an asymptotic result. For a bound with
an explicit starting index, see the companion project ../even-beta-counts/.

The determinant construction adapts the odd-zeta method of Fischler, Lai,
Sprang and Zudilin to the character modulo four. Only fractions with reduced
denominator divisible by four contribute. Restricting the numerator to
those fractions preserves the zeros used for cancellation and divides the
leading arithmetic cost by three.

Reproduction
------------
Run these commands from the asymptotic-even-beta-counts directory.

  python3 -B checks/verify_count_examples.py
  python3 -B tools/build_paper.py

The finite checks use Python 3 and SymPy. Building the PDF uses pdfLaTeX and
Poppler's pdfinfo. See checks/README.txt for the calculation and its scope.

Supporting records
------------------
  checks/count-examples.json                 Exact outputs and tail witnesses
  sources/reading-record.json                Editions, digests and reading scope
  review/independent-audit-2026-11-17.txt    Original proof and proposed improvement
  review/restricted-grid-correctness.txt     Review of the strengthened proof
  review/restricted-grid-exposition.txt      Review of its exposition
  review/publication-record.json             Audit versions and publication checks
  paper/build-record.json                    Source, PDF and build information
  SHA256SUMS.json                            Checksums of the publication files
