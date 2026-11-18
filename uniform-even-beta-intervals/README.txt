Uniform intervals for irrational even beta values
twinheterodox, 18 November 2026

Read output/pdf/uniform-even-beta-intervals.pdf. The complete LaTeX
source is paper/uniform-even-beta-intervals.tex.

Starting at any even index 2q+2 with q >= 1, the even beta values
through 2 floor(33q/2) include an irrational number. Farther out,
the upper endpoint can be shortened. For every positive odd m there
is an even s between m+1 and (299/20)m+C with beta(s) irrational,
where C is absolute. For odd m >= 4,000,000,001 one can take s < 15m.

The fixed-power construction uses a nonzero rational seed. A residue
representation gives bounded gaps between successful exponents for two
families of seeds. A congruence at prime exponents, explicit estimates
for primes, and two exact finite covers give the numerical endpoints.
Individual nonvanishing at every original exponent remains open.

Reproduction
------------
Run these commands from the uniform-even-beta-intervals directory.

  python3 -B tools/build_paper.py
  python3 -B checks/verify.py

The PDF build requires pdfLaTeX, its standard math packages, Latin Modern,
microtype, and Poppler's pdfinfo. The certificate suite requires Python 3
with python-flint. Both the five programs and the five reference tables
they reproduce are in checks/. See checks/README.txt for the commands
corresponding to individual stages of the proof.

The computations of the real phase and prime-density inputs are in the
earlier fifteen-beta-values paper, as cited in the manuscript. The new
programs reproduce the finite calculations specific to this paper.
sources/inputs.txt identifies the inputs. The reviews in review/ give
the scope of the proof checks and the separate exposition check.
