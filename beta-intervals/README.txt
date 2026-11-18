Two shorter intervals of even Dirichlet beta values
twinheterodox, 17 November 2026

Read output/pdf/beta-intervals.pdf. The complete LaTeX source is
paper/beta-intervals.tex.

Each of the following lists contains an irrational value.

  beta(6), beta(8), ..., beta(66)     31 values
  beta(8), beta(10), ..., beta(98)    46 values

One value in their overlap could satisfy both conclusions. The previously
proved upper endpoints for the same starting indices were 74 and 110.
The proof combines the isolated-pole denominator estimates with the
powered kernel, and proves nonvanishing for the two powers used here.

Reproduction
------------
Run these commands from the beta-intervals directory.

  python3 -B checks/interval_seeds.py
  python3 -B tools/build_paper.py

The seed recurrence uses the Python standard library. Building the PDF
requires pdfLaTeX and Poppler's pdfinfo. See checks/README.txt for the
second seed calculation and the independent audit programs.

Supporting records
------------------
  checks/recurrence-results.json             Seed residues from the recurrence
  checks/cross-check-results.json            Direct-product coefficient comparison
  checks/independent-audit/                  Independent code and recorded output
  sources/repository-inputs.json             Sources used for the construction and estimates
  review/independent-audit-2026-11-17.txt    Full audit of the original manuscript
  review/revision-review.txt                 Review of the revised decay and seed proof
  review/publication-record.json             Audit versions and publication checks
  paper/build-record.json                    Source, PDF and build information
  MANIFEST.json                             Files in the paper project
