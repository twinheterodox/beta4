AT LEAST ONE OF beta(4), beta(6), ..., beta(32) IS IRRATIONAL
twinheterodox; prepared with Codex; 11 November 2026

Canonical PDF: output/pdf/fifteen-beta-values.pdf
LaTeX source: paper/fifteen-beta-values.tex

The paper proves an irrational member of a set of fifteen even Dirichlet
beta values, improving Krisman's seventeen-value set through beta(36).
The new step is a denominator refinement: a pole alone in its residue
class has integral principal parts, and a zero identity avoids its
potential contribution to the rational constant's denominator. Only
repeated pole classes enter the large-prime cost.

The paper gives one complete proof with parameters r=3/4 and base index
N=96. It includes an isolated-pole example, the exact arithmetic and
phase bounds, and a nonzero seed transferred to infinitely many indices
by prime dilation. Appendix A supplies the small-prime estimates.

BUILD

From the repository root, with LaTeX and Poppler installed:

    python3 -B fifteen-beta-values/tools/build_paper.py

The build writes the canonical PDF and paper/build-record.json. It
rejects unresolved references and overfull boxes. Temporary files go
under build/ and are ignored by Git.

CHECK

    python3 -B fifteen-beta-values/checks/run_checks.py

See checks/README.txt for scope, optional dependencies and the command
that also regenerates the truncated-product seed calculation.

PREPARATION AND REVIEW

The denominator insight and fifteen-value conclusion were obtained
separately in two Codex research sessions using GPT-6-Astra max and
public sources, without access to each other's unpublished work. This
manuscript uses one parameter choice and incorporates the simpler
polynomial maximum argument developed in the second session.

An independent audit reported on 11 November 2026 found no mathematical
error or gap affecting the theorem and independently reproduced all
three computational inputs. Its report and the hashes of the audited
14-page version are in review/independent-audit-2026-11-11.txt and
review/independent-audit-2026-11-11.json. The auditor's unchanged programs
and numerical results are in checks/independent-audit/.

The current PDF adds an audit note and explanations of the pole-block
residue counts and the two odd-product lengths, and explains the integer
contradiction and the arithmetic integral before their machinery. These
later edits and their validation are recorded in
review/audit-response-2026-11-11.json and
review/exposition-update-2026-11-11.json.
Final source and PDF consistency checks are recorded in
review/publication-preparation-2026-11-11.json.
The original internal review remains in review/assembly-review.json
with its historical scope. The independent audit accepts the stated
classical inputs; it is not formal verification or a priority assessment.

Sources, versions, checksums and actual reading scopes are recorded in
sources/reading-record.json. Third-party PDFs are not redistributed.
