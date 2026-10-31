# Arithmetic obstructions in approximations to β(4)

This repository contains a research-paper draft on specific constructions that fail to prove the irrationality of

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

The paper studies two positive determinant families. In each, the real determinant becomes very small, but the **actual primitive integer polynomial**, after removing all common coefficient factors, has a value that grows. These are mathematical obstructions for the specified families. They do not exclude other proofs of irrationality.

**Status:** first complete main-paper draft, kept local and intended for `twinheterodox/beta4-obstructions`. No remote repository or release has been created. The broader technical supplement is organized but remains to be written; no new independent audit of the assembled manuscript is claimed.

Start with the **[main paper (PDF)](output/pdf/beta4-determinant-obstructions.pdf)**, its [LaTeX source](paper/beta4-determinant-obstructions.tex), and the [audit record](paper/audit.txt). The [reproduction instructions and results](checks/README.md) cover the core exact and interval checks. The supporting documents are:

- [Theorem inventory](planning/theorem-inventory.txt): 57 grouped claims, with precise scope, evidence, reading status and proposed placement.
- [Novelty audit](planning/novelty-audit.txt): established ingredients, closest primary comparisons, candidate contributions and unresolved priority questions.
- [Method catalog](planning/method-catalog.txt): coverage of all 12 method groups and 18 coordinated research rounds, including repaired and superseded conclusions.
- [Evidence guide](evidence/README.txt): 120 verbatim research snapshots, including proofs, review records and selected calculation scripts/results.
- [Source and reading record](sources/reading-record.json): versions, local PDF locations, checksums, reading scopes and unread dependencies.
- [Original paper plan](planning/paper-plan.txt) and [opening proposal](planning/draft-opening.txt): the selection-stage documents retained in the history.

The main results are a uniform primitive-growth obstruction for shifted consecutive principal determinants, and an exact large-prime valuation theorem leading to a primitive-growth obstruction for difference-variable Gram determinants. Other proven limitations and incomplete constructions belong in a technical supplement; finite searches are labeled separately. Priority of the specific statements remains subject to the recorded literature comparison and further scrutiny.

The work was prepared with Codex. Archived model cross-checks, exact rational calculations, interval calculations, source reading and formal verification are distinguished throughout. There has been no proof-assistant verification of these results. The current writing pass did not commission another independent model audit.

To check the packet's internal references and evidence hashes:

```sh
python3 tools/validate_packet.py
```

To regenerate the readable inventory and coverage catalog from the JSON inventory:

```sh
python3 tools/render_inventory.py
```

To build the main paper with LaTeX and Poppler installed:

```sh
python3 tools/build_paper.py
```

The complete original working archive remains in the sibling `beta4/` directory. It has not been moved into or replaced by this repository. Third-party PDFs remain in `beta4/literature/papers/`; their contents are not included here. The original material in this repository is licensed under Apache 2.0; cited works retain their own terms.
