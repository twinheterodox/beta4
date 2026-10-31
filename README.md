# Arithmetic obstructions in approximations to β(4)

This repository contains a focused paper and a companion investigation record on constructions for the irrationality of

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

The paper studies two positive determinant families. In each, the real determinant becomes very small, but the **actual primitive integer polynomial**, after removing all common coefficient factors, has a value that grows. These are mathematical obstructions for the specified families. They do not exclude other proofs of irrationality.

**Status:** complete drafts of the main paper and technical supplement, revised for readers unfamiliar with the constructions, kept local and intended for `twinheterodox/beta4-obstructions`. No remote repository or release has been created. The assembled drafts have not received a fresh independent audit.

Start with the focused paper; use the supplement for the broader investigation.

| Document | Contents | Source and review |
| --- | --- | --- |
| **[Main paper (16 pages)](output/pdf/beta4-determinant-obstructions.pdf)** | Complete proofs of the two determinant obstructions, with a shared proof outline and worked arithmetic example. | [LaTeX](paper/beta4-determinant-obstructions.tex) · [Assembly audit](paper/audit.txt) |
| **[Technical supplement (31 pages)](output/pdf/beta4-investigations.pdf)** | A guide to the constructions, selected complete proofs, exact reductions, unresolved arithmetic, and all 57 grouped records. | [Source guide](supplement/README.md) · [Assembly audit](supplement/audit.txt) |

Both abstracts and introductions now start from the irrationality criterion and explain the comparison between approximation error and arithmetic cost. The [accessibility revision record](planning/accessibility-review.txt) describes the added explanations and layout checks; the existing statements and proofs are unchanged.

The [reproduction instructions and results](checks/README.md) distinguish main-paper exact/interval checks from the supplement's certificate replays and exact reconstructions. The supporting documents are:

- [Theorem inventory](planning/theorem-inventory.txt): 57 grouped claims, with precise scope, evidence, reading status and proposed placement.
- [Novelty audit](planning/novelty-audit.txt): established ingredients, closest primary comparisons, candidate contributions and unresolved priority questions.
- [Method catalog](planning/method-catalog.txt): coverage of all 12 method groups and 18 coordinated research rounds, including repaired and superseded conclusions.
- [Evidence guide](evidence/README.txt): 120 original verbatim snapshots and seven supplemental source files, including proofs, review records and selected calculation scripts/results.
- [Source and reading record](sources/reading-record.json): versions, local PDF locations, checksums, reading scopes and unread dependencies.
- [Original paper plan](planning/paper-plan.txt) and [opening proposal](planning/draft-opening.txt): the selection-stage documents retained in the history.

The main results are a uniform primitive-growth obstruction for shifted consecutive principal determinants, and an exact large-prime valuation theorem leading to a primitive-growth obstruction for difference-variable Gram determinants. The supplement distinguishes further proved limitations from incomplete constructions and finite evidence. Its positive-source cone and ordinary-moment first variation also explain why particular-family failures do not settle the broader problem. Priority of the specific statements remains subject to the recorded literature comparison and further scrutiny.

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

To build the companion supplement:

```sh
python3 tools/render_supplement_register.py
python3 tools/build_supplement.py
```

The assembly audits describe their respective writing passes; the later accessibility revision is recorded separately. This checkout is under `beta4/obstructions/`, with its own Git history. The original working archive is in the enclosing `beta4/` directory, and third-party PDFs are in `../literature/papers/`. To compare the evidence with that local archive and check the literature files, run `python3 tools/validate_packet.py --archive-root .. --check-local-pdfs`. Earlier audit records describe the directory layout at the time of their review. The original material in this repository is licensed under Apache 2.0; cited works retain their own terms.
