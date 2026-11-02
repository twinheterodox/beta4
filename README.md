# Research on the irrationality of β(4)

This repository collects research on the irrationality of the Dirichlet beta value

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

Its scope includes constructions, proofs, computational checks and documented obstructions across approaches to this problem. The first contributions are a focused paper on two positive determinant families and a self-contained supplement covering the broader investigation.

In the two families studied in the main paper, the real determinant becomes very small, but the **actual primitive integer polynomial**, after removing all common coefficient factors, has a value that grows. These are mathematical obstructions for the specified families. They do not exclude other proofs of irrationality. The present papers do not establish the irrationality of β(4).

**Review status:** the main paper, technical supplement and accompanying code have passed an [independent model review](planning/independent-model-review-2026-11-02.json), reported on 2 November 2026.

Start with the focused paper; use the supplement for the broader investigation.

| Document | Contents | Source and review |
| --- | --- | --- |
| **[Main paper (21 pages)](output/pdf/beta4-determinant-obstructions.pdf)** | Explains the irrationality strategy, why these constructions are natural, their provenance and what their failure establishes, followed by complete proofs and a worked arithmetic example. | [LaTeX](paper/beta4-determinant-obstructions.tex) · [Final editorial reading](planning/final-reading-2026-11-02.txt) |
| **[Technical supplement (67 pages)](output/pdf/beta4-investigations.pdf)** | Self-contained case studies with explicit constructions, proofs, remaining questions and standalone certificates; a separate catalog covers all 57 grouped records. | [Source guide](supplement/README.md) · [Final editorial reading](planning/supplement-final-reading-2026-11-02.txt) |

The main paper explains the strategy before introducing the technical formulas: what would contradict rationality, why determinants might help, and what happens when their coefficients are made integral. Each main theorem is followed by a plain-language interpretation of its quantitative conclusion and scope. The supplement develops the substantial unfinished routes as reproducible case studies. Its paired hypergeometric identities have a standalone exact replay, and each unfinished case states the missing mathematical estimate.

The [structural review](planning/narrative-review.txt), [equation-by-equation review](planning/expository-review.txt), [accessibility review](planning/accessibility-review.txt) and [main-paper assembly audit](paper/audit.txt) record earlier revisions. Their preservation and checking statements apply to those historical versions.

The [reproduction instructions and results](checks/README.md) distinguish main-paper exact/interval checks from the supplement's certificate replays and exact reconstructions. The supporting documents are:

- [Theorem inventory](planning/theorem-inventory.txt): 57 grouped claims, with precise scope, evidence, reading status and proposed placement.
- [Novelty audit](planning/novelty-audit.txt): established ingredients, closest primary comparisons, candidate contributions and unresolved priority questions.
- [Method catalog](planning/method-catalog.txt): coverage of all 12 method groups and 18 coordinated research rounds, including repaired and superseded conclusions.
- [Evidence guide](evidence/README.txt): 120 original verbatim snapshots and seven supplemental source files, including proofs, review records and selected calculation scripts/results.
- [Source and reading record](sources/reading-record.json): versions, local PDF locations, checksums, reading scopes and unread dependencies.
- [Original paper plan](planning/paper-plan.txt) and [opening proposal](planning/draft-opening.txt): the selection-stage documents retained in the history.

The main results are a uniform primitive-growth obstruction for shifted consecutive principal determinants, and an exact large-prime valuation theorem leading to a primitive-growth obstruction for difference-variable Gram determinants. The supplement distinguishes further proved limitations from incomplete constructions and finite evidence. Its positive-source cone and ordinary-moment first variation also explain why particular-family failures do not settle the broader problem. Priority of the specific statements remains subject to the recorded literature comparison and further scrutiny.

The work was prepared with Codex. The independent model review, archived model cross-checks, exact rational calculations, interval calculations and source reading are recorded separately. The independent review is a model cross-check; it is not a human mathematical audit or proof-assistant verification.

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

The repository includes the evidence and certificates needed by these instructions; the original working archive is not required. Source ledgers retain URLs, hashes and local paths from the research environment. Third-party paper PDFs are not redistributed. The validator's optional `--archive-root` and `--check-local-pdfs` flags are for comparisons in that original environment. Historical audit records retain the directory layout and publication status at the time they were written.

The original material in this repository is licensed under [Apache 2.0](LICENSE); cited works retain their own terms.
