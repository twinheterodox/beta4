# Research on the irrationality of β(4)

This repository collects research on the irrationality of the Dirichlet beta value

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

Its scope includes constructions, proofs, computational checks and documented obstructions across approaches to this problem, together with results about the larger family of even beta values. The papers below develop two directions: a quantitative linear-independence bound for that family, and obstructions to particular determinant approaches to the individual value β(4). None of the papers establishes the irrationality of β(4).

## Papers

| Document | Result and contents | Source and checks |
| --- | --- | --- |
| **[Even beta values: linear independence (19 pages)](even-beta-independence/output/pdf/even-beta-independence.pdf)** | A lower bound of $0.21\sqrt{s/\log s}$ for the rational span of β(2), β(4), …, β(s), for every sufficiently large even s. Full proof of the half-integer adaptation of Fischler's method. | [Paper guide](even-beta-independence/README.md) · [LaTeX](even-beta-independence/paper/even-beta-independence.tex) · [Checks](even-beta-independence/checks/README.txt) |
| **[Determinant obstructions (21 pages)](output/pdf/beta4-determinant-obstructions.pdf)** | Two positive determinant families fail the intended irrationality test, even after all common coefficient factors are removed. Includes complete proofs and a worked arithmetic example. | [LaTeX](paper/beta4-determinant-obstructions.tex) · [Checks](checks/README.md) · [Final reading](planning/final-reading-2026-11-02.txt) |
| **[Determinant investigation supplement (67 pages)](output/pdf/beta4-investigations.pdf)** | Self-contained case studies with constructions, proofs, remaining questions and standalone certificates; a catalog covers all 57 grouped records. | [Source guide](supplement/README.md) · [Checks](checks/SUPPLEMENT.md) · [Final reading](planning/supplement-final-reading-2026-11-02.txt) |

**Review status:** an [independent model proof audit of the 19-page even-beta paper](even-beta-independence/review/independent-model-review-2026-11-03.json), reported on 3 November 2026, found no substantive mathematical gap. The [full report](even-beta-independence/review/independent-model-review-2026-11-03.txt) records the proof checks and additional computations. The [determinant papers and accompanying code](planning/independent-model-review-2026-11-02.json) passed a separately reported model cross-check. These reviews use the cited mathematical results and are not formal verification.

See the [changelog](CHANGELOG.md) for repository updates.

## Even beta values

The new paper proves, for every sufficiently large even integer $s$,

$$
\dim_{\mathbb Q}\operatorname{span}_{\mathbb Q}\{\beta(2),\beta(4),\ldots,\beta(s)\}
\geq 0.21\sqrt{\frac{s}{\log s}}.
$$

This concerns linear independence in a growing family, rather than the irrationality of an identified individual value. The proof uses half-integer sums to select the even arguments, controls denominators at an ordinary point of a differential system, and adapts Fischler's multiplicity argument through a quadratic change of variable. The [paper guide](even-beta-independence/README.md) includes the proof outline, source records and reproduction instructions. No explicit threshold in s is claimed, and priority remains subject to further scrutiny.

## Determinant investigations

In the two families studied in the determinant paper, the real determinant becomes very small, but the **actual primitive integer polynomial**, after removing all common coefficient factors, has a value that grows. These are mathematical obstructions for the specified infinite families. They do not exclude other proofs of irrationality. Start with that paper for the two main theorems and use its supplement for the broader investigation.

The determinant paper explains the strategy before introducing the technical formulas: what would contradict rationality, why determinants might help, and what happens when their coefficients are made integral. Each main theorem is followed by a plain-language interpretation of its quantitative conclusion and scope. The supplement develops the substantial unfinished routes as reproducible case studies. Its paired hypergeometric identities have a standalone exact replay, and each unfinished case states the missing mathematical estimate.

The [structural review](planning/narrative-review.txt), [equation-by-equation review](planning/expository-review.txt), [accessibility review](planning/accessibility-review.txt) and [main-paper assembly audit](paper/audit.txt) record earlier revisions. Their preservation and checking statements apply to those historical versions.

The [reproduction instructions and results](checks/README.md) distinguish main-paper exact/interval checks from the supplement's certificate replays and exact reconstructions. The supporting documents are:

- [Theorem inventory](planning/theorem-inventory.txt): 57 grouped claims, with precise scope, evidence, reading status and proposed placement.
- [Novelty audit](planning/novelty-audit.txt): established ingredients, closest primary comparisons, candidate contributions and unresolved priority questions.
- [Method catalog](planning/method-catalog.txt): coverage of all 12 method groups and 18 coordinated research rounds, including repaired and superseded conclusions.
- [Evidence guide](evidence/README.txt): 120 archived records and seven supplemental source files, including proofs, review records and selected calculation scripts/results. Personal filesystem paths are normalized in the published copies.
- [Source and reading record](sources/reading-record.json): versions, local PDF locations, checksums, reading scopes and unread dependencies.
- [Original paper plan](planning/paper-plan.txt) and [opening proposal](planning/draft-opening.txt): the selection-stage documents retained in the history.

The main determinant results are a uniform primitive-growth obstruction for shifted consecutive principal determinants, and an exact large-prime valuation theorem leading to a primitive-growth obstruction for difference-variable Gram determinants. The supplement distinguishes further proved limitations from incomplete constructions and finite evidence. Its positive-source cone and ordinary-moment first variation also explain why particular-family failures do not settle the broader problem. Priority of the specific statements remains subject to the recorded literature comparison and further scrutiny.

The work was prepared with Codex. The independent model review, archived model cross-checks, exact rational calculations, interval calculations and source reading are recorded separately. The independent review is a model cross-check; it is not a human mathematical audit or proof-assistant verification.

## Reproducing the results

For the even-beta paper, the [reproduction guide](even-beta-independence/README.md#reproduce-the-checks-and-build) gives the pinned dependency and build instructions. From the repository root:

```sh
python3 -m pip install -r even-beta-independence/checks/requirements.txt
python3 even-beta-independence/checks/check_local_system.py
python3 even-beta-independence/checks/check_linear_forms.py
python3 even-beta-independence/tools/build_paper.py
```

For the determinant investigation, check the packet's internal references and evidence hashes:

```sh
python3 tools/validate_packet.py
```

To regenerate the readable inventory and coverage catalog from the JSON inventory:

```sh
python3 tools/render_inventory.py
```

To build the determinant paper with LaTeX and Poppler installed:

```sh
python3 tools/build_paper.py
```

To build its companion supplement:

```sh
python3 tools/render_supplement_register.py
python3 tools/build_supplement.py
```

The repository includes the evidence and certificates needed by these instructions; the original working archive is not required. Source ledgers retain URLs, hashes and local paths with the home directory represented as `/home/research`. Third-party paper PDFs are not redistributed. The validator's optional `--archive-root` and `--check-local-pdfs` flags are for comparisons in the original environment; they account for this path normalization. Historical audit records retain the publication status and checking results at the time they were written; references to private prepublication commits are omitted.

The original material in this repository is licensed under [Apache 2.0](LICENSE); cited works retain their own terms.
