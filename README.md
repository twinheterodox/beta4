# Irrationality and linear independence of Dirichlet beta values

This repository collects research on even Dirichlet beta values. The project began with the irrationality question for

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

The papers below give explicit counts of irrational even beta values in initial segments and shifted intervals, a quantitative linear-independence bound for the growing family, and obstructions to particular determinant approaches to the individual value β(4). They are accompanied by LaTeX sources, computational checks and review records. None of the papers identifies β(4) as irrational.

## Papers

| Document | Result and contents | Source and checks |
| --- | --- | --- |
| **[Irrational even beta values in shifted intervals (v2, 20 pages)](shifted-even-beta/output/pdf/shifted-even-beta.pdf)** | For every positive odd m, an irrational even beta value with index from m+1 through 19m−1. More generally, any prescribed count with an upper endpoint linear in its starting index; improved finite counts from β(2). | [Paper guide](shifted-even-beta/README.txt) · [LaTeX](shifted-even-beta/paper/shifted-even-beta.tex) · [Checks](shifted-even-beta/checks/README.txt) |
| **[Explicit lower bounds for the number of irrational even beta values (20 pages)](even-beta-counts/output/pdf/even-beta-counts.pdf)** | At least $\lfloor\pi\sqrt{2/135}\sqrt{s/\log s}\rfloor$ irrational even beta values through β(s), for every even s ≥ 160, with sharper finite endpoints. | [Paper guide](even-beta-counts/README.txt) · [LaTeX](even-beta-counts/paper/even-beta-counts.tex) · [Checks](even-beta-counts/checks/README.txt) |
| **[At least three of β(2), β(4), …, β(364) are irrational. (16 pages)](finite-even-beta/output/pdf/finite-even-beta.pdf)** | At least three of β(2), β(4), …, β(364) are irrational. Complete proof and exact seed certificate. | [Paper guide](finite-even-beta/README.txt) · [LaTeX](finite-even-beta/paper/finite-even-beta.tex) · [Checks](finite-even-beta/checks/README.txt) |
| **[Even beta values: linear independence (19 pages)](even-beta-independence/output/pdf/even-beta-independence.pdf)** | A lower bound of $0.21\sqrt{s/\log s}$ for the rational span of β(2), β(4), …, β(s), for every sufficiently large even s. Full proof of the half-integer adaptation of Fischler's method. | [Paper guide](even-beta-independence/README.md) · [LaTeX](even-beta-independence/paper/even-beta-independence.tex) · [Checks](even-beta-independence/checks/README.txt) |
| **[Determinant obstructions (21 pages)](output/pdf/beta4-determinant-obstructions.pdf)** | Two positive determinant families fail the intended irrationality test, even after all common coefficient factors are removed. Includes complete proofs and a worked arithmetic example. | [LaTeX](paper/beta4-determinant-obstructions.tex) · [Checks](checks/README.md) · [Final reading](planning/final-reading-2026-11-02.txt) |
| **[Determinant investigation supplement (67 pages)](output/pdf/beta4-investigations.pdf)** | Self-contained case studies with constructions, proofs, remaining questions and standalone certificates; a catalog covers all 57 grouped records. | [Source guide](supplement/README.md) · [Checks](checks/SUPPLEMENT.md) · [Final reading](planning/supplement-final-reading-2026-11-02.txt) |

**Review status:** an [independent audit of the v2 shifted-interval proof](shifted-even-beta/review/independent-audit-2026-11-06.json), reported on 6 November 2026, found no substantive gap in the general argument and independently reproduced all 8,909 structural checks and 21 endpoint certificates. The PDF is a subsequent exposition of that audited plain-text proof; its assembly and additional worked-example check are [recorded separately](shifted-even-beta/review/assembly-review.json). A separate [independent audit of the count paper](even-beta-counts/review/independent-audit-2026-11-05.json), also reported on 5 November, found no substantive mathematical error or gap and independently reproduced the finite calculations. The [earlier linear-independence paper](even-beta-independence/review/independent-model-review-2026-11-03.json) received a separate model proof audit reported on 3 November 2026; its [full report](even-beta-independence/review/independent-model-review-2026-11-03.txt) records the proof checks and additional computations. The [determinant papers and accompanying code](planning/independent-model-review-2026-11-02.json) passed a separately reported model cross-check. Each report identifies its own scope; these audits are not formal verification.

See the [changelog](CHANGELOG.md) for repository updates.

## Irrational values in shifted intervals

**Versions:** v2 (current): [PDF](shifted-even-beta/output/pdf/shifted-even-beta.pdf) · [sources and checks](shifted-even-beta/) | v1: [PDF](https://github.com/twinheterodox/beta4/blob/shifted-even-beta-v1/shifted-even-beta/output/pdf/shifted-even-beta.pdf) · [sources and checks](https://github.com/twinheterodox/beta4/tree/shifted-even-beta-v1/shifted-even-beta).

Version 2 of the [shifted-interval paper](shifted-even-beta/output/pdf/shifted-even-beta.pdf) proves that, for every positive odd integer m, at least one of

$$
\beta(m+1),\beta(m+3),\ldots,\beta(19m-1)
$$

is irrational. Specific parameter choices give shorter intervals: for example, at least one of β(20), β(22), …, β(352) is irrational. The v1 endpoint was 67m−1, giving 1272 in this example.

More generally, for any integer B ≥ 1, the paper guarantees at least B irrational values among

$$
\beta(m+1),\beta(m+3),\ldots,\beta(C_Bm-1),
$$

where

$$
D_B=\sum_{\substack{1\le d\le2B-1\\d\ {\rm odd}}}\varphi(d),\qquad
C_B=1+2\left\lceil D_B\left(9+\frac{40}{9}\log(2B-1)\right)\right\rceil.
$$

Here φ is Euler's totient function and log is the natural logarithm. For each count B, the upper endpoint grows linearly with the starting index. Obtaining the same count by discarding initial values from the [earlier count paper's](even-beta-counts/output/pdf/even-beta-counts.pdf) square-root bound requires an endpoint of order m² log m. At the starting index 20, that paper's sharper finite count gives an endpoint of 3828, compared with 352 here.

The construction also gives individual parameter choices. The uniform column below uses the general formula above, which gives 19m−1 when B = 1.

| First even index | Last even index: specific parameters | Last even index: uniform bound | Irrational values (at least) |
| ---: | ---: | ---: | ---: |
| 4 | 54 | 56 | 1 |
| 6 | 92 | 94 | 1 |
| 8 | 128 | 132 | 1 |
| 20 | 352 | 360 | 1 |
| 20 | 1344 | 1614 | 2 |
| 20 | 3538 | 4350 | 3 |
| 100 | 1838 | 1880 | 1 |

Starting at β(2), the new construction gives counts 2, 3, 4, 5 and 10 through indices 70, 186, 360, 554 and 2794, respectively. At the fixed starting index 4, Krisman's earlier endpoint 36 remains shorter than 54.

The proof combines odd numerator grids and poles spaced by 2 with differentiation to remove the unwanted early beta values. A positive Gram matrix supplies the required rank for every starting index and count. The manuscript includes a worked two-grid example and explicit local coefficient scaling. Appendix B gives the sharper uniform multiplier 18.58 and explains its proximity to the limiting decay estimate.

## A uniform count from s = 160

The [count paper](even-beta-counts/output/pdf/even-beta-counts.pdf) proves

$$
I(s)\ge\left\lfloor\pi\sqrt{\frac{2}{135}}\sqrt{\frac{s}{\log s}}\right\rfloor
\qquad(s\ge160\text{ even}),
$$

where I(s) counts irrational values among β(2), β(4), …, β(s). The floor can be omitted for every even s ≥ 42,170.

The construction gives sharper counts at individual endpoints. The table compares these with the uniform bound above; the final column shows the expression inside the floor, rounded to two decimal places.

| Endpoint s | Construction: I(s) ≥ | Uniform bound: I(s) ≥ | Expression before flooring |
| ---: | ---: | ---: | ---: |
| 160 | 2 | 2 | 2.15 |
| 364 | 3 | 3 | 3.00 |
| 588 | 4 | 3 | 3.67 |
| 1030 | 5 | 4 | 4.66 |
| 1296 | 6 | 5 | 5.14 |
| 1992 | 7 | 6 | 6.19 |
| 2502 | 8 | 6 | 6.84 |
| 3260 | 9 | 7 | 7.68 |
| 3828 | 10 | 8 | 8.24 |

For example, the uniform bound guarantees three irrational values through β(588); the construction guarantees four.

A positive Gram matrix proves the required rank for every requested count; prime dilation transfers seed nonvanishing to infinitely many indices.

The [independent audit report](even-beta-counts/review/independent-audit-2026-11-05.txt) covers the general proof and independently written arithmetic checks, including all 253 finite interval comparisons. The [audit record](even-beta-counts/review/independent-audit-2026-11-05.json) identifies the 19-page PDF snapshot captured when the report was supplied. The current version adds the requested explanation of the uniform Stirling error and an audit note; these edits and their validation are [recorded separately](even-beta-counts/review/audit-response-2026-11-05.json). The audit was not formal verification or a fresh audit of all cited literature. The separate three-through-364 paper below is unchanged.

## A finite set of even beta values

The [finite-set paper](finite-even-beta/output/pdf/finite-even-beta.pdf) proves that at least three of β(2), β(4), …, β(364) are irrational. Three alternating sums allow any two chosen coefficients to be canceled; a finite seed and prime dilation ensure nonvanishing, and exact denominator and decay estimates finish the contradiction.

An [independent proof-audit report](finite-even-beta/review/independent-audit.json) for the preceding plain-text proof found no substantive mathematical gap and independently reconstructed its finite certificate. The typeset manuscript is a subsequent exposition revision and has not received a new independent audit. The reviewer identity and model were not specified for that report. The supplied exact checks, reported independent audit and source reading have separate scopes.

The count of three through 364 appears new among the sources checked. The introduction gives the relevant finite-subset precedents, including the block sizes and parameter conditions in the results of Sorokin and Huttner, and OpenAI's irrationality theorem for β(2).

Appendix A checks the finite estimates behind selected asymptotic results. The specified estimates and parameter families of Rivoal–Zudilin, Nishimoto, Fischler and the companion square-root paper do not yield the three-value count at 364. A separate exact checker reproduces the finite comparisons.

## Even beta values

The linear-independence paper proves, for every sufficiently large even integer $s$,

$$
\dim_{\mathbb Q}\mathrm{span}_{\mathbb Q}\{\beta(2),\beta(4),\ldots,\beta(s)\}
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

The work was prepared with Codex. Reported proof audits, model cross-checks, exact rational calculations, interval calculations and source reading are recorded separately with their respective scopes. No proof-assistant verification is claimed.

## Reproducing the results

For the shifted-interval paper, from the repository root:

```sh
python3 -B shifted-even-beta/checks/checks.py
python3 -B shifted-even-beta/checks/endpoints.py
python3 -B shifted-even-beta/checks/worked_example.py
python3 -B shifted-even-beta/tools/build_paper.py
```

For the count paper, from the repository root:

```sh
python3 -B even-beta-counts/checks/verify_construction.py
python3 -B even-beta-counts/checks/verify_endpoints.py
python3 -B even-beta-counts/tools/build_paper.py
```

For the finite-set paper, all certificate calculations use Python's standard library. From the repository root:

```sh
python3 finite-even-beta/checks/verify.py
python3 finite-even-beta/checks/check_asymptotic_comparisons.py
python3 finite-even-beta/tools/build_paper.py
```

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
