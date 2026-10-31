# Arithmetic obstructions in approximations to β(4)

This repository is being prepared for a research paper on specific constructions that fail to prove the irrationality of

β(4) = 1 − 1/3⁴ + 1/5⁴ − 1/7⁴ + ⋯.

The proposed paper studies two positive determinant families. In each, the real determinant becomes very small, but the **actual primitive integer polynomial**, after removing all common coefficient factors, has a value that grows. These are mathematical obstructions for the specified families. They do not exclude other proofs of irrationality.

**Status:** local writing-stage repository, intended for `twinheterodox/beta4-obstructions`. No remote repository or release has been created. This is a theorem inventory and paper plan, with preserved research evidence; it is not a finished manuscript.

Start with the [paper plan and proposed theorem statements](planning/paper-plan.txt). The [draft opening](planning/draft-opening.txt) shows the intended exposition. The supporting documents are:

- [Theorem inventory](planning/theorem-inventory.txt): 57 grouped claims, with precise scope, evidence, reading status and proposed placement.
- [Novelty audit](planning/novelty-audit.txt): established ingredients, closest primary comparisons, candidate contributions and unresolved priority questions.
- [Method catalog](planning/method-catalog.txt): coverage of all 12 method groups and 18 coordinated research rounds, including repaired and superseded conclusions.
- [Evidence guide](evidence/README.txt): 120 verbatim research snapshots, including proofs, review records and selected calculation scripts/results.
- [Source and reading record](sources/reading-record.json): versions, local PDF locations, checksums, reading scopes and unread dependencies.

The main candidates are a uniform primitive-growth obstruction for shifted consecutive principal determinants, and an exact large-prime valuation theorem leading to a primitive-growth obstruction for difference-variable Gram determinants. Other proven limitations and incomplete constructions belong in a technical supplement; finite searches are labeled separately.

The work was prepared with Codex. Archived model cross-checks, exact rational calculations, interval calculations, source reading and formal verification are distinguished throughout. There has been no proof-assistant verification of these results. The current writing pass did not commission another independent model audit.

To check the packet's internal references and evidence hashes:

```sh
python3 tools/validate_packet.py
```

To regenerate the readable inventory and coverage catalog from the JSON inventory:

```sh
python3 tools/render_inventory.py
```

The complete original working archive remains in the sibling `beta4/` directory. It has not been moved into or replaced by this repository. Third-party PDFs remain in `beta4/literature/papers/`; their contents are not included here. The original material in this repository is licensed under Apache 2.0; cited works retain their own terms.
