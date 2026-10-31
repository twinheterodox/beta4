EVIDENCE SNAPSHOT
31 October 2026

The research/ directory contains 120 byte-for-byte copies from the original
beta4 working archive. They total about 1.7 MB. manifest.json records the
original relative path, size and SHA256 of each copy. The original files
were not moved, rewritten or deleted.

The selection includes all eighteen coordinated reports, the archive map,
the underlying proofs referenced by the inventory, the existing model review
records for the two main determinant results, and selected exact/interval
check scripts and saved outputs. It is not yet a complete computational
release of every experiment in the original archive.

Start with these main-proof sources:

 research/collaboration/round-01/construction/findings.txt
   Sections 1-3: kernel, comparison and arbitrary consecutive offsets.
 research/collaboration/round-01/fresh/construction_review.txt
   Separate archived model reading of the arithmetic and branch arguments.
 research/collaboration/round-06/fresh/findings.txt
   Difference-variable construction, real bounds and exact prime law.
 research/collaboration/round-06/coordinator/gram_review.txt
   Archived model mathematical cross-check, with its stated limits.

The revised inventory is authoritative for present scope. The snapshots
deliberately preserve historical statements such as 'next target', 'paused',
and old conjectures or mismatches. Those are evidence of their original
stage, not current instructions or final mathematical conclusions.

The two selected mathematical proofs were read in the current writing pass.
Their archived calculations were not rerun in this pass. Checks of manifest
hashes and inventory references are checks of this writing packet, not new
mathematical verification. Existing model cross-checks and exact finite
calculations remain labeled as such. Nothing here is proof-assistant
verification or a priority certificate.

Selected scripts use python-flint/Arb and sometimes mpmath or SymPy. Several
write their output beside the script. To rerun one, copy it and any inputs
to a scratch directory first; do not overwrite these immutable snapshots.
The original environment is outside this repository. A pinned portable
reproduction harness is a drafting/release task, not a completed deliverable
of this selection stage.

Core check pairs preserved:

 research/collaboration/round-01/construction/checks.py and checks.json
 research/collaboration/round-01/coordinator/kernel_identity_checks.py
 research/collaboration/round-01/coordinator/kernel_identity_checks.json
 research/collaboration/round-06/fresh/gram_checks.py and gram_checks.json

The ordinary-variation correction is preserved in round-18/primitive-check.json
and primitive_check.py, together with its exact coefficient input. It must
not be replaced by values from an earlier incomplete content reduction.
