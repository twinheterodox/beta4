"""Validate writing-packet references/hashes, not mathematical statements."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--archive-root', type=Path,
                    help='Also compare snapshots with the original beta4 directory.')
parser.add_argument('--check-local-pdfs', action='store_true',
                    help='Verify the external PDFs at their recorded local paths.')
parser.add_argument('--report', type=Path, help='Write a JSON validation report.')
args = parser.parse_args()

def require(condition, message):
    if not condition:
        raise SystemExit('Validation failed: ' + message)

inventory = json.loads((ROOT / 'planning/theorem-inventory.json').read_text())
manifest = json.loads((ROOT / 'evidence/manifest.json').read_text())
original_file_count = len(manifest['files'])
extra_manifest = ROOT / 'evidence/supplement-manifest.json'
if extra_manifest.exists():
    manifest['files'] += json.loads(extra_manifest.read_text())['files']
reading = json.loads((ROOT / 'sources/reading-record.json').read_text())
claims = inventory['claims']
ids = [c['id'] for c in claims]
require(len(set(ids)) == len(ids), 'duplicate claim IDs')
source_ids = {s['id'] for s in reading['sources']}
paths = [row['path'] for row in manifest['files']]
require(len(set(paths)) == len(paths), 'duplicate manifest paths')
require(set(inventory['method_groups']) == {str(n) for n in range(1, 13)},
        'method group catalog must contain 1,...,12')
require({g for c in claims for g in c['method_groups']} == set(range(1, 13)),
        'uncovered method groups')
require({n for c in claims for n in c['rounds']} == set(range(1, 19)),
        'uncovered research rounds')

for c in claims:
    require(c['classification'] in inventory['classification_definitions'],
            f"unknown classification for {c['id']}")
    require(set(c['depends_on']) <= set(ids), f"unknown dependency for {c['id']}")
    require(set(c['literature_ids']) <= source_ids, f"unknown paper for {c['id']}")
    for p in c['proof_sources'] + c['check_sources']:
        require(p in paths, f"unmanifested evidence {p}")
    for field in ['statement', 'limitations', 'current_reading_scope', 'novelty_assessment']:
        require(bool(c[field].strip()), f"missing {field} for {c['id']}")

visited, visiting = set(), set()
by_id = {c['id']: c for c in claims}
def visit(cid):
    require(cid not in visiting, f'dependency cycle at {cid}')
    if cid in visited:
        return
    visiting.add(cid)
    for dep in by_id[cid]['depends_on']:
        visit(dep)
    visiting.remove(cid)
    visited.add(cid)
for cid in ids:
    visit(cid)

total_bytes = 0
for row in manifest['files']:
    p = ROOT / row['path']
    require(p.is_file() and not p.is_symlink(), f'missing or symlinked evidence {p}')
    require(p.resolve().is_relative_to(ROOT), f'evidence outside repository: {p}')
    raw = p.read_bytes()
    require(len(raw) == row['bytes'], f'size mismatch {p}')
    require(hashlib.sha256(raw).hexdigest() == row['sha256'], f'hash mismatch {p}')
    total_bytes += len(raw)
    if args.archive_root:
        original = args.archive_root / row['archive_relative_path']
        original_bytes = (original.read_bytes().replace(
            str(Path.home()).encode(), b'/home/research')
            if original.is_file() else None)
        require(original_bytes == raw,
                f'original differs or is missing: {original}')

pdfs_checked = 0
if args.check_local_pdfs:
    for s in reading['sources']:
        p = Path(s['local_path'])
        if p.is_relative_to('/home/research'):
            p = Path.home() / p.relative_to('/home/research')
        require(p.is_file(), f"missing PDF {s['id']}: {p}")
        require(hashlib.sha256(p.read_bytes()).hexdigest() == s['sha256'],
                f"PDF hash mismatch {s['id']}")
        pdfs_checked += 1

report = dict(status='passed', scope='Writing-packet consistency, references and byte hashes only; not mathematical verification.',
              claims=len(claims), evidence_files=len(paths), evidence_bytes=total_bytes,
              original_evidence_files=original_file_count,
              supplemental_evidence_files=len(paths)-original_file_count,
              method_groups=12, research_rounds=18, paper_sources=len(source_ids),
              original_archive_compared=bool(args.archive_root), local_pdfs_checked=pdfs_checked,
              classification_counts=dict(Counter(c['classification'] for c in claims)))
if args.report:
    args.report.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
