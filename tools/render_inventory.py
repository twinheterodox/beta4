"""Render the readable inventory and coverage catalog; standard library only."""
from collections import Counter
import json
from pathlib import Path
import textwrap

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / 'planning/theorem-inventory.json').read_text())
claims = data['claims']

def paragraph(label, text):
    return textwrap.fill(f'{label}: {text}', width=86, break_long_words=False,
                         break_on_hyphens=False, subsequent_indent='  ')

out = [
    'THEOREM AND CLAIM INVENTORY',
    f"Prepared {data['prepared']}; {len(claims)} grouped entries.",
    '', data['scope'], '',
    'The classifications describe the strength claimed in the archived proof.',
    'They do not certify that every underlying proof has been re-audited.',
    'Read the individual scope, reading status and novelty assessment.',
    'Paths below are relative to the repository root.', '',
    'CLASSIFICATIONS', '===============',
]
for key, definition in data['classification_definitions'].items():
    out.extend([paragraph(key, definition), ''])
for c in claims:
    title = f"{c['id']} — {c['title']}"
    out += [title, '=' * len(title),
            paragraph('Class', c['classification']),
            paragraph('Placement', c['placement']),
            paragraph('Statement', c['statement']),
            paragraph('Limits', c['limitations']),
            paragraph('Novelty', c['novelty_assessment']),
            paragraph('Reading in this pass', c['current_reading_scope']),
            paragraph('Verification status', c['verification'])]
    if c['depends_on']:
        out.append(paragraph('Dependencies', ', '.join(c['depends_on'])))
    if c['literature_ids']:
        out.append(paragraph('Literature', ', '.join(c['literature_ids'])))
    out.append('Proof/navigation sources:')
    out.extend('  ' + p for p in c['proof_sources'])
    if c['check_sources']:
        out.append('Preserved calculation sources/results:')
        out.extend('  ' + p for p in c['check_sources'])
    out += ['']
(ROOT / 'planning/theorem-inventory.txt').write_text('\n'.join(out).rstrip() + '\n')

out = ['METHOD CATALOG AND COVERAGE', f"Prepared {data['prepared']}", '',
       'This catalog covers all twelve method groups and eighteen coordinated',
       'rounds. It includes older exploratory work through P01 and O01; those',
       'are navigation entries, not collections of newly certified theorems.', '',
       'The original archive remains the complete research record. The present',
       'snapshot is selected evidence for writing, not an exhaustive copy.', '',
       'BY METHOD', '=========']
for gid, name in data['method_groups'].items():
    out += ['', f'{gid}. {name}']
    for c in claims:
        if int(gid) in c['method_groups']:
            out.append(textwrap.fill(f"  {c['id']} {c['title']} [{c['classification']}]",
                                     width=86, subsequent_indent='     ',
                                     break_long_words=False, break_on_hyphens=False))
out += ['', 'BY RESEARCH ROUND', '=================']
for n in range(1, 19):
    matched = [c['id'] for c in claims if n in c['rounds']]
    out += ['', f"Round {n:02d}: {', '.join(matched)}",
            f'  evidence/research/coordinated-round-{n:02d}.txt']
out += ['', 'SUPERSESSIONS AND SCOPE CORRECTIONS', '==================================']
for row in data['supersessions']:
    out += ['', paragraph('Historical issue', row['old']), paragraph('Current statement', row['current'])]
out += ['', 'COUNTS', '======']
for key, count in sorted(Counter(c['classification'] for c in claims).items()):
    out.append(f'{key}: {count}')
out += ['', 'Main paper: C01-C08. Context appendix: S01-S02 and P06-P07.',
        'Other items are supplemental candidates or catalog-only records.',
        'A class count is not a count of proved original publishable theorems.']
(ROOT / 'planning/method-catalog.txt').write_text('\n'.join(out) + '\n')
print(f'Rendered {len(claims)} claims, 12 method groups and 18 rounds.')
