"""Render stable claim/source navigation, preserving archived claim statuses."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
inventory = json.loads((ROOT / 'planning/theorem-inventory.json').read_text())
rewritten = {
    'S01': 'Unrestricted cone and degree law proved here; fixed-degree boundary formulas remain in the archive.',
    'S03': 'Repeated-quadratic obstruction proved here.',
    'S04': 'Uniform Euler--Boole obstruction proved here.',
    'S05': 'Fixed-endpoint even-depth obstruction proved here.',
    'S06': 'Prime-level sine-normalized obstruction proved here.',
    'S09': 'Rational-kernel inequality proved here with exact certificates; beta identification uses the declared Broadhurst input.',
    'E04': 'Scoped monodromy statement proved here.',
    'P07': 'Analytic sign and rate proved here; primitive arithmetic remains open.',
    'H01': 'Euler identity and one-endpoint reduction proved here; two-endpoint system retained as an archived claim.',
    'H02': 'Strict-index denominator word proved here; source-specific sharpness retained as an archived claim.',
}

def tex(s):
    return ''.join({'&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_',
                    '{':r'\{','}':r'\}','~':r'\textasciitilde{}',
                    '^':r'\textasciicircum{}','\\':r'\textbackslash{}'}.get(c,c)
                   for c in s).replace('—','---').replace('–','--')

parts = [r'\section{Claim and source register}\label{app:register}',
    'Each entry gives the archived classification, the scope of this writing pass,',
    'and one primary navigation source. Record labels elsewhere link here;',
    'the source links open files in the local checkout. All displayed paths',
    'are relative to the repository root.',
    'The complete proof and calculation path lists, limitations and literature links',
    r'are in \path{planning/theorem-inventory.json}. The frozen archive classification',
    'does not mean every claim has been independently re-audited.']
coverage = []
for c in inventory['claims']:
    cid = c['id']
    status = ('Companion-paper proof; not repeated here.' if cid.startswith('C') else
              rewritten.get(cid, 'Catalogued with its archived scope; not promoted to a newly audited theorem.'))
    sources = c['proof_sources']
    # Lead with the current corrected account; earlier records remain listed.
    if cid == 'M03':
        primary = sources[-1]
    elif cid == 'Q04':
        primary = next(p for p in sources if p.endswith('/complete-deformation.txt'))
    else:
        primary = sources[0]
    assert (ROOT / primary).is_file(), primary
    parts += [r'\paragraph{'+cid+' -- '+tex(c['title'])+r'}\label{record-'+cid+'}',
              tex(c['classification'].replace('_',' '))+'. '+tex(status),
              r'\par\smallskip\noindent{\footnotesize\href{run:../../'+primary+r'}{\nolinkurl{'+primary+r'}}}\par']
    coverage.append({'id':cid,'archived_classification':c['classification'],
                     'supplement_scope':status,'primary_path':primary,
                     'all_proof_paths':sources,'all_check_paths':c['check_sources']})
(ROOT/'supplement/claim-register.tex').write_text('\n'.join(parts)+'\n')
(ROOT/'supplement/coverage.json').write_text(json.dumps({
    'claims':len(coverage),'method_groups':len(inventory['method_groups']),
    'scope':'Coverage and review scope, not a certification of every archived claim.',
    'entries':coverage},indent=2)+'\n')
print(f'Rendered {len(coverage)} claim entries.')
