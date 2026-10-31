"""Replay the main-paper checks in scratch copies, preserving archival files.

Run with a Python interpreter containing checks/requirements.txt.
Outputs are written to checks/results/; no archived output is overwritten.
The scripts and the paper have model authorship; replay is not an independent
mathematical audit or a proof of the asymptotic theorems.
"""
import contextlib
import datetime
import hashlib
import importlib.metadata
import importlib.util
import io
import json
from pathlib import Path
import platform
import shutil
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'checks/results'
RESULTS.mkdir(parents=True, exist_ok=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def save(name, result):
    path = RESULTS / name
    path.write_text(json.dumps(result, indent=2) + '\n')
    return {'path': str(path.relative_to(ROOT)), 'sha256': digest(path)}


def main():
    started = time.monotonic()
    manifest = json.loads((ROOT / 'evidence/manifest.json').read_text())
    before = {r['path']: digest(ROOT/r['path']) for r in manifest['files']}
    assert all(before[r['path']] == r['sha256'] for r in manifest['files'])
    artifacts, reports = [], {}
    script_paths = [
        'evidence/research/collaboration/round-01/construction/checks.py',
        'evidence/research/collaboration/round-01/coordinator/kernel_identity_checks.py',
        'evidence/research/collaboration/round-06/fresh/gram_checks.py',
    ]
    with tempfile.TemporaryDirectory(prefix='beta4-manuscript-checks-') as tmp:
        tmp = Path(tmp)
        copied = []
        for i, source in enumerate(script_paths):
            folder = tmp / str(i)
            folder.mkdir()
            target = folder / Path(source).name
            shutil.copyfile(ROOT/source, target)
            copied.append(target)
        # Only the kernel/principal checks are relevant to this paper.
        # The unrelated slow-filter branch diagnostics are deliberately not run.
        kernel = module(copied[0], 'archived_kernel')
        identities = kernel.rank_one_checks()
        principal = kernel.principal_checks()
        old_kernel = json.loads((ROOT/Path(script_paths[0]).with_name('checks.json')).read_text())
        assert identities == old_kernel['rank_one_identities']
        assert [(r['nodes'], r['prime_layers']) for r in principal['cases']] == [
            (r['nodes'], r['prime_layers']) for r in old_kernel['principal']['cases']]
        for k in ['coefficient_prime_checks', 'entry_prime_checks']:
            assert principal[k] == old_kernel['principal'][k]
        artifacts.append(save('principal.json', {'rank_one_identities': identities,
                                                'principal': principal}))
        reports['principal'] = {'rank_one_identities': len(identities),
            'principal_cases': len(principal['cases']),
            'coefficient_prime_checks': principal['coefficient_prime_checks'],
            'entry_prime_checks': principal['entry_prime_checks'],
            'arb_precision_bits': principal['arb_precision_bits'],
            'exact_archived_records_match': True}
        # The appendix uses the same tau in both zero-node density bounds.
        from flint import arb
        theta = (arb(4).zeta(arb(1)/4)-arb(4).zeta(arb(3)/4))/256
        tau = 1/(arb(2).sqrt()*arb.pi())
        c0 = arb.pi()**3/192
        max_q = 1+2/(3*arb(6).sqrt())
        zero_checks = []
        for row in principal['cases']:
            nodes = row['nodes']
            if nodes[0] != 0:
                continue
            coeff = kernel.polynomial(nodes, kernel.rational_kernel(nodes))
            value = kernel.horner(coeff, theta)
            regularized = kernel.aq(kernel.cauchy_det(nodes[1:], nodes[1:])) / tau
            for a in nodes[1:]:
                regularized *= (a/(a+tau))**2
            lower = c0**len(nodes)*regularized
            upper = (c0*(max_q+5))**len(nodes)*regularized
            assert lower < value < upper
            zero_checks.append({'nodes':nodes, 'lower_log':str(lower.log()),
                                'value_log':str(value.log()), 'upper_log':str(upper.log())})
        artifacts.append(save('zero-node-bounds.json', {'cases':zero_checks}))
        reports['zero_node_appendix'] = {'interval_checks':len(zero_checks),
                                        'arb_precision_bits':kernel.ctx.prec}
        print('Principal kernel and coefficient checks passed.', flush=True)

        regularized = module(copied[1], 'archived_regularized')
        log = io.StringIO()
        with contextlib.redirect_stdout(log):
            regularized.run()
        new_identities = json.loads(copied[1].with_name('kernel_identity_checks.json').read_text())
        old_identities = json.loads((ROOT/Path(script_paths[1]).with_name('kernel_identity_checks.json')).read_text())
        assert new_identities == old_identities
        artifacts.append(save('regularized-cauchy.json', new_identities))
        reports['regularized_cauchy'] = {'exact_cases': new_identities['cases'],
                                         'exact_archived_records_match': True}
        print('Regularized Cauchy determinant identities passed.', flush=True)

        gram = module(copied[2], 'archived_gram')
        log = io.StringIO()
        with contextlib.redirect_stdout(log):
            gram.main()
        new_gram = json.loads(copied[2].with_name('gram_checks.json').read_text())
        old_gram = json.loads((ROOT/Path(script_paths[2]).with_name('gram_checks.json')).read_text())
        exact_fields = ['n', 'primitive_coefficients', 'kappa', 'prime_checks']
        assert [{k:r[k] for k in exact_fields} for r in new_gram['rows']] == [
            {k:r[k] for k in exact_fields} for r in old_gram['rows']]
        count_fields = ['moment_formula_checks', 'determinant_evaluation_checks',
                        'residue_checks', 'rank_checks']
        for k in count_fields:
            assert new_gram[k] == old_gram[k]
        artifacts.append(save('difference-gram.json', new_gram))
        reports['difference_gram'] = {k:new_gram[k] for k in count_fields}
        reports['difference_gram'].update({'matrix_sizes': list(range(1, 13)),
                                           'exact_archived_records_match': True})
        print('Difference Gram exact arithmetic and interval checks passed.', flush=True)

    extra = module(ROOT/'checks/manuscript_identities.py', 'manuscript_identities').run()
    artifacts.append(save('manuscript-identities.json', extra))
    reports['manuscript_identities'] = {
        'monic_polynomials': extra['monic_polynomials'],
        'orthogonality_equalities': extra['orthogonality_equalities'],
        'norm_checks': len(extra['exact_norms']),
        'gram_product_checks': len(extra['gram_product_checks']),
        'primitive_example': extra['primitive_example']}
    after = {r['path']: digest(ROOT/r['path']) for r in manifest['files']}
    assert before == after
    report = {
        'status': 'passed',
        'run_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'Finite exact and Arb interval checks of the main-paper formulas; not an asymptotic proof, independent model audit, human audit, priority check, or formal verification.',
        'python_version': platform.python_version(),
        'python_flint_version': importlib.metadata.version('python-flint'),
        'seconds': time.monotonic() - started,
        'source_scripts': [{'path': p, 'sha256': digest(ROOT/p)} for p in script_paths],
        'new_identity_script': {'path':'checks/manuscript_identities.py',
                                'sha256':digest(ROOT/'checks/manuscript_identities.py')},
        'all_120_archived_snapshots_unchanged': True,
        'reports': reports, 'artifacts': artifacts}
    (ROOT/'checks/reproduction.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status':'passed','reports':reports,
                      'archived_snapshots_unchanged':True}, indent=2))


if __name__ == '__main__':
    main()
