"""Verify all donor/type/quota three-draw diagnostics without new MuSiC fits.

Reconstruct cells from frozen global reference indices; use sparse count sums,
NumPy sample variance and SciPy pairwise distances as a separate code path.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import shutil
import traceback

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[1]
QUOTAS = [5, 10, 20, 25, 40, 50, 100, 250]
TYPE_TO_SOURCE = {'B cells': 'B', 'CD14+ Monocytes': 'cM', 'CD4 T cells': 'T4',
                  'CD8 T cells': 'T8', 'FCGR3A+ Monocytes': 'ncM', 'NK cells': 'NK'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 2**20), b''):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('selection', type=Path)
    parser.add_argument('run', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    selection_path = args.selection.resolve(); run_dir = args.run.resolve(); out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out / Path(__file__).name)
    record = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'selection': str(selection_path), 'run_directory': str(run_dir), 'new_fits': 0,
              'script_sha256': sha(__file__), 'scope': 'All 672 donor/type/quota profile diagnostics; no fitted gene weights or predictions rerun'}
    try:
        run_path = run_dir / 'run.json'; run = json.loads(run_path.read_bytes())
        assert run['status'] == 'completed bounded diagnostics; descriptive review pending'
        assert sha(selection_path) == run['selection_sha256']
        selection = json.loads(selection_path.read_bytes())
        inp_path = ROOT / selection['input_manifest']
        assert sha(inp_path) == selection['input_manifest_sha256']
        inp = json.loads(inp_path.read_bytes())
        frozen_inputs = {item['path']: item['sha256'] for item in inp['artifacts']}
        frozen_outputs = {item['path']: item['sha256'] for item in run['artifacts']}
        data = ROOT / inp['data_directory']
        sources = [selection_path, run_path, inp_path]
        for filename in ['counts.npz', 'cells.tsv', 'genes.tsv']:
            path = data / filename
            assert sha(path) == frozen_inputs[path.relative_to(ROOT).as_posix()]
            sources.append(path)
        for item in selection['artifacts']:
            assert sha(ROOT / item['path']) == item['sha256']
        ref_path = ROOT / selection['references_file']; sources.append(ref_path)
        refs = pd.read_csv(ref_path, sep='\t', keep_default_na=False)
        assert len(refs) == 504 and refs.reference_id.is_unique
        cells = pd.read_csv(data / 'cells.tsv', sep='\t', keep_default_na=False,
                            usecols=['donor', 'cell_type', 'total_counts', 'matrix_column_R'])
        counts = sparse.load_npz(data / 'counts.npz').tocsr()
        genes = (data / 'genes.tsv').read_text(encoding='utf-8').splitlines()
        assert counts.shape == (len(cells), len(genes)) == (93892, 30172)
        assert len(genes) == len(set(genes))
        assert np.array_equal(cells.matrix_column_R.to_numpy(), np.arange(1, len(cells) + 1))
        source_inventory_path = ROOT / selection['inventory_file']; sources.append(source_inventory_path)
        assert sha(source_inventory_path) == selection['inventory_sha256']
        inventory = pd.read_csv(source_inventory_path, keep_default_na=False).set_index('donor_id')
        donor = cells.donor.to_numpy(); cell_type = cells.cell_type.to_numpy()
        prefix_cells = {}; repeated_prefix_checks = 0
        for ref in refs.itertuples(index=False):
            index = np.asarray(ref.reference_columns_R.split('|'), dtype=int) - 1
            assert len(index) == len(np.unique(index)) == ref.budget * 6
            assert np.all((index >= 0) & (index < len(cells)))
            for d in ref.reference_donors.split('|'):
                for typ in selection['cell_types']:
                    sub = index[(donor[index] == d) & (cell_type[index] == typ)]
                    n = ref.major_count if d == ref.dominant_donor else ref.minor_count
                    assert len(sub) == n and n in QUOTAS
                    key = d, typ, ref.block, n
                    if key in prefix_cells:
                        assert np.array_equal(prefix_cells[key], sub), key
                        repeated_prefix_checks += 1
                    else:
                        prefix_cells[key] = sub
        assert len(prefix_cells) == 14 * 6 * 3 * 8
        nested_prefix_checks = 0
        for d in selection['donors']:
            for typ in selection['cell_types']:
                for block in range(3):
                    largest = prefix_cells[d, typ, block, 250]
                    for n in QUOTAS:
                        assert np.array_equal(prefix_cells[d, typ, block, n], largest[:n])
                        nested_prefix_checks += 1
        rows = []
        for d in selection['donors']:
            for typ in selection['cell_types']:
                N = int(inventory.loc[d, TYPE_TO_SOURCE[typ]])
                for n in QUOTAS:
                    profiles = []; mean_libraries = []
                    for block in range(3):
                        index = prefix_cells[d, typ, block, n]
                        gene_totals = np.asarray(counts[index].sum(axis=0, dtype=np.int64)).ravel()
                        total = int(gene_totals.sum())
                        assert total == int(cells.iloc[index].total_counts.sum()) and total > 0
                        profile = gene_totals.astype(float) / total
                        assert abs(profile.sum() - 1) <= 1e-12
                        profiles.append(profile); mean_libraries.append(total / n)
                    matrix = np.vstack(profiles)
                    lib = np.asarray(mean_libraries, dtype=float)
                    rows.append({'donor': d, 'cell_type': typ, 'reference_cells': n,
                        'source_inventory_N': N, 'sampling_fraction': n / N, 'draws': 3, 'fixed_profile_genes': len(genes),
                        'theta_variance_trace': float(np.var(matrix, axis=0, ddof=1).sum()),
                        'theta_pairwise_tv': float(pdist(matrix, metric='cityblock').mean() / 2),
                        'S_cv': float(np.std(lib, ddof=1) / np.mean(lib))})
        computed = pd.DataFrame(rows)
        assert len(computed) == 672
        keys = ['donor', 'cell_type', 'reference_cells']
        metrics = ['theta_variance_trace', 'theta_pairwise_tv', 'S_cv']
        metadata = ['source_inventory_N', 'sampling_fraction', 'draws', 'fixed_profile_genes']
        saved_parts = []
        for triple in selection['triple_keys']:
            path = run_dir / triple / 'profile_dispersion.csv'
            assert sha(path) == frozen_outputs[path.relative_to(ROOT).as_posix()]
            sources.append(path)
            frame = pd.read_csv(path, keep_default_na=False, float_precision='round_trip')
            assert len(frame) == 144
            frame['source_triple_key'] = triple
            saved_parts.append(frame)
        saved = pd.concat(saved_parts, ignore_index=True)
        assert len(saved) == 2016
        for key, group in saved.groupby(keys, sort=False):
            assert len(group) == 3, key
            values = group[metrics + metadata].to_numpy(float)
            assert np.array_equal(values, np.broadcast_to(values[0], values.shape)), key
        unique = saved.drop_duplicates(keys)
        joined = computed.merge(unique, on=keys, suffixes=('_independent', '_R'), validate='one_to_one')
        assert len(joined) == 672
        maxima = {}
        for metric in metrics + metadata:
            x = joined[metric + '_independent'].to_numpy(float)
            y = joined[metric + '_R'].to_numpy(float)
            assert np.isfinite(x).all() and np.isfinite(y).all()
            diff = float(np.max(np.abs(x - y)))
            assert diff <= 1e-10, (metric, diff)
            maxima[metric] = diff
        computed.to_csv(out / 'independent_all_profile_diagnostics.csv', index=False, float_format='%.17g')
        joined.to_csv(out / 'profile_diagnostic_comparison.csv', index=False, float_format='%.17g')
        record.update(status='passed independent all-profile diagnostic reconstruction',
            unique_donor_type_draw_quota_profiles=len(prefix_cells), profile_summary_rows=len(computed),
            saved_rows_checked_including_duplicates=len(saved), repeated_prefix_identity_checks=repeated_prefix_checks,
            nested_prefix_identity_checks=nested_prefix_checks, maximum_absolute_differences=maxima,
            metrics_checked=metrics, source_hashes=[{'path': str(path), 'sha256': sha(path)} for path in sources])
    except BaseException as error:
        record.update(status='failed; outputs retained', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save(out / 'verification.json', record)
    print(json.dumps({key: record[key] for key in ['status', 'profile_summary_rows', 'maximum_absolute_differences', 'new_fits']}, indent=2))


if __name__ == '__main__':
    main()
