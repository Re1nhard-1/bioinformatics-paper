"""Independently check diagnostic control bases and recorded positive gene weights.

No MuSiC fit is executed. Uses the frozen sparse count matrix directly rather
than the R adapter's per-donor matrix construction and profile cache.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import traceback

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-10


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 2**20), b''):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def compare_tables(a, b, keys, metrics):
    assert not a.duplicated(keys).any() and not b.duplicated(keys).any()
    merged = a[keys + metrics].merge(b[keys + metrics], on=keys, suffixes=('_check', '_saved'), validate='one_to_one')
    assert len(merged) == len(a), (len(merged), len(a))
    differences = {}
    for metric in metrics:
        x = merged[metric + '_check'].to_numpy(float)
        y = merged[metric + '_saved'].to_numpy(float)
        assert np.array_equal(np.isnan(x), np.isnan(y)), metric
        good = np.isfinite(x) & np.isfinite(y)
        assert np.all(good | (np.isnan(x) & np.isnan(y))), metric
        delta = float(np.max(np.abs(x[good] - y[good]))) if good.any() else 0.0
        assert delta <= TOL, (metric, delta)
        differences[metric] = delta
    return {'rows': len(merged), 'maximum_differences': differences}


def weight_metrics(values):
    norm = values / values.sum()
    top_n = int(np.ceil(len(norm) / 100))
    return {'n_weight_genes': len(norm),
            'normalized_weight_effective_fraction': float(1 / (norm @ norm) / len(norm)),
            'top1pct_weight_share': float(np.partition(norm, len(norm) - top_n)[-top_n:].sum())}


def pair_metrics(left, right):
    common = left.index.intersection(right.index)
    union = left.index.union(right.index)
    assert len(common) > 1
    a = left.loc[common].to_numpy(float)
    b = right.loc[common].to_numpy(float)
    ra = rankdata(a, method='average'); ra -= ra.mean()
    rb = rankdata(b, method='average'); rb -= rb.mean()
    denominator = np.sqrt((ra @ ra) * (rb @ rb))
    rho = float((ra @ rb) / denominator) if denominator else float('nan')
    return {'intersection_n': len(common), 'union_n': len(union),
            'Jaccard': len(common) / len(union), 'spearman_shared': rho,
            'normalized_weight_tv_intersection': float(np.abs(a / a.sum() - b / b.sum()).sum() / 2),
            'balanced_weight_fraction_on_intersection': float(a.sum() / left.sum()),
            'unequal_weight_fraction_on_intersection': float(b.sum() / right.sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('selection', type=Path)
    parser.add_argument('control', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    selection_path = args.selection.resolve(); control = args.control.resolve(); out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    record = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'selection': str(selection_path), 'control': str(control), 'new_fits': 0,
              'script_sha256': sha(__file__), 'tolerance': TOL}
    try:
        selection = json.loads(selection_path.read_bytes())
        inp_path = ROOT / selection['input_manifest']
        assert sha(inp_path) == selection['input_manifest_sha256']
        inp = json.loads(inp_path.read_bytes())
        expected_hash = {item['path']: item['sha256'] for item in inp['artifacts']}
        data_dir = ROOT / inp['data_directory']
        paths = [data_dir / name for name in ['counts.npz', 'cells.tsv', 'genes.tsv', 'target_positions.npy']]
        for path in paths:
            assert sha(path) == expected_hash[path.relative_to(ROOT).as_posix()]
        for item in selection['artifacts']:
            assert sha(ROOT / item['path']) == item['sha256']
        refs = pd.read_csv(ROOT / selection['references_file'], sep='\t', keep_default_na=False)
        first = selection['triple_keys'][0]
        first_refs = refs[refs.triple_key.eq(first)]
        dominant = sorted(first_refs.reference_donors.iloc[0].split('|'))[0]
        controls = first_refs[first_refs.block.eq(0) & (first_refs.level.eq('balanced') | first_refs.dominant_donor.eq(dominant))]
        assert len(controls) == 6 and controls.reference_id.is_unique
        cells = pd.read_csv(data_dir / 'cells.tsv', sep='\t', keep_default_na=False,
                            usecols=['cell_id', 'donor', 'cell_type', 'total_counts', 'matrix_column_R'])
        genes = pd.Index((data_dir / 'genes.tsv').read_text(encoding='utf-8').splitlines())
        counts = sparse.load_npz(data_dir / 'counts.npz').tocsr()
        assert counts.shape == (len(cells), len(genes)) == (93892, 30172)
        assert np.array_equal(cells.matrix_column_R.to_numpy(), np.arange(1, len(cells) + 1))
        assert genes.is_unique and cells.cell_id.is_unique
        donor = cells.donor.to_numpy(); cell_type = cells.cell_type.to_numpy()
        basis_rows = []; ref_support = {}
        for ref in controls.itertuples(index=False):
            positions = np.asarray(ref.reference_columns_R.split('|'), dtype=int) - 1
            assert len(positions) == len(np.unique(positions)) == ref.budget * 6
            assert np.all((positions >= 0) & (positions < len(cells)))
            ref_support[ref.reference_id] = np.asarray(counts[positions].sum(axis=0)).ravel() > 0
            for typ in inp['cell_types']:
                profiles = []; sizes = []
                for source_donor in ref.reference_donors.split('|'):
                    index = positions[(donor[positions] == source_donor) & (cell_type[positions] == typ)]
                    quota = ref.major_count if source_donor == ref.dominant_donor else ref.minor_count
                    assert len(index) == quota
                    summed = np.asarray(counts[index].sum(axis=0, dtype=np.int64)).ravel()
                    total = int(summed.sum())
                    assert total == int(cells.iloc[index].total_counts.sum()) and total > 0
                    profiles.append(summed.astype(float) / total); sizes.append(total / quota)
                matrix = np.vstack(profiles)
                # NumPy sample variance provides a separate implementation from R's row residual sum.
                mean_profile = np.mean(matrix, axis=0)
                variance = np.var(matrix, axis=0, ddof=1)
                assert np.all(mean_profile[~ref_support[ref.reference_id]] == 0)
                basis_rows.append({'reference_id': ref.reference_id, 'cell_type': typ,
                    'sum_sigma': float(variance.sum()), 'theta_squared_norm': float(mean_profile @ mean_profile),
                    'mean_cell_size': float(np.mean(sizes)), 'basis_genes': int(ref_support[ref.reference_id].sum()),
                    'fixed_profile_genes': len(genes)})
        basis = pd.DataFrame(basis_rows)
        saved_basis = pd.read_csv(control / 'basis.csv', keep_default_na=False, float_precision='round_trip')
        record['basis_comparison'] = compare_tables(basis, saved_basis, ['reference_id', 'cell_type'],
            ['sum_sigma', 'theta_squared_norm', 'mean_cell_size', 'basis_genes', 'fixed_profile_genes'])
        basis.to_csv(out / 'reconstructed_control_basis.csv', index=False)

        raw = pd.read_csv(control / 'control_weights.csv.gz', keep_default_na=False, float_precision='round_trip')
        assert not raw.duplicated(['reference_id', 'target_name', 'batch_mode', 'gene']).any()
        assert np.isfinite(raw.weight).all() and raw.weight.gt(0).all()
        assert set(raw.reference_id) == set(controls.reference_id)
        target_positions = np.load(data_dir / 'target_positions.npy', allow_pickle=False)
        target_support = {}
        for name in raw.target_name.unique():
            target_donor, target_id = name.rsplit('_', 1)
            positions = target_positions[inp['donors'].index(target_donor), int(target_id)]
            assert len(positions) == len(np.unique(positions)) == 300
            assert cells.iloc[positions].donor.eq(target_donor).all()
            target_support[name] = np.asarray(counts[positions].sum(axis=0)).ravel() > 0
        vectors = {}; metric_rows = []; support_checks = 0
        for (rid, target, style), group in raw.groupby(['reference_id', 'target_name', 'batch_mode'], sort=False):
            v = pd.Series(group.weight.to_numpy(float), index=group.gene)
            expected = genes[ref_support[rid] & target_support[target]]
            assert set(v.index) == set(expected), (rid, target, style)
            support_checks += 1; vectors[rid, target, style] = v
            metric_rows.append({'reference_id': rid, 'target_name': target, 'batch_mode': style, **weight_metrics(v.to_numpy(float))})
        weight_rows = pd.DataFrame(metric_rows)
        maximum_raw_difference = 0.0
        for (rid, target), group in raw.groupby(['reference_id', 'target_name'], sort=False):
            a = vectors[rid, target, 'pooled']; b = vectors[rid, target, 'separate']
            assert set(a.index) == set(b.index)
            maximum_raw_difference = max(maximum_raw_difference, float(np.max(np.abs(a.to_numpy() - b.loc[a.index].to_numpy()))))
        assert maximum_raw_difference <= TOL
        pair_rows = []
        for ref in controls[controls.level.eq('ratio10')].itertuples(index=False):
            balanced = controls[controls.budget.eq(ref.budget) & controls.level.eq('balanced')]
            assert len(balanced) == 1
            bid = balanced.reference_id.iloc[0]
            for target in raw.target_name.unique():
                for style in ['pooled', 'separate']:
                    pair_rows.append({'reference_id': ref.reference_id, 'balanced_reference_id': bid,
                        'target_name': target, 'batch_mode': style,
                        **pair_metrics(vectors[bid, target, style], vectors[ref.reference_id, target, style])})
        pair_rows = pd.DataFrame(pair_rows)
        weight_rows.to_csv(out / 'recomputed_control_weight_metrics.csv', index=False)
        pair_rows.to_csv(out / 'recomputed_control_weight_pairs.csv', index=False)
        record['raw_weights'] = {'positive_rows': len(raw), 'support_checks_against_raw_reference_and_target': support_checks,
            'maximum_pooled_separate_difference': maximum_raw_difference,
            'summary_rows': len(weight_rows), 'paired_rows': len(pair_rows)}
        production = control.parent / first
        if (production / 'weights.csv').exists() and (production / 'weight_pairs.csv').exists():
            wm = pd.read_csv(production / 'weights.csv', keep_default_na=False, float_precision='round_trip')
            pm = pd.read_csv(production / 'weight_pairs.csv', keep_default_na=False, float_precision='round_trip', na_values=['NA'])
            record['production_weight_summary_comparison'] = compare_tables(weight_rows[weight_rows.batch_mode.eq('pooled')], wm,
                ['reference_id', 'target_name'], ['n_weight_genes', 'normalized_weight_effective_fraction', 'top1pct_weight_share'])
            record['production_weight_pair_comparison'] = compare_tables(pair_rows[pair_rows.batch_mode.eq('pooled')], pm,
                ['reference_id', 'balanced_reference_id', 'target_name'], ['intersection_n', 'union_n', 'Jaccard', 'spearman_shared',
                'normalized_weight_tv_intersection', 'balanced_weight_fraction_on_intersection', 'unequal_weight_fraction_on_intersection'])
        else:
            record['production_metric_comparison'] = 'not performed; production summaries unavailable'
        record['source_hashes'] = [{'path': str(path), 'sha256': sha(path)} for path in
            [selection_path, control / 'basis.csv', control / 'control_weights.csv.gz', *paths]]
        record['status'] = 'passed independent control-basis, raw-weight support and batching verification'
    except BaseException as error:
        record.update(status='failed; outputs retained', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save(out / 'verification.json', record)
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
