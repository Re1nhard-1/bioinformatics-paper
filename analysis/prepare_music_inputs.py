"""Export raw-count PBMC inputs with the exact saved v3.1 cell identities."""
from __future__ import annotations
from pathlib import Path
from datetime import datetime, timezone
import gzip
import json
import warnings
import numpy as np
import pandas as pd
from scipy import io, sparse
from run_deconvolution_pilot import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'results/kang_input/20260916T220230568916Z/input_manifest.json'
SOURCE = ROOT / 'results/deconvolution_v3_1/20260916T221536546175Z'


def write_mtx(path, matrix):
    with gzip.GzipFile(filename=str(path), mode='wb', mtime=0) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            io.mmwrite(handle, matrix, field='integer', symmetry='general')


def main():
    manifest = json.loads(INPUT.read_bytes())
    source = json.loads((SOURCE / 'run.json').read_bytes())
    assert source['prediction_rows'] == 2419200
    for item in source['artifacts']:
        assert sha(ROOT / item['path']) == item['sha256'], item['path']
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = ROOT / 'data/processed/music' / run_id
    audit = ROOT / 'results/music_input' / run_id
    out.mkdir(parents=True); audit.mkdir(parents=True)
    donors, types = manifest['selected_donors'], manifest['selected_cell_types']
    genes = json.loads((ROOT / manifest['genes_path']).read_bytes())
    assert sha(ROOT / manifest['genes_path']) == manifest['genes_sha256']
    matrices, obs_by_donor, all_obs, all_x = {}, {}, [], []
    for rec in manifest['records']:
        donor = rec['donor']
        if donor not in donors:
            continue
        for name in ['matrix', 'cells']:
            assert sha(ROOT / rec[f'{name}_path']) == rec[f'{name}_sha256']
        x = sparse.load_npz(ROOT / rec['matrix_path']).tocsr()
        obs = pd.read_csv(ROOT / rec['cells_path'], dtype={'donor': str})
        assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), obs.umi_total)
        matrices[donor], obs_by_donor[donor] = x, obs
        keep = obs.cell_type.isin(types).to_numpy()
        all_x.append(x[keep]); all_obs.append(obs.loc[keep].copy())
    x = sparse.vstack(all_x, format='csr')
    obs = pd.concat(all_obs, ignore_index=True)
    assert x.shape == (11066, 35635) and not obs.cell_id.duplicated().any()
    obs.insert(0, 'matrix_column_R', np.arange(1, len(obs) + 1))
    obs.to_csv(out / 'cells.tsv', sep='\t', index=False, lineterminator='\n')
    pd.Series(genes).to_csv(out / 'genes.tsv', sep='\t', index=False, header=False, lineterminator='\n')
    write_mtx(out / 'counts.mtx.gz', x.T.tocoo())
    index = dict(zip(obs.cell_id, obs.matrix_column_R))
    cases, markers, target_info, target_blocks, checks = [], [], [], [], []
    for donor in donors:
        folder = SOURCE / f'held_out_{donor}'
        z = np.load(folder / 'targets.npz', allow_pickle=True)
        ids = z['cell_row_indices']
        old_obs, raw = obs_by_donor[donor], matrices[donor]
        assert np.array_equal(old_obs.cell_id.to_numpy()[ids], z['cell_ids'])
        assert old_obs.iloc[ids.ravel()].cell_type.isin(types).all()
        truth_counts = np.stack([[np.sum(old_obs.cell_type.to_numpy()[row] == ct) for ct in types] for row in ids])
        assert np.array_equal(truth_counts / 300, z['truth'])
        target = np.stack([np.asarray(raw[row].sum(axis=0)).ravel() for row in ids])
        assert np.equal(target, np.floor(target)).all() and target.min() >= 0
        library = np.asarray(raw.sum(axis=1)).ravel()
        reconstructed = np.stack([np.asarray((sparse.diags(10000 / library[row]) @ raw[row]).mean(axis=0)).ravel() for row in ids])
        difference = float(np.max(np.abs(reconstructed - z['expression'])))
        assert difference < 1e-8
        target_blocks.append(target)
        for m in range(60):
            row = {'target_name': f'{donor}_{m:02d}', 'held_out': donor, 'mixture_id': m, 'total_counts': int(target[m].sum())}
            row.update({f'true_{i}': float(z['truth'][m, i]) for i in range(len(types))})
            target_info.append(row)
        draws = pd.read_csv(folder / 'reference_samples.csv.gz', dtype={'held_out': str, 'reference_donor': str})
        chosen = draws[draws['repeat'].eq(0) & ~draws.scenario.str.startswith('ratio4_')]
        assert len(chosen) == 10 * 3 * 4 * 18
        for (triple, block, scenario), group in chosen.groupby(['triple_id', 'block', 'scenario'], sort=True):
            selected = [c for cell_ids in group.cell_ids for c in cell_ids.split('|')]
            assert len(selected) == len(set(selected)) == 360
            rows = np.array([index[c] - 1 for c in selected])
            selected_obs = obs.iloc[rows]
            assert donor not in set(selected_obs.donor) and selected_obs.donor.nunique() == 3
            assert selected_obs.groupby('cell_type').size().eq(60).all()
            actual = selected_obs.groupby(['donor', 'cell_type']).size().to_numpy()
            expected = [20] * 18 if scenario == 'balanced' else [5] * 12 + [50] * 6
            assert np.array_equal(np.sort(actual), np.sort(expected))
            cases.append({'case_id': f'{donor}_t{triple:02d}_b{block}_{scenario}', 'held_out': donor,
                'triple_id': int(triple), 'block': int(block), 'scenario': scenario,
                'reference_columns_R': '|'.join(map(str, rows + 1)),
                'source_repeat': 0, 'ratio': 1 if scenario == 'balanced' else 10})
        triples = json.loads((folder / 'reference_triples.json').read_bytes())
        for triple in triples:
            for gene in triple['feature_indices']:
                markers.append({'held_out': donor, 'triple_id': triple['triple_id'], 'gene_id': genes[gene]})
        checks.append({'held_out': donor, 'target_rows': 60, 'normalized_target_reconstruction_max_difference': difference,
                       'raw_pseudobulk_minimum_library': int(target.sum(axis=1).min()), 'raw_pseudobulk_maximum_library': int(target.sum(axis=1).max())})
    write_mtx(out / 'targets.mtx.gz', sparse.coo_matrix(np.vstack(target_blocks).T))
    for name, rows in [('targets', target_info), ('cases', cases), ('markers', markers)]:
        pd.DataFrame(rows).to_csv(out / f'{name}.tsv', sep='\t', index=False, lineterminator='\n')
    assert len(cases) == 720 and len(target_info) == 360
    record = {'run_id': run_id, 'status': 'raw count inputs prepared; no MuSiC predictions',
        'script': Path(__file__).relative_to(ROOT).as_posix(), 'script_sha256': sha(__file__),
        'input_manifest': INPUT.relative_to(ROOT).as_posix(), 'input_manifest_sha256': sha(INPUT),
        'source_run': (SOURCE / 'run.json').relative_to(ROOT).as_posix(), 'source_run_sha256': sha(SOURCE / 'run.json'),
        'data_directory': out.relative_to(ROOT).as_posix(), 'donors': donors, 'cell_types': types,
        'case_count': len(cases), 'targets_per_case': 60, 'reference_cells_per_case': 360,
        'selection': 'All six held-out donors and ten reference triples; all three source blocks, repeat 0 only; balanced and all three 50:5:5 scenarios; all original genes supplied with markers=NULL. Previous marker lists exported only for provenance, not used by this MuSiC comparison.',
        'checks': checks, 'raw_counts_shape': list(x.T.shape),
        'limitations': ['Raw-count sums preserve original sampled-cell identities but change the measurement model relative to v3 normalized-cell means.',
                       'Previous marker lists are retained for provenance but not used to restrict MuSiC fitting; the official implementation selects nonzero common support.',
                       'Three selected reference draws support only descriptive Monte Carlo variability; no population inference.'],
        'artifacts': [{'path': p.relative_to(ROOT).as_posix(), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(out.iterdir())]}
    save_json(audit / 'input_manifest.json', record)
    print(json.dumps({'manifest': (audit / 'input_manifest.json').relative_to(ROOT).as_posix(), 'data': str(out), 'cases': len(cases)}, indent=2))


if __name__ == '__main__':
    main()
