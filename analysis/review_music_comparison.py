"""Verify and summarize every frozen official MuSiC prediction without inference inflation."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import numpy as np
import pandas as pd
from run_deconvolution_pilot import sha, save_json

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    run = json.loads((source / 'run.json').read_bytes())
    assert run['status'] == 'completed' and run['prediction_rows'] == 86400
    for item in run['artifacts'] + run['source_files']:
        assert sha(ROOT / item['path']) == item['sha256'], item['path']
    files = sorted(source.glob('held_out_*/*_t*_b*.csv'))
    assert len(files) == 720
    data = pd.concat([pd.read_csv(p, dtype={'held_out': str}) for p in files], ignore_index=True)
    keys = ['held_out', 'triple_id', 'block', 'scenario', 'method', 'mixture_id']
    assert len(data) == 86400 and not data.duplicated(keys).any()
    assert data.groupby('held_out').size().eq(14400).all() and data.held_out.nunique() == 6
    assert data.groupby(['held_out', 'triple_id', 'block', 'scenario']).size().eq(120).all()
    pred_cols, truth_cols = [f'pred_{i}' for i in range(6)], [f'true_{i}' for i in range(6)]
    preds, truth = data[pred_cols].to_numpy(), data[truth_cols].to_numpy()
    assert np.isfinite(preds).all() and preds.min() >= -1e-10
    assert np.max(np.abs(preds.sum(axis=1) - 1)) < 1e-10
    assert data.groupby(['held_out', 'mixture_id'])[truth_cols].nunique().max().max() == 1
    errors = preds - truth
    mae_diff = float(np.max(np.abs(data.mae.to_numpy() - np.abs(errors).mean(axis=1))))
    rmse_diff = float(np.max(np.abs(data.rmse.to_numpy() - np.sqrt((errors**2).mean(axis=1)))))
    assert max(mae_diff, rmse_diff) < 1e-12
    complete = data.groupby(['held_out', 'triple_id', 'block', 'scenario', 'mixture_id']).method.nunique()
    assert complete.eq(2).all()
    diag = pd.concat([pd.read_csv(p) for p in sorted(source.glob('held_out_*/diagnostics.csv'))], ignore_index=True)
    assert len(diag) == 43200 and not diag.duplicated(['case_id', 'mixture_id']).any()
    maxiter = int(diag.convergence.eq('Reach Maxiter').sum())
    assert diag.convergence.str.startswith(('Converge at ', 'Reach Maxiter')).all()
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = ROOT / 'results/music_review' / run_id
    out.mkdir(parents=True)
    group = ['held_out', 'triple_id', 'scenario', 'ratio', 'method']
    scenarios = data.groupby(group)[['mae', 'rmse']].mean().mul(100).rename(columns={'mae': 'mae_pp', 'rmse': 'rmse_pp'})
    draw_sizes = data.groupby(group + ['mixture_id']).size()
    assert draw_sizes.eq(3).all()
    sd = data.groupby(group + ['mixture_id'])[pred_cols].std(ddof=1).mean(axis=1).groupby(group).mean().mul(100)
    scenarios = scenarios.join(sd.rename('conditional_sd_pp')).reset_index()
    baseline = scenarios[scenarios.ratio.eq(1)].drop(columns=['scenario', 'ratio']).rename(columns={
        'mae_pp': 'mae_balanced_pp', 'rmse_pp': 'rmse_balanced_pp', 'conditional_sd_pp': 'sd_balanced_pp'})
    strong = scenarios[scenarios.ratio.eq(10)].merge(baseline, on=['held_out', 'triple_id', 'method'], validate='many_to_one')
    strong = strong.rename(columns={'mae_pp': 'mae_strong_pp', 'rmse_pp': 'rmse_strong_pp', 'conditional_sd_pp': 'sd_strong_pp'})
    strong['mae_imbalance_change_pp'] = strong.mae_strong_pp - strong.mae_balanced_pp
    strong['sd_imbalance_change_pp'] = strong.sd_strong_pp - strong.sd_balanced_pp
    metrics = ['mae_balanced_pp', 'mae_strong_pp', 'mae_imbalance_change_pp', 'rmse_balanced_pp',
               'rmse_strong_pp', 'sd_balanced_pp', 'sd_strong_pp', 'sd_imbalance_change_pp']
    donor = strong.groupby(['held_out', 'method'])[metrics].mean().reset_index()
    aggregate = donor.groupby('method')[metrics].mean().reset_index()
    effects = pd.DataFrame(index=sorted(donor.held_out.unique()))
    effects.index.name = 'held_out'
    for method in ['music_nnls', 'music_weighted']:
        effects[f'mae_change_{method}_pp'] = donor[donor.method.eq(method)].set_index('held_out').mae_imbalance_change_pp
    effects['difference_in_differences_pp'] = effects.mae_change_music_weighted_pp - effects.mae_change_music_nnls_pp
    for allocation in ['balanced', 'strong']:
        pivot = donor.pivot(index='held_out', columns='method', values=f'mae_{allocation}_pp')
        effects[f'weighted_minus_nnls_{allocation}_pp'] = pivot.music_weighted - pivot.music_nnls
    effects = effects.reset_index()
    scenario_did = strong.pivot(index=['held_out', 'triple_id', 'scenario'], columns='method', values='mae_imbalance_change_pp').reset_index()
    scenario_did['difference_in_differences_pp'] = scenario_did.music_weighted - scenario_did.music_nnls
    blocks = data.groupby(['held_out', 'block', 'ratio', 'method']).mae.mean().mul(100).reset_index(name='mae_pp')
    block_pairs = blocks.pivot(index=['held_out', 'block', 'method'], columns='ratio', values='mae_pp').reset_index()
    block_pairs['mae_imbalance_change_pp'] = block_pairs[10] - block_pairs[1]
    block_effects = block_pairs.pivot(index=['held_out', 'block'], columns='method', values='mae_imbalance_change_pp').reset_index()
    block_effects['difference_in_differences_pp'] = block_effects.music_weighted - block_effects.music_nnls
    type_data = data[['held_out', 'method', 'ratio']].copy()
    for k in range(6):
        type_data[f'type_{k}_absolute_error_pp'] = np.abs(errors[:, k]) * 100
    type_summary = type_data.groupby(['held_out', 'method', 'ratio']).mean().reset_index()
    for name, table in [('scenario_metrics', scenarios), ('scenario_paired_effects', strong),
        ('donor_metrics', donor), ('aggregate_metrics', aggregate), ('donor_difference_in_differences', effects),
        ('scenario_difference_in_differences', scenario_did), ('seed_block_effects', block_effects),
        ('cell_type_errors', type_summary)]:
        table.to_csv(out / f'{name}.csv', index=False, lineterminator='\n')
    result = {'status': 'complete verified descriptive comparison' if maxiter == 0 else 'complete outputs; iteration-limit cases require interpretation',
        'source_run': (source / 'run.json').relative_to(ROOT).as_posix(), 'source_run_sha256': sha(source / 'run.json'),
        'script': Path(__file__).relative_to(ROOT).as_posix(), 'script_sha256': sha(__file__),
        'rows_verified': len(data), 'case_count': len(files), 'max_saved_mae_difference': mae_diff,
        'max_saved_rmse_difference': rmse_diff, 'complete_paired_grid': True,
        'diagnostics': {'maxiter_cases': maxiter, 'total_fits': len(diag),
            'nonfinite_software_variance': int((~diag.variance_finite).sum()),
            'effective_gene_range': [int(diag.n_features.min()), int(diag.n_features.max())]},
        'aggregate': aggregate.to_dict(orient='records'), 'donor_effects': effects.to_dict(orient='records'),
        'mean_difference_in_differences_pp': float(effects.difference_in_differences_pp.mean()),
        'donors_with_negative_did': int(effects.difference_in_differences_pp.lt(0).sum()),
        'specific_scenarios_with_negative_did': int(scenario_did.difference_in_differences_pp.lt(0).sum()),
        'specific_scenarios_count': len(scenario_did),
        'limitations': ['Six held-out source donors with overlapping reference donor subsets; no population inference.',
                       'Raw-count synthetic targets; not matched comparisons to normalized-cell v3 targets.',
                       'Both official branches use donor-equal factorized basis. Gene weights include residuals and cross-donor variance.',
                       'A negative difference-in-differences is not an absolute-performance recommendation.',
                       'Conditional SD estimated from only three reference draws.'],
        'artifacts': [{'path': p.relative_to(ROOT).as_posix(), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(out.iterdir())]}
    save_json(out / 'review.json', result)
    print(json.dumps({'output': str(out), 'aggregate': result['aggregate'], 'mean_did_pp': result['mean_difference_in_differences_pp'],
                      'diagnostics': result['diagnostics']}, indent=2))


if __name__ == '__main__':
    main()
