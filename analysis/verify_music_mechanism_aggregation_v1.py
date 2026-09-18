"""Verify diagnostic aggregation with direct equal-weight complete-grid means.

Independent of the reviewer's sequential aggregation implementation. Completeness
is checked before using direct means; all contrasts remain descriptive.
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

ROOT = Path(__file__).resolve().parents[1]
WEIGHT_METRICS = ['intersection_n', 'union_n', 'Jaccard', 'spearman_shared',
                 'normalized_weight_tv_intersection', 'balanced_weight_fraction_on_intersection',
                 'unequal_weight_fraction_on_intersection']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(a, b, keys, metrics):
    assert not a.duplicated(keys).any() and not b.duplicated(keys).any()
    x = a.set_index(keys).sort_index(); y = b.set_index(keys).sort_index()
    assert x.index.equals(y.index)
    result = {}
    for metric in metrics:
        aa = x[metric].to_numpy(float); bb = y[metric].to_numpy(float)
        assert np.array_equal(np.isnan(aa), np.isnan(bb))
        good = np.isfinite(aa) & np.isfinite(bb)
        assert np.all(good | (np.isnan(aa) & np.isnan(bb)))
        delta = float(np.max(abs(aa[good] - bb[good]))) if good.any() else 0.0
        assert delta <= 1e-10, (metric, delta)
        result[metric] = delta
    return {'rows': len(x), 'maximum_differences': result}


def counts(values):
    a = np.asarray(values, dtype=float)
    assert np.isfinite(a).all()
    return {'n': len(a), 'negative': int((a < 0).sum()), 'zero': int((a == 0).sum()), 'positive': int((a > 0).sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('review', type=Path)
    parser.add_argument('profile_check', type=Path)
    parser.add_argument('output', type=Path)
    a = parser.parse_args()
    run_dir = a.run.resolve(); review = a.review.resolve(); profile_check = a.profile_check.resolve(); out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out / Path(__file__).name)
    record = {'status': 'running', 'new_fits': 0, 'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Independent aggregation of complete diagnostic grid and contrary-direction checks; no biological inference',
              'script_sha256': sha(__file__)}
    try:
        run_path = run_dir / 'run.json'; run = json.loads(run_path.read_bytes())
        assert run['status'] == 'completed bounded diagnostics; descriptive review pending'
        selection_path = ROOT / run['selection']; selection = json.loads(selection_path.read_bytes())
        run_hashes = {x['path']: x['sha256'] for x in run['artifacts']}
        review_path = review / 'review.json'; report = json.loads(review_path.read_bytes())
        review_hashes = {x['path']: x['sha256'] for x in report['artifacts']}
        sources = [run_path, selection_path, review_path]
        def read_raw(triple, name):
            p = run_dir / triple / (name + '.csv')
            assert sha(p) == run_hashes[p.relative_to(ROOT).as_posix()]
            sources.append(p)
            return pd.read_csv(p, keep_default_na=False, na_values=['NA'], float_precision='round_trip')
        def read_review(name):
            p = review / (name + '.csv')
            assert sha(p) == review_hashes[p.relative_to(ROOT).as_posix()]
            sources.append(p)
            return pd.read_csv(p, keep_default_na=False, na_values=['NA'], float_precision='round_trip')
        basis = pd.concat([read_raw(t, 'basis') for t in selection['triple_keys']], ignore_index=True)
        pairs = pd.concat([read_raw(t, 'weight_pairs') for t in selection['triple_keys']], ignore_index=True)
        case_path = ROOT / selection['cases_file']; sources.append(case_path)
        cases = pd.read_csv(case_path, sep='\t', keep_default_na=False)
        bkeys = ['triple_key', 'block', 'budget', 'cell_type']
        base = basis[basis.level.eq('balanced')].set_index(bkeys)
        unequal = basis[basis.level.eq('ratio10')].copy()
        assert base.index.is_unique and len(unequal) == 2268
        unequal['sigma_trace_contrast'] = unequal.sum_sigma.to_numpy() - base.loc[pd.MultiIndex.from_frame(unequal[bkeys]), 'sum_sigma'].to_numpy()
        expanded = unequal.merge(cases[['reference_id', 'held_out', 'case_id']], on='reference_id', validate='many_to_many')
        assert len(expanded) == 9072 and not expanded.duplicated(['case_id', 'cell_type']).any()
        sigma_rows = []; weight_rows = []
        for donor in selection['donors']:
            for budget in [60, 120, 300]:
                s = expanded[expanded.held_out.eq(donor) & expanded.budget.eq(budget)]
                assert len(s) == 216 and s.cell_type.nunique() == 6 and s.triple_key.nunique() == 4
                assert s.groupby(['triple_key', 'block', 'dominant_donor']).size().eq(6).all()
                assert s.groupby(['triple_key', 'block']).size().eq(18).all()
                assert s.groupby('triple_key').size().eq(54).all()
                sigma_rows.append({'held_out': donor, 'budget': budget, 'sigma_trace_contrast': float(np.mean(s.sigma_trace_contrast.to_numpy()))})
                w = pairs[pairs.held_out.eq(donor) & pairs.budget.eq(budget)]
                assert len(w) == 288 and w.triple_key.nunique() == 4
                assert w.groupby(['triple_key', 'block', 'dominant_donor']).size().eq(8).all()
                assert w.groupby(['triple_key', 'block']).size().eq(24).all()
                assert w.groupby('triple_key').size().eq(72).all()
                assert w.groupby(['triple_key', 'block', 'dominant_donor']).mixture_id.nunique().eq(8).all()
                weight_rows.append({'held_out': donor, 'budget': budget,
                    **{m: float(np.mean(w[m].to_numpy(float))) for m in WEIGHT_METRICS}})
        sigma = pd.DataFrame(sigma_rows); weights = pd.DataFrame(weight_rows)
        record['sigma_donor_comparison'] = compare(sigma, read_review('sigma_donor_means'), ['held_out', 'budget'], ['sigma_trace_contrast'])
        record['weight_donor_comparison'] = compare(weights, read_review('weight_pairs_donor_means'), ['held_out', 'budget'], WEIGHT_METRICS)
        record['sigma_all_reference_type_directions'] = {str(b): counts(g.sigma_trace_contrast) for b, g in unequal.groupby('budget')}
        for key, frame, metric in [('sigma', sigma, 'sigma_trace_contrast'), ('weight_tv', weights, 'normalized_weight_tv_intersection')]:
            wide = frame.pivot(index='held_out', columns='budget', values=metric)
            record[key + '_donor_adjacent_directions'] = {'120_minus_60': counts(wide[120] - wide[60]), '300_minus_120': counts(wide[300] - wide[120])}
        profile_path = profile_check / 'independent_all_profile_diagnostics.csv'; sources.append(profile_path)
        pcheck = json.loads((profile_check / 'verification.json').read_bytes())
        assert pcheck['status'] == 'passed independent all-profile diagnostic reconstruction'
        profiles = pd.read_csv(profile_path, float_precision='round_trip')
        record['profile_trajectory_directions'] = {}
        for label, quotas in {'minor': [5, 10, 25], 'balanced': [20, 40, 100], 'dominant': [50, 100, 250]}.items():
            for metric in ['theta_variance_trace', 'theta_pairwise_tv', 'S_cv']:
                wide = profiles.pivot(index=['donor', 'cell_type'], columns='reference_cells', values=metric)
                first = wide[quotas[1]] - wide[quotas[0]]; second = wide[quotas[2]] - wide[quotas[1]]
                record['profile_trajectory_directions'][label + ':' + metric] = {
                    'endpoint': counts(wide[quotas[2]] - wide[quotas[0]]),
                    'first_adjacent': counts(first), 'second_adjacent': counts(second),
                    'any_adjacent_increase': int(((first > 0) | (second > 0)).sum())}
        pdonor = profiles.groupby(['donor', 'reference_cells'], sort=True)[['theta_variance_trace', 'theta_pairwise_tv', 'S_cv']].mean().reset_index()
        record['profile_donor_comparison'] = compare(pdonor, read_review('profile_dispersion_donor_means'), ['donor', 'reference_cells'], ['theta_variance_trace', 'theta_pairwise_tv', 'S_cv'])
        sigma.to_csv(out / 'direct_sigma_donor_means.csv', index=False, float_format='%.17g')
        weights.to_csv(out / 'direct_weight_donor_means.csv', index=False, float_format='%.17g')
        record.update(status='passed independent diagnostic aggregation and contrary-direction checks',
            source_hashes=[{'path': str(p), 'sha256': sha(p)} for p in sources])
    except BaseException as error:
        record.update(status='failed; outputs retained', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (out / 'verification.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: val for key, val in record.items() if key != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
