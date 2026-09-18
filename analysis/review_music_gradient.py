"""Audit all gradient predictions and report complete four-level paired effects."""
from pathlib import Path
from datetime import datetime, timezone
import argparse,json
import numpy as np
import pandas as pd
from run_deconvolution_pilot import sha,save_json
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/music_comparison/20260917T090859412332Z'
TOL=1e-10

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);a=p.parse_args()
    folder=a.run.resolve();run=json.loads((folder/'run.json').read_bytes())
    assert run['status']=='completed' and run['prediction_rows']==129600
    for item in run['source_files']+run['artifacts']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
    info=json.loads((ROOT/run['input_manifest']).read_bytes())
    cases=pd.read_csv(ROOT/info['case_file'],sep='\t',dtype={'held_out':str,'dominant_donor':str}).query("policy=='gradient'")
    new=pd.concat([pd.read_csv(p,dtype={'held_out':str,'dominant_donor':str}) for p in sorted(folder.glob('held_out_*/*_g*_d*.csv'))],ignore_index=True)
    assert len(new)==129600 and set(new.case_id)==set(cases.case_id)
    assert not new.duplicated(['case_id','method','mixture_id']).any()
    assert new.groupby(['case_id','method']).mixture_id.nunique().eq(60).all()
    assert new.groupby('case_id').size().eq(120).all()
    new['level']=new.major_count.astype(int);new['dominant']=new.dominant_donor
    assert set(new.level)=={30,40} and set(new.reference_cells)=={360} and set(new.policy)=={'gradient'}
    check=new[['case_id','held_out','triple_id','block','allocation','major_count','minor_count','dominant_donor']].drop_duplicates().sort_values('case_id').reset_index(drop=True)
    expected=cases[check.columns].sort_values('case_id').reset_index(drop=True)
    pd.testing.assert_frame_equal(check,expected,check_dtype=False)
    oldrun=json.loads((OLD/'run.json').read_bytes());oldhash={x['path']:x['sha256'] for x in oldrun['artifacts']}
    frames=[]
    for path in sorted(OLD.glob('held_out_*/*_t*_b*.csv')):
        assert sha(path)==oldhash[path.relative_to(ROOT).as_posix()]
        frames.append(pd.read_csv(path,dtype={'held_out':str}))
    old=pd.concat(frames,ignore_index=True);assert len(old)==86400
    old['level']=np.where(old.scenario.eq('balanced'),20,50)
    old['dominant']=old.scenario.map(lambda x:'balanced' if x=='balanced' else x.rsplit('_',1)[1])
    combined=pd.concat([old,new],ignore_index=True)
    keys=['held_out','triple_id','block','method','mixture_id']
    assert not combined.duplicated(keys+['level','dominant']).any()
    preds=[f'pred_{k}' for k in range(6)];truths=[f'true_{k}' for k in range(6)]
    x=combined[preds].to_numpy();y=combined[truths].to_numpy()
    assert np.isfinite(x).all() and np.isfinite(y).all() and x.min()>=-1e-10
    assert np.max(np.abs(x.sum(axis=1)-1))<1e-10 and np.max(np.abs(y.sum(axis=1)-1))<1e-12
    mae_diff=float(np.max(np.abs(combined.mae-np.abs(x-y).mean(axis=1))))
    rmse_diff=float(np.max(np.abs(combined.rmse-np.sqrt(((x-y)**2).mean(axis=1)))))
    assert max(mae_diff,rmse_diff)<1e-12
    balanced=combined[combined.level.eq(20)][keys+['mae','rmse']+truths]
    merged=combined[combined.level.ne(20)].merge(balanced,on=keys,suffixes=('','_balanced'),validate='many_to_one')
    assert len(merged)==194400
    assert np.array_equal(merged[truths].to_numpy(),merged[[c+'_balanced' for c in truths]].to_numpy())
    merged['mae_penalty_pp']=100*(merged.mae-merged.mae_balanced)
    merged['rmse_penalty_pp']=100*(merged.rmse-merged.rmse_balanced)
    metrics=['mae_penalty_pp','rmse_penalty_pp']
    effects=merged.groupby(['held_out','method','level'])[metrics].mean().reset_index()
    arrangements=merged.groupby(['held_out','triple_id','dominant','method','level'])[metrics].mean().reset_index()
    blocks=merged.groupby(['held_out','block','method','level'])[metrics].mean().reset_index()
    donor_group=['held_out','method','level']
    donor=combined.groupby(donor_group)[['mae','rmse']].mean().mul(100).rename(columns={'mae':'mae_pp','rmse':'rmse_pp'})
    draw_group=['held_out','method','level','triple_id','dominant','mixture_id']
    assert combined.groupby(draw_group).size().eq(3).all()
    sd=combined.groupby(draw_group)[preds].std(ddof=1).mean(axis=1).groupby(donor_group).mean().mul(100)
    donor=donor.join(sd.rename('conditional_sd_pp')).reset_index()
    absolute=donor.pivot(index=['held_out','method'],columns='level',values='mae_pp')
    for level in [30,40,50]:
        paired=effects[effects.level.eq(level)].set_index(['held_out','method']).mae_penalty_pp
        assert np.max(np.abs(absolute[level]-absolute[20]-paired))<1e-12
    aggregate=donor.groupby(['method','level'])[['mae_pp','rmse_pp','conditional_sd_pp']].mean().reset_index()
    penalty=effects.groupby(['method','level'])[metrics].mean().reset_index()
    did=effects.pivot(index=['held_out','level'],columns='method',values='mae_penalty_pp').reset_index()
    did['weighted_minus_nnls_penalty_pp']=did.music_weighted-did.music_nnls
    type_errors=combined[donor_group].copy()
    for k in range(6):type_errors[f'type_{k}_absolute_error_pp']=100*np.abs(combined[preds[k]]-combined[truths[k]])
    type_errors=type_errors.groupby(donor_group).mean().reset_index()
    trend_records=[];adjacent_tables=[];profile_tables=[]
    for scope,frame,index in [('donor',effects,['held_out','method']),('arrangement',arrangements,['held_out','triple_id','dominant','method']),('block',blocks,['held_out','block','method'])]:
        wide=frame.pivot(index=index,columns='level',values='mae_penalty_pp')
        wide[20]=0.;wide=wide[[20,30,40,50]]
        assert wide.notna().all().all()
        delta=wide.diff(axis=1)[[30,40,50]]
        profile=wide.reset_index().rename(columns={v:f'penalty_{v}_pp' for v in [20,30,40,50]})
        profile['nondecreasing']=delta.ge(-TOL).all(axis=1).to_numpy()
        profile['scope']=scope;profile_tables.append(profile)
        for level,start in [(30,20),(40,30),(50,40)]:
            rows=delta[[level]].reset_index().rename(columns={level:'adjacent_mae_change_pp'})
            rows['from_level']=start;rows['to_level']=level;rows['scope']=scope;adjacent_tables.append(rows)
            for method,part in rows.groupby('method'):
                val=part.adjacent_mae_change_pp
                trend_records.append({'scope':scope,'method':method,'from_level':start,'to_level':level,'n_profiles':len(val),
                    'positive':int(val.gt(TOL).sum()),'negative':int(val.lt(-TOL).sum()),'ties':int(val.abs().le(TOL).sum()),
                    'mean_change_pp':float(val.mean()),'min_change_pp':float(val.min()),'max_change_pp':float(val.max())})
    profiles=pd.concat(profile_tables,ignore_index=True);adjacent=pd.concat(adjacent_tables,ignore_index=True)
    profile_counts=profiles.groupby(['scope','method']).nondecreasing.agg(['sum','count']).reset_index().rename(columns={'sum':'nondecreasing','count':'total'})
    diagnostics=pd.concat([pd.read_csv(p) for p in folder.glob('held_out_*/diagnostics.csv')],ignore_index=True)
    assert len(diagnostics)==64800 and not diagnostics.duplicated(['case_id','mixture_id']).any()
    assert diagnostics.convergence.str.startswith(('Converge at ','Reach Maxiter')).all()
    diag={'weighted_fits':len(diagnostics),'iteration_limit':int(diagnostics.convergence.eq('Reach Maxiter').sum()),'nonfinite_variance':int((~diagnostics.variance_finite).sum()),'effective_gene_min':int(diagnostics.n_features.min()),'effective_gene_max':int(diagnostics.n_features.max())}
    out=ROOT/'results/music_gradient_review'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    tables={'donor_metrics':donor,'aggregate_metrics':aggregate,'donor_penalties':effects,'aggregate_penalties':penalty,'arrangement_penalties':arrangements,'block_penalties':blocks,'donor_difference_in_differences':did,'cell_type_errors':type_errors,'adjacent_differences':adjacent,'profile_monotonicity':profiles,'profile_counts':profile_counts,'trend_summary':pd.DataFrame(trend_records)}
    for name,frame in tables.items():frame.to_csv(out/f'{name}.csv',index=False,lineterminator='\n')
    record={'status':'complete descriptive gradient review','source_run':(folder/'run.json').relative_to(ROOT).as_posix(),'source_run_sha256':sha(folder/'run.json'),'script_sha256':sha(__file__),
        'new_prediction_rows':len(new),'reused_prediction_rows':len(old),'unique_combined_prediction_rows':len(combined),'paired_rows':len(merged),
        'saved_mae_max_difference':mae_diff,'saved_rmse_max_difference':rmse_diff,'numerical_tie_tolerance_pp':TOL,'diagnostics':diag,
        'aggregate':aggregate.to_dict('records'),'aggregate_penalties':penalty.to_dict('records'),
        'mean_difference_in_differences':did.groupby('level').weighted_minus_nnls_penalty_pp.mean().to_dict(),
        'trend_summary':trend_records,'profile_counts':profile_counts.to_dict('records'),
        'limitations':['Same selected six SLE donors and synthetic targets; not external validation.','Shared reference donors and compositions induce dependence; no population P/CI.','One budget60/type and three random draws; no general thresholds or assumed monotonic curve.','Official workflow support can change; not isolated causal mechanism.'],
        'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in out.iterdir()]}
    save_json(out/'review.json',record)
    print(json.dumps({'output':out.relative_to(ROOT).as_posix(),'aggregate':record['aggregate'],'profile_counts':record['profile_counts'],'diagnostics':diag},indent=2))
if __name__=='__main__':main()
