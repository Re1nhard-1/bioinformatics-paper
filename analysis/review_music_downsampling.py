"""Summarize the complete locked available-inventory policy contrast."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import numpy as np
import pandas as pd
from run_deconvolution_pilot import sha, save_json

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/music_comparison/20260917T090859412332Z'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    args=parser.parse_args()
    folder=args.run.resolve()
    run=json.loads((folder/'run.json').read_bytes())
    assert run['status']=='completed' and run['prediction_rows']==21600
    for rec in run['source_files']+run['artifacts']:
        assert sha(ROOT/rec['path'])==rec['sha256'],rec['path']
    new=pd.concat([pd.read_csv(p,dtype={'held_out':str}) for p in sorted(folder.glob('held_out_*/*_thin5.csv'))],ignore_index=True)
    keys=['held_out','triple_id','block','method','mixture_id']
    assert len(new)==21600 and not new.duplicated(keys).any()
    assert new.groupby(['held_out','triple_id','block']).size().eq(120).all()
    assert new.groupby('held_out').size().eq(3600).all()
    assert set(new.policy)=={'thin5'} and set(new.reference_cells)=={90}
    preds=[f'pred_{k}' for k in range(6)]; truths=[f'true_{k}' for k in range(6)]
    x=new[preds].to_numpy(); y=new[truths].to_numpy()
    assert np.isfinite(x).all() and x.min()>=-1e-10 and np.max(np.abs(x.sum(axis=1)-1))<1e-10
    mae_diff=float(np.max(np.abs(new.mae-np.abs(x-y).mean(axis=1))))
    rmse_diff=float(np.max(np.abs(new.rmse-np.sqrt(((x-y)**2).mean(axis=1)))))
    assert max(mae_diff,rmse_diff)<1e-12
    oldrun=json.loads((OLD/'run.json').read_bytes())
    oldhash={r['path']:r['sha256'] for r in oldrun['artifacts']}
    frames=[]
    for p in sorted(OLD.glob('held_out_*/*_t*_b*.csv')):
        assert sha(p)==oldhash[p.relative_to(ROOT).as_posix()]
        frames.append(pd.read_csv(p,dtype={'held_out':str}))
    old=pd.concat(frames,ignore_index=True)
    keep=old[old.ratio.eq(10)].copy()
    assert len(keep)==64800
    merge=keep.merge(new[keys+['mae','rmse']+truths],on=keys,suffixes=('_keep','_thin'),validate='many_to_one')
    assert len(merge)==64800
    assert np.array_equal(merge[[c+'_keep' for c in truths]].to_numpy(),merge[[c+'_thin' for c in truths]].to_numpy())
    merge['mae_thin_minus_keep_pp']=100*(merge.mae_thin-merge.mae_keep)
    merge['rmse_thin_minus_keep_pp']=100*(merge.rmse_thin-merge.rmse_keep)
    contrasts=['mae_thin_minus_keep_pp','rmse_thin_minus_keep_pp']
    paired=merge.groupby(['held_out','method'])[contrasts].mean().reset_index()
    arrangements=merge.groupby(['held_out','triple_id','scenario','method'])[contrasts].mean().reset_index()
    assert len(arrangements)==360
    blocks=merge.groupby(['held_out','block','method'])[contrasts].mean().reset_index()
    old['policy']=np.where(old.ratio.eq(10),'keep60','expanded60')
    new['policy']='thin15'
    combined=pd.concat([old,new],ignore_index=True)
    metric_group=['held_out','method','policy']
    donor=combined.groupby(metric_group)[['mae','rmse']].mean().mul(100).rename(columns={'mae':'mae_pp','rmse':'rmse_pp'})
    grouped=['held_out','method','policy','triple_id','scenario','mixture_id']
    assert combined.groupby(grouped).size().eq(3).all()
    conditional=combined.groupby(grouped)[preds].std(ddof=1).mean(axis=1).groupby(metric_group).mean().mul(100)
    donor=donor.join(conditional.rename('conditional_sd_pp')).reset_index()
    aggregate=donor.groupby(['method','policy'])[['mae_pp','rmse_pp','conditional_sd_pp']].mean().reset_index()
    check=donor.pivot(index=['held_out','method'],columns='policy',values='mae_pp')
    check['difference']=check.thin15-check.keep60
    direct=paired.set_index(['held_out','method']).mae_thin_minus_keep_pp
    assert float(np.max(np.abs(check['difference']-direct)))<1e-12
    type_errors=combined[metric_group].copy()
    for k in range(6):type_errors[f'type_{k}_absolute_error_pp']=100*np.abs(combined[preds[k]]-combined[truths[k]])
    type_errors=type_errors.groupby(metric_group).mean().reset_index()
    diagnostics=pd.concat([pd.read_csv(p) for p in folder.glob('held_out_*/diagnostics.csv')],ignore_index=True)
    assert len(diagnostics)==10800 and not diagnostics.duplicated(['case_id','mixture_id']).any()
    assert diagnostics.convergence.str.startswith(('Converge at ','Reach Maxiter')).all()
    diag={'weighted_fits':len(diagnostics),'iteration_limit':int(diagnostics.convergence.eq('Reach Maxiter').sum()),
          'nonfinite_variance':int((~diagnostics.variance_finite).sum()),
          'effective_gene_min':int(diagnostics.n_features.min()),'effective_gene_max':int(diagnostics.n_features.max())}
    details={}
    for method in ['music_weighted','music_nnls']:
        a=arrangements[arrangements.method.eq(method)].mae_thin_minus_keep_pp
        b=blocks[blocks.method.eq(method)].mae_thin_minus_keep_pp
        d=paired[paired.method.eq(method)].mae_thin_minus_keep_pp
        details[method]={'mean_mae_thin_minus_keep_pp':float(d.mean()),'donor_min_pp':float(d.min()),'donor_max_pp':float(d.max()),
            'donors_thinning_worse':int(d.gt(0).sum()),'donors_thinning_better':int(d.lt(0).sum()),
            'arrangements_thinning_worse':int(a.gt(0).sum()),'arrangements_thinning_better':int(a.lt(0).sum()),
            'arrangements':len(a),'blocks_thinning_worse':int(b.gt(0).sum()),'blocks_thinning_better':int(b.lt(0).sum()),'blocks':len(b)}
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results/music_downsampling_review'/run_id; out.mkdir(parents=True)
    for name,frame in [('donor_metrics',donor),('aggregate_metrics',aggregate),('donor_policy_effects',paired),
                       ('arrangement_policy_effects',arrangements),('block_policy_effects',blocks),('cell_type_errors',type_errors)]:
        frame.to_csv(out/f'{name}.csv',index=False,lineterminator='\n')
    report={'status':'complete descriptive review','source_run':(folder/'run.json').relative_to(ROOT).as_posix(),
            'source_run_sha256':sha(folder/'run.json'),'script_sha256':sha(__file__),
            'new_prediction_rows':len(new),'reused_prediction_rows':len(old),'paired_rows':len(merge),
            'saved_mae_max_difference':mae_diff,'saved_rmse_max_difference':rmse_diff,
            'diagnostics':diag,'aggregate':aggregate.to_dict(orient='records'),'policy_effects':details,
            'interpretation':'Positive thin-minus-keep means worse error from thinning; expanded60 cannot be obtained from the simulated keep60 inventory.',
            'limitations':['Same six SLE source donors and synthetic targets as earlier exploration; not independent validation.',
                          'The first-five reference is shared by three dominant-donor comparisons, not additional replicates.',
                          'Thinning discards75% cells and changes workflow gene support/precision; not isolated balance effect.',
                          'Three conditional reference draws; no population tests or confidence intervals.'],
            'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in out.iterdir()]}
    save_json(out/'review.json',report)
    print(json.dumps({'output':out.relative_to(ROOT).as_posix(),'policy_effects':details,'aggregate':report['aggregate'],'diagnostics':diag},indent=2))

if __name__=='__main__':main()
