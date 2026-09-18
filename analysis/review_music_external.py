"""Audit every external prediction and export complete descriptive summaries."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
LEVELS=['balanced','ratio2','ratio4','ratio10'];METHODS=['music_weighted','music_nnls'];BUDGETS=[60,120,300]
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
def signs(x):
    a=np.asarray(x);return {'mean_pp':float(a.mean()),'minimum_pp':float(a.min()),'maximum_pp':float(a.max()),'positive':int((a>1e-10).sum()),'negative':int((a< -1e-10).sum()),'numerically_zero':int((abs(a)<=1e-10).sum()),'n':len(a)}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);a=p.parse_args();rp=a.run.resolve();run=json.loads(rp.read_bytes())
    assert run['status']=='completed full external grid; descriptive review pending'
    for item in run['artifacts']+run['sources']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
    inp=json.loads((ROOT/run['input_manifest']).read_bytes());data=ROOT/inp['data_directory']
    cases=pd.read_csv(ROOT/inp['case_file'],sep='\t',keep_default_na=False).set_index('case_id');targets=pd.read_csv(data/'targets.tsv',sep='\t').set_index('target_name')
    frames=[];diagnostics=[]
    for fold in sorted(run['folds'],key=lambda x:x['triple_key']):
        d=rp.parent/fold['triple_key'];f=pd.read_csv(d/'predictions.csv',keep_default_na=False);diag=pd.read_csv(d/'diagnostics.csv',keep_default_na=False)
        assert len(f)==43200 and len(diag)==21600 and (f.triple_key==fold['triple_key']).all()
        frames.append(f);diagnostics.append(diag)
    df=pd.concat(frames,ignore_index=True);dg=pd.concat(diagnostics,ignore_index=True);del frames,diagnostics
    key=['case_id','method','mixture_id'];assert len(df)==604800 and not df.duplicated(key).any() and set(df.case_id)==set(cases.index)
    assert set(df.method)==set(METHODS) and set(df.mixture_id)==set(range(60))
    assert (df.groupby(['case_id','method'],observed=True).size()==60).all()
    meta=['reference_id','held_out','triple_id','triple_key','block','budget','level','dominant_donor']
    expected=cases.loc[df.case_id,meta].reset_index(drop=True)
    assert df[meta].equals(expected)
    truth=df[[f'true_{k}' for k in range(6)]].to_numpy();pred=df[[f'pred_{k}' for k in range(6)]].to_numpy()
    assert np.isfinite(pred).all() and pred.min()>=-1e-10 and np.max(abs(pred.sum(1)-1))<1e-10
    assert np.max(abs(truth-targets.loc[df.target_name,[f'true_{k}' for k in range(6)]].to_numpy()))<1e-14
    assert np.array_equal(df.held_out.to_numpy(),targets.loc[df.target_name,'held_out'].to_numpy())
    assert np.array_equal(df.mixture_id.to_numpy(),targets.loc[df.target_name,'mixture_id'].to_numpy())
    err=pred-truth;max_metric_delta=max(float(np.max(abs(abs(err).mean(1)-df.mae))),float(np.max(abs(np.sqrt((err**2).mean(1))-df.rmse))))
    assert max_metric_delta<1e-12
    assert len(dg)==302400 and not dg.duplicated(['case_id','mixture_id']).any()
    assert set(map(tuple,dg[['case_id','mixture_id']].to_numpy()))==set(map(tuple,df[df.method==METHODS[0]][['case_id','mixture_id']].to_numpy()))
    assert (dg.reference_id.to_numpy()==cases.loc[dg.case_id,'reference_id'].to_numpy()).all()
    assert (dg.held_out.to_numpy()==cases.loc[dg.case_id,'held_out'].to_numpy()).all()
    assert (dg.target_name.to_numpy()==(dg.held_out+'_'+dg.mixture_id.map(lambda x:f'{x:02d}')).to_numpy()).all()
    assert np.isfinite(dg[['n_features','nnls_coefficient_sum','weighted_coefficient_sum']].to_numpy()).all()
    out=ROOT/'results/music_external_review'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    tables={}
    def export(name,frame):
        frame.to_csv(out/(name+'.csv'),index=False,lineterminator='\n');tables[name]=len(frame);return frame
    df['mae_pp']=df.mae*100;df['rmse_pp']=df.rmse*100
    for k in range(6):df[f'absolute_error_{k}_pp']=abs(err[:,k])*100
    metrics=['mae_pp','rmse_pp']+[f'absolute_error_{k}_pp' for k in range(6)]
    grouping=['held_out','method','budget','level']
    donor=export('donor_metrics',df.groupby(grouping,observed=True)[metrics].mean().reset_index())
    assert len(donor)==336
    overall=export('aggregate_metrics',donor.groupby(['method','budget','level'],observed=True)[metrics].mean().reset_index())
    pairkey=['held_out','triple_id','block','budget','method','mixture_id']
    baseline=df[df.level=='balanced'][pairkey+metrics].rename(columns={x:'baseline_'+x for x in metrics})
    assert not baseline.duplicated(pairkey).any()
    paired=df[df.level!='balanced'][pairkey+['dominant_donor','level']+metrics].merge(baseline,on=pairkey,validate='many_to_one')
    assert len(paired)==544320
    penalties=[]
    for m in metrics:
        c=m.replace('_pp','_penalty_pp');paired[c]=paired[m]-paired['baseline_'+m];penalties.append(c)
    dp=export('donor_penalties',paired.groupby(grouping,observed=True)[penalties].mean().reset_index())
    assert len(dp)==252
    ap=export('aggregate_penalties',dp.groupby(['method','budget','level'],observed=True)[penalties].mean().reset_index())
    # Save each dominant choice, not merely the donor averages.
    arrangement=export('reference_arrangement_penalties',paired.groupby(['held_out','triple_id','dominant_donor','method','budget','level'],observed=True)[penalties].mean().reset_index())
    block=export('donor_block_penalties',paired.groupby(['held_out','block','method','budget','level'],observed=True)[penalties].mean().reset_index())
    # Two separately formed aggregate identities verify sharing/averaging the baseline.
    wide=donor.pivot(index=['held_out','method','budget'],columns='level',values='mae_pp')
    recalculated=dp.set_index(grouping).mae_penalty_pp
    max_pair_delta=0.0
    for level in LEVELS[1:]:
        expected_diff=wide[level]-wide.balanced
        observed=recalculated.xs(level,level='level').reindex(expected_diff.index)
        max_pair_delta=max(max_pair_delta,float(np.max(abs(expected_diff-observed))))
    assert max_pair_delta<1e-10
    contrasts=[]
    for (held,method,level),sub in dp.groupby(['held_out','method','level'],observed=True):
        sub=sub.set_index('budget')
        for low,high in [(60,120),(120,300),(60,300)]:
            row={'held_out':held,'method':method,'level':level,'low_budget':low,'high_budget':high}
            for metric in penalties:row[metric.replace('_penalty_pp','_penalty_change_pp')]=float(sub.loc[high,metric]-sub.loc[low,metric])
            contrasts.append(row)
    bc=export('donor_budget_penalty_contrasts',pd.DataFrame(contrasts));assert len(bc)==252
    primary=bc[(bc.method=='music_weighted')&(bc.level=='ratio10')&(bc.low_budget==60)&(bc.high_budget==300)].copy()
    primary=primary.set_index('held_out').loc[inp['donors']].reset_index();export('primary_endpoint_all_donors',primary)
    method_penalty=dp.pivot(index=['held_out','budget','level'],columns='method',values='mae_penalty_pp').reset_index()
    method_penalty['weighted_minus_nnls_penalty_pp']=method_penalty.music_weighted-method_penalty.music_nnls
    export('donor_method_penalty_differences',method_penalty)
    absolute=donor.pivot(index=['held_out','budget','level'],columns='method',values='mae_pp').reset_index()
    absolute['weighted_minus_nnls_mae_pp']=absolute.music_weighted-absolute.music_nnls;export('donor_method_absolute_differences',absolute)
    # Conditional variation across the three frozen reference-cell draws only.
    sdkeys=['held_out','triple_id','dominant_donor','method','budget','level','mixture_id']
    groups=df.groupby(sdkeys,observed=True);assert (groups.size()==3).all()
    sd=groups[[f'pred_{k}' for k in range(6)]].std(ddof=1)*100
    sd['mean_type_sd_pp']=sd.mean(axis=1)
    ds=export('donor_conditional_sd',sd.reset_index().groupby(grouping,observed=True)[['mean_type_sd_pp']].mean().reset_index())
    export('aggregate_conditional_sd',ds.groupby(['method','budget','level'],observed=True).mean(numeric_only=True).reset_index())
    # Complete adjacent transitions and non-monotonicity at three descriptive scopes.
    transitions=[];monotonic=[]
    for scope,frame,ids,value in [
        ('donor',dp,['held_out','method','budget'],'mae_penalty_pp'),
        ('reference_arrangement',arrangement,['held_out','triple_id','dominant_donor','method','budget'],'mae_penalty_pp'),
        ('donor_block',block,['held_out','block','method','budget'],'mae_penalty_pp')]:
        w=frame.pivot(index=ids,columns='level',values=value);w['balanced']=0.0
        diff=np.diff(w[LEVELS].to_numpy(),axis=1)
        for idx,values in zip(w.index,diff):
            base=dict(zip(ids,idx));monotonic.append({'scope':scope,**base,'all_adjacent_nondecreasing':bool((values>=-1e-10).all()),'any_decrease':bool((values< -1e-10).any())})
            for j,val in enumerate(values):transitions.append({'scope':scope,**base,'from_level':LEVELS[j],'to_level':LEVELS[j+1],'mae_change_pp':float(val)})
    export('all_adjacent_allocation_changes',pd.DataFrame(transitions));mo=export('all_monotonicity_profiles',pd.DataFrame(monotonic).fillna(''))
    summaries=[]
    for scope,frame,col,group in [
        ('donor_penalty',dp,'mae_penalty_pp',['method','budget','level']),
        ('arrangement_penalty',arrangement,'mae_penalty_pp',['method','budget','level']),
        ('donor_budget_penalty_change',bc,'mae_penalty_change_pp',['method','level','low_budget','high_budget']),
        ('method_penalty_difference',method_penalty,'weighted_minus_nnls_penalty_pp',['budget','level']),
        ('method_absolute_difference',absolute,'weighted_minus_nnls_mae_pp',['budget','level'])]:
        for idx,sub in frame.groupby(group,observed=True):summaries.append({'scope':scope,**dict(zip(group,idx if isinstance(idx,tuple) else (idx,))),**signs(sub[col])})
    export('complete_direction_summaries',pd.DataFrame(summaries).fillna(''))
    summary={'status':'passed complete external descriptive audit','run':rp.relative_to(ROOT).as_posix(),'run_sha256':sha(rp),
        'producer':{'path':Path(__file__).relative_to(ROOT).as_posix(),'sha256':sha(__file__)},'input_manifest_sha256':sha(ROOT/run['input_manifest']),
        'prediction_rows':len(df),'logical_cases':len(cases),'unique_references':int(df.reference_id.nunique()),'held_out_donors':len(inp['donors']),
        'maximum_proportion_sum_error':float(np.max(abs(pred.sum(1)-1))),'maximum_metric_recalculation_difference':max_metric_delta,
        'maximum_paired_aggregation_difference_pp':max_pair_delta,'paired_unbalanced_rows':len(paired),
        'weighted_diagnostics':{'rows':len(dg),'convergence_counts':{str(k):int(v) for k,v in dg.convergence.value_counts().items()},
            'nonfinite_variance_count':int((~dg.variance_finite.astype(bool)).sum()),'effective_gene_minimum':int(dg.n_features.min()),'effective_gene_maximum':int(dg.n_features.max())},
        'primary_weighted_I300_minus_I60_pp':signs(primary.mae_penalty_change_pp),'table_rows':tables,
        'scope':'Descriptive fixed 14 high-inventory normal/na donors, shared references and synthetic single-cell sums; no population inference, measured bulk validation, or causal between-study contrast.'}
    summary['artifacts']=[{'path':f.relative_to(ROOT).as_posix(),'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(out.glob('*.csv'))]
    save(out/'review.json',summary);print('OUTPUT',out.relative_to(ROOT).as_posix());print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
