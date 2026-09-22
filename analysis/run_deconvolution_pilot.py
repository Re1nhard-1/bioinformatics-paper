"""Exploratory donor-imbalance experiment on normalized-cell synthetic mixtures."""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.optimize import nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
SEED=20260917
N_MIX=60
N_CELLS=300
N_REPEATS=8
MARKERS_PER_TYPE=100

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')

def fit(reference,mixtures):
    if np.linalg.matrix_rank(reference)!=reference.shape[1]: raise ValueError('Rank deficient signature')
    guesses=[]; sums=[]
    for y in mixtures:
        coef,_=nnls(reference,y,maxiter=100)
        total=coef.sum()
        if not np.isfinite(coef).all() or total<=0: raise ValueError('Invalid NNLS coefficients')
        guesses.append(coef/total); sums.append(total)
    return np.asarray(guesses),np.asarray(sums)

def record_predictions(rows,reference,mixtures,truth,types,info,checks):
    test_p=np.arange(1,len(types)+1,dtype=float); test_p/=test_p.sum()
    recovered,_=fit(reference,(reference@test_p)[None,:])
    recovery=float(np.max(np.abs(recovered[0]-test_p)))
    checks['max_exact_recovery_error']=max(checks['max_exact_recovery_error'],recovery)
    if recovery>1e-6: raise ValueError('Exact-mixture numerical control failed')
    checks['max_signature_condition']=max(checks['max_signature_condition'],float(np.linalg.cond(reference)))
    pred,sums=fit(reference,mixtures)
    error=pred-truth
    for m in range(len(truth)):
        row={**info,'mixture_id':m,'mae':float(np.abs(error[m]).mean()),
             'rmse':float(np.sqrt(np.square(error[m]).mean())),
             'coefficient_sum':float(sums[m])}
        for k,ct in enumerate(types):
            row[f'true_{ct}']=float(truth[m,k]); row[f'pred_{ct}']=float(pred[m,k])
        rows.append(row)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_manifest',type=Path)
    args=parser.parse_args()
    started=time.perf_counter()
    manifest_path=args.input_manifest.resolve()
    manifest=json.loads(manifest_path.read_bytes())
    run=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results/deconvolution_pilot'/run
    out.mkdir(parents=True,exist_ok=False)
    donors=[r['donor'] for r in manifest['records']]
    genes=json.loads((ROOT/manifest['genes_path']).read_bytes())
    assert sha(ROOT/manifest['genes_path'])==manifest['genes_sha256']
    normalized={}; obs={}; pools={}; means={}; detections={}
    types=sorted(set.intersection(*[
        {ct for ct,n in r['cell_type_counts'].items() if n>=50} for r in manifest['records']]))
    if len(types)<3: raise ValueError('Fewer than three eligible types')
    for rec in manifest['records']:
        donor=rec['donor']
        for key,hashkey in [('raw_path','sha256'),('matrix_path','matrix_sha256'),('cells_path','cells_sha256')]:
            if sha(ROOT/rec[key])!=rec[hashkey]: raise ValueError('Input hash mismatch: '+key)
        matrix=sparse.load_npz(ROOT/rec['matrix_path']).astype(np.float64)
        frame=pd.read_csv(ROOT/rec['cells_path'])
        if matrix.shape!=(len(frame),len(genes)): raise ValueError('Matrix annotation mismatch')
        if frame['cell_id'].duplicated().any(): raise ValueError('Duplicate donor cell IDs')
        library=np.asarray(matrix.sum(axis=1)).ravel()
        if not np.allclose(library,frame['umi_total']): raise ValueError('UMI totals do not match')
        matrix=(sparse.diags(10000/library)@matrix).tocsr()
        if not np.allclose(np.asarray(matrix.sum(axis=1)).ravel(),10000): raise ValueError('Normalization failed')
        normalized[donor]=matrix; obs[donor]=frame
        pools[donor]={ct:np.flatnonzero(frame['cell_type'].to_numpy()==ct) for ct in types}
        means[donor]=np.stack([np.asarray(matrix[pools[donor][ct]].mean(axis=0)).ravel() for ct in types])
        detections[donor]=np.stack([np.asarray((matrix[pools[donor][ct]]>0).mean(axis=0)).ravel() for ct in types])
    # Only sample-count metadata determined type eligibility; test expression never selects features.
    rng=np.random.default_rng(SEED)
    composition=np.vstack([np.ones(len(types))/len(types),rng.dirichlet(np.ones(len(types)),N_MIX-1)])
    cell_counts=np.stack([rng.multinomial(N_CELLS,p) for p in composition])
    cell_counts[0]=N_CELLS//len(types)
    cell_counts[0,-1]+=N_CELLS-cell_counts[0].sum()
    truth=cell_counts/N_CELLS
    pd.DataFrame(cell_counts,columns=types).to_csv(out/'mixture_cell_counts.csv',index_label='mixture_id',lineterminator='\n')
    checks={'max_exact_recovery_error':0.,'max_balanced_reference_difference':0.,
            'max_signature_condition':0.,'donor_disjoint':True,'normalization_passed':True}
    rows=[]; feature_records=[]; sampling_records=[]
    for fold,test in enumerate(donors):
        train=[d for d in donors if d!=test]
        assert test not in train and len(train)==3
        train_ids={(d,str(c)) for d in train for c in obs[d]['cell_id']}
        test_ids={(test,str(c)) for c in obs[test]['cell_id']}
        if train_ids & test_ids: raise ValueError('Train/test donor-cell leakage')
        average=np.mean([means[d] for d in train],axis=0)
        detection=np.mean([detections[d] for d in train],axis=0)
        chosen=[]
        for k,ct in enumerate(types):
            other=(average.sum(axis=0)-average[k])/(len(types)-1)
            score=np.log2((average[k]+0.1)/(other+0.1))
            valid=np.flatnonzero((detection[k]>=.1)&(average[k]>=.5)&(score>0))
            ranked=valid[np.lexsort((valid,-score[valid]))][:MARKERS_PER_TYPE]
            if len(ranked)<5: raise ValueError('Insufficient eligible markers')
            chosen.extend(ranked.tolist())
            feature_records.extend([{'held_out':test,'cell_type':ct,'gene':genes[g],
                'specificity_log2':float(score[g]),'rank':i+1} for i,g in enumerate(ranked)])
        features=np.array(sorted(set(chosen)))
        projected={d:normalized[d][:,features].toarray() for d in donors}
        target_rng=np.random.default_rng(SEED+1000+fold)
        mixture=[]
        for m,counts in enumerate(cell_counts):
            indices=np.concatenate([target_rng.choice(pools[test][ct],int(counts[k]),replace=True)
                                    for k,ct in enumerate(types) if counts[k]>0])
            if len(indices)!=N_CELLS: raise ValueError('Wrong mixture cell count')
            mixture.append(projected[test][indices].mean(axis=0))
        mixture=np.asarray(mixture)
        np.savez_compressed(out/f'{test}_targets.npz',expression=mixture,truth=truth,feature_indices=features)
        scenarios=[('balanced','none',1,np.array([20,20,20]))]
        for ratio,major,minor in [(4,40,10),(10,50,5)]:
            for d,dominant in enumerate(train):
                allocations=np.full(3,minor); allocations[d]=major
                scenarios.append((f'ratio{ratio}_{dominant}',dominant,ratio,allocations))
        full_means=np.stack([means[d][:,features] for d in train]) # donor x type x gene
        for arm in ['fixed_means','sampled_budget']:
            repeats=[-1] if arm=='fixed_means' else range(N_REPEATS)
            for repeat in repeats:
                if arm=='sampled_budget':
                    ref_rng=np.random.default_rng(SEED+2000+fold*100+repeat)
                    shuffled={d:{ct:ref_rng.permutation(pools[d][ct]) for ct in types} for d in train}
                for name,dominant,ratio,allocation in scenarios:
                    if arm=='fixed_means': donor_means=full_means
                    else:
                        donor_means=[]
                        for d,n in zip(train,allocation):
                            per_type=[]
                            for ct in types:
                                indices=shuffled[d][ct][:n]
                                if len(indices)!=n or len(np.unique(indices))!=n: raise ValueError('Invalid reference sample')
                                per_type.append(projected[d][indices].mean(axis=0))
                                sampling_records.append({'held_out':test,'reference_repeat':repeat,'scenario':name,
                                    'reference_donor':d,'cell_type':ct,'n_cells':int(n),
                                    'cell_ids':'|'.join(obs[d].iloc[indices]['cell_id'].astype(str))})
                            donor_means.append(per_type)
                        donor_means=np.asarray(donor_means)
                    pooled=np.einsum('d,dkg->kg',allocation/allocation.sum(),donor_means).T
                    equal=donor_means.mean(axis=0).T
                    if ratio==1:
                        diff=float(np.max(np.abs(pooled-equal)))
                        checks['max_balanced_reference_difference']=max(checks['max_balanced_reference_difference'],diff)
                        if diff>1e-10: raise ValueError('Balanced methods differ')
                    for method,reference in [('pooled',pooled),('donor_equal',equal)]:
                        info={'arm':arm,'held_out':test,'reference_repeat':repeat,'scenario':name,
                              'dominant':dominant,'ratio':ratio,'method':method,'n_features':len(features)}
                        record_predictions(rows,reference,mixture,truth,types,info,checks)
        print(f'{test}: {len(features)} markers; completed both arms',flush=True)
    results=pd.DataFrame(rows)
    expected=len(donors)*7*(N_REPEATS+1)*N_MIX*2
    if len(results)!=expected: raise ValueError(f'Incomplete experiment {len(results)} != {expected}')
    keys=['arm','held_out','reference_repeat','scenario','dominant','ratio','mixture_id']
    comparisons=results.pivot(index=keys,columns='method',values='mae').reset_index()
    comparisons.columns.name=None
    comparisons['equal_minus_pooled_pp']=100*(comparisons['donor_equal']-comparisons['pooled'])
    paired_truth=results.groupby(keys)[[f'true_{ct}' for ct in types]].nunique().max().max()
    if paired_truth!=1: raise ValueError('Unpaired target proportions')
    results.to_csv(out/'predictions.csv.gz',index=False,float_format='%.12g',compression={'method':'gzip','mtime':0},lineterminator='\n')
    comparisons.to_csv(out/'paired_comparisons.csv.gz',index=False,float_format='%.12g',compression={'method':'gzip','mtime':0},lineterminator='\n')
    pd.DataFrame(feature_records).to_csv(out/'marker_panels.csv',index=False,lineterminator='\n')
    pd.DataFrame(sampling_records).to_csv(out/'reference_samples.csv.gz',index=False,compression={'method':'gzip','mtime':0},lineterminator='\n')
    summary=results.groupby(['arm','held_out','scenario','dominant','ratio','method'],as_index=False).agg(
        mean_mae=('mae','mean'),mean_rmse=('rmse','mean'),mean_coefficient_sum=('coefficient_sum','mean'))
    summary['mean_mae_pp']=100*summary['mean_mae']
    summary.to_csv(out/'scenario_summary.csv',index=False,lineterminator='\n')
    paired_summary=comparisons.groupby(['arm','held_out','scenario','dominant','ratio'],as_index=False).agg(
        equal_minus_pooled_pp=('equal_minus_pooled_pp','mean'))
    paired_summary.to_csv(out/'paired_summary.csv',index=False,lineterminator='\n')
    by_type=[]
    for ct in types:
        frame=results[['arm','held_out','scenario','ratio','method']].copy()
        frame['cell_type']=ct
        frame['bias_pp']=100*(results[f'pred_{ct}']-results[f'true_{ct}'])
        frame['absolute_error_pp']=frame['bias_pp'].abs()
        by_type.append(frame.groupby(['arm','held_out','scenario','ratio','method','cell_type'],as_index=False).mean())
    pd.concat(by_type).to_csv(out/'cell_type_errors.csv',index=False,lineterminator='\n')
    ratio_summary=summary.groupby(['arm','held_out','ratio','method'],as_index=False)['mean_mae_pp'].mean()
    ratio_summary.to_csv(out/'donor_ratio_summary.csv',index=False,lineterminator='\n')
    aggregate=ratio_summary.groupby(['arm','ratio','method'],as_index=False)['mean_mae_pp'].mean()
    aggregate.to_csv(out/'aggregate_summary.csv',index=False,lineterminator='\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors={'pooled':'#225a92','donor_equal':'#d06b36'}
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    sample=ratio_summary.query("arm == 'sampled_budget'")
    for ax,donor in zip(axes.flat,donors):
        for method in ['pooled','donor_equal']:
            part=sample[(sample.held_out==donor)&(sample.method==method)].sort_values('ratio')
            ax.plot(range(3),part.mean_mae_pp,'o-',color=colors[method],label=method.replace('_',' '),lw=2)
        ax.set_title(f'Held-out {donor}')
        ax.set_xticks(range(3),['20:20:20','40:10:10','50:5:5'])
        ax.set_xlabel('Reference cells per type (dominant donor rotated)')
        ax.set_ylabel('Mean absolute error (percentage points)')
        ax.grid(axis='y',alpha=.18)
    axes[0,0].legend(frameon=False)
    fig.suptitle('Donor imbalance: fixed 60-cell signature budget per type',fontweight='bold',fontsize=14)
    fig.text(.5,.015,'Exploratory normalized-cell synthetic mixtures; 4 source donors; no population confidence intervals.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.95))
    fig.savefig(out/'pilot_donor_errors.png',dpi=170)
    fig.savefig(out/'pilot_donor_errors.pdf')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.7))
    strong=paired_summary.query("arm == 'sampled_budget' and ratio == 10")
    for i,donor in enumerate(donors):
        part=strong[strong.held_out==donor]
        ax.scatter(np.full(len(part),i)+np.linspace(-.12,.12,len(part)),part.equal_minus_pooled_pp,s=55,color='#225a92')
        for x,(_,row) in zip(np.full(len(part),i)+np.linspace(-.12,.12,len(part)),part.iterrows()):
            ax.annotate(row.dominant.replace('human','ref'),(x,row.equal_minus_pooled_pp),xytext=(4,5),textcoords='offset points',fontsize=8)
    ax.axhline(0,color='#777777',lw=1)
    ax.set_xticks(range(4),donors)
    ax.set_xlabel('Held-out donor')
    ax.set_ylabel('Equal-donor MAE minus pooled MAE\n(percentage points; positive = equal averaging worse)')
    ax.set_title('Strong reference imbalance: all dominant-donor choices',fontweight='bold')
    ax.grid(axis='y',alpha=.15)
    fig.tight_layout()
    fig.savefig(out/'pilot_dominant_donor_tradeoffs.png',dpi=170)
    plt.close(fig)
    fixed=results.query("arm == 'fixed_means' and method == 'donor_equal'")
    checks['fixed_equal_invariance_max_mae_range']=float(fixed.groupby(['held_out','mixture_id']).mae.agg(lambda x:x.max()-x.min()).max())
    if checks['fixed_equal_invariance_max_mae_range']>1e-12: raise ValueError('Fixed equal reference was not invariant')
    checks['prediction_rows']=len(results); checks['expected_prediction_rows']=expected
    run_record={'run_id':run,'status':'completed exploratory pilot','topic':'donor reference imbalance',
        'input_manifest':manifest_path.relative_to(ROOT).as_posix(),'input_manifest_sha256':sha(manifest_path),
        'script':'analysis/run_deconvolution_pilot.py','script_sha256':sha(__file__),
        'protocol':'analysis/DECONVOLUTION_PILOT_PROTOCOL.md','protocol_sha256':sha(ROOT/'analysis/DECONVOLUTION_PILOT_PROTOCOL.md'),
        'seed':SEED,'n_mixtures_per_donor':N_MIX,'cells_per_mixture':N_CELLS,'reference_repeats':N_REPEATS,
        'eligible_cell_types':types,'source_donors':donors,'elapsed_seconds':round(time.perf_counter()-started,3),
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
        'checks':checks,'aggregate':aggregate.to_dict(orient='records'),
        'limitations':['Synthetic mixtures of normalized single cells, not real bulk validation.',
            'Only four source donors and five eligible types; conditional proportions only.',
            'Feature selection uses additional training cells beyond signature budget.',
            'No population p-values, independent validation or established-method superiority claim.',
            'Novelty and publication readiness not established.'],
        'artifacts':[{ 'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.iterdir())]}
    save_json(out/'run.json',run_record)
    print(json.dumps({'run':str(out),'checks':checks,'aggregate':aggregate.to_dict(orient='records')},indent=2))

if __name__=='__main__': main()
