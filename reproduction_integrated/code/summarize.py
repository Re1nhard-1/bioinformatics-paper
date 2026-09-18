"""Recalculate primary MAE/RMSE results from five complete saved prediction stages."""
import argparse,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from common import PACKAGE,read,save,verify,sha,verify_manifest

METRICS=['mae_pp','rmse_pp']
def compare_table(actual,expected,keys,values):
 a=actual.set_index(keys).sort_index();e=expected.set_index(keys).sort_index()
 assert a.index.equals(e.index) and not a.index.duplicated().any(),keys
 delta=float(np.max(np.abs(a[values].to_numpy()-e[values].to_numpy())))
 assert delta<1e-10,delta
 return dict(rows=len(a),values=len(a)*len(values),maximum_difference_pp=delta)
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for stage in ['endpoints','gradient','thinning','external','alternative']:p.add_argument('--'+stage,type=Path,required=True)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--recomputed',action='store_true',help='For new full-route output: skip historical byte equality, retain file/key/error/result checks')
 a=p.parse_args();verify_manifest();plan=read(PACKAGE/'STUDY.json');out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 start=time.monotonic();record=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(__file__),historical_prediction_bytes_checked=not a.recomputed,package_manifest_sha256=sha(PACKAGE/'MANIFEST.json'),prediction_rows={},stage_error_checks={},tables={})
 all_data={};input_hashes=[]
 def export(name,df):df.to_csv(out/(name+'.csv'),index=False,lineterminator='\n');return df
 def expected(name,file):return pd.read_csv(PACKAGE/'expected/summaries'/name/(file+'.csv'),dtype={'held_out':str},keep_default_na=False)
 def compare(name,actual,exp,keys,vals=METRICS):record['tables'][name]=compare_table(actual,exp,keys,vals)
 try:
  for label,spec in plan['stages'].items():
   folder=getattr(a,label);frames=[];maxmae=0.;maxrmse=0.
   for item in plan['prediction_files'][label]:
    path=folder/item['path'];h=sha(path)
    if not a.recomputed:assert h==item['sha256'],item['path']
    input_hashes.append(dict(stage=label,path=item['path'],sha256=h))
    f=pd.read_csv(path,dtype={'held_out':str,'dominant_donor':str},keep_default_na=False)
    truth=f[[f'true_{i}' for i in range(6)]].to_numpy();pred=f[[f'pred_{i}' for i in range(6)]].to_numpy()
    assert np.isfinite(truth).all() and np.isfinite(pred).all() and pred.min()>=-1e-10 and truth.min()>=0
    assert np.max(abs(truth.sum(axis=1)-1))<1e-10 and np.max(abs(pred.sum(axis=1)-1))<1e-10
    err=pred-truth;mae=abs(err).mean(axis=1);rmse=np.sqrt((err*err).mean(axis=1))
    maxmae=max(maxmae,float(np.max(abs(mae-f.mae.to_numpy()))));maxrmse=max(maxrmse,float(np.max(abs(rmse-f.rmse.to_numpy()))))
    assert maxmae<1e-10 and maxrmse<1e-10
    f['mae_pp']=mae*100;f['rmse_pp']=rmse*100
    if label=='endpoints':f['level']=np.where(f.ratio==1,20,50)
    if label=='gradient':f['level']=f.major_count.astype(int)
    frames.append(f)
   df=pd.concat(frames,ignore_index=True);assert len(df)==spec['full_prediction_rows']
   assert not df.duplicated(['case_id','method','mixture_id']).any()
   cases=pd.read_csv(PACKAGE/spec['case_file'],sep='\t',dtype={'held_out':str},keep_default_na=False)
   if label in ['gradient','thinning']:cases=cases[cases.policy==('gradient' if label=='gradient' else 'thin5')]
   assert set(df.case_id)==set(cases.case_id)
   assert set(df.method)=={'music_weighted','music_nnls'} and df.groupby(['case_id','method']).size().eq(60).all()
   assert df.groupby(['case_id','method']).mixture_id.agg(lambda x:set(x)==set(range(60))).all()
   index=cases.set_index('case_id');observed=df[['case_id','held_out','triple_id','block']].drop_duplicates().set_index('case_id').sort_index()
   assert len(observed)==len(cases) and observed.equals(index.loc[observed.index,['held_out','triple_id','block']])
   record['prediction_rows'][label]=len(df);record['stage_error_checks'][label]=dict(mae_maximum_difference=maxmae,rmse_maximum_difference=maxrmse)
   all_data[label]=df;print('Read and checked',label,len(df),'rows',flush=True)
  original=pd.concat([all_data['endpoints'],all_data['gradient']],ignore_index=True)
  donor=export('original_donor_metrics',original.groupby(['held_out','method','level'])[METRICS].mean().reset_index())
  compare('original_donor_metrics',donor,expected('original','donor_metrics'),['held_out','method','level'])
  agg=export('original_aggregate_metrics',donor.groupby(['method','level'])[METRICS].mean().reset_index())
  compare('original_aggregate_metrics',agg,expected('original','aggregate_metrics'),['method','level'])
  # The thinning baseline reuses the three strong-allocation parents; they are not new samples.
  keep=all_data['endpoints'].copy();keep['policy']=np.where(keep.ratio==10,'keep60','expanded60')
  thin=all_data['thinning'].copy();thin['policy']='thin15'
  td=export('thinning_donor_metrics',pd.concat([keep,thin],ignore_index=True).groupby(['held_out','method','policy'])[METRICS].mean().reset_index())
  compare('thinning_donor_metrics',td,expected('thinning','donor_metrics'),['held_out','method','policy'])
  ta=export('thinning_aggregate_metrics',td.groupby(['method','policy'])[METRICS].mean().reset_index())
  compare('thinning_aggregate_metrics',ta,expected('thinning','aggregate_metrics'),['method','policy'])
  tw=td.pivot(index=['held_out','method'],columns='policy',values=METRICS)
  te=pd.DataFrame({m.replace('_pp','_thin_minus_keep_pp'):tw[(m,'thin15')]-tw[(m,'keep60')] for m in METRICS}).reset_index()
  export('thinning_donor_policy_effects',te);compare('thinning_donor_policy_effects',te,expected('thinning','donor_policy_effects'),['held_out','method'],['mae_thin_minus_keep_pp','rmse_thin_minus_keep_pp'])
  primaries={}
  for label in ['external','alternative']:
   df=all_data[label];keys=['held_out','method','budget','level']
   dm=export(label+'_donor_metrics',df.groupby(keys)[METRICS].mean().reset_index())
   compare(label+'_donor_metrics',dm,expected(label,'donor_metrics'),keys)
   am=export(label+'_aggregate_metrics',dm.groupby(keys[1:])[METRICS].mean().reset_index())
   compare(label+'_aggregate_metrics',am,expected(label,'aggregate_metrics'),keys[1:])
   # Independently join paired target/reference conditions before averaging the contrast.
   pk=['held_out','triple_id','block','budget','method','mixture_id']
   base=df[df.level=='balanced'][pk+METRICS]
   paired=df[df.level=='ratio10'].merge(base,on=pk,validate='many_to_one',suffixes=('','_baseline'))
   assert len(paired)==len(df[df.level=='ratio10'])
   for m in METRICS:paired[m.replace('_pp','_penalty_pp')]=paired[m]-paired[m+'_baseline']
   pc=['mae_penalty_pp','rmse_penalty_pp']
   dp=export(label+'_strongest_donor_penalties',paired.groupby(['held_out','method','budget'])[pc].mean().reset_index())
   ex=expected(label,'donor_penalties');ex=ex[ex.level=='ratio10'];compare(label+'_strongest_donor_penalties',dp,ex,['held_out','method','budget'],pc)
   w=dp[dp.method=='music_weighted'].pivot(index='held_out',columns='budget',values=pc)
   pr=pd.DataFrame({m.replace('_penalty_pp','_penalty_change_pp'):w[(m,300)]-w[(m,60)] for m in pc}).reset_index()
   export(label+'_primary_endpoint',pr);compare(label+'_primary_endpoint',pr,expected(label,'primary_endpoint_all_donors'),['held_out'],['mae_penalty_change_pp','rmse_penalty_change_pp'])
   primaries[label]=pr
   record[label+'_primary_summary']=dict(mean_pp=float(pr.mae_penalty_change_pp.mean()),negative=int((pr.mae_penalty_change_pp < -1e-10).sum()),positive=int((pr.mae_penalty_change_pp > 1e-10).sum()),n=len(pr))
  paired=primaries['external'].merge(primaries['alternative'],on='held_out',suffixes=('_original','_alternative'),validate='one_to_one')
  paired['alternative_minus_original_pp']=paired.mae_penalty_change_pp_alternative-paired.mae_penalty_change_pp_original
  export('grouping_paired_primary',paired)
  # Saved two-design table contains both methods; primary check here is weighted only.
  exp=expected('alternative','paired_design_comparison');exp=exp[exp.method=='music_weighted']
  compare('grouping_paired_primary',paired,exp.rename(columns={'followup_minus_original_pp':'alternative_minus_original_pp'}),['held_out'],['alternative_minus_original_pp'])
  record.update(status='passed',total_prediction_rows=sum(record['prediction_rows'].values()),grouping_mean_change_pp=float(paired.alternative_minus_original_pp.mean()),checked_numeric_values=sum(x['values'] for x in record['tables'].values()),maximum_summary_difference_pp=max(x['maximum_difference_pp'] for x in record['tables'].values()))
 except Exception as e:record.update(status='failed; partial outputs retained',exception_type=type(e).__name__);raise
 finally:
  record.update(elapsed_seconds=time.monotonic()-start,scope='Saved-prediction error and selected primary-summary replay; no new fitting, population inference or complete historical-analysis regeneration.')
  save(out/'input_hashes.json',input_hashes);save(out/'verification.json',record)
 print(record['status'],record.get('total_prediction_rows'),record.get('maximum_summary_difference_pp'),flush=True)
if __name__=='__main__':main()
