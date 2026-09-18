"""Replay five unchanged R executors from frozen processed inputs. Default: technical cases."""
import argparse,os,shutil,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from common import PACKAGE,sha,read,save,verify,verify_manifest

def compare(actual,expected,keys):
 a=pd.read_csv(actual,keep_default_na=False);e=pd.read_csv(expected,keep_default_na=False)
 assert len(a)==len(e) and not a.duplicated(keys).any() and set(a.columns)==set(e.columns)
 a=a.sort_values(keys).reset_index(drop=True);e=e.sort_values(keys).reset_index(drop=True)
 assert a[keys].equals(e[keys]),'Prediction/diagnostic keys differ'
 maximum=0.
 for c in e.columns:
  if pd.api.types.is_numeric_dtype(e[c]) and e[c].dtype!=bool:
   aa=pd.to_numeric(a[c]).to_numpy();ee=e[c].to_numpy();assert np.isfinite(aa).all() and np.isfinite(ee).all()
   d=float(np.max(np.abs(aa-ee)));assert d<=1e-10,(c,d);maximum=max(maximum,d)
  else:assert a[c].equals(e[c]),'Discrete column differs: '+c
 return dict(rows=len(a),maximum_numeric_difference=maximum,all_columns_checked=len(e.columns))

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--inputs',type=Path,required=True,help='Directory with original/ and external/ input subdirectories')
 p.add_argument('--rscript',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True);p.add_argument('--official-source',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['all','endpoints','gradient','thinning','external','alternative'],default='all')
 p.add_argument('--full',action='store_true',help='Explicitly enable the complete selected scientific grids')
 a=p.parse_args();verify_manifest();plan=read(PACKAGE/'STUDY.json');source=read(PACKAGE/'environment/official_source.json')
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 record=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),mode='full' if a.full else 'technical',stages={},script_sha256=sha(__file__),rscript_sha256=sha(a.rscript),package_manifest_sha256=sha(PACKAGE/'MANIFEST.json'),official_commit=source['commit'],workers=1)
 try:
  for f in source['files']:verify(a.official_source/f['path'],f['sha256'])
  labels=list(plan['stages']) if a.stage=='all' else [a.stage]
  for cohort in sorted({plan['stages'][k]['cohort'] for k in labels}):
   for item in plan['cohorts'][cohort]['files']:verify(a.inputs/cohort/item['filename'],item['sha256'])
  work=out/'control_workspace';(work/'environment/music').mkdir(parents=True)
  shutil.copy2(PACKAGE/'environment/music/runtime.R',work/'environment/music/runtime.R')
  shutil.copy2(PACKAGE/'environment/r_packages.tsv',work/'environment/r_packages.tsv')
  env=dict(os.environ,LC_ALL='C',LANG='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MUSIC_LIBRARY=str(a.r_library.resolve()),MUSIC_SOURCE_ROOT=str(a.official_source.resolve()))
  for label in labels:
   spec=plan['stages'][label];folder=out/label;folder.mkdir();data=(a.inputs/spec['cohort']).resolve()
   inp=dict(donors=spec['donors'],cell_types=spec['cell_types'],data_directory=str(data),case_file=str(PACKAGE/spec['case_file']))
   if spec['cohort']=='external':
    inp.update(reference_file=str(PACKAGE/spec['reference_file']),selection_manifest=str(PACKAGE/spec['selection_manifest']),donor_matrices=[dict(donor=x['donor'],path=str(data/x['filename'])) for x in spec['donor_matrices']])
   control=work/(label+'_input.json');save(control,inp)
   verify(PACKAGE/spec['script'],spec['script_sha256'])
   units=spec['units'] if a.full else [spec['technical_unit']];stage=dict(units=[],prediction_rows=0);record['stages'][label]=stage
   for unit in units:
    dest=folder/('held_out_'+unit if spec['cohort']=='original' else unit);dest.mkdir()
    mode='full' if a.full else spec['technical_mode']
    cmd=[str(a.rscript.resolve()),'--vanilla',str(PACKAGE/spec['script']),str(control),str(dest),unit,mode]
    if spec['cohort']=='external':cmd.append('pooled')
    with (dest/'console.log').open('w',encoding='utf-8') as log:
     subprocess.run(cmd,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=7200 if a.full else 900)
    entry=dict(unit=unit,mode=mode)
    if a.full:
     fold=read(dest/'fold.json');assert fold['prediction_rows']==spec['full_rows_per_unit'];entry['prediction_rows']=fold['prediction_rows']
    else:
     expected=PACKAGE/'expected/technical'/label;entry['comparisons']={}
     for name in spec['expected_predictions']:entry['comparisons'][name]=compare(dest/name,expected/name,['case_id','method','mixture_id'])
     entry['prediction_rows']=sum(x['rows'] for x in entry['comparisons'].values());assert entry['prediction_rows']==spec['technical_rows']
     entry['diagnostics']={}
     for f in expected.glob('*diagnostics.csv'):entry['diagnostics'][f.name]=compare(dest/f.name,f,['case_id','mixture_id'])
     if spec['cohort']=='original':assert read(dest/'fold.json')['traced_vs_plain_maximum_difference']==0
     else:assert read(dest/'controls.json')['batch_equivalence_passed']
    stage['units'].append(entry);stage['prediction_rows']+=entry['prediction_rows']
   assert stage['prediction_rows']==spec['full_prediction_rows' if a.full else 'technical_rows']
   print(label,stage['prediction_rows'],'rows passed',flush=True)
  record.update(status='passed',prediction_rows=sum(x['prediction_rows'] for x in record['stages'].values()))
 except Exception as e:
  record.update(status='failed; partial outputs retained',exception_type=type(e).__name__);raise
 finally:
  record.update(elapsed_seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat())
  record['artifacts']=[dict(path=f.relative_to(out).as_posix(),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(out.rglob('*')) if f.is_file() and f.name!='run.json']
  save(out/'run.json',record)
 print(record['status'],record.get('prediction_rows'),flush=True)
if __name__=='__main__':main()
