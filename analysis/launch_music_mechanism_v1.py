"""Run identifier-frozen diagnostics only after technical equivalence gates."""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,os,shutil,subprocess,threading,time,traceback
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
R=Path('C:/Program Files/R/R-4.5.2/bin/Rscript.exe')
ENV=dict(os.environ,LC_ALL='C',LANG='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*2**20),b''):h.update(b)
    return h.hexdigest()
def rel(p):return Path(p).resolve().relative_to(ROOT).as_posix()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')

def prediction_check(new,old):
    key=['reference_id','target_name','method']
    assert not new.duplicated(key).any() and not old.duplicated(key).any()
    cols=['mae','rmse']+[f'{prefix}_{i}' for prefix in ['true','pred'] for i in range(6)]
    join=new[key+cols].merge(old[key+cols],on=key,suffixes=('_new','_old'),validate='one_to_one')
    assert len(join)==len(new)
    a=join[[c+'_new' for c in cols]].to_numpy(float);b=join[[c+'_old' for c in cols]].to_numpy(float)
    diff=float(np.max(np.abs(a-b)));assert np.isfinite(a).all() and diff<=1e-10,diff
    err=new[[f'pred_{i}' for i in range(6)]].to_numpy(float)-new[[f'true_{i}' for i in range(6)]].to_numpy(float)
    em=float(np.max(abs(np.mean(abs(err),axis=1)-new.mae.to_numpy(float))))
    er=float(np.max(abs(np.sqrt(np.mean(err**2,axis=1))-new.rmse.to_numpy(float))))
    assert max(em,er)<=1e-10
    return {'rows':len(new),'maximum_original_difference':diff,'maximum_mae_recalculation_difference':em,'maximum_rmse_recalculation_difference':er}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('selection',type=Path);ap.add_argument('--adapter',type=Path,default=Path('analysis/run_music_mechanism_v1.R'));a=ap.parse_args()
    sp=a.selection.resolve();adapter=a.adapter.resolve();sel=json.loads(sp.read_bytes())
    assert adapter.is_file()
    for z in sel['sources']+sel['artifacts']:assert sha(ROOT/z['path'])==z['sha256'],z['path']
    ip=ROOT/sel['input_manifest'];inp=json.loads(ip.read_bytes())
    assert sha(ip)==sel['input_manifest_sha256']
    for z in inp['artifacts']:assert sha(ROOT/z['path'])==z['sha256'],z['path']
    oldroot=ROOT/sel['original_run'];old=json.loads((oldroot/'run.json').read_bytes())
    assert sha(oldroot/'run.json')==sel['original_run_sha256']
    oldhash={z['path']:z['sha256'] for z in old['artifacts']}
    for t in sel['triple_keys']:
        p=oldroot/t/'predictions.csv';assert sha(p)==oldhash[rel(p)]
    out=ROOT/'results/music_mechanism'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    for p in [Path(__file__),adapter,ROOT/'analysis/MUSIC_MECHANISM_PROTOCOL_V1.md']:shutil.copyfile(p,out/p.name)
    frozen_adapter=out/adapter.name
    sourcepaths=[sp,adapter,Path(__file__),ROOT/'analysis/MUSIC_MECHANISM_PROTOCOL_V1.md',ROOT/'environment/music/runtime.R',R]
    record={'status':'running controls','started_utc':datetime.now(timezone.utc).isoformat(),'selection':rel(sp),'selection_sha256':sha(sp),
        'input_manifest':sel['input_manifest'],'input_manifest_sha256':sha(ip),'original_run':sel['original_run'],
        'sources':[{'path':rel(p) if p.is_relative_to(ROOT) else str(p),'sha256':sha(p)} for p in sourcepaths],
        'planned_weighted_fits':16128,'planned_prediction_rows':32256,'workers':sel['workers'],'fitting_limit_seconds':sel['fitting_limit_seconds']}
    save(out/'started.json',record);print('OUTPUT',rel(out),flush=True)
    start=time.monotonic();active=set();lock=threading.Lock();stop=threading.Event()
    def run(t,mode):
        if stop.is_set():raise RuntimeError('Diagnostic dispatch stopped')
        folder=out/('controls' if mode=='control' else t);folder.mkdir(exist_ok=False)
        cmd=[str(R),'--vanilla',str(frozen_adapter),str(sp),str(folder),t,mode]
        with (folder/'console.log').open('w',encoding='utf-8') as log:
            proc=subprocess.Popen(cmd,cwd=ROOT,env=ENV,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            with lock:active.add(proc)
            try:
                remaining=sel['fitting_limit_seconds']-(time.monotonic()-start)
                if remaining<=0:raise TimeoutError('Fitting limit exceeded')
                code=proc.wait(timeout=remaining)
                if code:raise RuntimeError(f'{t}/{mode}: R exit {code}; see {rel(folder)}/console.log')
            finally:
                if proc.poll() is None:proc.kill();proc.wait()
                with lock:active.discard(proc)
        return folder
    try:
        first=sel['triple_keys'][0];control=run(first,'control')
        ctl=json.loads((control/'controls.json').read_bytes())
        assert ctl['batch_equivalence_passed'] and ctl['basis_checks_passed']
        assert ctl['trace_equivalence_passed'] and ctl['maxiter_count']==0
        assert all(x['equivalence_passed'] for x in ctl['checks'])
        nc=pd.read_csv(control/'control_predictions.csv',keep_default_na=False)
        oldc=pd.read_csv(oldroot/first/'predictions.csv',keep_default_na=False)
        assert len(nc)==768 and set(nc.batch_mode)=={'pooled','separate'}
        control_checks={k:prediction_check(g,oldc) for k,g in nc.groupby('batch_mode',sort=True)}
        save(out/'control_gate.json',{'status':'passed R subset/batching/basis controls and original-prediction comparison','r_controls':ctl,'original_comparisons':control_checks})
        print('CONTROL GATE PASSED',flush=True)
        record['status']='running full bounded diagnostics';save(out/'progress.json',record)
        completed=[]
        with ThreadPoolExecutor(max_workers=sel['workers']) as pool:
            futures={pool.submit(run,t,'full'):t for t in sel['triple_keys']}
            try:
                for fut in as_completed(futures):
                    folder=fut.result();fold=json.loads((folder/'fold.json').read_bytes());completed.append(fold)
                    assert fold['status']=='completed' and fold['basis_checks_passed'] and fold['maxiter_count']==0
                    save(out/'progress.json',{'status':'running','completed_triples':len(completed),'total_triples':14,'elapsed_seconds':time.monotonic()-start,'last_triple':folder.name})
                    print('COMPLETED',len(completed),'/14',folder.name,round(time.monotonic()-start,1),'s',flush=True)
                    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<=sel['output_limit_bytes']
            except BaseException:
                stop.set()
                with lock:
                    for proc in list(active):
                        if proc.poll() is None:proc.kill()
                for fut in futures:fut.cancel()
                raise
        checks=[]
        for t in sel['triple_keys']:
            new=pd.read_csv(out/t/'predictions.csv',keep_default_na=False);expected=pd.read_csv(oldroot/t/'predictions.csv',keep_default_na=False)
            assert len(new)==2304
            checks.append({'triple_key':t,**prediction_check(new,expected)})
        assert sum(x['rows'] for x in checks)==32256
        save(out/'original_prediction_checks.json',{'status':'passed all diagnostic predictions against original outputs','checks':checks})
        for s in record['sources']:assert sha(ROOT/s['path'])==s['sha256'] if not Path(s['path']).is_absolute() else sha(s['path'])==s['sha256']
        record.update(status='completed bounded diagnostics; descriptive review pending',finished_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,folds=completed,prediction_checks=checks)
    except BaseException as e:
        record.update(status='failed; outputs retained',error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start)
        save(out/'failed.json',record);print(record['error'],flush=True);raise
    finally:
        record['artifacts']=[{'path':rel(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['run.json','failed.json']]
        save(out/'run.json',record)
    print('COMPLETED ALL DIAGNOSTICS',round(record['elapsed_seconds'],1),'s',flush=True)

if __name__=='__main__':main()
