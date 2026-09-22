"""Run the two frozen intermediate allocations after mandatory replay/trace controls."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import json
import os
import subprocess
import time
import numpy as np
import pandas as pd
from run_deconvolution_pilot import sha, save_json

ROOT=Path(__file__).resolve().parents[1]
RSCRIPT=Path('C:/Program Files/R/R-4.5.2/bin/Rscript.exe')
ORIGINAL=ROOT/'results/music_comparison/20260917T090859412332Z'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest',type=Path)
    args=parser.parse_args()
    manifest=args.manifest.resolve()
    info=json.loads(manifest.read_bytes())
    for item in info['sources']+info['unchanged_inputs']+info['artifacts']:
        assert sha(ROOT/item['path'])==item['sha256'],item['path']
    original=json.loads((ORIGINAL/'run.json').read_bytes())
    assert original['status']=='completed' and original['prediction_rows']==86400
    for item in original['source_files']:
        assert sha(ROOT/item['path'])==item['sha256'],item['path']
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results/music_gradient'/run_id; out.mkdir(parents=True)
    source_paths=[Path(__file__),ROOT/'analysis/run_music_gradient.R',ROOT/'analysis/MUSIC_GRADIENT_PROTOCOL.md']
    record={'status':'running','run_id':run_id,'started_utc':datetime.now(timezone.utc).isoformat(),
            'input_manifest':manifest.relative_to(ROOT).as_posix(),'input_manifest_sha256':sha(manifest),
            'source_files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in source_paths],
            'original_run':(ORIGINAL/'run.json').relative_to(ROOT).as_posix(),'original_run_sha256':sha(ORIGINAL/'run.json'),
            'official_commit':original['official_commit'],'rscript_sha256':sha(RSCRIPT),
            'total_budget_seconds':1800,'parallel_workers':3,'controls':[],'folds':[]}
    save_json(out/'started.json',record)
    started=time.perf_counter()
    env=dict(os.environ,LC_ALL='C',LANG='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    def call(donor,mode,folder):
        folder.mkdir(parents=True)
        command=[str(RSCRIPT),'--vanilla','analysis/run_music_gradient.R',str(manifest),str(folder),donor,mode]
        with (folder/'console.log').open('w',encoding='utf-8') as handle:
            subprocess.run(command,cwd=ROOT,env=env,stdout=handle,stderr=subprocess.STDOUT,check=True,
                           timeout=max(1,1800-(time.perf_counter()-started)))
        fold=json.loads((folder/'fold.json').read_bytes())
        return fold
    name='run.json'
    try:
        replay=out/'controls/replay'
        record['controls'].append(call('101','replay',replay))
        oldpath=ORIGINAL/'held_out_101/101_t00_b0_balanced.csv'
        oldhash=next(x['sha256'] for x in original['artifacts'] if x['path']==oldpath.relative_to(ROOT).as_posix())
        assert sha(oldpath)==oldhash
        key=['case_id','method','mixture_id']
        actual=pd.read_csv(replay/'101_t00_b0_balanced.csv').sort_values(key)
        expected=pd.read_csv(oldpath).sort_values(key)
        assert actual[key].reset_index(drop=True).equals(expected[key].reset_index(drop=True))
        numeric=['mae','rmse']+[f'{prefix}_{k}' for prefix in ['true','pred'] for k in range(6)]
        delta=float(np.max(np.abs(actual[numeric].to_numpy()-expected[numeric].to_numpy())))
        assert delta<1e-10
        record['original_60_target_replay_max_abs_difference']=delta
        tech=call('101','technical',out/'controls/gradient_technical')
        assert tech['traced_vs_plain_maximum_difference']==0 and tech['prediction_rows']==6
        record['controls'].append(tech)
        print('Replay and gradient trace controls passed; now executing all 1080 locked cases.',flush=True)
        save_json(out/'controls_passed.json',record)
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures={pool.submit(call,d,'full',out/f'held_out_{d}'):d for d in info['donors']}
            for future in as_completed(futures):
                fold=future.result()
                assert fold['cases']==180 and fold['prediction_rows']==21600
                record['folds'].append(fold)
                print('Completed held-out donor',futures[future],flush=True)
        assert sum(f['prediction_rows'] for f in record['folds'])==129600
        record.update(status='completed',prediction_rows=129600,case_count=1080,controls_prediction_rows=126)
    except Exception as error:
        name='failed.json'
        record.update(status='failed; partial outputs retained',error=repr(error))
    record['elapsed_seconds']=time.perf_counter()-started
    record['artifacts']=[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size}
                         for p in sorted(out.rglob('*')) if p.is_file()]
    save_json(out/name,record)
    print(json.dumps({'output':(out/name).relative_to(ROOT).as_posix(),'status':record['status'],'elapsed_seconds':record['elapsed_seconds']},indent=2))
    if name=='failed.json':raise SystemExit(1)

if __name__=='__main__':main()
