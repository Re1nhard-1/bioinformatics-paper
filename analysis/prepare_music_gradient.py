"""Freeze paired intermediate allocations using original donor/type permutations."""
from pathlib import Path
from datetime import datetime, timezone
import json
import pandas as pd
from run_deconvolution_pilot import sha, save_json
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/'results/music_input/20260917T090438982999Z/input_manifest.json'
HISTORY=ROOT/'results/deconvolution_v3_1/20260916T221536546175Z'

def main():
    old=json.loads(ORIGINAL.read_bytes())
    for x in old['artifacts']: assert sha(ROOT/x['path'])==x['sha256']
    data=ROOT/old['data_directory']
    cells=pd.read_csv(data/'cells.tsv',sep='\t',dtype={'donor':str}).set_index('matrix_column_R',drop=False)
    original_cases=pd.read_csv(data/'cases.tsv',sep='\t',dtype={'held_out':str})
    history_hash={x['path']:x['sha256'] for x in json.loads((HISTORY/'run.json').read_bytes())['artifacts']}
    new=[]; sources=[]; old_prefix_groups=0
    for donor in old['donors']:
        p=HISTORY/f'held_out_{donor}/reference_samples.csv.gz'
        assert sha(p)==history_hash[p.relative_to(ROOT).as_posix()]
        sources.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)})
        hist=pd.read_csv(p,dtype={'held_out':str,'reference_donor':str})
        hist=hist[hist['repeat'].eq(0)]
        for (triple,block),group in original_cases[original_cases.held_out.eq(donor)].groupby(['triple_id','block']):
            permutations={}; originals={}
            for row in group.itertuples():
                ids=list(map(int,row.reference_columns_R.split('|')))
                obs=cells.loc[ids]
                saved=hist[hist.triple_id.eq(triple)&hist.block.eq(block)&hist.scenario.eq(row.scenario)].set_index(['reference_donor','cell_type'])
                for (d,t),part in obs.groupby(['donor','cell_type'],sort=True):
                    sequence=part.matrix_column_R.astype(int).tolist()
                    assert part.cell_id.tolist()==saved.loc[(d,t),'cell_ids'].split('|')
                    if len(sequence)==50:
                        assert (d,t) not in permutations
                        permutations[d,t]=sequence
                originals[row.scenario]=row
            assert len(permutations)==18
            for row in group.itertuples():
                ids=list(map(int,row.reference_columns_R.split('|')))
                for key,part in cells.loc[ids].groupby(['donor','cell_type'],sort=True):
                    assert part.matrix_column_R.astype(int).tolist()==permutations[key][:len(part)]
                    old_prefix_groups+=1
            refs=sorted({d for d,t in permutations})
            assert len(refs)==3 and donor not in refs
            for major,minor in [(30,15),(40,10)]:
                for dominant in refs:
                    chosen=[v for (d,t),permutation in sorted(permutations.items()) for v in permutation[:major if d==dominant else minor]]
                    obs=cells.loc[chosen]
                    assert len(chosen)==len(set(chosen))==360
                    assert obs.groupby('cell_type').size().eq(60).all()
                    for (d,t),part in obs.groupby(['donor','cell_type']):assert len(part)==(major if d==dominant else minor)
                    new.append({'case_id':f'{donor}_t{triple:02d}_b{block}_g{major}_d{dominant}',
                        'held_out':donor,'triple_id':int(triple),'block':int(block),
                        'scenario':f'gradient{major}_{dominant}','policy':'gradient',
                        'allocation':f'{major}:{minor}:{minor}','major_count':major,'minor_count':minor,
                        'dominant_donor':dominant,'reference_cells':360,'cells_per_type':60,
                        'reference_columns_R':'|'.join(map(str,chosen)),
                        'paired_balanced_case':originals['balanced'].case_id,
                        'paired_strong_case':originals[f'ratio10_{dominant}'].case_id})
    assert len(new)==1080
    control=original_cases.iloc[0].to_dict()
    assert control['case_id']=='101_t00_b0_balanced'
    control.update(policy='replay',allocation='20:20:20',major_count=20,minor_count=20,dominant_donor='',reference_cells=360,cells_per_type=60)
    new.append(control)
    out=ROOT/'results/music_gradient_input'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    frame=pd.DataFrame(new)
    assert not frame.case_id.duplicated().any()
    frame.to_csv(out/'cases.tsv',sep='\t',index=False,lineterminator='\n')
    for p in [Path(__file__),ROOT/'analysis/MUSIC_GRADIENT_PROTOCOL.md']:
        sources.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)})
    record={'status':'prepared; no new predictions','data_directory':old['data_directory'],'case_file':(out/'cases.tsv').relative_to(ROOT).as_posix(),
        'donors':old['donors'],'cell_types':old['cell_types'],'case_count':1080,'prediction_rows':129600,
        'source_manifest':ORIGINAL.relative_to(ROOT).as_posix(),'source_manifest_sha256':sha(ORIGINAL),
        'old_donor_type_prefix_groups_checked':old_prefix_groups,'sources':sources,
        'unchanged_inputs':[x for x in old['artifacts'] if Path(x['path']).name in ['counts.mtx.gz','targets.mtx.gz','cells.tsv','genes.tsv','targets.tsv']],
        'artifacts':[{'path':(out/'cases.tsv').relative_to(ROOT).as_posix(),'sha256':sha(out/'cases.tsv')}]}
    save_json(out/'input_manifest.json',record)
    print((out/'input_manifest.json').relative_to(ROOT).as_posix())
if __name__=='__main__':main()
