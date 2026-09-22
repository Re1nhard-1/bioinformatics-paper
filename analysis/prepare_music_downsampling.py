"""Create one outcome-blind, paired donor/type thinning design from saved draws."""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import pandas as pd
from run_deconvolution_pilot import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT/'results/music_input/20260917T090438982999Z/input_manifest.json'
HISTORY = ROOT/'results/deconvolution_v3_1/20260916T221536546175Z'

def main():
    old = json.loads(ORIGINAL.read_bytes())
    for item in old['artifacts']:
        assert sha(ROOT/item['path']) == item['sha256']
    data = ROOT/old['data_directory']
    cells = pd.read_csv(data/'cells.tsv',sep='\t',dtype={'donor':str}).set_index('matrix_column_R',drop=False)
    cases = pd.read_csv(data/'cases.tsv',sep='\t',dtype={'held_out':str})
    new, parents, source_files = [], [], []
    for donor in old['donors']:
        path=HISTORY/f'held_out_{donor}/reference_samples.csv.gz'
        historic=pd.read_csv(path,dtype={'held_out':str,'reference_donor':str})
        expected=next(x for x in json.loads((HISTORY/'run.json').read_bytes())['artifacts'] if x['path']==path.relative_to(ROOT).as_posix())
        assert sha(path)==expected['sha256']
        source_files.append({'path':path.relative_to(ROOT).as_posix(),'sha256':sha(path)})
        for (triple,block),group in cases[cases.held_out.eq(donor)].groupby(['triple_id','block']):
            chosen=[]
            for row in group.itertuples():
                ids=[int(x) for x in row.reference_columns_R.split('|')]
                obs=cells.loc[ids]
                by=obs.groupby(['donor','cell_type'],sort=True)
                subset=[int(x) for _,part in by for x in part.matrix_column_R.iloc[:5]]
                assert len(subset)==len(set(subset))==90 and set(subset)<=set(ids)
                assert donor not in set(cells.loc[subset].donor)
                saved=historic[historic.triple_id.eq(triple)&historic.block.eq(block)&historic['repeat'].eq(0)&historic.scenario.eq(row.scenario)]
                saved=saved.sort_values(['reference_donor','cell_type'])
                id_strings=[c for s in saved.cell_ids for c in s.split('|')[:5]]
                assert list(cells.loc[subset].cell_id)==id_strings
                chosen.append(tuple(subset))
                if row.scenario != 'balanced':
                    parents.append({'thin_case':f'{donor}_t{triple:02d}_b{block}_thin5','parent_case':row.case_id,
                                    'held_out':donor,'triple_id':int(triple),'block':int(block),'parent_scenario':row.scenario})
            assert len(chosen)==4 and len(set(chosen))==1
            new.append({'case_id':f'{donor}_t{triple:02d}_b{block}_thin5','held_out':donor,
                        'triple_id':int(triple),'block':int(block),'scenario':'thin5','policy':'thin5',
                        'reference_columns_R':'|'.join(map(str,chosen[0])), 'cells_per_donor_type':5,
                        'reference_cells':90,'cells_per_type':15})
    control=cases.iloc[0].to_dict()
    assert control['case_id']=='101_t00_b0_balanced'
    control.update(policy='replay',cells_per_donor_type=20,reference_cells=360,cells_per_type=60)
    new.append(control)
    assert len(new)==181 and len(parents)==540
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'results/music_downsampling_input'/run_id; out.mkdir(parents=True)
    pd.DataFrame(new).to_csv(out/'cases.tsv',sep='\t',index=False,lineterminator='\n')
    pd.DataFrame(parents).to_csv(out/'parent_pairs.tsv',sep='\t',index=False,lineterminator='\n')
    paths=[Path(__file__),ROOT/'analysis/MUSIC_DOWNSAMPLING_PROTOCOL.md']
    record={'status':'prepared; no new predictions','data_directory':old['data_directory'],
            'case_file':(out/'cases.tsv').relative_to(ROOT).as_posix(),'donors':old['donors'],'cell_types':old['cell_types'],
            'case_count':180,'prediction_rows':21600,'reference_cells':90,'discarded_fraction':0.75,
            'parent_pairs':540,'all_four_original_allocations_share_same_first_five_subset':True,
            'source_manifest':ORIGINAL.relative_to(ROOT).as_posix(),'source_manifest_sha256':sha(ORIGINAL),
            'sources':source_files+[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in paths],
            'unchanged_inputs':[x for x in old['artifacts'] if Path(x['path']).name in ['counts.mtx.gz','targets.mtx.gz','cells.tsv','genes.tsv','targets.tsv']],
            'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in out.iterdir()]}
    save_json(out/'input_manifest.json',record)
    print((out/'input_manifest.json').relative_to(ROOT).as_posix())

if __name__=='__main__':main()
