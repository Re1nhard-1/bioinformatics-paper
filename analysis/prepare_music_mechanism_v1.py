"""Freeze identifier-selected post-hoc diagnostics; no fits or outcome selection."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,shutil
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'results/music_external_input/20260917T122709190802Z/input_manifest.json'
OLD=ROOT/'results/music_external/20260917T125429748635Z'
PROTOCOL=ROOT/'analysis/MUSIC_MECHANISM_PROTOCOL_V1.md'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*2**20),b''):h.update(b)
    return h.hexdigest()
def rel(p):return Path(p).resolve().relative_to(ROOT).as_posix()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')

def main():
    m=json.loads(INPUT.read_bytes());old=json.loads((OLD/'run.json').read_bytes())
    assert sha(INPUT)==old['input_manifest_sha256']
    for item in m['artifacts']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
    refs=pd.read_csv(ROOT/m['reference_file'],sep='\t',keep_default_na=False)
    refs=refs[refs.level.isin(['balanced','ratio10'])].copy()
    refs=refs.sort_values(['triple_key','block','budget','level','dominant_donor'])
    cases=pd.read_csv(ROOT/m['case_file'],sep='\t',keep_default_na=False)
    cases=cases[cases.reference_id.isin(refs.reference_id)].copy()
    mids=sorted([0]+sorted(range(1,60),key=lambda i:hashlib.sha256(f'music-mechanism-v1|mixture|{i}'.encode()).hexdigest())[:7])
    targets=pd.read_csv(ROOT/m['data_directory']/'targets.tsv',sep='\t',keep_default_na=False)
    targets=targets[targets.mixture_id.isin(mids)].copy()
    assert len(refs)==504 and len(cases)==2016 and len(targets)==112
    assert refs.groupby('triple_key').size().eq(36).all() and cases.reference_id.value_counts().eq(4).all()
    assert targets.groupby('held_out').size().eq(8).all()
    inventory=ROOT/'results/music_external_selection/20260917T122427453355Z/cohort_inventory.csv'
    inv=pd.read_csv(inventory,keep_default_na=False);cells=pd.read_csv(ROOT/m['data_directory']/'cells.tsv',sep='\t',usecols=['donor','author_cell_type'])
    names={'B':'B cells','cM':'CD14+ Monocytes','T4':'CD4 T cells','T8':'CD8 T cells','ncM':'FCGR3A+ Monocytes','NK':'NK cells'}
    cover=inv.melt(id_vars=['donor_id'],value_vars=list(names),var_name='author_type',value_name='source_inventory_N')
    exported=cells.groupby(['donor','author_cell_type']).size().rename('exported_union_N').reset_index()
    cover=cover.merge(exported,left_on=['donor_id','author_type'],right_on=['donor','author_cell_type'],validate='one_to_one')
    cover['cell_type']=cover.author_type.map(names);cover['exported_fraction']=cover.exported_union_N/cover.source_inventory_N
    assert len(cover)==84 and cover.exported_union_N.sum()==93892 and cover.source_inventory_N.sum()==114524
    out=ROOT/'results/music_mechanism_selection'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    for name,df in [('references',refs),('cases',cases),('targets',targets)]:df.to_csv(out/f'{name}.tsv',sep='\t',index=False,lineterminator='\n')
    cover.to_csv(out/'inventory_coverage.csv',index=False,lineterminator='\n')
    shutil.copyfile(PROTOCOL,out/PROTOCOL.name);shutil.copyfile(__file__,out/Path(__file__).name)
    record={'status':'frozen post-hoc diagnostic selection; no new fits','created_utc':datetime.now(timezone.utc).isoformat(),
      'input_manifest':rel(INPUT),'input_manifest_sha256':sha(INPUT),'original_run':rel(OLD),'original_run_sha256':sha(OLD/'run.json'),
      'references_file':rel(out/'references.tsv'),'cases_file':rel(out/'cases.tsv'),'targets_file':rel(out/'targets.tsv'),
      'inventory_file':rel(inventory),'inventory_sha256':sha(inventory),'donors':m['donors'],'cell_types':m['cell_types'],
      'mixture_ids':mids,'unique_references':504,'logical_cases':2016,'targets':112,'weighted_fits':16128,'prediction_rows':32256,
      'triple_keys':sorted(refs.triple_key.unique()),'workers':3,'fitting_limit_seconds':1800,'output_limit_bytes':2**30,
      'sources':[{'path':rel(p),'sha256':sha(p)} for p in [INPUT,PROTOCOL,Path(__file__),inventory,OLD/'run.json']],
      'artifacts':[{'path':rel(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(out.iterdir())]}
    save(out/'selection.json',record)
    print('OUTPUT',rel(out/'selection.json'),flush=True)
    print(json.dumps({k:record[k] for k in ['mixture_ids','unique_references','targets','weighted_fits']}))

if __name__=='__main__':main()
