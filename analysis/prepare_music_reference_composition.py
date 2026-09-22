"""Freeze a single new reference ring while retaining all original cells/targets."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OLD_INPUT = ROOT / 'results/music_external_input/20260917T122709190802Z/input_manifest.json'
PROTOCOL = ROOT / 'analysis/MUSIC_REFERENCE_COMPOSITION_PROTOCOL.md'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8*2**20), b''): h.update(b)
    return h.hexdigest()

def rel(p): return Path(p).relative_to(ROOT).as_posix()
def save(p, x): p.write_text(json.dumps(x, ensure_ascii=False, indent=2)+'\n', encoding='utf-8', newline='\n')

def main():
    inp = json.loads(OLD_INPUT.read_bytes())
    sp = ROOT / inp['selection_manifest']; sel = json.loads(sp.read_bytes())
    assert sha(sp) == inp['selection_manifest_sha256']
    for item in inp['artifacts'] + sel['artifacts']:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    original = inp['donors']; assert len(original) == len(set(original)) == 14
    ring = sorted(original, key=lambda d: hashlib.sha256(f'music-external-reference-composition-v1|ring|{d}'.encode()).hexdigest())
    assert set(ring) == set(original)
    prefix_file = ROOT / sel['data_directory'] / 'reference_prefix_rows.npy'
    prefixes = np.load(prefix_file, allow_pickle=False); assert prefixes.shape == (3,14,6,250)
    cells = pd.read_csv(ROOT/inp['data_directory']/'cells.tsv', sep='\t', keep_default_na=False)
    row_to_col = dict(zip(cells.obs_row_0.astype(int), range(1, len(cells)+1)))
    refs = {}; cases = []; incidence = []; omitted = []
    for i, held in enumerate(ring):
        others = [ring[(i+j)%14] for j in range(1,14)]
        omitted.append({'held_out':held,'omitted_reference_donor':others[-1]})
        for ti in range(4):
            donors = sorted(others[3*ti:3*ti+3])
            triple = '_'.join(f'd{original.index(d):02d}' for d in donors)
            incidence.extend({'held_out':held,'triple_id':ti,'triple_key':triple,'reference_donor':d} for d in donors)
            for block in range(3):
                for budget in [60,300]:
                    allocations = [('balanced','',budget//3,budget//3)] + [('ratio10',d,5*budget//6,budget//12) for d in donors]
                    for level, dominant, major, minor in allocations:
                        rid = f'rc1_{triple}_b{block}_n{budget}_{level}' + (f'_d{original.index(dominant):02d}' if dominant else '')
                        if rid not in refs:
                            rows = np.concatenate([prefixes[block,original.index(d),k,:major if d==dominant else minor] for d in donors for k in range(6)])
                            cols = [row_to_col[int(r)] for r in rows]
                            refs[rid] = dict(reference_id=rid,triple_key=triple,reference_donors='|'.join(donors),block=block,budget=budget,level=level,dominant_donor=dominant,major_count=major,minor_count=minor,reference_cells=6*budget,reference_columns_R='|'.join(map(str,cols)))
                        cases.append(dict(case_id=f'rc1_d{original.index(held):02d}_t{ti}_b{block}_n{budget}_{level}'+(f'_d{original.index(dominant):02d}' if dominant else ''),held_out=held,triple_id=ti,reference_id=rid,triple_key=triple,block=block,budget=budget,level=level,dominant_donor=dominant))
    rf = pd.DataFrame(refs.values()); cf = pd.DataFrame(cases); inc = pd.DataFrame(incidence); om = pd.DataFrame(omitted)
    assert len(rf) == 336 and len(cf) == 1344 and rf.reference_id.is_unique and cf.case_id.is_unique
    assert rf.groupby('triple_key').size().eq(24).all() and rf.triple_key.nunique() == 14
    assert cf.reference_id.value_counts().eq(4).all() and inc.reference_donor.value_counts().eq(12).all()
    assert om.omitted_reference_donor.value_counts().eq(1).all()
    # Independently validate the exported membership/quotas against cell metadata and original prefix axes.
    sets = {}
    for ref in rf.itertuples(index=False):
        cols = np.array(list(map(int,ref.reference_columns_R.split('|')))); cm = cells.iloc[cols-1]
        assert len(cols) == len(set(cols)) == 6*ref.budget
        ds = ref.reference_donors.split('|'); assert set(cm.donor) == set(ds)
        helds = set(cf.loc[cf.reference_id==ref.reference_id,'held_out']); assert not helds & set(ds)
        assert len(helds)==4
        for donor in ds:
            quota = ref.budget//3 if ref.level=='balanced' else (5*ref.budget//6 if donor==ref.dominant_donor else ref.budget//12)
            for k, typ in enumerate(inp['cell_types']):
                subset = cm[(cm.donor==donor)&(cm.cell_type==typ)]
                assert len(subset)==quota
                assert np.array_equal(subset.obs_row_0.to_numpy(),prefixes[ref.block,original.index(donor),k,:quota])
        key=(ref.triple_key,ref.block,ref.level,ref.dominant_donor)
        if ref.budget==60: sets[key]=set(cols)
        else: assert sets[key].issubset(set(cols))
    for i, held in enumerate(ring):
        local = inc[inc.held_out==held]
        assert local.reference_donor.nunique()==12 and held not in set(local.reference_donor)
        for ti in range(4):
            expected={ring[(i+1+3*ti+j)%14] for j in range(3)}
            assert set(local.loc[local.triple_id==ti,'reference_donor'])==expected
    replay=[]
    for i,held in enumerate(ring):
        budget=[60,300][(i//2)%2];level=['balanced','ratio10'][i%2]
        chosen=cf[(cf.held_out==held)&(cf.triple_id==i%4)&(cf.block==(i//4)%3)&(cf.budget==budget)&(cf.level==level)]
        if level=='ratio10': chosen=chosen[chosen.dominant_donor==sorted(rf.set_index('reference_id').loc[chosen.iloc[0].reference_id,'reference_donors'].split('|'))[0]]
        assert len(chosen)==1;replay.append(chosen.iloc[0].to_dict())
    assert len({(r['budget'],r['level']) for r in replay})==4
    old_refs=pd.read_csv(ROOT/inp['reference_file'],sep='\t',keep_default_na=False)
    old_om=pd.read_csv(sp.parent/'omitted_reference_donors.csv',keep_default_na=False).set_index('held_out')
    om['original_omitted_reference_donor']=om.held_out.map(old_om.omitted_reference_donor)
    overlap=sorted(set(rf.triple_key)&set(old_refs.triple_key))
    out=ROOT/'results/music_reference_composition_selection'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    rf.to_csv(out/'references.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0},lineterminator='\n')
    cf.to_csv(out/'cases.tsv',sep='\t',index=False,lineterminator='\n')
    inc.to_csv(out/'reference_incidence.csv',index=False,lineterminator='\n');om.to_csv(out/'omitted_donors.csv',index=False,lineterminator='\n')
    save(out/'selection.json',dict(status='frozen plan; follow-up predictions not run',donors=ring,original_prefix_donor_order=original,cell_types=inp['cell_types'],replay_cases=replay,original_triple_overlap=overlap,changed_omissions=int((om.omitted_reference_donor!=om.original_omitted_reference_donor).sum()),checks='All 336 reference quotas/prefix axes/nesting, 1344 case identities and fold incidences passed',source_input=rel(OLD_INPUT),source_input_sha256=sha(OLD_INPUT),sources=[{'path':rel(p),'sha256':sha(p)} for p in [PROTOCOL,Path(__file__),prefix_file,sp]],artifacts=[{'path':rel(p),'sha256':sha(p)} for p in out.iterdir()]))
    new=dict(inp)
    new.update(status='original verified matrices/targets with frozen alternative reference grouping',donors=ring,selection_manifest=rel(out/'selection.json'),selection_manifest_sha256=sha(out/'selection.json'),reference_file=rel(out/'references.tsv.gz'),case_file=rel(out/'cases.tsv'),logical_cases=1344,unique_references=336,scientific_prediction_rows=161280,source_input_manifest=rel(OLD_INPUT),source_input_sha256=sha(OLD_INPUT))
    save(out/'input_manifest.json',new)
    print('OUTPUT',rel(out));print(json.dumps({'ring':ring,'overlapping_unique_triples':len(overlap),'changed_omissions':int((om.omitted_reference_donor!=om.original_omitted_reference_donor).sum()),'references':336,'cases':1344,'predictions_planned':161280},indent=2))

if __name__=='__main__': main()
