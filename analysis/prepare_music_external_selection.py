"""Freeze outcome-free external donor, reference and target cell selections."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / 'results/cellxgene_metadata/20260917T111912321314Z/audit.json'
INVENTORY = ROOT / 'results/cellxgene_inventory/20260917T113340800025Z/audit.json'
OLD = ROOT / 'results/music_input/20260917T090438982999Z/input_manifest.json'
TYPE_MAP = {'B': 'B cells', 'cM': 'CD14+ Monocytes', 'T4': 'CD4 T cells', 'T8': 'CD8 T cells', 'ncM': 'FCGR3A+ Monocytes', 'NK': 'NK cells'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 2**20), b''): h.update(b)
    return h.hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')


def rng(key):
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(key.encode('utf-8')).digest()[:16], 'big')))


def main():
    meta = json.loads(META.read_bytes()); inventory = json.loads(INVENTORY.read_bytes()); old = json.loads(OLD.read_bytes())
    for obj in [meta, inventory]:
        for item in obj['artifacts']: assert sha(ROOT / item['path']) == item['sha256'], item['path']
    types = old['cell_types']; authors = list(TYPE_MAP)
    assert types == list(TYPE_MAP.values()), types
    donor_table = pd.read_csv(INVENTORY.parent / 'donor_six_type_inventory.csv', dtype={'donor_id': str}, keep_default_na=False)
    cohort = donor_table[(donor_table.disease == 'normal') & (donor_table.disease_state == 'na') & (donor_table.minimum_six_type_cells >= 250)].copy()
    ring = sorted(cohort.donor_id, key=lambda d: hashlib.sha256(f'music-external-v1|ring|{d}'.encode()).hexdigest())
    assert len(ring) == len(set(ring)) == 14
    assert not set(ring) & {'101','1015','1016','1244','1256','1488','107','1039'}
    data_source = ROOT / meta['derived_data_directory']
    cats = json.loads((data_source / 'obs_categories.json').read_bytes()); arrays = np.load(data_source / 'obs_columns.npz', allow_pickle=False)
    indices = pd.read_csv(data_source / 'cell_index.tsv.gz', sep='\t', keep_default_na=False)
    assert np.array_equal(indices.obs_row_0, np.arange(meta['cells'])) and indices.cell_id.is_unique
    original_targets = ROOT / old['data_directory'] / 'targets.tsv'
    expected_hash = next(x['sha256'] for x in old['artifacts'] if x['path'] == original_targets.relative_to(ROOT).as_posix())
    assert sha(original_targets) == expected_hash
    old_targets = pd.read_csv(original_targets, sep='\t', dtype={'held_out': str})
    composition = None
    for donor in old['donors']:
        truth = old_targets[old_targets.held_out == donor].sort_values('mixture_id')[[f'true_{i}' for i in range(6)]].to_numpy()
        counts = np.rint(truth * 300).astype(np.int32)
        assert counts.shape == (60,6) and np.max(abs(counts / 300 - truth)) < 1e-14 and np.all(counts.sum(1) == 300)
        if composition is None: composition = counts
        else: np.testing.assert_array_equal(composition, counts)
    assert np.array_equal(composition.max(0), [179,159,201,185,201,200])
    prefixes = np.empty((3,14,6,250), dtype=np.int32)
    target_rows = np.empty((14,60,300), dtype=np.int32)
    for di, donor in enumerate(ring):
        donor_mask = arrays['donor_id'] == cats['donor_id'].index(donor)
        assert np.all(arrays['disease'][donor_mask] == cats['disease'].index('normal'))
        assert np.all(arrays['disease_state'][donor_mask] == cats['disease_state'].index('na'))
        assert len(np.unique(arrays['sample_uuid'][donor_mask])) == len(np.unique(arrays['suspension_uuid'][donor_mask])) == 1
        for ti, author in enumerate(authors):
            pool = np.flatnonzero(donor_mask & (arrays['author_cell_type'] == cats['author_cell_type'].index(author)))
            assert len(pool) == int(cohort.set_index('donor_id').loc[donor, author]) and len(pool) >= 250
            for block in range(3):
                prefixes[block,di,ti] = rng(f'music-external-v1|reference|{donor}|{author}|{block}').permutation(pool)[:250]
            for m in range(60):
                offset = int(composition[m,:ti].sum()); count = int(composition[m,ti])
                assert count <= len(pool)
                target_rows[di,m,offset:offset+count] = rng(f'music-external-v1|target|{donor}|{author}|{m}').permutation(pool)[:count]
        assert all(len(set(row)) == 300 for row in target_rows[di])
    selected_rows = np.unique(np.concatenate([prefixes.ravel(), target_rows.ravel()]))
    selected = indices.iloc[selected_rows].copy().reset_index(drop=True)
    selected.insert(0,'matrix_column_R',np.arange(1,len(selected)+1))
    for key in ['donor_id','author_cell_type','Processing_Cohort','sample_uuid','suspension_uuid','library_uuid','disease','disease_state']:
        selected[key] = np.asarray(cats[key],dtype=str)[arrays[key][selected_rows]]
    selected['donor'] = selected.donor_id; selected['cell_type'] = selected.author_cell_type.map(TYPE_MAP)
    assert selected.cell_type.notna().all()
    cases=[]; memberships=[]; refs={}; omitted=[]
    for di, held in enumerate(ring):
        others=[ring[(di+j)%14] for j in range(1,14)]; omitted.append({'held_out':held,'omitted_reference_donor':others[-1]})
        for ti in range(4):
            triple=others[ti*3:ti*3+3]; triple_key='_'.join(f'd{ring.index(d):02d}' for d in sorted(triple))
            memberships.extend({'held_out':held,'triple_id':ti,'reference_donor':d,'triple_key':triple_key} for d in triple)
            for block in range(3):
                for budget in [60,120,300]:
                    allocations=[('balanced','',budget//3,budget//3)]+[(level,d,major,minor) for level,major,minor in [('ratio2',budget//2,budget//4),('ratio4',2*budget//3,budget//6),('ratio10',5*budget//6,budget//12)] for d in sorted(triple)]
                    for level,dominant,major,minor in allocations:
                        ref_id=f'{triple_key}_b{block}_n{budget}_{level}' + (f'_d{ring.index(dominant):02d}' if dominant else '')
                        if ref_id not in refs:
                            rows=np.concatenate([prefixes[block,ring.index(d),t,:major if d==dominant else minor] for d in sorted(triple) for t in range(6)])
                            assert len(rows)==len(set(rows))==6*budget
                            pos=np.searchsorted(selected_rows,rows);assert np.array_equal(selected_rows[pos],rows)
                            refs[ref_id]={'reference_id':ref_id,'triple_key':triple_key,'reference_donors':'|'.join(sorted(triple)),'block':block,'budget':budget,'level':level,'dominant_donor':dominant,'major_count':major,'minor_count':minor,'reference_cells':6*budget,'reference_columns_R':'|'.join(map(str,pos+1))}
                        cases.append({'case_id':f'd{di:02d}_t{ti}_b{block}_n{budget}_{level}'+(f'_d{ring.index(dominant):02d}' if dominant else ''),'held_out':held,'triple_id':ti,'reference_id':ref_id,'triple_key':triple_key,'block':block,'budget':budget,'level':level,'dominant_donor':dominant})
    assert len(cases)==5040 and len(refs)==1260
    logical=pd.DataFrame(cases); incidence=pd.DataFrame(memberships)
    assert logical.reference_id.value_counts().eq(4).all()
    assert incidence.groupby('reference_donor').size().eq(12).all()
    assert all(d not in set(g.reference_donor) for d,g in incidence.groupby('held_out'))
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'); out=ROOT/'results/music_external_selection'/run_id; data=ROOT/'data/processed/music_external_selection'/run_id
    out.mkdir(parents=True);data.mkdir(parents=True)
    np.save(data/'reference_prefix_rows.npy',prefixes,allow_pickle=False);np.save(data/'target_rows.npy',target_rows,allow_pickle=False);np.save(data/'selected_rows.npy',selected_rows,allow_pickle=False)
    selected.to_csv(data/'cells.tsv',sep='\t',index=False,lineterminator='\n')
    pd.DataFrame(refs.values()).to_csv(data/'references.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0},lineterminator='\n')
    pd.DataFrame({'mixture_id':range(60),**{f'count_{i}':composition[:,i] for i in range(6)}}).to_csv(out/'target_compositions.csv',index=False,lineterminator='\n')
    logical.to_csv(out/'cases.tsv',sep='\t',index=False,lineterminator='\n');incidence.to_csv(out/'reference_incidence.csv',index=False,lineterminator='\n');pd.DataFrame(omitted).to_csv(out/'omitted_reference_donors.csv',index=False,lineterminator='\n')
    cohort.set_index('donor_id').loc[ring].reset_index().to_csv(out/'cohort_inventory.csv',index=False,lineterminator='\n')
    summary=selected.groupby(['donor','author_cell_type','Processing_Cohort'],observed=True).size().rename('selected_cells').reset_index();summary.to_csv(out/'selected_cell_covariates.csv',index=False,lineterminator='\n')
    meta_checks=pd.read_csv(INVENTORY.parent/'donor_field_consistency.csv',dtype={'donor_id':str},keep_default_na=False);meta_checks[meta_checks.donor_id.isin(ring)].to_csv(out/'cohort_covariates.csv',index=False,lineterminator='\n')
    record={'status':'frozen selections; no external predictions','run_id':run_id,'donors':ring,'cell_types':types,'author_types':authors,'data_directory':data.relative_to(ROOT).as_posix(),'source_metadata_audit':META.relative_to(ROOT).as_posix(),'source_inventory_audit':INVENTORY.relative_to(ROOT).as_posix(),'logical_cases':5040,'unique_reference_matrices':1260,'targets':840,'scientific_prediction_rows':604800,'selected_cells':len(selected_rows),'selected_stored_entries':int(np.diff(np.load(data_source/'raw_indptr.npy',allow_pickle=False))[selected_rows].sum()),'seed_algorithm':'SHA256 first16 bytes unsigned big endian -> PCG64; see frozen protocol','sources':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in [META,INVENTORY,OLD,original_targets,Path(__file__),ROOT/'analysis/MUSIC_EXTERNAL_PROTOCOL.md']],'artifacts':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for folder in [out,data] for p in sorted(folder.iterdir())]}
    save(out/'selection.json',record)
    print('OUTPUT',out.relative_to(ROOT).as_posix());print(json.dumps({k:v for k,v in record.items() if k not in ['sources','artifacts']},indent=2))


if __name__=='__main__': main()
