"""Reconstruct frozen MC memberships independently from source metadata and seed keys."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, shutil
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*2**20),b''): h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('selection',type=Path);a=ap.parse_args()
    sp=a.selection.resolve();s=json.loads(sp.read_bytes());d=ROOT/s['data_directory']
    out=ROOT/'results/music_mc_selection_check'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True)
    report={'status':'running','selection':sp.relative_to(ROOT).as_posix(),'selection_sha256':sha(sp),'producer_sha256':sha(__file__)}
    shutil.copyfile(__file__,out/Path(__file__).name)
    try:
        for x in s['sources']+s['artifacts']: assert sha(ROOT/x['path'])==x['sha256'],x['path']
        assert s['blocks']==list(range(3,15)) and s['protocol_sha256']=='2cb1aee483eccecfae85e29fbd4c396ac710052889d63a6998ef3c40d1184988'
        meta=json.loads((ROOT/s['source_metadata_audit']).read_bytes());md=ROOT/meta['derived_data_directory']
        cats=json.loads((md/'obs_categories.json').read_bytes());arrays=np.load(md/'obs_columns.npz',allow_pickle=False)
        donors=s['donors'];authors=s['author_types'];types=s['cell_types']
        prefixes=np.load(d/'reference_prefix_rows.npy',allow_pickle=False);assert prefixes.shape==(12,14,6,250)
        all_pools=[]
        for i,donor in enumerate(donors):
            for j,typ in enumerate(authors):
                pool=np.where((arrays['donor_id']==cats['donor_id'].index(donor)) & (arrays['author_cell_type']==cats['author_cell_type'].index(typ)))[0]
                all_pools.append(pool)
                for bi,b in enumerate(range(3,15)):
                    seed=int(hashlib.sha256(f'music-external-v1|reference|{donor}|{typ}|{b}'.encode()).hexdigest()[:32],16)
                    expected=np.random.Generator(np.random.PCG64(seed)).permutation(pool)[:250]
                    np.testing.assert_array_equal(expected,prefixes[bi,i,j])
        assert sum(map(len,all_pools))==114524
        oldsel=json.loads((ROOT/s['source_selection_manifest']).read_bytes());od=ROOT/oldsel['data_directory']
        assert (d/'target_rows.npy').read_bytes()==(od/'target_rows.npy').read_bytes()
        rows=np.load(d/'selected_rows.npy',allow_pickle=False);oldrows=np.load(od/'selected_rows.npy',allow_pickle=False)
        np.testing.assert_array_equal(rows,np.union1d(oldrows,prefixes.ravel()))
        cells=pd.read_csv(d/'cells.tsv',sep='\t',keep_default_na=False);np.testing.assert_array_equal(cells.obs_row_0,rows)
        source_index=pd.read_csv(md/'cell_index.tsv.gz',sep='\t',keep_default_na=False)
        np.testing.assert_array_equal(cells.cell_id,source_index.iloc[rows].cell_id)
        refs=pd.read_csv(d/'references.tsv.gz',sep='\t',keep_default_na=False);cases=pd.read_csv(sp.parent/'cases.tsv',sep='\t',keep_default_na=False)
        assert len(refs)==1344 and refs.reference_id.is_unique and len(cases)==5376 and cases.case_id.is_unique
        assert refs.groupby('triple_key').size().eq(96).all() and refs.triple_key.nunique()==14
        assert cases.groupby('reference_id').size().eq(4).all()
        rowsets={}
        for r in refs.itertuples(index=False):
            names=r.reference_donors.split('|');assert len(set(names))==3
            cols=np.fromstring(r.reference_columns_R,sep='|',dtype=int)-1
            assert len(cols)==r.budget*6 and len(set(cols))==len(cols)
            expected=[]
            for donor in sorted(names):
                n=r.budget//3 if r.level=='balanced' else (r.budget*5//6 if donor==r.dominant_donor else r.budget//12)
                for typ in range(6):expected.extend(prefixes[r.block-3,donors.index(donor),typ,:n])
            np.testing.assert_array_equal(rows[cols],expected)
            cm=cells.iloc[cols];counts=cm.groupby(['donor','cell_type']).size()
            assert len(counts)==18 and set(cm.donor)==set(names) and set(cm.cell_type)==set(types)
            rowsets[(r.triple_key,r.block,r.level,r.dominant_donor,r.budget)]=set(expected)
        for key,small in rowsets.items():
            if key[-1]==60:assert small.issubset(rowsets[key[:-1]+(300,)])
        rbyid=refs.set_index('reference_id')
        for c in cases.itertuples(index=False):
            i=donors.index(c.held_out);others=[donors[(i+j)%14] for j in range(1,14)]
            names=others[3*c.triple_id:3*c.triple_id+3];r=rbyid.loc[c.reference_id]
            assert set(names)==set(r.reference_donors.split('|')) and c.held_out not in names
            for k in ['triple_key','block','budget','level','dominant_donor']: assert getattr(c,k)==r[k]
        report.update(status='passed independent metadata and membership reconstruction',eligible_cells=114524,selected_cells=len(rows),new_cells=len(rows)-len(oldrows),prefixes_reconstructed=1008,reference_matrices=1344,logical_cases=5376,targets_unchanged=True,held_out_exclusion=True,quota_and_nesting_verified=True)
    except Exception as e:
        report.update(status='failed',error=repr(e));raise
    finally:
        (out/'check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report,indent=2));print('OUTPUT',out)
if __name__=='__main__':main()
