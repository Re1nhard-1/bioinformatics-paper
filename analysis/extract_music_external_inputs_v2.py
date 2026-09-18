"""Validate every selected raw count and export frozen external MuSiC inputs."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, gzip, hashlib, json, shutil, time, warnings
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, io
from cellxgene_range_reader import RangeReader

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*2**20),b''):h.update(b)
    return h.hexdigest()


def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')


def write_mtx(path,matrix):
    with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',DeprecationWarning)
            io.mmwrite(handle,matrix,field='integer',symmetry='general')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('selection',type=Path);args=parser.parse_args()
    selection_path=args.selection.resolve();sel=json.loads(selection_path.read_bytes());source=ROOT/sel['data_directory']
    for item in sel['sources']+sel['artifacts']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
    meta_path=ROOT/sel['source_metadata_audit'];meta=json.loads(meta_path.read_bytes());mdata=ROOT/meta['derived_data_directory']
    access=json.loads((ROOT/'results/cellxgene_access/20260917T111504127905Z/access.json').read_bytes())
    assert shutil.disk_usage(ROOT).free>8*2**30
    rows=np.load(source/'selected_rows.npy',allow_pickle=False);target_rows=np.load(source/'target_rows.npy',allow_pickle=False)
    cells=pd.read_csv(source/'cells.tsv',sep='\t',keep_default_na=False);genes=pd.read_csv(mdata/'raw_genes.tsv',sep='\t',keep_default_na=False)
    ptr=np.load(mdata/'raw_indptr.npy',allow_pickle=False)
    assert len(rows)==sel['selected_cells'] and np.array_equal(rows,cells.obs_row_0) and np.all(np.diff(rows)>0)
    assert genes.gene_id.is_unique and len(genes)==30172
    lengths=np.diff(ptr)[rows];nnz=int(lengths.sum());assert nnz==sel['selected_stored_entries']
    outptr=np.r_[0,np.cumsum(lengths,dtype=np.int64)];vals=np.empty(nnz,dtype=np.int32);cols=np.empty(nnz,dtype=np.int32)
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out=ROOT/'results/music_external_input'/run_id;data=ROOT/'data/processed/music_external'/run_id
    out.mkdir(parents=True);data.mkdir(parents=True)
    record={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),'selection_manifest':selection_path.relative_to(ROOT).as_posix(),'selection_manifest_sha256':sha(selection_path),'script_sha256':sha(__file__),'range_reader_sha256':sha(ROOT/'analysis/cellxgene_range_reader.py'),'data_directory':data.relative_to(ROOT).as_posix()}
    save(out/'started.json',record);print('OUTPUT',out.relative_to(ROOT).as_posix(),flush=True)
    started=time.monotonic();checks=[];r=RangeReader(access['file_url'],ROOT/'data/raw/cellxgene_ranges/c55dc602-d168-4d15-acc1-5de4f2f5d551',expected_size=access['head_headers']['Content-Length'],expected_etag=access['head_headers']['ETag'],max_network_bytes=256*2**20,deadline_seconds=2700)
    try:
        with h5py.File(r,'r') as f:
            assert tuple(f['raw/X'].attrs['shape'])==(meta['cells'],len(genes))
            np.testing.assert_array_equal(f['raw/X/indptr'][:],ptr)
            np.testing.assert_array_equal(f['raw/var/_index'].asstr()[:],genes.gene_id.to_numpy())
            ids=f['obs/index'].asstr()[:];np.testing.assert_array_equal(ids[rows],cells.cell_id.to_numpy());del ids
            for wi,window in enumerate(np.unique(rows//10000)):
                if time.monotonic()-started>2700:raise TimeoutError('Selected-count extraction deadline')
                positions=np.flatnonzero(rows//10000==window);first=int(rows[positions[0]]);last=int(rows[positions[-1]])+1
                start=int(ptr[first]);end=int(ptr[last]);raw_vals=f['raw/X/data'][start:end];raw_cols=f['raw/X/indices'][start:end]
                span=sparse.csr_matrix((raw_vals,raw_cols,ptr[first:last+1]-start),shape=(last-first,len(genes)))
                part=span[rows[positions]-first].tocsr();del span,raw_vals,raw_cols
                assert np.isfinite(part.data).all() and (part.data>=0).all() and np.equal(part.data,np.floor(part.data)).all()
                assert ((part.indices>=0)&(part.indices<len(genes))).all()
                assert part.data.max()<=np.iinfo(np.int32).max
                had_sorted=part.has_sorted_indices;part.sort_indices();assert part.has_canonical_format,'Duplicate gene coordinate within a selected cell'
                np.testing.assert_array_equal(np.diff(part.indptr),lengths[positions])
                assert np.all(np.asarray(part.sum(axis=1)).ravel()>0)
                begin=int(outptr[positions[0]]);finish=int(outptr[positions[-1]+1]);assert finish-begin==part.nnz
                vals[begin:finish]=part.data.astype(np.int32);cols[begin:finish]=part.indices
                checks.append({'source_window':int(window),'selected_cells':len(positions),'stored_entries':int(part.nnz),'minimum':float(part.data.min()),'maximum':float(part.data.max()),'source_indices_already_sorted':bool(had_sorted),'valid_finite_nonnegative_integer':True,'unique_gene_coordinates':True})
                if wi%10==0:print('Validated windows',wi+1,'selected cells',int(positions[-1])+1,flush=True)
            counts=sparse.csr_matrix((vals,cols,outptr),shape=(len(rows),len(genes)))
            assert counts.has_canonical_format and counts.nnz==nnz
            # One deterministic source row per donor plus first/last selected rows.
            spot=sorted(set([0,len(rows)-1]+[int(np.flatnonzero(cells.donor.to_numpy()==d)[0]) for d in sel['donors']]))
            for pos in spot:
                row=int(rows[pos]);start=int(ptr[row]);end=int(ptr[row+1])
                direct=sparse.csr_matrix((f['raw/X/data'][start:end],f['raw/X/indices'][start:end],np.array([0,end-start])),shape=(1,len(genes)))
                assert (direct-counts[pos]).nnz==0
        cells['total_counts']=np.asarray(counts.sum(axis=1,dtype=np.int64)).ravel()
        cells.to_csv(data/'cells.tsv',sep='\t',index=False,lineterminator='\n');pd.Series(genes.gene_id).to_csv(data/'genes.tsv',sep='\t',index=False,header=False,lineterminator='\n')
        sparse.save_npz(data/'counts.npz',counts,compressed=True)
        donor_matrices=[]
        for donor in sel['donors']:
            positions=np.flatnonzero(cells.donor.to_numpy()==donor);path=data/f'counts_{donor}.mtx.gz'
            write_mtx(path,counts[positions].T.tocoo());donor_matrices.append({'donor':donor,'path':path.relative_to(ROOT).as_posix(),'cells':len(positions),'genes':len(genes),'stored_entries':int(counts[positions].nnz),'sha256':sha(path)})
            print('Exported donor',donor,flush=True)
        compositions=pd.read_csv(selection_path.parent/'target_compositions.csv')[[f'count_{i}' for i in range(6)]].to_numpy(dtype=np.int32)
        matrix=np.empty((len(sel['donors'])*60,len(genes)),dtype=np.int64);targets=[]
        target_positions=np.searchsorted(rows,target_rows);assert np.array_equal(rows[target_positions],target_rows)
        type_codes=pd.Categorical(cells.cell_type,categories=sel['cell_types']).codes
        for di,donor in enumerate(sel['donors']):
            for m in range(60):
                idx=target_positions[di,m];assert len(set(idx))==300 and cells.iloc[idx].donor.eq(donor).all()
                actual=np.bincount(type_codes[idx],minlength=6);np.testing.assert_array_equal(actual,compositions[m])
                total=np.asarray(counts[idx].sum(axis=0,dtype=np.int64)).ravel();assert int(total.sum())==int(cells.iloc[idx].total_counts.sum())
                matrix[di*60+m]=total
                target={'target_name':f'{donor}_{m:02d}','held_out':donor,'mixture_id':m,'total_counts':int(total.sum())}
                target.update({f'true_{i}':float(actual[i]/300) for i in range(6)});targets.append(target)
        write_mtx(data/'targets.mtx.gz',sparse.coo_matrix(matrix.T));np.save(data/'target_positions.npy',target_positions,allow_pickle=False)
        pd.DataFrame(targets).to_csv(data/'targets.tsv',sep='\t',index=False,lineterminator='\n')
        references=pd.read_csv(source/'references.tsv.gz',sep='\t',keep_default_na=False);references.to_csv(data/'references.tsv',sep='\t',index=False,lineterminator='\n')
        shutil.copyfile(selection_path.parent/'cases.tsv',data/'cases.tsv')
        pd.DataFrame(checks).to_csv(out/'selected_count_validation.csv',index=False,lineterminator='\n')
        pd.DataFrame({'matrix_column_R':np.array(spot)+1,'source_row_0':rows[spot],'direct_hdf5_matches':True}).to_csv(out/'direct_row_checks.csv',index=False,lineterminator='\n')
        record.update(status='complete selected raw counts and targets; no predictions',donors=sel['donors'],cell_types=sel['cell_types'],selected_cells=len(rows),genes=len(genes),selected_stored_entries=nnz,all_selected_entries_validated=True,whole_source_matrix_validated=False,raw_count_minimum=int(vals.min()),raw_count_maximum=int(vals.max()),positive_cell_libraries=True,cell_library_minimum=int(cells.total_counts.min()),cell_library_maximum=int(cells.total_counts.max()),targets=840,logical_cases=5040,unique_references=1260,scientific_prediction_rows=604800,direct_rows_checked=len(spot),donor_matrices=donor_matrices,case_file=(data/'cases.tsv').relative_to(ROOT).as_posix(),reference_file=(data/'references.tsv').relative_to(ROOT).as_posix())
    except Exception as error:
        record.update(status='failed; partial inputs retained',error=repr(error));raise
    finally:
        save(out/'range_access.json',r.report());r.close()
        record.update(finished_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started)
        record['artifacts']=[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for folder in [out,data] for p in sorted(folder.iterdir()) if p.is_file()]
        save(out/'input_manifest.json',record)
        print(record['status'],'elapsed',record['elapsed_seconds'],flush=True)


if __name__=='__main__':main()
