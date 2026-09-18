"""Extend verified external counts only for missing reference cells; preserve targets."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import json
import shutil
import time
import warnings
import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse
from cellxgene_range_reader import RangeReader

ROOT = Path(__file__).resolve().parents[1]
ACCESS = ROOT / "results/cellxgene_access/20260917T111504127905Z/access.json"
CACHE = ROOT / "data/raw/cellxgene_ranges/c55dc602-d168-4d15-acc1-5de4f2f5d551"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def verify(items):
    for item in items:
        assert sha(ROOT / item["path"]) == item["sha256"], item["path"]


def write_mtx(path, matrix):
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            io.mmwrite(handle, matrix, field="integer", symmetry="general")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", type=Path)
    parser.add_argument("--max-network-mib", type=int, default=256)
    args = parser.parse_args()
    assert 0 <= args.max_network_mib <= 256
    selection_path = args.selection.resolve()
    sel = json.loads(selection_path.read_bytes())
    verify(sel["sources"] + sel["artifacts"])
    assert sel["blocks"] == list(range(3, 15)) and sel["logical_cases"] == 5376
    old_path = ROOT / sel["source_input_manifest"]
    assert sha(old_path) == sel["source_input_sha256"]
    old = json.loads(old_path.read_bytes())
    verify(old["artifacts"])
    old_data = ROOT / old["data_directory"]
    source = ROOT / sel["data_directory"]
    meta = json.loads((ROOT / sel["source_metadata_audit"]).read_bytes())
    metadata_dir = ROOT / meta["derived_data_directory"]
    access = json.loads(ACCESS.read_bytes())
    assert access["file_url"] == "https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad"
    assert int(access["head_headers"]["Content-Length"]) == 12218105530
    assert shutil.disk_usage(ROOT).free > 8 * 2**30
    rows = np.load(source / "selected_rows.npy", allow_pickle=False)
    target_rows = np.load(source / "target_rows.npy", allow_pickle=False)
    cells = pd.read_csv(source / "cells.tsv", sep="\t", keep_default_na=False)
    old_cells = pd.read_csv(old_data / "cells.tsv", sep="\t", keep_default_na=False)
    old_rows = old_cells.obs_row_0.to_numpy()
    old_positions = np.searchsorted(rows, old_rows)
    np.testing.assert_array_equal(rows[old_positions], old_rows)
    np.testing.assert_array_equal(cells.iloc[old_positions].cell_id.to_numpy(), old_cells.cell_id.to_numpy())
    np.testing.assert_array_equal(cells.iloc[old_positions].donor.to_numpy(), old_cells.donor.to_numpy())
    np.testing.assert_array_equal(cells.iloc[old_positions].cell_type.to_numpy(), old_cells.cell_type.to_numpy())
    genes = (old_data / "genes.tsv").read_text(encoding="utf-8").splitlines()
    assert len(genes) == len(set(genes)) == 30172
    assert len(rows) == sel["selected_cells"] and np.array_equal(rows, cells.obs_row_0)
    assert np.all(np.diff(rows) > 0) and cells.cell_id.is_unique
    ptr = np.load(metadata_dir / "raw_indptr.npy", allow_pickle=False)
    lengths = np.diff(ptr)[rows]
    nnz = int(lengths.sum())
    assert nnz == sel["selected_stored_entries"]
    present = np.zeros(len(rows), dtype=bool)
    present[old_positions] = True
    missing_positions = np.flatnonzero(~present)
    assert len(missing_positions) == sel["new_cells_requiring_counts"]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = ROOT / "results/music_mc_extension_input" / stamp
    data = ROOT / "data/processed/music_mc_extension" / stamp
    out.mkdir(parents=True); data.mkdir(parents=True)
    shutil.copyfile(Path(__file__), out / Path(__file__).name)
    started = time.monotonic()
    record = {"status": "running count extension; no predictions", "started_utc": datetime.now(timezone.utc).isoformat(), "selection_manifest": rel(selection_path), "selection_manifest_sha256": sha(selection_path), "source_input_manifest": rel(old_path), "source_input_sha256": sha(old_path), "data_directory": rel(data), "script_sha256": sha(__file__), "range_reader_sha256": sha(ROOT / "analysis/cellxgene_range_reader.py"), "source_access_sha256": sha(ACCESS), "blocks": sel["blocks"]}
    save(out / "started.json", record)
    print("OUTPUT", rel(out), flush=True)
    cache_identity = json.loads((CACHE / "identity.json").read_bytes())
    assert cache_identity["url"] == access["file_url"] and cache_identity["size"] == int(access["head_headers"]["Content-Length"]) and cache_identity["etag"] == access["head_headers"]["ETag"]
    cache_before = {"identity": cache_identity, "cached_blocks": len(list(CACHE.glob("block_*.bin"))), "cached_bytes": sum(p.stat().st_size for p in CACHE.glob("block_*.bin"))}
    save(out / "cache_before.json", cache_before)
    print("Cache", cache_before["cached_blocks"], "blocks", cache_before["cached_bytes"], "bytes; missing cells", len(missing_positions), flush=True)
    reader = RangeReader(access["file_url"], CACHE, expected_size=cache_identity["size"], expected_etag=cache_identity["etag"], max_network_bytes=args.max_network_mib * 2**20, deadline_seconds=2700)
    try:
        original = sparse.load_npz(old_data / "counts.npz").tocsr()
        assert original.shape == (93892, 30172) and original.nnz == 75749702 and original.has_canonical_format
        assert np.issubdtype(original.dtype, np.integer)
        np.testing.assert_array_equal(np.diff(original.indptr), lengths[old_positions])
        outptr = np.r_[0, np.cumsum(lengths, dtype=np.int64)]
        vals = np.empty(nnz, dtype=np.int32)
        cols = np.empty(nnz, dtype=np.int32)
        for old_i, new_i in enumerate(old_positions):
            old_start, old_end = original.indptr[old_i:old_i + 2]
            start, end = outptr[new_i:new_i + 2]
            vals[start:end] = original.data[old_start:old_end]
            cols[start:end] = original.indices[old_start:old_end]
        checks = []
        with h5py.File(reader, "r") as h5:
            assert tuple(h5["raw/X"].attrs["shape"]) == (meta["cells"], 30172)
            np.testing.assert_array_equal(h5["raw/X/indptr"][:], ptr)
            np.testing.assert_array_equal(h5["raw/var/_index"].asstr()[:], np.asarray(genes))
            source_ids = h5["obs/index"].asstr()[:]
            np.testing.assert_array_equal(source_ids[rows], cells.cell_id.to_numpy())
            del source_ids
            windows = np.unique(rows[missing_positions] // 10000)
            for wi, window in enumerate(windows):
                if time.monotonic() - started > 2700:
                    raise TimeoutError("Count extraction exceeded technical time bound")
                positions = missing_positions[rows[missing_positions] // 10000 == window]
                first, last = int(rows[positions[0]]), int(rows[positions[-1]]) + 1
                raw_start, raw_end = int(ptr[first]), int(ptr[last])
                raw_vals, raw_cols = h5["raw/X/data"][raw_start:raw_end], h5["raw/X/indices"][raw_start:raw_end]
                span = sparse.csr_matrix((raw_vals, raw_cols, ptr[first:last + 1] - raw_start), shape=(last - first, len(genes)))
                part = span[rows[positions] - first].tocsr()
                del span, raw_vals, raw_cols
                assert np.isfinite(part.data).all() and (part.data >= 0).all() and np.equal(part.data, np.floor(part.data)).all()
                assert np.all((part.indices >= 0) & (part.indices < len(genes)))
                assert part.data.max() <= np.iinfo(np.int32).max
                part.sort_indices()
                assert part.has_canonical_format, "Duplicate gene coordinate in new source cell"
                np.testing.assert_array_equal(np.diff(part.indptr), lengths[positions])
                assert np.all(np.asarray(part.sum(axis=1)).ravel() > 0)
                for local_i, new_i in enumerate(positions):
                    begin, finish = part.indptr[local_i:local_i + 2]
                    start, end = outptr[new_i:new_i + 2]
                    vals[start:end] = part.data[begin:finish].astype(np.int32)
                    cols[start:end] = part.indices[begin:finish]
                present[positions] = True
                checks.append({"source_window": int(window), "new_cells": len(positions), "stored_entries": int(part.nnz), "minimum": float(part.data.min()), "maximum": float(part.data.max()), "finite_nonnegative_integer": True, "unique_gene_coordinates": True, "rowptr_lengths_match": True})
                if wi % 10 == 0 or wi == len(windows) - 1:
                    print("Validated missing-cell windows", wi + 1, "/", len(windows), "; network bytes", reader.network_bytes, flush=True)
            # Direct source checks cover both newly extracted and reused rows.
            spots = sorted(set([int(p) for p in missing_positions[[0, -1]]] + [int(old_positions[np.flatnonzero(old_cells.donor.to_numpy() == d)[0]]) for d in sel["donors"]]))
            for position in spots:
                row = int(rows[position]); begin, finish = int(ptr[row]), int(ptr[row + 1])
                direct = sparse.csr_matrix((h5["raw/X/data"][begin:finish], h5["raw/X/indices"][begin:finish], np.array([0, finish - begin])), shape=(1, len(genes)))
                start, end = outptr[position:position + 2]
                observed = sparse.csr_matrix((vals[start:end], cols[start:end], np.array([0, end - start])), shape=(1, len(genes)))
                assert (direct - observed).nnz == 0
        assert present.all()
        counts = sparse.csr_matrix((vals, cols, outptr), shape=(len(rows), len(genes)))
        assert counts.has_canonical_format and counts.nnz == nnz
        overlap_checks = []
        for donor in sel["donors"]:
            old_ids = np.flatnonzero(old_cells.donor.to_numpy() == donor)
            delta = counts[old_positions[old_ids]] - original[old_ids]
            assert delta.nnz == 0, "Old-cell count changed"
            overlap_checks.append({"donor": donor, "old_cells": len(old_ids), "count_difference_nnz": int(delta.nnz)})
        del original
        cells["total_counts"] = np.asarray(counts.sum(axis=1, dtype=np.int64)).ravel()
        assert cells.total_counts.gt(0).all()
        np.testing.assert_array_equal(cells.iloc[old_positions].total_counts.to_numpy(), old_cells.total_counts.to_numpy())
        cells.to_csv(data / "cells.tsv", sep="\t", index=False, lineterminator="\n")
        sparse.save_npz(data / "counts.npz", counts, compressed=True)
        donor_matrices = []
        for donor in sel["donors"]:
            positions = np.flatnonzero(cells.donor.to_numpy() == donor)
            matrix = counts[positions]
            path = data / f"counts_{donor}.mtx.gz"
            write_mtx(path, matrix.T.tocoo())
            donor_matrices.append({"donor": donor, "path": rel(path), "cells": len(positions), "genes": len(genes), "stored_entries": int(matrix.nnz), "sha256": sha(path)})
            print("Exported donor", donor, flush=True)
        preserved = []
        for name in ["genes.tsv", "targets.tsv", "targets.mtx.gz"]:
            shutil.copyfile(old_data / name, data / name)
            assert sha(data / name) == sha(old_data / name)
            preserved.append({"file": name, "old_path": rel(old_data / name), "new_path": rel(data / name), "sha256": sha(data / name), "byte_identical": True})
        shutil.copyfile(source / "target_rows.npy", data / "target_rows.npy")
        target_positions = np.searchsorted(rows, target_rows)
        np.testing.assert_array_equal(rows[target_positions], target_rows)
        np.save(data / "target_positions.npy", target_positions, allow_pickle=False)
        targets = pd.read_csv(data / "targets.tsv", sep="\t", keep_default_na=False)
        for di, donor in enumerate(sel["donors"]):
            for m in range(60):
                positions = target_positions[di, m]
                assert len(np.unique(positions)) == 300 and cells.iloc[positions].donor.eq(donor).all()
                observed = cells.iloc[positions].cell_type.value_counts().reindex(sel["cell_types"], fill_value=0).to_numpy()
                truth = targets[(targets.held_out == donor) & (targets.mixture_id == m)][[f"true_{i}" for i in range(6)]].to_numpy()[0]
                np.testing.assert_array_equal(observed, np.rint(truth * 300).astype(int))
        refs = pd.read_csv(source / "references.tsv.gz", sep="\t", keep_default_na=False)
        refs.to_csv(data / "references.tsv", sep="\t", index=False, lineterminator="\n")
        shutil.copyfile(selection_path.parent / "cases.tsv", data / "cases.tsv")
        pd.DataFrame(checks).to_csv(out / "new_count_checks.csv", index=False, lineterminator="\n")
        pd.DataFrame(overlap_checks).to_csv(out / "old_count_identity.csv", index=False, lineterminator="\n")
        save(out / "preserved_targets_and_genes.json", {"files": preserved, "target_rows_byte_identical": sha(data / "target_rows.npy") == sha(source / "target_rows.npy"), "target_positions_remapped_only": True})
        record.update(status="complete MC extension inputs; no predictions", donors=sel["donors"], cell_types=sel["cell_types"], selected_cells=len(rows), complete_eligible_cells=114524, genes=len(genes), selected_stored_entries=nnz, new_cells=len(missing_positions), new_stored_entries=int(lengths[missing_positions].sum()), reused_old_cells=len(old_rows), all_old_counts_exact=True, all_new_entries_validated=True, whole_source_matrix_validated=False, raw_count_minimum=int(vals.min()), raw_count_maximum=int(vals.max()), positive_cell_libraries=True, targets=840, logical_cases=5376, unique_references=1344, scientific_prediction_rows=645120, weighted_fits=322560, direct_rows_checked=len(spots), donor_matrices=donor_matrices, reference_file=rel(data / "references.tsv"), case_file=rel(data / "cases.tsv"), preserved_targets_and_genes=preserved)
    except Exception as error:
        record.update(status="failed; partial inputs retained", error=repr(error))
        raise
    finally:
        save(out / "range_access.json", reader.report())
        reader.close()
        record.update(finished_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - started)
        record["artifacts"] = [{"path": rel(p), "sha256": sha(p), "bytes": p.stat().st_size} for folder in [out, data] for p in sorted(folder.iterdir()) if p.is_file()]
        save(out / "input_manifest.json", record)
        print(record["status"], "elapsed", record["elapsed_seconds"], flush=True)


if __name__ == "__main__":
    main()
