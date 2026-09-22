"""Freeze twelve new reference draws from complete eligible source-cell pools."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OLD_INPUT = ROOT / "results/music_external_input/20260917T122709190802Z/input_manifest.json"
PROTOCOL = ROOT / "analysis/MUSIC_MC_EXTENSION_PROTOCOL_V1.md"
BLOCKS = list(range(3, 15))


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 2**20), b""):
            value.update(chunk)
    return value.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def verify_records(records):
    for item in records:
        assert sha(ROOT / item["path"]) == item["sha256"], item["path"]


def rng(key):
    seed = int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:16], "big")
    return np.random.Generator(np.random.PCG64(seed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    assert len(args.protocol_sha256) == 64 and sha(PROTOCOL) == args.protocol_sha256, "Protocol is not the frozen authorized version"
    old = json.loads(OLD_INPUT.read_bytes())
    verify_records(old["artifacts"])
    old_selection_path = ROOT / old["selection_manifest"]
    assert sha(old_selection_path) == old["selection_manifest_sha256"]
    old_selection = json.loads(old_selection_path.read_bytes())
    verify_records(old_selection["sources"] + old_selection["artifacts"])
    old_data = ROOT / old["data_directory"]
    old_selection_data = ROOT / old_selection["data_directory"]
    metadata_path = ROOT / old_selection["source_metadata_audit"]
    inventory_path = ROOT / old_selection["source_inventory_audit"]
    metadata = json.loads(metadata_path.read_bytes())
    inventory = json.loads(inventory_path.read_bytes())
    verify_records(metadata["artifacts"])
    verify_records(inventory["artifacts"])
    source = ROOT / metadata["derived_data_directory"]
    categories = json.loads((source / "obs_categories.json").read_bytes())
    arrays = np.load(source / "obs_columns.npz", allow_pickle=False)
    indices = pd.read_csv(source / "cell_index.tsv.gz", sep="\t", keep_default_na=False)
    assert indices.cell_id.is_unique and np.array_equal(indices.obs_row_0, np.arange(metadata["cells"]))
    donors, types, authors = old_selection["donors"], old_selection["cell_types"], old_selection["author_types"]
    assert len(donors) == 14 and len(types) == len(authors) == 6
    type_map = dict(zip(authors, types))
    cohort = pd.read_csv(old_selection_path.parent / "cohort_inventory.csv", keep_default_na=False).set_index("donor_id").loc[donors]
    assert int(cohort[authors].to_numpy().sum()) == 114524
    assert list(donors) == sorted(donors, key=lambda d: hashlib.sha256(f"music-external-v1|ring|{d}".encode()).hexdigest())
    old_prefixes = np.load(old_selection_data / "reference_prefix_rows.npy", allow_pickle=False)
    prefixes = np.empty((12, 14, 6, 250), dtype=np.int32)
    source_pool_checks = []
    for di, donor in enumerate(donors):
        donor_mask = arrays["donor_id"] == categories["donor_id"].index(donor)
        assert np.all(arrays["disease"][donor_mask] == categories["disease"].index("normal"))
        assert np.all(arrays["disease_state"][donor_mask] == categories["disease_state"].index("na"))
        for ti, author in enumerate(authors):
            pool = np.flatnonzero(donor_mask & (arrays["author_cell_type"] == categories["author_cell_type"].index(author)))
            assert len(pool) == int(cohort.loc[donor, author]) and len(pool) >= 250
            for block in range(3):
                reproduced = rng(f"music-external-v1|reference|{donor}|{author}|{block}").permutation(pool)[:250]
                np.testing.assert_array_equal(reproduced, old_prefixes[block, di, ti])
            for bi, block in enumerate(BLOCKS):
                prefixes[bi, di, ti] = rng(f"music-external-v1|reference|{donor}|{author}|{block}").permutation(pool)[:250]
                assert len(np.unique(prefixes[bi, di, ti])) == 250
            source_pool_checks.append({"donor": donor, "author_cell_type": author, "eligible_cells": len(pool), "old_three_prefixes_exact": True})
    old_rows = np.load(old_selection_data / "selected_rows.npy", allow_pickle=False)
    target_rows = np.load(old_selection_data / "target_rows.npy", allow_pickle=False)
    selected_rows = np.union1d(old_rows, prefixes.ravel())
    assert np.isin(target_rows, selected_rows).all() and len(old_rows) == 93892
    selected = indices.iloc[selected_rows].copy().reset_index(drop=True)
    selected.insert(0, "matrix_column_R", np.arange(1, len(selected) + 1))
    for key in ["donor_id", "author_cell_type", "Processing_Cohort", "sample_uuid", "suspension_uuid", "library_uuid", "disease", "disease_state"]:
        selected[key] = np.asarray(categories[key], dtype=str)[arrays[key][selected_rows]]
    selected["donor"] = selected.donor_id
    selected["cell_type"] = selected.author_cell_type.map(type_map)
    assert selected.cell_type.notna().all()
    references, cases, memberships = {}, [], []
    for di, held in enumerate(donors):
        others = [donors[(di + j) % 14] for j in range(1, 14)]
        for ti in range(4):
            triple = others[ti * 3:ti * 3 + 3]
            triple_key = "_".join(f"d{donors.index(d):02d}" for d in sorted(triple))
            assert held not in triple
            memberships.extend({"held_out": held, "triple_id": ti, "reference_donor": d, "triple_key": triple_key} for d in triple)
            for bi, block in enumerate(BLOCKS):
                for budget in [60, 300]:
                    allocations = [("balanced", "", budget // 3, budget // 3)] + [("ratio10", d, 5 * budget // 6, budget // 12) for d in sorted(triple)]
                    for level, dominant, major, minor in allocations:
                        suffix = f"_d{donors.index(dominant):02d}" if dominant else ""
                        ref_id = f"{triple_key}_b{block}_n{budget}_{level}" + suffix
                        if ref_id not in references:
                            rows = np.concatenate([prefixes[bi, donors.index(d), t, :major if d == dominant else minor] for d in sorted(triple) for t in range(6)])
                            assert len(rows) == len(np.unique(rows)) == budget * 6
                            positions = np.searchsorted(selected_rows, rows)
                            np.testing.assert_array_equal(selected_rows[positions], rows)
                            references[ref_id] = {"reference_id": ref_id, "triple_key": triple_key, "reference_donors": "|".join(sorted(triple)), "block": block, "budget": budget, "level": level, "dominant_donor": dominant, "major_count": major, "minor_count": minor, "reference_cells": 6 * budget, "reference_columns_R": "|".join(map(str, positions + 1))}
                        cases.append({"case_id": f"d{di:02d}_t{ti}_b{block}_n{budget}_{level}" + suffix, "held_out": held, "triple_id": ti, "reference_id": ref_id, "triple_key": triple_key, "block": block, "budget": budget, "level": level, "dominant_donor": dominant})
    logical = pd.DataFrame(cases)
    incidence = pd.DataFrame(memberships)
    old_incidence = pd.read_csv(old_selection_path.parent / "reference_incidence.csv", keep_default_na=False)
    pd.testing.assert_frame_equal(incidence, old_incidence)
    assert len(logical) == 5376 and len(references) == 1344 and logical.case_id.is_unique
    assert logical.reference_id.value_counts().eq(4).all() and logical.triple_key.nunique() == 14
    assert logical.groupby("triple_key").size().eq(384).all()
    assert set(logical.block) == set(BLOCKS)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = ROOT / "results/music_mc_extension_selection" / stamp
    data = ROOT / "data/processed/music_mc_extension_selection" / stamp
    out.mkdir(parents=True); data.mkdir(parents=True)
    for name, value in [("reference_prefix_rows.npy", prefixes), ("selected_rows.npy", selected_rows)]:
        np.save(data / name, value, allow_pickle=False)
    shutil.copyfile(old_selection_data / "target_rows.npy", data / "target_rows.npy")
    selected.to_csv(data / "cells.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(references.values()).to_csv(data / "references.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0}, lineterminator="\n")
    logical.to_csv(out / "cases.tsv", sep="\t", index=False, lineterminator="\n")
    incidence.to_csv(out / "reference_incidence.csv", index=False, lineterminator="\n")
    pd.DataFrame(source_pool_checks).to_csv(out / "source_pool_checks.csv", index=False, lineterminator="\n")
    for name in ["target_compositions.csv", "cohort_inventory.csv", "omitted_reference_donors.csv"]:
        shutil.copyfile(old_selection_path.parent / name, out / name)
    shutil.copyfile(Path(__file__), out / Path(__file__).name)
    source_paths = [PROTOCOL, Path(__file__), OLD_INPUT, old_selection_path, metadata_path, inventory_path, old_data / "targets.tsv", old_data / "targets.mtx.gz", old_selection_data / "target_rows.npy"]
    record = {"status": "frozen MC extension selections; no new predictions", "run_id": stamp, "donors": donors, "cell_types": types, "author_types": authors, "blocks": BLOCKS, "budgets": [60, 300], "levels": ["balanced", "ratio10"], "data_directory": rel(data), "source_metadata_audit": rel(metadata_path), "source_inventory_audit": rel(inventory_path), "source_input_manifest": rel(OLD_INPUT), "source_input_sha256": sha(OLD_INPUT), "source_selection_manifest": rel(old_selection_path), "source_selection_sha256": sha(old_selection_path), "protocol": rel(PROTOCOL), "protocol_sha256": sha(PROTOCOL), "complete_eligible_cells": 114524, "selected_cells": len(selected_rows), "old_input_cells_reused": 93892, "new_cells_requiring_counts": len(selected_rows) - 93892, "selected_stored_entries": int(np.diff(np.load(source / "raw_indptr.npy", allow_pickle=False))[selected_rows].sum()), "logical_cases": 5376, "unique_reference_matrices": 1344, "targets": 840, "scientific_prediction_rows": 645120, "weighted_fits": 322560, "seed_algorithm": "SHA256 first16 bytes unsigned big endian -> PCG64; original key family, block 3..14", "old_three_prefixes_reconstructed_exactly": True, "target_rows_byte_identical": sha(data / "target_rows.npy") == sha(old_selection_data / "target_rows.npy"), "sources": [{"path": rel(p), "sha256": sha(p)} for p in source_paths], "artifacts": [{"path": rel(p), "sha256": sha(p), "bytes": p.stat().st_size} for folder in [out, data] for p in sorted(folder.iterdir()) if p.is_file()]}
    save(out / "selection.json", record)
    print("SELECTION", rel(out / "selection.json"), flush=True)
    print(json.dumps({k: record[k] for k in ["selected_cells", "new_cells_requiring_counts", "selected_stored_entries", "logical_cases", "unique_reference_matrices", "scientific_prediction_rows"]}), flush=True)


if __name__ == "__main__":
    main()
