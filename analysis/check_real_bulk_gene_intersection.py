from pathlib import Path
import gzip
import json

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk_scale/gene_intersection_diagnostic.json"
assert not OUT.exists(), "Reuse saved diagnosis."
manifest = json.loads((ROOT / "results/music_mc_extension_input/20260918T154429497808Z/input_manifest.json").read_text())
data = ROOT / manifest["data_directory"]
genes = (data / "genes.tsv").read_text().splitlines()
bulk = pd.read_csv(ROOT / "data/processed/real_bulk/bulk_counts.tsv.gz", sep="\t", usecols=["gene_id"])
bulk_genes = set(bulk.gene_id)
missing = [g for g in genes if g not in bulk_genes]
gene_rows = [genes.index(g) for g in missing]
meta = pd.read_csv(ROOT / "data/processed/cellxgene_metadata/20260917T111912321314Z/raw_genes.tsv", sep="\t").set_index("gene_id")
cells = pd.read_csv(data / "cells.tsv", sep="\t")
refs = pd.read_csv(ROOT / manifest["reference_file"], sep="\t")
ref_id = "d00_d01_d02_b3_n300_balanced"
ref = refs.loc[refs.reference_id.eq(ref_id)].iloc[0]
selected = cells.iloc[np.array(ref.reference_columns_R.split("|"), dtype=int) - 1].copy()
selected["omitted_umi"] = 0.0
for donor in selected.donor.unique():
    entry = next(x for x in manifest["donor_matrices"] if x["donor"] == donor)
    donor_cells = cells.loc[cells.donor.eq(donor)].reset_index(drop=True)
    chosen = selected.loc[selected.donor.eq(donor), "cell_id"]
    positions = pd.Index(donor_cells.cell_id).get_indexer(chosen)
    assert (positions >= 0).all()
    with gzip.open(ROOT / entry["path"], "rb") as f:
        matrix = mmread(f).tocsr()
    omitted = np.asarray(matrix[gene_rows, :][:, positions].sum(axis=0)).ravel()
    selected.loc[selected.donor.eq(donor), "omitted_umi"] = omitted
    del matrix
selected["omitted_fraction"] = selected.omitted_umi / selected.total_counts
grouped = selected.groupby("cell_type")[["omitted_umi", "total_counts"]].sum()
grouped["omitted_fraction"] = grouped.omitted_umi / grouped.total_counts
result = {
    "reference_id": ref_id,
    "reference_genes": len(genes),
    "common_genes": len(genes) - len(missing),
    "missing_genes": [{"gene_id": g, "symbol": meta.loc[g, "feature_name"]} for g in missing],
    "cells": len(selected),
    "pooled_omitted_fraction": float(selected.omitted_umi.sum() / selected.total_counts.sum()),
    "maximum_cell_omitted_fraction": float(selected.omitted_fraction.max()),
    "by_type": grouped.reset_index().to_dict("records"),
    "scope": "Input-only diagnosis of the effect of intersecting the single-cell matrix before computing library sizes; no deconvolution refit; not a test of Salmon counts or a reconstruction of Zhao's unpublished gene set."
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
