from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk_scale"
assert not (OUT / "residual_summary.json").exists(), "Reuse saved diagnosis."
gtf = pd.read_csv(ROOT / "data/raw/real_bulk/human.gene_sums.G026.gtf.gz", sep="\t", comment="#", header=None)
gtf["gene_id"] = gtf[8].str.extract('gene_id "([^"]+)"')[0].str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
gtf["symbol"] = gtf[8].str.extract('gene_name "([^"]+)"')[0]
meta = gtf[["gene_id", "symbol"]]
profile = pd.read_csv(OUT / "diagnostic_reference_profiles.csv").merge(meta, on="gene_id", validate="one_to_one")
rp = profile.symbol.str.match(r"^RP[LS]")
excluded = set(profile.loc[rp, "gene_id"])
(OUT / "excluded_rpl_rps_genes.txt").write_text("\n".join(sorted(excluded)) + "\n")
residual = pd.read_csv(OUT / "residuals.csv.gz").merge(meta, on="gene_id", validate="many_to_one")
result = {"reference_id": "d00_d01_d02_b3_n300_balanced", "selection_rule": "GENCODE v26 gene_name matches ^RP[LS]", "excluded_genes": len(excluded)}
result["reference_rpl_rps_fraction"] = {str(k): float(profile.loc[rp, f"theta_{k}"].sum() / profile[f"theta_{k}"].sum()) for k in range(6)}
fractions = []
for arm, filename in [("direct_counts", "real_bulk/bulk_counts.tsv.gz"), ("length_scaled", "real_bulk_scale/length_scaled.tsv.gz"), ("author_tpm", "real_bulk_scale/author_tpm.tsv.gz")]:
    bulk = pd.read_csv(ROOT / "data/processed" / filename, sep="\t", index_col=0).reindex(profile.gene_id).fillna(0)
    mass = bulk.loc[bulk.index.isin(excluded)].sum() / bulk.sum()
    fractions.extend({"arm": arm, "sample_id": sample, "rpl_rps_fraction": float(value)} for sample, value in mass.items())
result["bulk_rpl_rps_fraction"] = fractions
residual["squared_residual"] = residual.residual.pow(2)
result["largest_unweighted_residuals"] = {}
for method, frame in residual.groupby("method"):
    top = frame.groupby(["gene_id", "symbol"])[["observed", "fitted", "squared_residual"]].mean().sort_values("squared_residual", ascending=False).head(20).reset_index()
    result["largest_unweighted_residuals"][method] = top.to_dict("records")
nnls = residual.loc[residual.method.eq("nnls")].merge(profile[["gene_id", *[f"D_{k}" for k in range(6)]]], on="gene_id", validate="many_to_one")
nnls["group"] = np.where(nnls.gene_id.isin(excluded), "RPL_RPS_prefix", "other")
gradient = []
for k in [0, 2, 3]:
    nnls["D_times_residual"] = nnls[f"D_{k}"] * nnls.residual
    grouped = nnls.groupby(["sample_id", "group"]).D_times_residual.sum().unstack()
    for sample, row in grouped.iterrows():
        gradient.append({"sample_id": sample, "type_index": k, "RPL_RPS_prefix": float(row.RPL_RPS_prefix), "other": float(row.other), "total": float(row.sum())})
result["nnls_D_times_residual"] = gradient
result["gradient_interpretation"] = "At a zero NNLS coefficient, positive D-transpose residual would favor entry; negative discourages entry. This is the ordinary objective, not MuSiC's weighted objective."
(OUT / "residual_summary.json").write_text(json.dumps(result, indent=2) + "\n")
print("Reference RPL/RPS fractions:", result["reference_rpl_rps_fraction"])
print(pd.DataFrame(fractions).groupby("arm").rpl_rps_fraction.agg(["min", "max"]).to_string())
print(pd.DataFrame(gradient).groupby("type_index")[["RPL_RPS_prefix", "other", "total"]].mean().to_string())
