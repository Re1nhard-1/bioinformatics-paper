from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk"
DATA = ROOT / "data/processed/real_bulk"
assert not (OUT / "summary.csv").exists(), "Reuse the saved summaries."
files = sorted(OUT.glob("predictions_*.csv"))
assert len(files) == 14
pred = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
diag = pd.concat([pd.read_csv(p) for p in sorted(OUT.glob("diagnostics_*.csv"))], ignore_index=True)
truth = pd.read_csv(DATA / "flow_truth.csv")
sources = json.loads((DATA / "sources.json").read_text())
ref_input = json.loads(Path(sources["reference_input_manifest"]).read_text())
refs = pd.read_csv(ROOT / ref_input["reference_file"], sep="\t")
assert len(pred) == 32256 and len(diag) == 16128
assert pred.groupby("reference_id").size().eq(24).all()
assert diag.groupby("reference_id").size().eq(12).all()
assert not pred.duplicated(["reference_id", "sample_id", "method"]).any()
assert not diag.duplicated(["reference_id", "sample_id"]).any()
assert set(pred.reference_id) == set(refs.reference_id)
assert set(pred.sample_id) == set(truth.sample_id)
six = pred[[f"pred_{i}" for i in range(6)]].to_numpy()
assert np.isfinite(six).all() and (six >= -1e-10).all()
assert np.allclose(six.sum(axis=1), 1, rtol=0, atol=1e-10)
lineages = ["B", "T", "NK", "Monocytes"]
pred["B"] = pred.pred_0
pred["T"] = pred.pred_2 + pred.pred_3
pred["NK"] = pred.pred_5
pred["Monocytes"] = pred.pred_1 + pred.pred_4
scored = pred.merge(truth[["sample_id", *lineages]], on="sample_id", suffixes=("_pred", "_truth"), validate="many_to_one")
errors = scored[[f"{k}_pred" for k in lineages]].to_numpy() - scored[[f"{k}_truth" for k in lineages]].to_numpy()
scored["mae_pp"] = np.abs(errors).mean(axis=1) * 100
keys = ["method", "sample_id", "block", "triple_key", "budget"]
balanced = scored.loc[scored.level.eq("balanced"), [*keys, "mae_pp"]].rename(columns={"mae_pp": "balanced_mae_pp"})
contrast = scored.loc[scored.level.eq("ratio10"), [*keys, "dominant_donor", "mae_pp"]].merge(balanced, on=keys, validate="many_to_one")
contrast = contrast.rename(columns={"mae_pp": "ratio10_mae_pp"})
contrast["I_pp"] = contrast.ratio10_mae_pp - contrast.balanced_mae_pp
assert len(contrast) == 24192
assert contrast.groupby(["method", "sample_id", "block", "triple_key", "budget"]).size().eq(3).all()
metrics = ["balanced_mae_pp", "ratio10_mae_pp", "I_pp"]
donors = contrast.groupby(["method", "sample_id", "budget"], as_index=False)[metrics].mean()
blocks = contrast.groupby(["method", "block", "budget"], as_index=False)[metrics].mean()
arrangements = contrast.groupby(["method", "triple_key", "dominant_donor", "budget"], as_index=False)[metrics].mean()
summaries = []
endpoints = []
for method in ["weighted", "nnls"]:
    dm = donors.loc[donors.method.eq(method)]
    bm = blocks.loc[blocks.method.eq(method)]
    for budget in [60, 300]:
        a = arrangements.loc[arrangements.method.eq(method) & arrangements.budget.eq(budget), "I_pp"]
        b = bm.loc[bm.budget.eq(budget), "I_pp"]
        average = dm.loc[dm.budget.eq(budget), metrics].mean()
        summaries.append({"method": method, "budget": budget, **average.to_dict(),
            "I_mcse_pp": b.std(ddof=1) / np.sqrt(12), "arrangement_mean_abs_I_pp": a.abs().mean(),
            "arrangement_min_I_pp": a.min(), "arrangement_max_I_pp": a.max(),
            "arrangement_positive_count": int(a.gt(0).sum()), "arrangement_count": len(a)})
    de = dm.pivot(index="sample_id", columns="budget", values="I_pp")
    be = bm.pivot(index="block", columns="budget", values="I_pp")
    delta = de[300] - de[60]
    bdelta = be[300] - be[60]
    endpoints.append({"method": method, "endpoint_change_pp": delta.mean(),
        "endpoint_mcse_pp": bdelta.std(ddof=1) / np.sqrt(12),
        "donor_min_endpoint_pp": delta.min(), "donor_max_endpoint_pp": delta.max(),
        "donors_negative_endpoint": int(delta.lt(-1e-12).sum()), "donors_positive_endpoint": int(delta.gt(1e-12).sum()),
        "donors_numerically_zero_endpoint": int(delta.abs().le(1e-12).sum()),
        "donor_count": len(delta)})
summary = pd.DataFrame(summaries)
endpoint = pd.DataFrame(endpoints)
agreement = scored.loc[scored.level.eq("balanced")].groupby(["method", "sample_id", "budget"], as_index=False)[
    [f"{k}_pred" for k in lineages] + [f"{k}_truth" for k in lineages]
].mean()
type_errors = []
for group, frame in scored.groupby(["method", "budget", "level"]):
    for lineage in lineages:
        error = (frame[f"{lineage}_pred"] - frame[f"{lineage}_truth"]) * 100
        type_errors.append(dict(zip(["method", "budget", "level"], group)) | {
            "lineage": lineage, "mean_error_pp": error.mean(), "mae_pp": error.abs().mean()})
scored[["reference_id", "sample_id", "method", "mae_pp"]].to_csv(OUT / "prediction_errors.csv", index=False)
contrast.to_csv(OUT / "paired_contrasts.csv", index=False)
donors.to_csv(OUT / "donor_summary.csv", index=False)
blocks.to_csv(OUT / "block_summary.csv", index=False)
arrangements.to_csv(OUT / "arrangement_summary.csv", index=False)
endpoint.to_csv(OUT / "endpoint_summary.csv", index=False)
agreement.to_csv(OUT / "balanced_agreement.csv", index=False)
pd.DataFrame(type_errors).to_csv(OUT / "lineage_errors.csv", index=False)
quality = {"reference_configurations": int(pred.reference_id.nunique()), "bulk_samples": int(pred.sample_id.nunique()),
    "weighted_fits": len(diag), "prediction_rows_both_outputs": len(pred),
    "convergence_counts": {str(k): int(v) for k, v in diag.convergence.value_counts(dropna=False).items()},
    "nonfinite_variance_count": int((~diag.variance_finite).sum()),
    "effective_gene_range": [int(diag.n_features.min()), int(diag.n_features.max())],
    "represented_PBMC_percent_range": sources["represented_PBMC_percent_range"],
    "protocol_sha256": __import__("hashlib").sha256((ROOT / "analysis/REAL_BULK_PROTOCOL.md").read_bytes()).hexdigest()}
(OUT / "fit_quality.json").write_text(json.dumps(quality, indent=2) + "\n", encoding="utf-8")
summary.to_csv(OUT / "summary.csv", index=False)
print(summary.to_string(index=False))
print(endpoint.to_string(index=False))
print(json.dumps(quality, indent=2))
