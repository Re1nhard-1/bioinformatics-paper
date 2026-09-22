from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk_scale"
assert not (OUT / "summary.csv").exists(), "Reuse saved summaries."
refs = pd.read_csv(ROOT / "data/processed/real_bulk_scale/references.tsv", sep="\t")
truth = pd.read_csv(ROOT / "data/processed/real_bulk/flow_truth.csv")
new = pd.concat([pd.read_csv(p) for p in sorted(OUT.glob("predictions_*.csv"))], ignore_index=True)
diag = pd.concat([pd.read_csv(p) for p in sorted(OUT.glob("diagnostics_*.csv"))], ignore_index=True)
assert len(new) == 672 and len(diag) == 336
assert new.groupby(["arm", "reference_id"]).size().eq(24).all()
assert diag.groupby(["arm", "reference_id"]).size().eq(12).all()
assert not new.duplicated(["arm", "reference_id", "sample_id", "method"]).any()
assert set(new.reference_id) == set(refs.reference_id)
raw = pd.concat([pd.read_csv(p) for p in sorted((ROOT / "results/real_bulk").glob("predictions_*.csv"))], ignore_index=True)
raw = raw.loc[raw.reference_id.isin(refs.reference_id)].copy()
assert len(raw) == 336
raw["arm"] = "direct_counts"
pred = pd.concat([raw, new], ignore_index=True)
lineages = ["B", "T", "NK", "Monocytes"]
pred["B"] = pred.pred_0
pred["T"] = pred.pred_2 + pred.pred_3
pred["NK"] = pred.pred_5
pred["Monocytes"] = pred.pred_1 + pred.pred_4
scored = pred.merge(truth[["sample_id", *lineages]], on="sample_id", suffixes=("_pred", "_truth"), validate="many_to_one")
errors = scored[[f"{k}_pred" for k in lineages]].to_numpy() - scored[[f"{k}_truth" for k in lineages]].to_numpy()
scored["mae_pp"] = np.nan
is_new = scored.arm.ne("direct_counts")
scored.loc[is_new, "mae_pp"] = np.abs(errors[is_new]).mean(axis=1) * 100
old_errors = pd.read_csv(ROOT / "results/real_bulk/prediction_errors.csv").set_index(["reference_id", "sample_id", "method"]).mae_pp
raw_keys = pd.MultiIndex.from_frame(scored.loc[~is_new, ["reference_id", "sample_id", "method"]])
scored.loc[~is_new, "mae_pp"] = old_errors.loc[raw_keys].to_numpy()
assert scored.mae_pp.notna().all()
donors = scored.groupby(["arm", "method", "sample_id"], as_index=False).mae_pp.mean()
baseline = donors.loc[donors.arm.eq("direct_counts"), ["method", "sample_id", "mae_pp"]].rename(columns={"mae_pp": "direct_mae_pp"})
donors = donors.merge(baseline, on=["method", "sample_id"], validate="many_to_one")
donors["change_from_direct_pp"] = donors.mae_pp - donors.direct_mae_pp
rows = []
for (arm, method), frame in donors.groupby(["arm", "method"]):
    fits = scored.loc[scored.arm.eq(arm) & scored.method.eq(method)]
    rows.append({"arm": arm, "method": method, "mae_pp": frame.mae_pp.mean(), "donor_min_mae_pp": frame.mae_pp.min(),
                 "donor_max_mae_pp": frame.mae_pp.max(), "change_from_direct_pp": frame.change_from_direct_pp.mean(),
                 "donors_improved": int(frame.change_from_direct_pp.lt(-1e-12).sum()), "donors": len(frame),
                 "B_zero_fits": int(fits.B_pred.le(1e-10).sum()), "T_zero_fits": int(fits.T_pred.le(1e-10).sum()), "fits": len(fits)})
summary = pd.DataFrame(rows)
lineage_rows = []
for (arm, method), frame in scored.groupby(["arm", "method"]):
    for lineage in lineages:
        difference = (frame[f"{lineage}_pred"] - frame[f"{lineage}_truth"]) * 100
        lineage_rows.append({"arm": arm, "method": method, "lineage": lineage, "mean_error_pp": difference.mean(), "mae_pp": difference.abs().mean()})
agreement = scored.groupby(["arm", "method", "sample_id"], as_index=False)[[f"{k}_{kind}" for k in lineages for kind in ["pred", "truth"]]].mean()
quality = {}
for arm, frame in diag.groupby("arm"):
    iterations = pd.to_numeric(frame.convergence.str.extract(r"^Converge at (\d+)$")[0], errors="coerce")
    quality[arm] = {"fits": len(frame), "all_converged": bool(iterations.notna().all()), "convergence_min": float(iterations.min()),
                    "convergence_max": float(iterations.max()), "variance_finite": bool(frame.variance_finite.all()),
                    "features_min": int(frame.n_features.min()), "features_max": int(frame.n_features.max())}
scored[["arm", "reference_id", "sample_id", "method", "mae_pp"]].to_csv(OUT / "prediction_errors.csv", index=False)
donors.to_csv(OUT / "donor_summary.csv", index=False)
summary.to_csv(OUT / "summary.csv", index=False)
pd.DataFrame(lineage_rows).to_csv(OUT / "lineage_errors.csv", index=False)
agreement.to_csv(OUT / "agreement.csv", index=False)
(OUT / "fit_quality.json").write_text(json.dumps(quality, indent=2) + "\n")
print(summary.to_string(index=False))
print(json.dumps(quality))
