from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk_scale"
assert not (OUT / "rpl_control_summary.csv").exists(), "Reuse saved control results."
pred = pd.read_csv(OUT / "rpl_control_predictions.csv")
diag = pd.read_csv(OUT / "rpl_control_diagnostics.csv")
assert len(pred) == 24 and len(diag) == 12
assert not pred.duplicated(["sample_id", "method"]).any()
truth = pd.read_csv(ROOT / "data/processed/real_bulk/flow_truth.csv").set_index("sample_id")
pred["B"] = pred.pred_0
pred["T"] = pred.pred_2 + pred.pred_3
pred["NK"] = pred.pred_5
pred["Monocytes"] = pred.pred_1 + pred.pred_4
lineages = ["B", "T", "NK", "Monocytes"]
pred["mae_pp"] = np.abs(pred[lineages].to_numpy() - truth.loc[pred.sample_id, lineages].to_numpy()).mean(axis=1) * 100
old = pd.read_csv(ROOT / "results/real_bulk/prediction_errors.csv")
old = old.loc[old.reference_id.eq(pred.reference_id.iloc[0])].rename(columns={"mae_pp": "direct_mae_pp"})
scored = pred.merge(old, on=["reference_id", "sample_id", "method"], validate="one_to_one")
scored["change_pp"] = scored.mae_pp - scored.direct_mae_pp
scored.to_csv(OUT / "rpl_control_scored.csv", index=False)
rows = []
for method, frame in scored.groupby("method"):
    rows.append({"method": method, "direct_mae_pp": frame.direct_mae_pp.mean(), "rpl_excluded_mae_pp": frame.mae_pp.mean(),
                 "change_pp": frame.change_pp.mean(), "participants_improved": int(frame.change_pp.lt(-1e-12).sum()),
                 "B_zero_fits": int(frame["B"].le(1e-10).sum()), "T_zero_fits": int(frame["T"].le(1e-10).sum()), "participants": len(frame)})
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "rpl_control_summary.csv", index=False)
quality = {"all_converged": bool(diag.convergence.str.match(r"^Converge at \d+$").all()),
           "variance_finite": bool(diag.variance_finite.all()), "features_min": int(diag.features.min()), "features_max": int(diag.features.max())}
(OUT / "rpl_control_quality.json").write_text(json.dumps(quality, indent=2) + "\n")
print(summary.to_string(index=False))
print(json.dumps(quality))
