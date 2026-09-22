from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scientific_style import COLORS, configure_style, figure_mm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/real_bulk_salmon"
if "--plot-only" in sys.argv:
    pred = pd.read_csv(OUT / "comparison.csv")
    truth = pd.read_csv(OUT / "flow_truth.csv").set_index("lineage").flow_truth
    lineages = ["B", "T", "NK", "Monocytes"]
    columns = pred.columns
    decision = json.loads((OUT / "pilot_decision.json").read_text())
else:
    assert not (OUT / "comparison.csv").exists(), "Reuse the saved pilot summary."
    new = pd.read_csv(OUT / "predictions.csv")
    reference = "d00_d01_d02_b3_n300_balanced"
    old = pd.read_csv(ROOT / "results/real_bulk/predictions_d00_d01_d02.csv")
    old = old.loc[old.reference_id.eq(reference) & old.sample_id.eq("453W")].copy()
    assert len(new) == 4 and len(old) == 2
    old["arm"] = "direct_counts"
    pred = pd.concat([old, new], ignore_index=True)
    pred["B"] = pred.pred_0
    pred["T"] = pred.pred_2 + pred.pred_3
    pred["NK"] = pred.pred_5
    pred["Monocytes"] = pred.pred_1 + pred.pred_4
    lineages = ["B", "T", "NK", "Monocytes"]
    truth = pd.read_csv(ROOT / "data/processed/real_bulk/flow_truth.csv").set_index("sample_id").loc["453W", lineages].astype(float)
    old_errors = pd.read_csv(ROOT / "results/real_bulk/prediction_errors.csv")
    old_errors = old_errors.loc[old_errors.reference_id.eq(reference) & old_errors.sample_id.eq("453W")].set_index("method")
    pred["mae_pp"] = np.nan
    for index, row in pred.iterrows():
        pred.loc[index, "mae_pp"] = (old_errors.loc[row.method, "mae_pp"] if row.arm == "direct_counts"
                                     else np.abs(row[lineages].astype(float).to_numpy() - truth.to_numpy()).mean() * 100)
    baseline = pred.loc[pred.arm.eq("direct_counts")].set_index("method").mae_pp
    pred["change_from_direct_pp"] = pred.mae_pp - pred.method.map(baseline)
    pred["B_zero"] = pred.B.le(1e-10)
    pred["T_zero"] = pred["T"].le(1e-10)
    columns = ["arm", "method", *lineages, "mae_pp", "change_from_direct_pp", "B_zero", "T_zero"]
    pred[columns].to_csv(OUT / "comparison.csv", index=False)
    pd.DataFrame({"lineage": lineages, "flow_truth": truth.to_numpy()}).to_csv(OUT / "flow_truth.csv", index=False)
    qc = json.loads((OUT / "input_quality.json").read_text())
    diag = pd.read_csv(OUT / "diagnostics.csv")
    primary = pred.loc[pred.arm.eq("salmon_saved_basis") & pred.method.eq("weighted")].iloc[0]
    primary_diag = diag.set_index("arm").loc["salmon_saved_basis"]
    gates = {"read_retention_ge_70pct": qc["read_retention"] >= .70,
             "mapping_ge_70pct": qc["mapping_percent"] >= 70,
             "weighted_converged": str(primary_diag.convergence).startswith("Converge at "),
             "variance_finite": bool(primary_diag.variance_finite),
             "B_recovered": bool(primary.B > .0001), "T_recovered": bool(primary["T"] > .0001),
             "weighted_mae_improved": bool(primary.change_from_direct_pp < 0)}
    decision = {"sample_id": "453W", "reference_id": reference, "gates": gates,
                "expand_to_remaining_participants": all(gates.values()),
                "interpretation_limit": "One selected post-hoc diagnostic sample and one fixed reference; not cohort validation or exact reproduction of Zhao's pipeline."}
    (OUT / "pilot_decision.json").write_text(json.dumps(decision, indent=2) + "\n")
configure_style()
fig = figure_mm(183, 88)
axes = fig.subplots(1, 2, sharey=True)
arms = ["direct_counts", "salmon_saved_basis", "salmon_shared_prefilter"]
labels = ["Saved direct counts", "Salmon: saved basis", "Salmon: shared genes first"]
colors = [COLORS["neutral"], COLORS["blue"], COLORS["orange"]]
for ax, method, letter in zip(axes, ["weighted", "nnls"], ["a", "b"]):
    positions = np.arange(4)
    rows = pred.loc[pred.method.eq(method)].set_index("arm")
    for shift, arm, label, color, hatch in zip([-.24, 0, .24], arms, labels, colors, ["", "//", ".."]):
        ax.bar(positions + shift, rows.loc[arm, lineages].astype(float).to_numpy() * 100,
               width=.22, color=color, hatch=hatch, edgecolor="white", linewidth=.25, label=label, zorder=2)
        zero = rows.loc[arm, lineages].astype(float).to_numpy() == 0
        ax.scatter(positions[zero] + shift, np.zeros(zero.sum()), s=8, marker="s", color=color, zorder=3, clip_on=False)
    ax.scatter(positions, truth * 100, marker="_", s=230, color="black", linewidths=1.3, label="Flow cytometry", zorder=4)
    ax.set_xticks(positions, ["B", "T", "NK", "Monocytes"])
    ax.set_ylim(0, 100)
    ax.set_title("Weighted MuSiC" if method == "weighted" else "Ordinary NNLS", loc="left")
    ax.text(-.12, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=8)
    error_text = " / ".join(f"{rows.loc[arm, 'mae_pp']:.2f}" for arm in arms)
    ax.text(0, -.23, f"MAE (pp), left to right: {error_text}", transform=ax.transAxes, fontsize=6.5)
axes[0].set_ylabel("Fraction within four represented lineages (%)")
handles, legend_labels = axes[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="upper center", bbox_to_anchor=(.54, .98), ncol=2, fontsize=6.5)
fig.subplots_adjust(left=.09, right=.985, top=.76, bottom=.24, wspace=.2)
stem = ROOT.parent / "figures/real_bulk_salmon_pilot"
stem.parent.mkdir(exist_ok=True)
for suffix in ("png", "pdf", "svg"):
    fig.savefig(stem.with_suffix('.' + suffix), dpi=450)
plt.close(fig)
caption = ("Quantification diagnostic for Monaco sample 453W, using the same balanced three-donor reference "
           "(300 cells per cell type). Bars show cell fractions from saved direct counts, Salmon counts with the saved "
           "reference basis, and Salmon counts with shared genes selected before reference construction. Black marks "
           "show independent flow cytometry, renormalized to the four represented lineages (90.911% of PBMCs). Small squares on the baseline mark zero estimates. "
           "MAE is the mean absolute error across these four fractions in percentage points. Results describe one sample and one reference. "
           "Ordinary NNLS is MuSiC's unweighted comparator.\n")
(OUT / "figure_caption.txt").write_text(caption, encoding="utf-8")
print(pred[columns].to_string(index=False))
print(json.dumps(decision, indent=2))
