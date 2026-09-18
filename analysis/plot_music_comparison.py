"""Plot absolute error, imbalance changes, and their interaction from a fixed review.

This is presentation-only code: no new estimators, confidence intervals, tests,
parameter choices, exclusions, or reanalysis of individual predictions.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scientific_style import COLORS, configure_style, figure_mm, save_figure_bundle, sha256

ROOT = Path(__file__).resolve().parents[1]
DONORS = ["101", "1015", "1016", "1244", "1256", "1488"]
METHODS = [
    ("music_nnls", "Official MuSiC: ordinary NNLS", COLORS["blue"], "o", -.12),
    ("music_weighted", "Official MuSiC: weighted fit", COLORS["vermillion"], "s", .12),
]


def check_recorded_hash(review: dict, path: Path):
    matches = [item for item in review.get("artifacts", [])
               if Path(item.get("path", item.get("file", ""))).name == path.name]
    if len(matches) != 1 or sha256(path) != matches[0]["sha256"]:
        raise ValueError(f"Missing or mismatched frozen input hash: {path.name}")


def fixed_ticks(ax, low: float, high: float):
    """Choose readable ticks without placing hidden off-range tick text."""
    ax.set_xlim(low, high)
    ticks = MaxNLocator(nbins=3).tick_values(low, high)
    ax.set_xticks(ticks[(ticks >= low - 1e-12) & (ticks <= high + 1e-12)])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}".replace("-", "−")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review", type=Path, help="Completed results/music_review/<UTC> directory")
    args = parser.parse_args()
    source = args.review.resolve()
    review_path = source / "review.json"
    review = json.loads(review_path.read_bytes())
    if review.get("rows_verified") != 86400 or not review.get("complete_paired_grid", False):
        raise ValueError("Plotting requires the completed, verified 86,400-row paired review.")
    metric_path = source / "donor_metrics.csv"
    did_path = source / "donor_difference_in_differences.csv"
    for path in (metric_path, did_path):
        check_recorded_hash(review, path)
    data = pd.read_csv(metric_path, dtype={"held_out": str})
    did = pd.read_csv(did_path, dtype={"held_out": str}).set_index("held_out")
    if (len(data) != 12 or set(data.held_out) != set(DONORS)
            or set(data.method) != {method[0] for method in METHODS}
            or data.duplicated(["held_out", "method"]).any()
            or len(did) != 6 or set(did.index) != set(DONORS) or did.index.duplicated().any()):
        raise ValueError("Expected exactly six donors and two official MuSiC outputs.")
    numeric = data.select_dtypes(include="number")
    if not np.isfinite(numeric).all().all() or not np.isfinite(did.select_dtypes(include="number")).all().all():
        raise ValueError("Non-finite source data.")
    for donor in DONORS:
        rows = data.query("held_out == @donor").set_index("method")
        actual = rows.mae_strong_pp - rows.mae_balanced_pp
        if not np.allclose(actual, rows.mae_imbalance_change_pp, atol=1e-10, rtol=0):
            raise ValueError("Inconsistent strong-minus-balanced differences.")
        for method, *_ in METHODS:
            if abs(actual.loc[method] - did.loc[donor, f"mae_change_{method}_pp"]) > 1e-10:
                raise ValueError("Donor effect source files disagree.")
        expected = actual.loc["music_weighted"] - actual.loc["music_nnls"]
        if abs(expected - did.loc[donor, "difference_in_differences_pp"]) > 1e-10:
            raise ValueError("Difference-in-differences mismatch.")
        for allocation in ("balanced", "strong"):
            expected = rows.loc["music_weighted", f"mae_{allocation}_pp"] - rows.loc["music_nnls", f"mae_{allocation}_pp"]
            if abs(expected - did.loc[donor, f"weighted_minus_nnls_{allocation}_pp"]) > 1e-10:
                raise ValueError("Absolute method comparison mismatch.")
    style = configure_style()
    fig = figure_mm(183, 113)
    axes = fig.subplots(1, 3, sharey=True, gridspec_kw={"width_ratios": [1.15, 1, 1]})
    fig.subplots_adjust(left=.10, right=.98, bottom=.225, top=.735, wspace=.35)
    records, handles = [], []
    for method, label, color, marker, offset in METHODS:
        part = data.query("method == @method").set_index("held_out").loc[DONORS]
        ys = np.arange(6) + offset
        for i, donor in enumerate(DONORS):
            balanced = float(part.loc[donor, "mae_balanced_pp"])
            strong = float(part.loc[donor, "mae_strong_pp"])
            change = float(part.loc[donor, "mae_imbalance_change_pp"])
            axes[0].plot([balanced, strong], [ys[i], ys[i]], color=color, lw=.8, zorder=1)
            for allocation, value in (("balanced", balanced), ("strong", strong)):
                face = "white" if allocation == "balanced" else color
                axes[0].scatter(value, ys[i], marker=marker, s=19, facecolor=face, edgecolor=color, linewidth=.8, zorder=3)
                records.append({"panel": "a", "held_out": donor, "method": method,
                                "quantity": "absolute_mae", "allocation": allocation,
                                "value_pp": value, "unit": "percentage points"})
            axes[1].scatter(change, ys[i], marker=marker, s=19, facecolor=color, edgecolor=color, linewidth=.8, zorder=3)
            records.append({"panel": "b", "held_out": donor, "method": method,
                            "quantity": "strong_minus_balanced_mae", "allocation": "paired_difference",
                            "value_pp": change, "unit": "percentage points"})
        handles.append(Line2D([], [], color=color, marker=marker, markerfacecolor=color,
                              markeredgewidth=.8, linestyle="none", markersize=4, label=label))
    for i, donor in enumerate(DONORS):
        value = float(did.loc[donor, "difference_in_differences_pp"])
        axes[2].scatter(value, i, marker="D", s=17, color=COLORS["black"], zorder=3)
        records.append({"panel": "c", "held_out": donor, "method": "weighted_minus_nnls",
                        "quantity": "difference_in_differences", "allocation": "paired_interaction",
                        "value_pp": value, "unit": "percentage points"})
    for allocation, filled, label in [("balanced", False, "Balanced: 20:20:20"),
                                      ("strong", True, "Strong imbalance: 50:5:5")]:
        handles.append(Line2D([], [], color="black", marker="o", markerfacecolor="black" if filled else "white",
                              markeredgewidth=.8, linestyle="none", markersize=4, label=label))
    titles = ["Absolute error", "Effect of imbalance", "Difference in differences"]
    for index, (ax, title, letter) in enumerate(zip(axes, titles, ["a", "b", "c"])):
        ax.set_yticks(np.arange(6), DONORS)
        ax.set_ylim(5.45, -.5)
        ax.tick_params(axis="y", length=0, pad=5, labelleft=(index == 0))
        ax.spines["left"].set_visible(False)
        ax.set_title(title, loc="left", pad=9)
        ax.text(-.085, 1.04, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        if index > 0:
            ax.axvline(0, color=COLORS["neutral"], linestyle=(0, (3, 2)), lw=.65, zorder=0)
    axes[0].set_ylabel("Held-out donor")
    max_mae = data[["mae_balanced_pp", "mae_strong_pp"]].to_numpy().max()
    fixed_ticks(axes[0], 0, np.ceil(max_mae * 1.13))
    for ax, values in [(axes[1], data.mae_imbalance_change_pp),
                       (axes[2], did.difference_in_differences_pp)]:
        bound = max(.1, np.ceil(np.abs(values).max() * 1.15 * 10) / 10)
        fixed_ticks(ax, -bound, bound)
    axes[0].set_xlabel("Mean absolute error\n(percentage points)", labelpad=6)
    axes[1].set_xlabel("Strong − balanced MAE\n(percentage points)", labelpad=6)
    axes[2].set_xlabel("Weighted − NNLS change\n(percentage points)", labelpad=6)
    fig.text(.10, .958, "Official MuSiC workflow under reference donor imbalance", fontsize=7, fontweight="bold")
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.10, .91), ncol=2,
               columnspacing=2.4, handletextpad=.6)
    fig.text(.10, .075, "Positive values in c indicate a larger imbalance penalty; compare absolute error in a.", fontsize=6.5)
    fig.text(.10, .035, "Six held-out donors; raw-count synthetic mixtures; three reference draws; descriptive results.", fontsize=6.5)
    caption = f"""# Figure. Absolute error and response to reference donor imbalance in the official MuSiC workflow

**a**, Mean absolute error (MAE) for balanced (20:20:20; open symbols) and strongly
imbalanced (50:5:5; filled symbols) reference cell allocations, with 60 reference
cells per type from three donors. Each horizontal segment connects the two
allocations for the same held-out donor and method only; it is not an error bar.
**b**, Strong-minus-balanced MAE for each method. **c**, The difference between
those changes: weighted MuSiC minus the official ordinary NNLS output. The
vertical dashed lines in b and c denote zero. All quantities are in percentage
points. A positive value in c means a larger imbalance penalty for the weighted
estimator; it does not mean higher absolute error overall. A negative value
would mean a more favourable response to imbalance for the weighted estimator,
but would not imply lower absolute MAE. Panel a must be read alongside b and c.
Donor order is fixed by numeric ID, not by observed effect.

Blue circles represent `Est.prop.allgene` (ordinary NNLS) and vermilion squares
represent `Est.prop.weighted`, both returned by the same fixed official
`music_prop` implementation. Both share its donor-equal reference basis and
cell-size construction. They are not the previous cell-pooled versus equal-donor
NNLS comparison. Black diamonds in c denote the within-donor interaction.

The source cohort comprises six eligible GSE96583 PBMC donors with SLE under
the control-culture condition. All genes enter as raw counts, with the official
MuSiC common-support and nonzero-gene filtering (`markers=NULL`). Targets are
raw-count sums of 300 cells, using the previously frozen cell identities and
six-type compositions; these are synthetic mixtures, not measured bulk samples.
Each donor's MAE averages the 60 fixed targets, three reference draws, all ten
three-donor training subsets, six cell types, and (for strong imbalance) all
three dominant-donor choices. Each strong setting shares its paired balanced
baseline; duplicated use does not create additional independent observations.

Six source donors define the biological scope, while training references overlap
between held-out folds. No population confidence intervals, hypothesis tests,
significance stars or inferential error bars are plotted. Three reference draws
provide limited Monte Carlo precision. Conditional sampling SD is reported in
the adjacent review tables, not shown here and not treated as a population CI.
The raw-count mixture model differs from the earlier normalized-cell v3.1
mixtures; cross-run absolute MAEs cannot isolate method effects.

Source: `{source.relative_to(ROOT).as_posix()}/donor_metrics.csv` and
`{source.relative_to(ROOT).as_posix()}/donor_difference_in_differences.csv`.
The adjacent source-data file records all 42 plotted values. Source hashes and
paired-difference arithmetic are checked before plotting. No extra analyses or
post-result subset selections are performed by this script.
"""
    output = ROOT / "results/music_figures" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output.mkdir(parents=True)
    manifest = save_figure_bundle(fig, output / "music_absolute_and_imbalance_effects", pd.DataFrame(records),
                                 caption, [metric_path, did_path, review_path], Path(__file__), {
        "title": "Official MuSiC absolute error and donor-imbalance effects",
        "font": style, "biological_donors": 6, "selected_donor_ids": DONORS,
        "methods": [item[0] for item in METHODS], "reference_draws": 3,
        "statistical_claim": "Descriptive paired effects; folds share training donors; no inferential intervals/tests.",
        "plotted_values": 42,
        "contrast_definition": "(weighted strong − weighted balanced) − (NNLS strong − NNLS balanced)",
        "absolute_error_caveat": "A negative difference in differences does not imply lower absolute MAE.",
    })
    plt.close(fig)
    print(json.dumps({"output": str(output), "checks": manifest["checks"], "plotted_values": len(records)}, indent=2))


if __name__ == "__main__":
    main()
