from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/real_bulk"
OUTPUT = ROOT.parent / "figures"
ORANGE = "#D55E00"
BLUE = "#0072B2"
GREY = "#B6BCC2"
COLORS = {"B": "#0072B2", "T": "#D55E00", "NK": "#009E73", "Monocytes": "#CC79A7"}
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": .6, "lines.linewidth": .8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white",
})
agreement = pd.read_csv(DATA / "balanced_agreement.csv")
donors = pd.read_csv(DATA / "donor_summary.csv")
summary = pd.read_csv(DATA / "summary.csv")
value_columns = [f"{lineage}_{kind}" for lineage in COLORS for kind in ["pred", "truth"]]
upper = max(60, math.ceil(agreement.loc[agreement.method.eq("weighted"), value_columns].to_numpy().max() * 10) * 10)
fig, axes = plt.subplots(2, 2, figsize=(170 / 25.4, 153 / 25.4))
fig.subplots_adjust(left=.11, right=.965, bottom=.11, top=.855, hspace=.72, wspace=.48)

for ax, budget, letter in zip(axes[0], [60, 300], "ab"):
    frame = agreement.loc[agreement.method.eq("weighted") & agreement.budget.eq(budget)]
    ax.plot([0, upper], [0, upper], color="#777777", linestyle=(0, (3, 3)), linewidth=.7, zorder=0)
    for lineage, color in COLORS.items():
        ax.scatter(frame[f"{lineage}_truth"] * 100, frame[f"{lineage}_pred"] * 100,
                   s=19, color=color, edgecolor="white", linewidth=.35, alpha=.9)
    ax.set(xlim=(-2, upper), ylim=(-2, upper), xticks=list(range(0, upper + 1, 20)), yticks=list(range(0, upper + 1, 20)),
           xlabel="Flow proportion (%)", ylabel="Estimated proportion (%)")
    ax.set_title(f"Balanced reference · {budget} cells per type", loc="left", pad=9, fontsize=8.5)
    ax.text(-.22, 1.07, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)
    mae = summary.loc[summary.method.eq("weighted") & summary.budget.eq(budget), "balanced_mae_pp"].iloc[0]
    ax.text(.97, .96, f"Mean MAE: {mae:.3f} pp", ha="right", va="top", transform=ax.transAxes,
            fontsize=7.5, bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})

span = donors.I_pp.max() - donors.I_pp.min()
limits = (min(0, donors.I_pp.min()) - .12 * span, max(0, donors.I_pp.max()) + .12 * span)
for ax, method, letter, title, color, marker in zip(axes[1], ["weighted", "nnls"], "cd",
        ["Weighted MuSiC", "Ordinary NNLS"], [ORANGE, BLUE], ["s", "o"]):
    frame = donors.loc[donors.method.eq(method)]
    for _, group in frame.groupby("sample_id"):
        ax.plot([0, 1], group.set_index("budget").loc[[60, 300], "I_pp"], color=GREY, linewidth=.7)
    means = summary.loc[summary.method.eq(method)].set_index("budget").loc[[60, 300], "I_pp"]
    ax.plot([0, 1], means, color=color, marker=marker, linewidth=1.6, markersize=4, zorder=3)
    ax.axhline(0, color="#444444", linestyle=(0, (3, 3)), linewidth=.7)
    ax.set(xlim=(-.15, 1.15), ylim=limits, xticks=[0, 1], xticklabels=[60, 300],
           xlabel="Reference cells per type", ylabel="MAE difference\n(percentage points)")
    ax.set_title(title, loc="left", pad=9)
    ax.text(-.22, 1.07, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)

fig.text(.11, .982, "Measured PBMC RNA-seq and independent flow cytometry", fontsize=10, va="top")
fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=color,
                           markeredgecolor="none", markersize=4, label=label) for label, color in COLORS.items()],
           loc="upper center", bbox_to_anchor=(.54, .945), frameon=False, ncol=4, handletextpad=.25,
           columnspacing=1.5, fontsize=8)
fig.legend(handles=[Line2D([0], [0], color=GREY, label="Bulk participant"),
                    Line2D([0], [0], color="#333333", linewidth=1.6, label="Mean of 12 participants")],
           loc="lower center", bbox_to_anchor=(.54, .012), frameon=False, ncol=2, fontsize=8)
OUTPUT.mkdir(exist_ok=True)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUTPUT / f"Figure_S5.{ext}", dpi=450)
plt.close(fig)
print(OUTPUT / "Figure_S5.png")
