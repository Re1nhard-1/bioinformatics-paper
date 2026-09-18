"""Summarize frozen arrangement-level MAE contrasts; no prediction fitting.

Each source value already averages 60 targets and three reference draws within
one held-out donor, reference triple and dominant-donor arrangement. Quantiles
describe these finite, dependent arrangements, not population uncertainty.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scientific_style import COLORS, configure_style, figure_mm, save_figure_bundle

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "original_grouping": ROOT / "results/music_external_review/20260917T133518483857Z",
    "alternative_grouping": ROOT / "results/music_reference_composition_review/20260917T143411436461Z",
}
METHODS = ["music_weighted", "music_nnls"]
SIGN_TOLERANCE_PP = 1e-10
KEYS = ["source", "method", "budget"]
STATISTICS = ["signed_mean_pp", "signed_median_pp", "mean_absolute_pp",
              "median_absolute_pp", "signed_q10_pp", "signed_q50_pp", "signed_q90_pp",
              "signed_minimum_pp", "signed_maximum_pp", "signed_q90_minus_q10_pp"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def stats(values: np.ndarray) -> dict:
    q10, q50, q90 = np.quantile(values, [0.1, 0.5, 0.9], method="linear")
    return {
        "arrangements": len(values),
        "signed_mean_pp": float(values.mean()),
        "signed_median_pp": float(np.median(values)),
        "mean_absolute_pp": float(np.abs(values).mean()),
        "median_absolute_pp": float(np.median(np.abs(values))),
        "signed_q10_pp": float(q10), "signed_q50_pp": float(q50), "signed_q90_pp": float(q90),
        "signed_minimum_pp": float(values.min()), "signed_maximum_pp": float(values.max()),
        "signed_q90_minus_q10_pp": float(q90 - q10),
        "positive": int((values > SIGN_TOLERANCE_PP).sum()),
        "negative": int((values < -SIGN_TOLERANCE_PP).sum()),
        "ties": int((np.abs(values) <= SIGN_TOLERANCE_PP).sum()),
    }


def independently_summarize(values: list[float]) -> dict:
    ordered = sorted(values)
    def quantile(q: float) -> float:
        index = (len(ordered) - 1) * q
        lo, hi = math.floor(index), math.ceil(index)
        return ordered[lo] + (ordered[hi] - ordered[lo]) * (index - lo)
    return {
        "arrangements": len(values),
        "signed_mean_pp": math.fsum(values) / len(values),
        "signed_median_pp": statistics.median(values),
        "mean_absolute_pp": math.fsum(map(abs, values)) / len(values),
        "median_absolute_pp": statistics.median(map(abs, values)),
        "signed_q10_pp": quantile(0.1), "signed_q50_pp": quantile(0.5),
        "signed_q90_pp": quantile(0.9),
        "signed_minimum_pp": min(values), "signed_maximum_pp": max(values),
        "signed_q90_minus_q10_pp": quantile(0.9) - quantile(0.1),
        "positive": sum(x > SIGN_TOLERANCE_PP for x in values),
        "negative": sum(x < -SIGN_TOLERANCE_PP for x in values),
        "ties": sum(abs(x) <= SIGN_TOLERANCE_PP for x in values),
    }


def main() -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = ROOT / "results/music_arrangement_distribution_v1" / stamp
    out.mkdir(parents=True, exist_ok=False)
    source_code = Path(__file__).resolve()
    assert b"\r" not in source_code.read_bytes(), "Freeze code as LF before execution"
    shutil.copyfile(source_code, out / source_code.name)
    started = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "started; no new fitting", "script_sha256": sha(source_code),
        "unit": "MAE percentage-point difference: ratio10 minus balanced",
        "sign_tolerance_pp": SIGN_TOLERANCE_PP,
        "quantile_rule": "linear interpolation at (n-1)*q, NumPy method=linear",
        "sources": [],
    }
    save_json(out / "started.json", started)
    try:
        frames = []
        python_groups: dict[tuple, list[float]] = {}
        python_donor_groups: dict[tuple, list[float]] = {}
        inputs = []
        for label, folder in SOURCES.items():
            path = folder / "reference_arrangement_penalties.csv"
            review_path = folder / "review.json"
            review = json.loads(review_path.read_bytes())
            listed = next(x for x in review["artifacts"] if x["path"] == path.relative_to(ROOT).as_posix())
            assert sha(path) == listed["sha256"]
            assert sha(ROOT / review["producer"]["path"]) == review["producer"]["sha256"]
            frame = pd.read_csv(path)
            frame["source_csv_line"] = np.arange(len(frame)) + 2
            frame = frame.loc[frame.level.eq("ratio10")].copy()
            expected_budgets = {60, 120, 300} if label == "original_grouping" else {60, 300}
            assert set(frame.budget) == expected_budgets and set(frame.method) == set(METHODS)
            assert frame.held_out.nunique() == 14
            assert not frame.duplicated(["held_out", "triple_id", "dominant_donor", "method", "budget"]).any()
            assert frame.groupby(["held_out", "method", "budget"]).size().eq(12).all()
            assert frame.groupby(["method", "budget"]).size().eq(168).all()
            assert np.isfinite(frame.mae_penalty_pp).all()
            frame.insert(0, "source", label)
            frame["source_csv"] = path.relative_to(ROOT).as_posix()
            frame["source_csv_sha256"] = sha(path)
            frame["source_review"] = review_path.relative_to(ROOT).as_posix()
            frame["source_review_sha256"] = sha(review_path)
            frame["unit"] = "percentage points"
            frame["contrast"] = "MAE(ratio10) - MAE(balanced)"
            frame["targets_averaged"] = 60
            frame["reference_draws_averaged"] = 3
            frames.append(frame)
            # This path rereads original text with stdlib CSV rather than using
            # NumPy/Pandas values or summary functions.
            with path.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    if row["level"] != "ratio10":
                        continue
                    key = (label, row["method"], int(row["budget"]))
                    value = float(row["mae_penalty_pp"])
                    python_groups.setdefault(key, []).append(value)
                    python_donor_groups.setdefault(key + (row["held_out"],), []).append(value)
            started["sources"].append({"source": label, "path": path.relative_to(ROOT).as_posix(),
                "sha256": sha(path), "review_sha256": sha(review_path),
                "source_producer": review["producer"], "selected_rows": len(frame)})
            inputs.extend([path, review_path, ROOT / review["producer"]["path"]])
        full = pd.concat(frames, ignore_index=True)
        assert len(full) == 1680
        full.to_csv(out / "all_arrangement_source_rows.csv", index=False, lineterminator="\n")
        rows = []
        donor_rows = []
        for key, group in full.groupby(KEYS, sort=True):
            rows.append(dict(zip(KEYS, key), **stats(group.mae_penalty_pp.to_numpy()),
                             held_out_donors=14, arrangements_per_donor=12,
                             unit="percentage points"))
        for key, group in full.groupby(KEYS + ["held_out"], sort=True):
            donor_rows.append(dict(zip(KEYS + ["held_out"], key),
                                   **stats(group.mae_penalty_pp.to_numpy()), unit="percentage points"))
        summary = pd.DataFrame(rows)
        donor_summary = pd.DataFrame(donor_rows)
        summary.to_csv(out / "table_s13_distribution_summary.csv", index=False, lineterminator="\n")
        donor_summary.to_csv(out / "donor_distribution_summary.csv", index=False, lineterminator="\n")
        differences = []
        for table, independent_groups, keys in [(summary, python_groups, KEYS),
                                                (donor_summary, python_donor_groups, KEYS + ["held_out"])]:
            for row in table.to_dict("records"):
                independently = independently_summarize(independent_groups[tuple(row[k] for k in keys)])
                for column in STATISTICS + ["arrangements", "positive", "negative", "ties"]:
                    differences.append(abs(row[column] - independently[column]))
                assert row["positive"] + row["negative"] + row["ties"] == row["arrangements"]
        assert max(differences) < 1e-12
        pooled = donor_summary.groupby(KEYS).signed_mean_pp.mean().sort_index()
        direct = summary.set_index(KEYS).signed_mean_pp.sort_index()
        assert np.max(np.abs(pooled.to_numpy() - direct.to_numpy())) < 1e-12

        configure_style(6.5)
        fig = figure_mm(170, 151)
        axes = fig.subplots(2, 3, sharex=True, sharey=True)
        fig.subplots_adjust(left=0.148, right=0.982, bottom=0.105, top=0.885, wspace=0.17, hspace=0.26)
        fig.text(0.148, 0.973, "Original reference grouping: 10:1:1 versus balanced", fontsize=7, va="top")
        fig.text(0.148, 0.944, "Each point averages 60 targets and 3 reference draws; 12 arrangements per held-out donor.", fontsize=6, va="top")
        original = full.loc[full.source.eq("original_grouping")].copy()
        donors = sorted(original.held_out.unique())
        plotted = []
        for row_index, method in enumerate(METHODS):
            color = COLORS["vermillion"] if method == "music_weighted" else COLORS["blue"]
            marker = "s" if method == "music_weighted" else "o"
            method_label = "Weighted" if method == "music_weighted" else "Ordinary NNLS"
            for col_index, budget in enumerate([60, 120, 300]):
                ax = axes[row_index, col_index]
                subset = original.loc[original.method.eq(method) & original.budget.eq(budget)]
                for donor_index, donor in enumerate(donors):
                    if donor_index % 2 == 0:
                        ax.axhspan(donor_index - 0.5, donor_index + 0.5, facecolor="#F1F1F1", zorder=0)
                    group = subset.loc[subset.held_out.eq(donor)].sort_values(["triple_id", "dominant_donor"]).copy()
                    group["plot_y"] = donor_index + np.linspace(-0.28, 0.28, 12)
                    group["panel"] = chr(ord("a") + row_index * 3 + col_index)
                    ax.scatter(group.mae_penalty_pp, group.plot_y, s=5.5, color=color, marker=marker,
                               linewidths=0.25, edgecolors="white", zorder=3, alpha=0.87)
                    plotted.append(group)
                ax.axvline(0, color="black", linewidth=0.6, linestyle=(0, (2, 2)), zorder=2)
                ax.set_yticks(range(len(donors)), labels=donors, fontsize=5.8)
                ax.set_ylim(len(donors) - 0.45, -0.55)
                ax.set_xlim(-2.5, 3)
                ax.set_xticks([-2, -1, 0, 1, 2, 3])
                ax.tick_params(axis="y", length=0, pad=2)
                ax.tick_params(axis="x", labelbottom=True)
                ax.spines["left"].set_visible(False)
                ax.text(0, 1.055, chr(ord("a") + row_index * 3 + col_index), transform=ax.transAxes,
                        fontsize=8, fontweight="bold", ha="left", va="bottom")
                ax.text(0.11, 1.055, f"{method_label}, B = {budget}", transform=ax.transAxes,
                        fontsize=6.5, ha="left", va="bottom")
        fig.text(0.565, 0.047, "Arrangement-averaged MAE difference (percentage points)", ha="center", fontsize=7)
        fig.text(0.565, 0.018, "Positive: higher error with 10:1:1 allocation. Points are dependent configurations, not independent samples.",
                 ha="center", fontsize=5.8)
        figure_source = pd.concat(plotted, ignore_index=True)
        assert len(figure_source) == 1008
        caption = (
            "Figure S4. Arrangement-level MAE contrasts in the original external reference grouping. "
            "Panels a-c show the weighted MuSiC output and panels d-f the ordinary NNLS output at "
            "per-cell-type budgets B = 60, 120 and 300. Each point is MAE under 10:1:1 within-type "
            "donor allocation minus MAE under balanced allocation, in percentage points, averaged "
            "over the same 60 target mixtures and three reference draws for a fixed held-out donor, "
            "reference triple and dominant donor. Each panel contains 168 arrangements, displayed "
            "as 12 arrangements within each of 14 held-out-donor rows. Small deterministic vertical "
            "offsets separate points within donor rows and have no scientific meaning. The dashed "
            "line is zero; positive values indicate higher MAE under 10:1:1 allocation. All panels "
            "share the same horizontal scale. No population error bars or confidence intervals are "
            "shown: arrangements share target mixtures and reference sources and are not independent "
            "biological replicates. Mean absolute contrast decreases across these budgets for the "
            "weighted output, but the signed 10th-to-90th percentile span does not decrease monotonically."
        )
        bundle = save_figure_bundle(fig, out / "figure_s4_arrangement_distribution", figure_source, caption,
            inputs, source_code, {"title": "Arrangement-level MAE contrasts by held-out donor",
             "scope": "Original external grouping only; alternative grouping retained in tables",
             "rows_per_panel": 168, "donor_rows_per_panel": 14, "arrangements_per_donor": 12,
             "sample_unit": "dependent reference arrangement after target/draw averaging",
             "error_bars": "none", "new_predictions": 0, "sign_tolerance_pp": SIGN_TOLERANCE_PP})
        plt.close(fig)
        save_json(out / "independent_summary_check.json", {
            "status": "passed", "method": "stdlib csv reread, math.fsum, statistics.median, explicit linear quantile",
            "overall_groups_checked": len(summary), "donor_groups_checked": len(donor_summary),
            "scalar_comparisons": len(differences), "maximum_absolute_difference": max(differences),
            "checks": ["Frozen source hashes", "168 arrangements per source/method/budget",
                       "14 donors and 12 arrangements per donor", "signed counts sum to n",
                       "equal-donor mean equals arrangement mean", "no new fitting or predictions"],
            "new_fits": 0, "new_predictions": 0})
        result = {**started, "status": "completed; visual review pending", "output_directory": out.relative_to(ROOT).as_posix(),
                  "summary_rows": len(summary), "donor_summary_rows": len(donor_summary),
                  "retained_source_rows": len(full), "plotted_rows": len(figure_source),
                  "figure_dimensions_mm": bundle["figure_dimensions_mm"],
                  "independent_maximum_difference": max(differences),
                  "artifacts": [{"file": p.name, "sha256": sha(p), "bytes": p.stat().st_size}
                                for p in sorted(out.iterdir()) if p.is_file()]}
        save_json(out / "run.json", result)
        print(out.relative_to(ROOT).as_posix())
        print(summary[["source", "method", "budget", "signed_mean_pp", "mean_absolute_pp", "signed_q10_pp", "signed_q90_pp", "positive", "negative"]].to_string(index=False))
    except Exception:
        save_json(out / "failure.json", {"status": "failed; preserve this run", "traceback": traceback.format_exc(),
                    "producer_sha256": sha(source_code), "utc": datetime.now(timezone.utc).isoformat()})
        raise


if __name__ == "__main__":
    main()
