"""Audit frozen MuSiC diagnostic outputs and export descriptive summaries/figure.

This review does not replace the 60-target scientific endpoint. No population
test, interval, causal mechanism claim or new fitted condition is introduced.
v2 replaces the unsupported Arial superscript-minus glyph with math text only.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, NullLocator
import numpy as np
import pandas as pd

from scientific_style import COLORS, configure_style, figure_mm, save_figure_bundle, sha256

ROOT = Path(__file__).resolve().parents[1]
BUDGETS = [60, 120, 300]
QUOTAS = [5, 10, 20, 25, 40, 50, 100, 250]
TRAJECTORIES = {"minor": [5, 10, 25], "balanced": [20, 40, 100], "dominant": [50, 100, 250]}
PROFILE_METRICS = ["theta_variance_trace", "theta_pairwise_tv", "S_cv"]
WEIGHT_METRICS = ["intersection_n", "union_n", "Jaccard", "spearman_shared",
                  "normalized_weight_tv_intersection", "balanced_weight_fraction_on_intersection",
                  "unequal_weight_fraction_on_intersection"]


def resolve(path):
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def relative(path):
    p = Path(path).resolve()
    return p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else str(p)


def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def finite(frame, columns):
    assert np.isfinite(frame[columns].to_numpy(float)).all(), columns


def equal_group(frame, keys, metrics):
    # Undefined rank correlations remain undefined, rather than silently giving
    # different targets or arrangements different weight after NA exclusion.
    aggregation = {m: (lambda s: s.mean(skipna=False)) if m == "spearman_shared" else "mean" for m in metrics}
    return frame.groupby(keys, observed=True, sort=True, dropna=False).agg(aggregation).reset_index()


def check_group_size(frame, keys, expected):
    sizes = frame.groupby(keys, observed=True, sort=True, dropna=False).size()
    assert sizes.eq(expected).all(), (keys,expected,sizes.value_counts().to_dict())


def directions(values):
    v = np.asarray(values, dtype=float)
    assert np.isfinite(v).all()
    return {"n": len(v), "negative": int((v < 0).sum()), "zero": int((v == 0).sum()),
            "positive": int((v > 0).sum()), "mean": float(v.mean()),
            "minimum": float(v.min()), "maximum": float(v.max())}


def make_trajectories(frame, ids):
    rows = []
    for name, quotas in TRAJECTORIES.items():
        for metric in PROFILE_METRICS:
            wide = frame.pivot(index=ids, columns="reference_cells", values=metric)
            for index, row in wide.iterrows():
                index = index if isinstance(index, tuple) else (index,)
                low, middle, high = (float(row[q]) for q in quotas)
                d1, d2 = middle - low, high - middle
                rows.append({**dict(zip(ids, index)), "trajectory": name, "metric": metric,
                    "low_quota": quotas[0], "middle_quota": quotas[1], "high_quota": quotas[2],
                    "low_value": low, "middle_value": middle, "high_value": high,
                    "endpoint_change": high - low,
                    "endpoint_ratio": high / low if low > 0 else np.nan,
                    "ratio_defined": low > 0, "endpoint_direction": "lower" if high < low else "higher" if high > low else "equal",
                    "first_adjacent_change": d1, "second_adjacent_change": d2,
                    "both_adjacent_nonincreasing": d1 <= 0 and d2 <= 0,
                    "both_adjacent_strictly_decreasing": d1 < 0 and d2 < 0,
                    "any_adjacent_increase": d1 > 0 or d2 > 0})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Completed diagnostic directory (or its run.json)")
    args = parser.parse_args()
    rp = args.run.resolve()
    rp = rp / "run.json" if rp.is_dir() else rp
    run = json.loads(rp.read_bytes()); directory = rp.parent
    assert run["status"] == "completed bounded diagnostics; descriptive review pending"
    checked_sources = []
    for item in run["artifacts"] + run["sources"]:
        path = resolve(item["path"])
        assert sha256(path) == item["sha256"], item["path"]
        checked_sources.append({"path": relative(path), "sha256": item["sha256"]})
    sp = resolve(run["selection"])
    assert sha256(sp) == run["selection_sha256"]
    selection = json.loads(sp.read_bytes())
    for item in selection["sources"] + selection["artifacts"]:
        assert sha256(resolve(item["path"])) == item["sha256"], item["path"]
    donors, types = selection["donors"], selection["cell_types"]
    assert len(donors) == 14 and len(types) == 6 and len(selection["mixture_ids"]) == 8
    cases = pd.read_csv(resolve(selection["cases_file"]), sep="\t", keep_default_na=False)
    references = pd.read_csv(resolve(selection["references_file"]), sep="\t", keep_default_na=False)
    targets = pd.read_csv(resolve(selection["targets_file"]), sep="\t", keep_default_na=False)
    assert len(cases) == 2016 and len(references) == 504 and len(targets) == 112
    assert not cases.case_id.duplicated().any() and not references.reference_id.duplicated().any()
    assert not cases.duplicated(["reference_id", "held_out"]).any()
    assert cases.groupby("reference_id").size().eq(4).all()
    expected_counts = {"predictions": 2304, "basis": 216, "weights": 1152,
                       "weight_pairs": 864, "profile_dispersion": 144}
    frames = {name: [] for name in expected_counts}
    assert len(run["folds"]) == 14
    assert {f["triple_key"] for f in run["folds"]} == set(selection["triple_keys"])
    for triple in selection["triple_keys"]:
        fold = json.loads((directory / triple / "fold.json").read_bytes())
        assert fold["status"] == "completed" and fold["basis_checks_passed"] and fold["maxiter_count"] == 0
        assert fold["traced_fits"] == 1152 and fold["reference_matrices"] == 36 and fold["logical_cases"] == 144
        assert fold["prediction_rows"] == 2304 and fold["profile_dispersion_rows"] == 144
        assert max(fold["basis_checks_maximum_difference"].values()) <= 1e-10
        for name, count in expected_counts.items():
            # Missing Spearman coefficients use the explicit R NA marker. Empty
            # dominant-donor identifiers are legitimate balanced metadata.
            frame = pd.read_csv(directory / triple / f"{name}.csv", keep_default_na=False, na_values=["NA"])
            assert len(frame) == count, (triple, name, len(frame))
            if name != "profile_dispersion":
                assert frame.triple_key.eq(triple).all()
            else:
                frame["source_triple_key"] = triple
            frames[name].append(frame)
    data = {name: pd.concat(parts, ignore_index=True) for name, parts in frames.items()}
    predictions, basis, weights, pairs, profiles = (data[n] for n in expected_counts)
    ref_index = references.set_index("reference_id")
    ref_fields = ["triple_key","block","budget","level","dominant_donor"]
    target_index = targets.set_index("target_name")
    permitted_reference_targets = set(map(tuple,cases[["reference_id","held_out"]].to_numpy()))
    for frame in [basis,weights,pairs]:
        assert frame[ref_fields].equals(ref_index.loc[frame.reference_id,ref_fields].reset_index(drop=True))
    for frame in [weights,pairs]:
        assert np.array_equal(frame.held_out.to_numpy(),target_index.loc[frame.target_name,"held_out"].to_numpy())
        assert np.array_equal(frame.mixture_id.to_numpy(),target_index.loc[frame.target_name,"mixture_id"].to_numpy())
        assert set(map(tuple,frame[["reference_id","held_out"]].to_numpy())) <= permitted_reference_targets
    matched_bases = ref_index.loc[pairs.balanced_reference_id]
    assert matched_bases.level.eq("balanced").all()
    assert np.array_equal(matched_bases[["triple_key","block","budget"]].to_numpy(),
                          pairs[["triple_key","block","budget"]].to_numpy())
    assert len(predictions) == 32256 and not predictions.duplicated(["case_id", "method", "mixture_id"]).any()
    assert set(predictions.method) == {"music_weighted", "music_nnls"}
    assert set(predictions.mixture_id) == set(selection["mixture_ids"])
    assert predictions.groupby(["case_id", "method"]).size().eq(8).all()
    meta = ["reference_id", "held_out", "triple_id", "triple_key", "block", "budget", "level", "dominant_donor"]
    expected = cases.set_index("case_id").loc[predictions.case_id, meta].reset_index(drop=True)
    assert predictions[meta].equals(expected)
    truth = predictions[[f"true_{k}" for k in range(6)]].to_numpy(float)
    pred = predictions[[f"pred_{k}" for k in range(6)]].to_numpy(float)
    assert np.isfinite(truth).all() and np.isfinite(pred).all()
    assert pred.min() >= -1e-10 and np.max(abs(pred.sum(axis=1)-1)) <= 1e-10
    assert np.max(abs(truth-target_index.loc[predictions.target_name,[f"true_{k}" for k in range(6)]].to_numpy(float))) <= 1e-12
    error = pred-truth
    metric_delta = max(float(np.max(abs(np.mean(abs(error),axis=1)-predictions.mae))),
                       float(np.max(abs(np.sqrt(np.mean(error**2,axis=1))-predictions.rmse))))
    assert metric_delta <= 1e-10
    assert len(basis) == 3024 and not basis.duplicated(["reference_id", "cell_type"]).any()
    assert basis.groupby("reference_id").size().eq(6).all() and set(basis.cell_type) == set(types)
    basis_metrics = ["sum_sigma", "theta_squared_norm", "mean_cell_size", "basis_genes", "fixed_profile_genes"]
    finite(basis, basis_metrics)
    assert (basis.sum_sigma >= 0).all() and (basis.theta_squared_norm > 0).all() and (basis.mean_cell_size > 0).all()
    assert basis.fixed_profile_genes.eq(30172).all() and basis.basis_genes.between(1,30172).all()
    weight_keys = ["reference_id", "target_name"]
    assert len(weights) == 16128 and not weights.duplicated(weight_keys).any()
    assert weights.groupby("reference_id").size().eq(32).all()
    finite(weights, ["n_weight_genes", "normalized_weight_effective_fraction", "top1pct_weight_share"])
    assert weights.n_weight_genes.between(2,30172).all()
    assert weights.normalized_weight_effective_fraction.between(0,1+1e-12).all()
    assert weights.top1pct_weight_share.between(0,1+1e-12).all()
    assert len(pairs) == 12096 and not pairs.duplicated(weight_keys).any()
    assert pairs.level.eq("ratio10").all() and pairs.groupby("reference_id").size().eq(32).all()
    finite(pairs, [m for m in WEIGHT_METRICS if m != "spearman_shared"])
    defined = pairs.spearman_defined.astype(str).str.lower().eq("true")
    assert np.array_equal(defined.to_numpy(), pairs.spearman_shared.notna().to_numpy())
    assert pairs.loc[defined,"spearman_shared"].between(-1-1e-12,1+1e-12).all()
    for metric in ["Jaccard", "normalized_weight_tv_intersection", "balanced_weight_fraction_on_intersection", "unequal_weight_fraction_on_intersection"]:
        assert pairs[metric].between(0,1+1e-12).all()
    assert (pairs.intersection_n > 1).all() and (pairs.union_n >= pairs.intersection_n).all()
    assert np.max(abs(pairs.Jaccard-pairs.intersection_n/pairs.union_n)) <= 1e-12
    weight_index = weights.set_index(weight_keys)
    bn = weight_index.loc[pd.MultiIndex.from_frame(pairs[["balanced_reference_id","target_name"]]), "n_weight_genes"].to_numpy()
    un = weight_index.loc[pd.MultiIndex.from_frame(pairs[weight_keys]), "n_weight_genes"].to_numpy()
    assert np.array_equal(pairs.union_n.to_numpy(),bn+un-pairs.intersection_n.to_numpy())
    pairs["spearman_defined"] = defined
    # Verify all duplicate profile statistics exactly before deduplicating.
    pkeys = ["donor", "cell_type", "reference_cells"]
    numeric_profile = ["source_inventory_N", "sampling_fraction", "draws", "fixed_profile_genes"] + PROFILE_METRICS
    finite(profiles,numeric_profile)
    assert profiles.reference_cells.isin(QUOTAS).all() and profiles.draws.eq(3).all()
    assert profiles.fixed_profile_genes.eq(30172).all() and profiles.sampling_fraction.between(0,1).all()
    assert np.max(abs(profiles.sampling_fraction-profiles.reference_cells/profiles.source_inventory_N)) <= 1e-14
    duplicate_audit = []
    for key, group in profiles.groupby(pkeys, observed=True, sort=True):
        assert len(group) == 3, key
        values = group[numeric_profile].to_numpy(float)
        assert np.array_equal(values,np.broadcast_to(values[0],values.shape)), key
        duplicate_audit.append({**dict(zip(pkeys,key)),"repeated_triple_records":len(group),
                                "maximum_numeric_difference":0.0})
    unique_profiles = profiles.drop_duplicates(pkeys).drop(columns="source_triple_key").sort_values(pkeys).reset_index(drop=True)
    assert len(unique_profiles) == 672 and len(duplicate_audit) == 672
    assert unique_profiles.groupby(["donor","cell_type"]).size().eq(8).all()
    assert unique_profiles.groupby(["donor","reference_cells"]).size().eq(6).all()
    assert set(unique_profiles.donor) == set(donors) and set(unique_profiles.cell_type) == set(types)
    assert (unique_profiles[PROFILE_METRICS] >= 0).all().all()
    inventory = pd.read_csv(resolve(selection["inventory_file"])).set_index("donor_id")
    inventory_types = dict(zip(types,["B","cM","T4","T8","ncM","NK"]))
    for row in unique_profiles.itertuples():
        assert row.source_inventory_N == inventory.loc[row.donor,inventory_types[row.cell_type]]

    out = ROOT / "results/music_mechanism_review" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name)
    tables = {}
    def export(name, frame):
        frame.to_csv(out/f"{name}.csv",index=False,lineterminator="\n",float_format="%.17g")
        tables[name] = len(frame)
        return frame

    export("profile_dispersion_all_triple_records",profiles)
    export("profile_duplicate_audit",pd.DataFrame(duplicate_audit))
    export("profile_dispersion_unique",unique_profiles)
    profile_donor = export("profile_dispersion_donor_means",equal_group(unique_profiles,["donor","reference_cells"],PROFILE_METRICS))
    profile_overall = export("profile_dispersion_equal_donor_means",equal_group(profile_donor,["reference_cells"],PROFILE_METRICS))
    type_trajectories = export("profile_donor_type_trajectories",make_trajectories(unique_profiles,["donor","cell_type"]))
    donor_trajectories = export("profile_donor_mean_trajectories",make_trajectories(profile_donor,["donor"]))
    assert len(type_trajectories) == 756 and len(donor_trajectories) == 126
    trajectory_summaries = []
    for scope, frame in [("donor_type",type_trajectories),("donor_mean",donor_trajectories)]:
        for (trajectory,metric), group in frame.groupby(["trajectory","metric"],observed=True):
            valid = group.ratio_defined.astype(bool)
            trajectory_summaries.append({"scope":scope,"trajectory":trajectory,"metric":metric,
                **directions(group.endpoint_change),"defined_ratios":int(valid.sum()),
                "mean_endpoint_ratio":float(group.loc[valid,"endpoint_ratio"].mean()) if valid.any() else np.nan,
                "median_endpoint_ratio":float(group.loc[valid,"endpoint_ratio"].median()) if valid.any() else np.nan,
                "both_adjacent_nonincreasing":int(group.both_adjacent_nonincreasing.sum()),
                "both_adjacent_strictly_decreasing":int(group.both_adjacent_strictly_decreasing.sum()),
                "any_adjacent_increase":int(group.any_adjacent_increase.sum()),
                "mean_low_value":group.low_value.mean(),"mean_middle_value":group.middle_value.mean(),
                "mean_high_value":group.high_value.mean(),
                "ratio_of_mean_endpoint_values":group.high_value.mean()/group.low_value.mean() if group.low_value.mean()>0 else np.nan})
    export("profile_trajectory_direction_summaries",pd.DataFrame(trajectory_summaries))

    # Reference summaries become held-out configurations only through this exact
    # four-case mapping, then each hierarchical dimension receives equal weight.
    export("basis_all_reference_types",basis)
    bkey = ["triple_key","block","budget","cell_type"]
    base = basis[basis.level.eq("balanced")][bkey+["reference_id","sum_sigma"]].rename(columns={"reference_id":"balanced_reference_id","sum_sigma":"balanced_sum_sigma"})
    assert not base.duplicated(bkey).any()
    paired_basis = basis[basis.level.eq("ratio10")].merge(base,on=bkey,validate="many_to_one")
    paired_basis["sigma_trace_contrast"] = paired_basis.sum_sigma-paired_basis.balanced_sum_sigma
    assert len(paired_basis) == 2268
    export("sigma_reference_type_pairs",paired_basis)
    case_map = cases[["reference_id","case_id","held_out","triple_id"]]
    expanded_basis = export("sigma_held_out_type_pairs",paired_basis.merge(case_map,on="reference_id",validate="many_to_many"))
    assert len(expanded_basis) == 9072 and not expanded_basis.duplicated(["case_id","cell_type"]).any()
    keys = ["held_out","triple_key","block","budget","dominant_donor"]
    check_group_size(expanded_basis,keys,6)
    sigma_type_mean = export("sigma_equal_type_means",equal_group(expanded_basis,keys,["sigma_trace_contrast"]))
    assert len(sigma_type_mean) == 1512
    check_group_size(sigma_type_mean,[k for k in keys if k!="block"],3)
    sigma_draw_mean = export("sigma_equal_draw_means",equal_group(sigma_type_mean,[k for k in keys if k!="block"],["sigma_trace_contrast"]))
    check_group_size(sigma_draw_mean,["held_out","triple_key","budget"],3)
    sigma_dominant_mean = export("sigma_equal_dominant_means",equal_group(sigma_draw_mean,["held_out","triple_key","budget"],["sigma_trace_contrast"]))
    check_group_size(sigma_dominant_mean,["held_out","budget"],4)
    sigma_donor = export("sigma_donor_means",equal_group(sigma_dominant_mean,["held_out","budget"],["sigma_trace_contrast"]))
    check_group_size(sigma_donor,["budget"],14)
    sigma_overall = export("sigma_equal_donor_means",equal_group(sigma_donor,["budget"],["sigma_trace_contrast"]))
    assert len(sigma_donor) == 42

    export("weights_all_targets",weights)
    export("weight_pairs_all_targets",pairs)
    check_group_size(pairs,keys,8)
    weight_target_mean = export("weight_pairs_equal_target_means",equal_group(pairs,keys,WEIGHT_METRICS))
    check_group_size(weight_target_mean,[k for k in keys if k!="block"],3)
    weight_draw_mean = export("weight_pairs_equal_draw_means",equal_group(weight_target_mean,[k for k in keys if k!="block"],WEIGHT_METRICS))
    check_group_size(weight_draw_mean,["held_out","triple_key","budget"],3)
    weight_dominant_mean = export("weight_pairs_equal_dominant_means",equal_group(weight_draw_mean,["held_out","triple_key","budget"],WEIGHT_METRICS))
    check_group_size(weight_dominant_mean,["held_out","budget"],4)
    weight_donor = export("weight_pairs_donor_means",equal_group(weight_dominant_mean,["held_out","budget"],WEIGHT_METRICS))
    check_group_size(weight_donor,["budget"],14)
    weight_overall = export("weight_pairs_equal_donor_means",equal_group(weight_donor,["budget"],WEIGHT_METRICS))
    assert len(weight_donor) == 42
    availability = pairs.assign(undefined_spearman=(~pairs.spearman_defined).astype(int)).groupby(["held_out","budget"]).agg(
        comparison_rows=("target_name","size"),undefined_spearman=("undefined_spearman","sum")).reset_index()
    assert availability.comparison_rows.eq(288).all()
    export("weight_rank_correlation_availability",availability)

    # MAE summaries use the eight frozen target IDs only and are diagnostic,
    # leaving the original 60-target endpoint and reported manuscript intact.
    predictions["mae_pp"] = predictions.mae*100
    predictions["rmse_pp"] = predictions.rmse*100
    export("subset_predictions",predictions)
    absolute = export("subset_donor_absolute_errors",equal_group(predictions,["held_out","method","budget","level"],["mae_pp","rmse_pp"]))
    export("subset_equal_donor_absolute_errors",equal_group(absolute,["method","budget","level"],["mae_pp","rmse_pp"]))
    pkey = ["held_out","triple_key","block","budget","target_name","method"]
    pbase = predictions[predictions.level.eq("balanced")][pkey+["reference_id","mae_pp","rmse_pp"]].rename(columns={
        "reference_id":"balanced_reference_id","mae_pp":"balanced_mae_pp","rmse_pp":"balanced_rmse_pp"})
    assert not pbase.duplicated(pkey).any()
    penalty = predictions[predictions.level.eq("ratio10")].merge(pbase,on=pkey,validate="many_to_one")
    penalty["mae_penalty_pp"] = penalty.mae_pp-penalty.balanced_mae_pp
    penalty["rmse_penalty_pp"] = penalty.rmse_pp-penalty.balanced_rmse_pp
    assert len(penalty) == 24192
    export("subset_target_allocation_contrasts",penalty)
    penalty_metrics = ["mae_penalty_pp","rmse_penalty_pp"]
    stages = equal_group(penalty,keys+["method"],penalty_metrics)
    stages = equal_group(stages,[k for k in keys if k!="block"]+["method"],penalty_metrics)
    stages = equal_group(stages,["held_out","triple_key","budget","method"],penalty_metrics)
    penalty_donor = export("subset_donor_allocation_contrasts",equal_group(stages,["held_out","budget","method"],penalty_metrics))
    penalty_overall = export("subset_equal_donor_allocation_contrasts",equal_group(penalty_donor,["budget","method"],penalty_metrics))
    endpoint = penalty_donor.pivot(index=["held_out","method"],columns="budget",values="mae_penalty_pp")
    endpoint["I300_minus_I60_pp"] = endpoint[300]-endpoint[60]
    export("subset_donor_budget_contrasts",endpoint.rename(columns={b:f"I{b}_pp" for b in BUDGETS}).reset_index())
    absolute_wide = absolute.pivot(index=["held_out","method","budget"],columns="level",values="mae_pp")
    penalty_identity = penalty_donor.set_index(["held_out","method","budget"]).mae_penalty_pp.reindex(absolute_wide.index)
    max_pair_difference = float(np.max(abs(absolute_wide.ratio10-absolute_wide.balanced-penalty_identity)))
    assert max_pair_difference <= 1e-10

    # Long-form figure data contain both the actual raw statistic and displayed
    # scaling, without NaN filler columns that obscure export validation.
    figure_rows = []
    panels = [("a",profile_donor[profile_donor.reference_cells.isin([5,10,25])],"donor","reference_cells","theta_variance_trace",1.0),
              ("b",sigma_donor,"held_out","budget","sigma_trace_contrast",1e6),
              ("c",weight_donor,"held_out","budget","normalized_weight_tv_intersection",1.0)]
    for panel, frame, donor_col, xcol, metric, scale in panels:
        for row in frame.to_dict("records"):
            figure_rows.append({"panel":panel,"series":"individual_donor","donor":row[donor_col],
                "x":row[xcol],"metric":metric,"value_raw":row[metric],"display_scale":scale,"value_plotted":row[metric]*scale})
        for x, group in frame.groupby(xcol,observed=True):
            assert len(group)==14
            figure_rows.append({"panel":panel,"series":"equal_donor_mean","donor":"all_14",
                "x":x,"metric":metric,"value_raw":group[metric].mean(),"display_scale":scale,"value_plotted":group[metric].mean()*scale})
    figure_data = pd.DataFrame(figure_rows)
    assert len(figure_data)==135
    style = configure_style(7)
    fig = figure_mm(183,95)
    axes = [fig.add_axes([left,.34,.21,.54]) for left in [.105,.435,.765]]
    titles = ["Reference-profile\ndraw dispersion","Between-donor\nvariance contrast","Fitted gene-weight\nallocation difference"]
    labels = ["Across-draw variance trace\n(relative-expression units²)",
              "Σ trace: 10:1:1 − balanced\n(relative-expression units²; " + r"$\times10^{-6}$" + ")",
              "Normalized-weight distance\n(shared-support TV)"]
    for panel,ax,title,label in zip("abc",axes,titles,labels):
        part = figure_data[figure_data.panel.eq(panel)];xs = [5,10,25] if panel=="a" else BUDGETS
        for donor in donors:
            y = part[(part.series.eq("individual_donor")) & part.donor.eq(donor)].set_index("x").loc[xs,"value_plotted"]
            ax.plot(range(3),y,color="#A7A7A7",linewidth=.55,marker="o",markersize=1.7,zorder=2)
        mean = part[part.series.eq("equal_donor_mean")].set_index("x").loc[xs,"value_plotted"]
        ax.plot(range(3),mean,color=COLORS["vermillion"],linewidth=1.0,marker="s",markersize=3.3,zorder=4)
        ax.set_xticks(range(3),[str(x) for x in xs]);ax.set_xlim(-.13,2.13)
        ax.set_xlabel("Cells per donor and type" if panel=="a" else "Total reference cells per type",labelpad=5)
        ax.set_ylabel(label,labelpad=5);ax.set_title(title,pad=7)
        ax.text(-.28,1.10,panel,transform=ax.transAxes,ha="left",va="bottom",fontsize=8,fontweight="bold")
        if panel=="a":
            assert (part.value_plotted>0).all()
            ax.set_yscale("log")
            lo,hi = ax.get_ylim()
            exponents = range(int(np.floor(np.log10(lo))),int(np.ceil(np.log10(hi)))+1)
            ticks = [10.0**e for e in exponents if lo<=10.0**e<=hi]
            if len(ticks)<2:
                ticks = np.geomspace(lo*1.05,hi/1.05,3).tolist()
            ax.set_yticks(ticks,[f"{value:.1e}" for value in ticks])
            ax.set_ylim(lo,hi)
            ax.yaxis.set_minor_locator(NullLocator());ax.yaxis.set_minor_formatter(NullFormatter())
        elif panel=="b":
            ax.axhline(0,color="black",linewidth=.55,linestyle=":",zorder=1)
        else:
            ax.set_ylim(bottom=0)
    fig.legend(handles=[Line2D([0],[0],color="#A7A7A7",marker="o",markersize=2,label="All 14 donor means"),
                        Line2D([0],[0],color=COLORS["vermillion"],marker="s",markersize=3.3,label="Equal-donor mean")],
               loc="center",bbox_to_anchor=(.53,.16),ncol=2,columnspacing=2)
    fig.text(.105,.035,"Conditional descriptive diagnostics; three fixed reference draws.\nShared simulations; no confidence intervals or causal interpretation.",ha="left",va="bottom",fontsize=7)
    caption = """# Conditional reference-profile and gene-weight diagnostics

**a**, Trace of the sample covariance of relative-expression profiles across three fixed reference-cell draws, summed over all 30,172 genes and averaged equally over six types within each reference donor. Quotas of 5, 10 and 25 cells per donor/type correspond to the minor-donor counts under 10:1:1 allocation at budgets of 60, 120 and 300 cells per type. The vertical axis is logarithmic. **b**, Paired difference in official MuSiC between-donor variance trace (10:1:1 minus balanced), averaged equally over types, draws, dominant choices and reference triples within each held-out donor. Displayed values are multiplied by 10⁶; the underlying unit is squared relative expression. **c**, Total-variation distance between fitted gene weights from 10:1:1 and balanced references, separately normalized on their shared positive finite gene support. Comparisons average eight fixed targets, three draws, three dominant choices and four triples within each held-out donor. Support overlap and retained weight fractions are reported in the accompanying complete tables; missing weights are not treated as zero.

Grey lines show all 14 donor summaries; vermilion squares show their equal-weight mean. Donors in panel a are reference-profile donors; those in panels b and c label overlapping held-out evaluation configurations. Horizontal positions are categorical. Panel b preserves positive and negative contrasts. Three draws provide imprecise conditional dispersion estimates, and these summaries do not establish biological-variance decomposition or mediation of prediction error. No confidence intervals, significance tests or causal claims are shown. The eight-target subset diagnoses fitted quantities and does not replace the original 60-target-per-donor performance endpoint. This is a research diagnostic figure, not a replacement manuscript figure.
"""
    figure_manifest = save_figure_bundle(fig,out/"mechanism_diagnostics",figure_data,caption,
        inputs=[rp,out/"profile_dispersion_donor_means.csv",out/"sigma_donor_means.csv",out/"weight_pairs_donor_means.csv"],
        producer=Path(__file__),details={"title":"Conditional MuSiC mechanism diagnostics","style":style,
            "scope":"Post-hoc descriptive diagnostics; three draws; no causal attribution",
            "donors":14,"types":6,"targets_per_held_out_donor":8,"manual_visual_review_completed":False})
    plt.close(fig)
    review = {"status":"passed bounded diagnostic numerical review; manual figure review pending",
        "run":relative(rp),"run_sha256":sha256(rp),"selection":relative(sp),"selection_sha256":sha256(sp),
        "producer":{"path":relative(__file__),"sha256":sha256(__file__)},"checked_source_and_artifact_hashes":checked_sources,
        "complete_triples":14,"prediction_rows":len(predictions),"unique_profile_rows":len(unique_profiles),
        "duplicate_profile_rows":len(profiles),"profile_duplicate_maximum_difference":0.0,
        "maximum_prediction_metric_recalculation_difference":metric_delta,
        "maximum_subset_paired_aggregation_difference_pp":max_pair_difference,
        "undefined_weight_rank_correlations":int((~defined).sum()),
        "source_sampling_fraction_minimum":float(unique_profiles.sampling_fraction.min()),
        "source_sampling_fraction_maximum":float(unique_profiles.sampling_fraction.max()),
        "sigma_donor_directions_by_budget":{str(b):directions(g.sigma_trace_contrast) for b,g in sigma_donor.groupby("budget")},
        "weight_tv_equal_donor_means":{str(r.budget):float(r.normalized_weight_tv_intersection) for r in weight_overall.itertuples()},
        "subset_weighted_budget_contrast_pp":directions(endpoint.xs("music_weighted",level="method").I300_minus_I60_pp),
        "table_rows":tables,"figure_checks":figure_manifest["checks"],"manual_visual_review_completed":False,
        "scope":"All diagnostic results are descriptive. Reference folds, targets and cells overlap. Three draws do not identify population uncertainty; Sigma includes between-donor differences and sampling variation. No causal mechanism or 25-cell threshold is established; no measured-bulk validation is added.",
        "main_endpoint_replaced":False,"p_values_or_confidence_intervals_computed":False}
    review["artifacts"] = [{"path":relative(p),"sha256":sha256(p),"bytes":p.stat().st_size}
                           for p in sorted(out.iterdir()) if p.is_file()]
    save(out/"review.json",review)
    print("OUTPUT",relative(out),flush=True)
    print(json.dumps({k:review[k] for k in ["status","unique_profile_rows","undefined_weight_rank_correlations",
        "sigma_donor_directions_by_budget","weight_tv_equal_donor_means","subset_weighted_budget_contrast_pp"]},indent=2),flush=True)


if __name__ == "__main__":
    main()
