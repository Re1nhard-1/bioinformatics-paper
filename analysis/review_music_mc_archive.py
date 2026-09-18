"""Audit D041 predictions by triple; summarize conditional global-block precision."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METHODS = ["music_weighted", "music_nnls"]
BUDGETS = [60, 300]
LEVELS = ["balanced", "ratio10"]
METRICS = ["mae_pp", "rmse_pp"] + [f"absolute_error_{k}_pp" for k in range(6)]
PENALTIES = [name.replace("_pp", "_penalty_pp") for name in METRICS]
TOL = 1e-10


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8", newline="\n")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def direction(values):
    x = np.asarray(values, dtype=float)
    require(len(x) > 0 and np.isfinite(x).all(), "Invalid direction/distribution values")
    quantiles = np.quantile(x, [0.1, 0.5, 0.9], method="linear")
    return dict(n=len(x), mean_pp=float(x.mean()), mean_absolute_pp=float(abs(x).mean()),
                median_absolute_pp=float(np.median(abs(x))), minimum_pp=float(x.min()),
                maximum_pp=float(x.max()), p10_pp=float(quantiles[0]), median_pp=float(quantiles[1]),
                p90_pp=float(quantiles[2]), p90_minus_p10_pp=float(quantiles[2] - quantiles[0]),
                positive=int((x > TOL).sum()), negative=int((x < -TOL).sum()),
                numerically_zero=int((abs(x) <= TOL).sum()))


def attach_set(frame, name):
    return frame.assign(block_set=name)


def combine_3_12(old, new, keys, values):
    require(not old.duplicated(keys).any() and not new.duplicated(keys).any(), "Duplicate combination keys")
    joined = old[keys + values].merge(new[keys + values], on=keys, how="outer", suffixes=("_old", "_new"),
                                     validate="one_to_one", indicator=True)
    require(joined._merge.eq("both").all(), "Old/new combination keys differ")
    output = joined[keys].copy()
    for name in values:
        output[name] = (3 * joined[name + "_old"] + 12 * joined[name + "_new"]) / 15
    return output


def make_budget_contrast(frame, keys, values):
    low = frame[frame.budget == 60][keys + values]
    high = frame[frame.budget == 300][keys + values]
    joined = low.merge(high, on=keys, how="outer", suffixes=("_60", "_300"), validate="one_to_one",
                       indicator=True)
    require(joined._merge.eq("both").all(), "Incomplete paired budget contrast")
    output = joined[keys].copy()
    for name in values:
        output[name.replace("_penalty_pp", "_budget_contrast_pp")] = joined[name + "_300"] - joined[name + "_60"]
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Completed new run.json")
    parser.add_argument("original_review", type=Path, help="Frozen original external review.json")
    args = parser.parse_args()
    run_path = args.run.resolve()
    original_path = args.original_review.resolve()
    run = json.loads(run_path.read_bytes())
    original_review = json.loads(original_path.read_bytes())
    require(run["status"] == "completed endpoint extension; descriptive review pending", "Run incomplete")
    require(original_review["status"] == "passed complete external descriptive audit", "Original audit incomplete")
    require(len(run["folds"]) == 14 and not run["failures"] and run["prediction_rows"] == 645120,
            "Run count or failures invalid")
    initial_hashes = {rel(run_path): sha(run_path), rel(original_path): sha(original_path),
                      rel(Path(__file__)): sha(__file__)}
    omitted_runtime = []
    for item in run["artifacts"] + run["sources"] + run["verified_input_artifacts"]:
        if item["path"].startswith((".tools/R-library/MuSiC/", "data/raw/MuSiC_runtime/")):
            omitted_runtime.append(item)
            continue
        require(sha(ROOT / item["path"]) == item["sha256"], f"Frozen file hash mismatch: {item['path']}")
    input_path = ROOT / run["input_manifest"]
    require(sha(input_path) == run["input_manifest_sha256"], "Input manifest hash mismatch")
    inp = json.loads(input_path.read_bytes())
    require(original_review["input_manifest_sha256"] == inp["source_input_sha256"],
            "Original review and new input are not from the same fixed source input")
    data_dir = ROOT / inp["data_directory"]
    cases = pd.read_csv(ROOT / inp["case_file"], sep="\t", keep_default_na=False).set_index("case_id")
    targets = pd.read_csv(data_dir / "targets.tsv", sep="\t", keep_default_na=False).set_index("target_name")
    refs = pd.read_csv(ROOT / inp["reference_file"], sep="\t", keep_default_na=False)
    require(cases.index.is_unique and len(cases) == 5376 and len(refs) == 1344, "Expected case/reference grid")
    require(len(targets) == 840 and targets.index.is_unique, "Expected original fixed targets")
    require(set(cases.block) == set(range(3, 15)) and set(cases.budget) == set(BUDGETS)
            and set(cases.level) == set(LEVELS), "Unexpected extension design")
    require(set(cases.held_out) == set(inp["donors"]) and len(inp["donors"]) == 14, "Donor set mismatch")
    require(set(cases.triple_key) == {fold["triple_key"] for fold in run["folds"]}, "Triple set mismatch")
    artifact_hashes = {item["path"]: item["sha256"] for item in original_review["artifacts"]}

    def load_original(name):
        path = original_path.parent / (name + ".csv")
        require(rel(path) in artifact_hashes and sha(path) == artifact_hashes[rel(path)],
                f"Original table hash mismatch: {name}")
        initial_hashes[rel(path)] = sha(path)
        frame = pd.read_csv(path, keep_default_na=False)
        selected = frame[frame.budget.isin(BUDGETS)].copy()
        if name != "donor_metrics":
            selected = selected[selected.level == "ratio10"].copy()
        else:
            selected = selected[selected.level.isin(LEVELS)].copy()
        return selected

    old_arrangement = load_original("reference_arrangement_penalties")
    old_block = load_original("donor_block_penalties")
    old_donor_metrics = load_original("donor_metrics")
    require(set(old_block.block) == {0, 1, 2}, "Original blocks changed")
    require(len(old_arrangement) == 672 and len(old_block) == 168 and len(old_donor_metrics) == 112,
            "Original endpoint subset counts changed")
    out = ROOT / "results/music_mc_extension_review" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True, exist_ok=False)
    source_dir = out / "source"
    source_dir.mkdir()
    shutil.copyfile(__file__, source_dir / Path(__file__).name)
    record = dict(status="running complete endpoint audit", run=rel(run_path), run_sha256=sha(run_path),
                  input_manifest=rel(input_path), input_manifest_sha256=sha(input_path),
                  original_review=rel(original_path), original_review_sha256=sha(original_path),
                  producer=dict(path=rel(__file__), sha256=sha(__file__)),
                  original_tables=[dict(path=path, sha256=digest) for path, digest in initial_hashes.items()
                                   if path.endswith(".csv")],
                  started_utc=datetime.now(timezone.utc).isoformat())
    save(out / "started.json", record)
    print("OUTPUT", rel(out), flush=True)
    tables = {}

    def export(name, frame):
        path = out / (name + ".csv")
        frame.to_csv(path, index=False, lineterminator="\n")
        tables[name] = len(frame)
        return frame

    try:
        chunks = []
        diagnostics_counts = Counter()
        validated_rows = validated_diagnostics = 0
        max_metric_error = max_sum_error = max_truth_error = 0.0
        variance_false = maxiter = 0
        feature_min, feature_max = np.inf, -np.inf
        per_triple = []
        meta = ["reference_id", "held_out", "triple_id", "triple_key", "block", "budget", "level", "dominant_donor"]
        numeric_truth = [f"true_{k}" for k in range(6)]
        numeric_pred = [f"pred_{k}" for k in range(6)]
        for fold in sorted(run["folds"], key=lambda item: item["triple_key"]):
            triple = fold["triple_key"]
            folder = run_path.parent / triple
            frame = pd.read_csv(folder / "predictions.csv", keep_default_na=False)
            diag = pd.read_csv(folder / "diagnostics.csv", keep_default_na=False)
            local_cases = cases[cases.triple_key == triple]
            key = ["case_id", "method", "mixture_id"]
            require(len(frame) == 46080 and len(diag) == 23040 and len(local_cases) == 384,
                    f"Wrong row counts: {triple}")
            require(not frame.duplicated(key).any() and set(frame.case_id) == set(local_cases.index),
                    f"Prediction keys invalid: {triple}")
            require(set(frame.method) == set(METHODS) and set(frame.mixture_id) == set(range(60))
                    and frame.groupby(["case_id", "method"]).size().eq(60).all(), "Missing case targets")
            expected = local_cases.loc[frame.case_id, meta].reset_index(drop=True)
            require(frame[meta].equals(expected), f"Prediction metadata mismatch: {triple}")
            truth = frame[numeric_truth].to_numpy(dtype=float)
            pred = frame[numeric_pred].to_numpy(dtype=float)
            require(np.isfinite(pred).all() and np.isfinite(truth).all()
                    and pred.min() >= -TOL and pred.max() <= 1 + TOL, "Invalid proportion values")
            sum_error = float(np.max(abs(pred.sum(axis=1) - 1)))
            require(sum_error < TOL and np.max(abs(truth.sum(axis=1) - 1)) < TOL, "Proportions do not sum to one")
            truth_error = float(np.max(abs(truth - targets.loc[frame.target_name, numeric_truth].to_numpy())))
            require(truth_error < 1e-14, "Target truth changed")
            for name in ("held_out", "mixture_id"):
                require(np.array_equal(frame[name].to_numpy(), targets.loc[frame.target_name, name].to_numpy()),
                        f"Target metadata changed: {name}")
            errors = pred - truth
            recalculated_mae = abs(errors).mean(axis=1)
            recalculated_rmse = np.sqrt((errors**2).mean(axis=1))
            require(np.isfinite(frame[["mae", "rmse"]].to_numpy()).all(), "Nonfinite saved metric")
            metric_error = max(float(np.max(abs(recalculated_mae - frame.mae))),
                               float(np.max(abs(recalculated_rmse - frame.rmse))))
            require(metric_error < 1e-12, "Saved MAE/RMSE do not match proportions")
            dkeys = ["case_id", "mixture_id"]
            weighted = frame[frame.method == "music_weighted"]
            require(not diag.duplicated(dkeys).any() and diag[dkeys].sort_values(dkeys).reset_index(drop=True).equals(
                weighted[dkeys].sort_values(dkeys).reset_index(drop=True)), "Diagnostic keys do not cover every weighted fit")
            for name in ("reference_id", "held_out"):
                require(np.array_equal(diag[name].to_numpy(), cases.loc[diag.case_id, name].to_numpy()),
                        f"Diagnostic {name} mismatch")
            for name in ("held_out", "mixture_id"):
                require(np.array_equal(diag[name].to_numpy(), targets.loc[diag.target_name, name].to_numpy()),
                        f"Diagnostic target {name} mismatch")
            require(np.isfinite(diag[["n_features", "nnls_coefficient_sum", "weighted_coefficient_sum"]].to_numpy()).all(),
                    "Nonfinite diagnostic coefficients/support")
            require((diag.n_features > 0).all() and (diag.n_features <= inp["genes"]).all(), "Invalid gene support size")
            finite_flags = diag.variance_finite.astype(str).str.lower()
            require(finite_flags.isin(["true", "false"]).all(), "Malformed finite-variance flags")
            require(diag.convergence.astype(str).str.fullmatch(r"Converge at [0-9]+|Reach Maxiter").all(),
                    "Unrecognized convergence status")
            false_count = int(finite_flags.eq("false").sum())
            limit_count = int(diag.convergence.eq("Reach Maxiter").sum())
            require(false_count == fold["nonfinite_variance_count"] and limit_count == fold["maxiter_count"],
                    "Diagnostic counts differ from fold manifest")
            require(int(diag.n_features.min()) == fold["effective_gene_range"][0]
                    and int(diag.n_features.max()) == fold["effective_gene_range"][1], "Support range mismatch")
            diagnostics_counts.update(diag.convergence.value_counts().to_dict())
            variance_false += false_count
            maxiter += limit_count
            feature_min = min(feature_min, diag.n_features.min())
            feature_max = max(feature_max, diag.n_features.max())
            frame["mae_pp"] = recalculated_mae * 100
            frame["rmse_pp"] = recalculated_rmse * 100
            for k in range(6):
                frame[f"absolute_error_{k}_pp"] = abs(errors[:, k]) * 100
            reduced = frame.groupby(["case_id"] + meta + ["method"], sort=True)[METRICS].mean().reset_index()
            require(len(reduced) == 768, "Invalid case aggregation count")
            chunks.append(reduced)
            validated_rows += len(frame)
            validated_diagnostics += len(diag)
            max_metric_error = max(max_metric_error, metric_error)
            max_sum_error = max(max_sum_error, sum_error)
            max_truth_error = max(max_truth_error, truth_error)
            per_triple.append(dict(triple_key=triple, prediction_rows=len(frame), weighted_fits=len(diag),
                                   metric_maximum_difference=metric_error, maxiter_count=limit_count,
                                   nonfinite_variance_count=false_count))
            print("VALIDATED", len(per_triple), "/14", triple, flush=True)
        require(validated_rows == 645120 and validated_diagnostics == 322560, "Incomplete full row audit")
        cm = pd.concat(chunks, ignore_index=True)
        require(len(cm) == 10752 and not cm.duplicated(["case_id", "method"]).any(), "Duplicate case means")
        export("new_case_metrics", cm)
        pair_keys = ["held_out", "triple_id", "triple_key", "block", "budget", "method"]
        baseline = cm[cm.level == "balanced"][pair_keys + METRICS].rename(columns={name: "baseline_" + name for name in METRICS})
        require(not baseline.duplicated(pair_keys).any(), "Duplicated shared baseline")
        paired = cm[cm.level == "ratio10"].merge(baseline, on=pair_keys, validate="many_to_one", how="left")
        require(len(paired) == 8064 and paired[["baseline_" + name for name in METRICS]].notna().all().all(),
                "Missing paired baseline")
        for metric, penalty in zip(METRICS, PENALTIES):
            paired[penalty] = paired[metric] - paired["baseline_" + metric]
        arrangement_keys = ["held_out", "triple_id", "dominant_donor", "method", "budget", "level"]
        export("new_arrangement_block_penalties", paired[pair_keys + ["dominant_donor", "level"] + PENALTIES])
        new_arrangement = paired.groupby(arrangement_keys)[PENALTIES].mean().reset_index()
        require(len(new_arrangement) == 672, "New arrangement count mismatch")
        block_keys = ["held_out", "block", "method", "budget", "level"]
        new_block = paired.groupby(block_keys)[PENALTIES].mean().reset_index()
        require(len(new_block) == 672, "New donor-block count mismatch")
        all_block = pd.concat([old_block[block_keys + PENALTIES], new_block], ignore_index=True)
        require(len(all_block) == 840 and not all_block.duplicated(block_keys).any(), "Combined block keys invalid")
        export("donor_block_penalties_all_15", all_block.sort_values(block_keys))
        donor_keys = ["held_out", "method", "budget", "level"]
        new_donor_metrics = cm.groupby(donor_keys)[METRICS].mean().reset_index()
        export("new_donor_block_metrics", cm.groupby(["held_out", "block", "method", "budget", "level"])[METRICS].mean().reset_index())
        all_donor_metrics = pd.concat([
            attach_set(old_donor_metrics[donor_keys + METRICS], "original_3"),
            attach_set(new_donor_metrics, "new_12"),
            attach_set(combine_3_12(old_donor_metrics, new_donor_metrics, donor_keys, METRICS), "combined_15")], ignore_index=True)
        export("donor_metrics_by_block_set", all_donor_metrics)
        export("aggregate_metrics_by_block_set", all_donor_metrics.groupby(["block_set", "method", "budget", "level"])[METRICS].mean().reset_index())
        combined_arrangement = combine_3_12(old_arrangement, new_arrangement, arrangement_keys, PENALTIES)
        all_arrangement = pd.concat([attach_set(old_arrangement[arrangement_keys + PENALTIES], "original_3"),
                                     attach_set(new_arrangement, "new_12"),
                                     attach_set(combined_arrangement, "combined_15")], ignore_index=True)
        incidence = cases.reset_index()[["held_out", "triple_id", "triple_key"]].drop_duplicates()
        all_arrangement = all_arrangement.merge(incidence, on=["held_out", "triple_id"], validate="many_to_one")
        export("reference_arrangement_penalties_by_block_set", all_arrangement)
        donor_penalties = []
        global_blocks = all_block.groupby(["block", "method", "budget"])[PENALTIES].mean().reset_index()
        export("global_block_allocation_contrasts", global_blocks)
        global_budget = make_budget_contrast(global_blocks, ["block", "method"], PENALTIES)
        export("global_block_budget_contrasts", global_budget)
        donor_budget = make_budget_contrast(all_block, ["held_out", "block", "method", "level"], PENALTIES)
        export("donor_block_budget_contrasts", donor_budget)
        precision = []
        trajectory = []
        directions = []
        max_identity = 0.0
        sets = [("original_3", list(range(3))), ("new_12", list(range(3, 15))), ("combined_15", list(range(15)))]
        for block_set, blocks in sets:
            donor = all_block[all_block.block.isin(blocks)].groupby(donor_keys)[PENALTIES].mean().reset_index()
            donor_penalties.append(attach_set(donor, block_set))
            absolutes = all_donor_metrics[all_donor_metrics.block_set == block_set]
            recalculated = absolutes.pivot(index=["held_out", "method", "budget"], columns="level", values="mae_pp")
            difference = recalculated.ratio10 - recalculated.balanced
            observed = donor.set_index(["held_out", "method", "budget"]).mae_penalty_pp.reindex(difference.index)
            max_identity = max(max_identity, float(np.max(abs(difference - observed))))
            arr_mean = all_arrangement[all_arrangement.block_set == block_set].groupby(donor_keys)[PENALTIES].mean()
            block_mean = donor.set_index(donor_keys)[PENALTIES].reindex(arr_mean.index)
            max_identity = max(max_identity, float(np.max(abs(arr_mean.to_numpy() - block_mean.to_numpy()))))
            for method in METHODS:
                for metric, penalty in zip(METRICS, PENALTIES):
                    budget_name = penalty.replace("_penalty_pp", "_budget_contrast_pp")
                    quantities = [
                        ("I60", global_blocks[(global_blocks.block.isin(blocks)) & (global_blocks.method == method)
                                              & (global_blocks.budget == 60)][penalty]),
                        ("I300", global_blocks[(global_blocks.block.isin(blocks)) & (global_blocks.method == method)
                                               & (global_blocks.budget == 300)][penalty]),
                        ("C300_minus_60", global_budget[(global_budget.block.isin(blocks)) & (global_budget.method == method)][budget_name])]
                    for quantity, values in quantities:
                        x = values.to_numpy()
                        require(len(x) == len(blocks), "Global-block sampling unit incomplete")
                        sd = float(x.std(ddof=1))
                        mcse = float(sd / np.sqrt(len(x)))
                        eligible = method == "music_weighted" and metric == "mae_pp"
                        precision.append(dict(block_set=block_set, method=method, metric=metric, estimand=quantity,
                                              n_blocks=len(x), mean_pp=float(x.mean()), block_sd_pp=sd,
                                              conditional_mcse_pp=mcse,
                                              precision_target_pp=0.05 if eligible else "",
                                              meets_precision_target=(mcse <= 0.05) if eligible else ""))
                        directions.append(dict(scope="global_block", block_set=block_set, method=method,
                                               metric=metric, estimand=quantity, **direction(x)))
        require(max_identity < 1e-10, "Paired/absolute/arrangement averaging identities disagree")
        donor_penalties = pd.concat(donor_penalties, ignore_index=True)
        export("donor_penalties_by_block_set", donor_penalties)
        export("aggregate_penalties_by_block_set", donor_penalties.groupby(["block_set", "method", "budget", "level"])[PENALTIES].mean().reset_index())
        donor_changes = make_budget_contrast(donor_penalties, ["block_set", "held_out", "method", "level"], PENALTIES)
        export("donor_budget_contrasts_by_block_set", donor_changes)
        arrangement_changes = make_budget_contrast(all_arrangement,
            ["block_set", "held_out", "triple_id", "triple_key", "dominant_donor", "method", "level"], PENALTIES)
        export("arrangement_budget_contrasts_by_block_set", arrangement_changes)
        distribution = []
        for (block_set, method, budget), group in all_arrangement.groupby(["block_set", "method", "budget"]):
            require(len(group) == 168, "Arrangement distribution must contain 168 averages")
            for penalty in PENALTIES:
                distribution.append(dict(block_set=block_set, method=method, budget=int(budget), metric=penalty,
                                         averaging="60 targets then specified blocks within arrangement; abs applied afterwards",
                                         **direction(group[penalty])))
        export("arrangement_effect_distributions", pd.DataFrame(distribution))
        for scope, frame, groups, values in [
            ("donor_allocation", donor_penalties, ["block_set", "method", "budget"], PENALTIES),
            ("donor_budget", donor_changes, ["block_set", "method"], [x.replace("_penalty_pp", "_budget_contrast_pp") for x in PENALTIES]),
            ("arrangement_budget", arrangement_changes, ["block_set", "method"], [x.replace("_penalty_pp", "_budget_contrast_pp") for x in PENALTIES]),
        ]:
            for identity, group in frame.groupby(groups):
                info = dict(zip(groups, identity))
                for value in values:
                    directions.append(dict(scope=scope, **info, metric=value, **direction(group[value])))
        for n in (3, 6, 9, 12):
            blocks = list(range(3, 3 + n))
            for method in METHODS:
                for quantity, frame, col in [
                    ("I60", global_blocks[global_blocks.budget == 60], "mae_penalty_pp"),
                    ("I300", global_blocks[global_blocks.budget == 300], "mae_penalty_pp"),
                    ("C300_minus_60", global_budget, "mae_budget_contrast_pp")]:
                    x = frame[(frame.method == method) & frame.block.isin(blocks)][col].to_numpy()
                    require(len(x) == n, "Trajectory is missing a global block")
                    trajectory.append(dict(method=method, estimand=quantity, new_blocks=n,
                                           first_block=3, last_block=3 + n - 1, mean_pp=float(x.mean()),
                                           conditional_mcse_pp=float(x.std(ddof=1) / np.sqrt(n)),
                                           use="diagnostic only; no stopping or seed selection"))
        absolute_method = all_donor_metrics.pivot(index=["block_set", "held_out", "budget", "level"], columns="method", values="mae_pp").reset_index()
        absolute_method["weighted_minus_nnls_mae_pp"] = absolute_method.music_weighted - absolute_method.music_nnls
        export("donor_method_absolute_differences", absolute_method)
        penalty_method = donor_penalties.pivot(index=["block_set", "held_out", "budget", "level"], columns="method", values="mae_penalty_pp").reset_index()
        penalty_method["weighted_minus_nnls_penalty_pp"] = penalty_method.music_weighted - penalty_method.music_nnls
        export("donor_method_penalty_differences", penalty_method)
        for scope, frame, value in [("method_absolute", absolute_method, "weighted_minus_nnls_mae_pp"),
                                    ("method_allocation", penalty_method, "weighted_minus_nnls_penalty_pp")]:
            for identity, group in frame.groupby(["block_set", "budget", "level"]):
                directions.append(dict(scope=scope, **dict(zip(["block_set", "budget", "level"], identity)),
                                       metric=value, **direction(group[value])))
        precision_frame = export("conditional_global_block_precision", pd.DataFrame(precision))
        export("new_block_precision_trajectory", pd.DataFrame(trajectory))
        export("complete_direction_summaries", pd.DataFrame(directions).fillna(""))
        export("per_triple_checks", pd.DataFrame(per_triple))
        export("cell_type_metric_key", pd.DataFrame(dict(type_index=range(6), cell_type=inp["cell_types"])))
        for path, digest in initial_hashes.items():
            require(sha(ROOT / path) == digest, f"Source changed during review: {path}")
        target_rows = precision_frame[(precision_frame.block_set == "new_12")
                                      & (precision_frame.method == "music_weighted")
                                      & (precision_frame.metric == "mae_pp")]
        require(len(target_rows) == 3, "Three prespecified precision targets required")
        record.update(archived_runtime_provenance_not_revalidated=omitted_runtime,
                      saved_prediction_recalculation_only=True,
                      status="passed complete endpoint Monte Carlo audit",
                      prediction_rows=validated_rows, weighted_diagnostic_rows=validated_diagnostics,
                      logical_cases=5376, unique_references=1344,
                      maximum_metric_recalculation_difference=max_metric_error,
                      maximum_proportion_sum_error=max_sum_error, maximum_truth_difference=max_truth_error,
                      maximum_aggregation_identity_difference_pp=max_identity,
                      weighted_diagnostics=dict(convergence_counts={str(k): int(v) for k, v in diagnostics_counts.items()},
                                                maxiter_count=maxiter, nonfinite_variance_count=variance_false,
                                                effective_gene_minimum=int(feature_min), effective_gene_maximum=int(feature_max)),
                      prespecified_new_weighted_precision=target_rows.to_dict("records"),
                      all_three_precision_targets_met=bool(target_rows.meets_precision_target.all()),
                      scope="Conditional reference-cell randomization at fixed finite inventory, donors, original grouping and targets. Whole global block is the MC unit; no donor-population P values or confidence intervals. Original 3, new 12 and combined 15 are reported separately; all 12 new blocks retained regardless of precision or direction.")
    except BaseException as error:
        record.update(status="failed endpoint review; partial outputs retained", error=repr(error))
        raise
    finally:
        record.update(finished_utc=datetime.now(timezone.utc).isoformat(), table_rows=tables,
                      artifacts=[dict(path=rel(path), sha256=sha(path), bytes=path.stat().st_size)
                                 for path in sorted(out.rglob("*")) if path.is_file() and path.name != "review.json"])
        save(out / "review.json", record)
        print(record["status"], flush=True)
        if "prespecified_new_weighted_precision" in record:
            print(json.dumps(record["prespecified_new_weighted_precision"], indent=2))


if __name__ == "__main__":
    main()
