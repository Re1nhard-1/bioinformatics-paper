"""Reconstruct selected MC summaries using only standard-library arithmetic."""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("music_weighted", "music_nnls")
BUDGETS = (60, 300)
SIGN_TOL = 1e-10
COMPARE_TOL = 1e-10


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for data in iter(lambda: handle.read(8 * 2**20), b""):
            digest.update(data)
    return digest.hexdigest()


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def linear_quantile(values, probability):
    values = sorted(values)
    position = probability * (len(values) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    return values[low] + (position - low) * (values[high] - values[low])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review", type=Path)
    args = parser.parse_args()
    review_path = args.review.resolve()
    review = json.loads(review_path.read_bytes())
    require(review["status"] == "passed complete endpoint Monte Carlo audit", "Primary review is incomplete")
    original_path = ROOT / review["original_review"]
    require(sha(original_path) == review["original_review_sha256"], "Original audit hash differs")
    original = json.loads(original_path.read_bytes())
    recorded = {entry["path"]: entry["sha256"] for entry in review["artifacts"]}
    original_recorded = {entry["path"]: entry["sha256"] for entry in original["artifacts"]}
    source_hashes = {rel(review_path): sha(review_path), rel(original_path): sha(original_path),
                     rel(Path(__file__)): sha(__file__)}

    def load(directory, name, hashes):
        path = directory / (name + ".csv")
        require(rel(path) in hashes and sha(path) == hashes[rel(path)], f"Table hash mismatch: {name}")
        source_hashes[rel(path)] = sha(path)
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    new = load(review_path.parent, "new_case_metrics", recorded)
    old_block = load(original_path.parent, "donor_block_penalties", original_recorded)
    old_arrangements = load(original_path.parent, "reference_arrangement_penalties", original_recorded)
    reported_precision = load(review_path.parent, "conditional_global_block_precision", recorded)
    reported_blocks = load(review_path.parent, "global_block_allocation_contrasts", recorded)
    reported_changes = load(review_path.parent, "global_block_budget_contrasts", recorded)
    reported_arrangements = load(review_path.parent, "reference_arrangement_penalties_by_block_set", recorded)
    reported_distributions = load(review_path.parent, "arrangement_effect_distributions", recorded)
    reported_trajectory = load(review_path.parent, "new_block_precision_trajectory", recorded)
    require(len(new) == 10752, "Expected all 10752 new case means")
    seen = set()
    balanced = {}
    direct_values = defaultdict(list)
    for row in new:
        unique = row["case_id"], row["method"]
        require(unique not in seen, "Duplicate case mean")
        seen.add(unique)
        block, budget = int(row["block"]), int(row["budget"])
        require(block in range(3, 15) and budget in BUDGETS and row["method"] in METHODS,
                "Unexpected case-mean design")
        value = float(row["mae_pp"])
        require(math.isfinite(value), "Nonfinite case mean")
        direct_values[(block, row["method"], budget, row["level"])].append(value)
        if row["level"] == "balanced":
            key = row["held_out"], int(row["triple_id"]), block, budget, row["method"]
            require(key not in balanced, "Duplicate balanced baseline")
            balanced[key] = value
    require(len(balanced) == 2688, "Expected 2688 shared balanced means")
    new_arrangement_values = defaultdict(list)
    paired_block_values = defaultdict(list)
    for row in new:
        if row["level"] == "balanced":
            continue
        require(row["level"] == "ratio10", "Unexpected unequal allocation")
        block, budget = int(row["block"]), int(row["budget"])
        base_key = row["held_out"], int(row["triple_id"]), block, budget, row["method"]
        require(base_key in balanced, "Missing shared baseline")
        contrast = float(row["mae_pp"]) - balanced[base_key]
        arrangement = row["held_out"], int(row["triple_id"]), row["dominant_donor"], row["method"], budget
        new_arrangement_values[arrangement].append(contrast)
        paired_block_values[(block, row["method"], budget)].append(contrast)
    require(len(new_arrangement_values) == 672 and all(len(x) == 12 for x in new_arrangement_values.values()),
            "Incomplete new arrangements")
    block_estimates = {}
    for key, values in paired_block_values.items():
        require(len(values) == 168, "Expected 168 paired arrangements per global block/method/budget")
        block_estimates[key] = statistics.fmean(values)
        unequal = direct_values[key + ("ratio10",)]
        equal = direct_values[key + ("balanced",)]
        require(len(unequal) == 168 and len(equal) == 56, "Unbalanced case-count structure")
        require(abs(block_estimates[key] - (statistics.fmean(unequal) - statistics.fmean(equal))) < COMPARE_TOL,
                "Baseline averaging identity failed")
    old_values = defaultdict(list)
    for row in old_block:
        budget = int(row["budget"])
        if budget in BUDGETS and row["level"] == "ratio10":
            key = int(row["block"]), row["method"], budget
            old_values[key].append(float(row["mae_penalty_pp"]))
    require(len(old_values) == 12 and all(len(x) == 14 for x in old_values.values()), "Expected original 3 block means")
    for key, values in old_values.items():
        require(key not in block_estimates, "Old and new blocks overlap")
        block_estimates[key] = statistics.fmean(values)
    require(len(block_estimates) == 60, "Expected 15 blocks x 2 budgets x 2 methods")
    comparisons = []

    def compare(name, expected, observed):
        error = abs(expected - float(observed))
        require(math.isfinite(error) and error < COMPARE_TOL, f"{name}: {expected} versus {observed}")
        comparisons.append((name, error))

    require(len(reported_blocks) == 60 and len(reported_changes) == 30, "Reported global-block table incomplete")
    for row in reported_blocks:
        key = int(row["block"]), row["method"], int(row["budget"])
        compare("global block I", block_estimates[key], row["mae_penalty_pp"])
    for row in reported_changes:
        block, method = int(row["block"]), row["method"]
        compare("global block C", block_estimates[block, method, 300] - block_estimates[block, method, 60],
                row["mae_budget_contrast_pp"])
    block_sets = {"original_3": list(range(3)), "new_12": list(range(3, 15)), "combined_15": list(range(15))}

    def values_for(blocks, method, estimand):
        if estimand == "I60":
            return [block_estimates[b, method, 60] for b in blocks]
        if estimand == "I300":
            return [block_estimates[b, method, 300] for b in blocks]
        require(estimand == "C300_minus_60", "Unknown estimand")
        return [block_estimates[b, method, 300] - block_estimates[b, method, 60] for b in blocks]

    checks = []
    primary_rows = [row for row in reported_precision if row["metric"] == "mae_pp"]
    require(len(primary_rows) == 18, "Expected three estimands x two methods x three block sets")
    for row in primary_rows:
        blocks = block_sets[row["block_set"]]
        values = values_for(blocks, row["method"], row["estimand"])
        mean = statistics.fmean(values)
        sd = statistics.stdev(values)
        mcse = sd / math.sqrt(len(values))
        require(int(row["n_blocks"]) == len(values), "MCSE denominator is not global blocks")
        compare("global block mean", mean, row["mean_pp"])
        compare("global block SD", sd, row["block_sd_pp"])
        compare("global block MCSE", mcse, row["conditional_mcse_pp"])
        if row["method"] == "music_weighted":
            require(float(row["precision_target_pp"]) == .05, "Precision target changed")
            require(row["meets_precision_target"].lower() == str(mcse <= .05).lower(), "Precision label incorrect")
        if row["block_set"] == "new_12" and row["method"] == "music_weighted":
            checks.append(dict(estimand=row["estimand"], mean_pp=mean, conditional_mcse_pp=mcse,
                               meets_precision_target=mcse <= .05))
    for row in reported_trajectory:
        n = int(row["new_blocks"])
        require(n in (3, 6, 9, 12), "Unexpected running diagnostic checkpoint")
        values = values_for(range(3, 3 + n), row["method"], row["estimand"])
        compare("trajectory mean", statistics.fmean(values), row["mean_pp"])
        compare("trajectory MCSE", statistics.stdev(values) / math.sqrt(n), row["conditional_mcse_pp"])
    old_arr = {}
    for row in old_arrangements:
        budget = int(row["budget"])
        if budget in BUDGETS and row["level"] == "ratio10":
            key = row["held_out"], int(row["triple_id"]), row["dominant_donor"], row["method"], budget
            require(key not in old_arr, "Duplicate old arrangement")
            old_arr[key] = float(row["mae_penalty_pp"])
    require(set(old_arr) == set(new_arrangement_values), "Old/new arrangement identities differ")
    arrangements = {}
    for key, values in new_arrangement_values.items():
        arrangements["original_3", key] = old_arr[key]
        arrangements["new_12", key] = statistics.fmean(values)
        arrangements["combined_15", key] = (old_arr[key] * 3 + math.fsum(values)) / 15
    require(len(reported_arrangements) == len(arrangements) == 2016, "Arrangement table count mismatch")
    seen_arrangements = set()
    for row in reported_arrangements:
        key = (row["block_set"], (row["held_out"], int(row["triple_id"]), row["dominant_donor"], row["method"], int(row["budget"])))
        require(key not in seen_arrangements, "Duplicate reported arrangement")
        seen_arrangements.add(key)
        compare("arrangement averaged signed contrast", arrangements[key], row["mae_penalty_pp"])
    distributions = defaultdict(list)
    for (block_set, key), value in arrangements.items():
        distributions[block_set, key[3], key[4]].append(value)
    checked_distribution_rows = 0
    for row in reported_distributions:
        if row["metric"] != "mae_penalty_pp":
            continue
        values = distributions[row["block_set"], row["method"], int(row["budget"])]
        require(len(values) == int(row["n"]) == 168, "Wrong arrangement distribution denominator")
        values_abs = [abs(x) for x in values]
        expected = dict(mean_pp=statistics.fmean(values), mean_absolute_pp=statistics.fmean(values_abs),
                        median_absolute_pp=statistics.median(values_abs), minimum_pp=min(values), maximum_pp=max(values),
                        p10_pp=linear_quantile(values, .1), median_pp=statistics.median(values),
                        p90_pp=linear_quantile(values, .9))
        expected["p90_minus_p10_pp"] = expected["p90_pp"] - expected["p10_pp"]
        for name, value in expected.items():
            compare("arrangement " + name, value, row[name])
        require(sum(x > SIGN_TOL for x in values) == int(row["positive"])
                and sum(x < -SIGN_TOL for x in values) == int(row["negative"])
                and sum(abs(x) <= SIGN_TOL for x in values) == int(row["numerically_zero"]), "Arrangement sign counts differ")
        checked_distribution_rows += 1
    require(checked_distribution_rows == 12, "Expected 12 MAE arrangement distributions")
    for path, digest in source_hashes.items():
        require(sha(ROOT / path) == digest, f"Source changed during check: {path}")
    out = ROOT / "results/music_mc_summary_check" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True, exist_ok=False)
    result = dict(status="passed standard-library MC summary check", review=rel(review_path),
                  review_sha256=sha(review_path), producer=dict(path=rel(__file__), sha256=sha(__file__)),
                  sources=[dict(path=path, sha256=digest) for path, digest in sorted(source_hashes.items())],
                  arithmetic_comparisons=len(comparisons), maximum_absolute_difference=max(error for _, error in comparisons),
                  global_block_allocation_values=60, global_budget_values=30, precision_rows=18,
                  arrangement_values=2016, arrangement_distribution_rows=12, primary_new_weighted=checks,
                  scope="Separate standard-library aggregation from new case means and hash-verified original CSVs. Checks MAE summaries, whole-block MCSE, trajectory, paired baseline weighting, draw-weighted combination, and arrangement magnitudes/quantiles/signs. Does not repeat fitting or recheck raw prediction columns; does not constitute independent human or biological review.")
    (out / "check.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                                    encoding="utf-8", newline="\n")
    print("OUTPUT", rel(out), flush=True)
    print(json.dumps(dict(status=result["status"], maximum_absolute_difference=result["maximum_absolute_difference"],
                          primary_new_weighted=checks), indent=2))


if __name__ == "__main__":
    main()
