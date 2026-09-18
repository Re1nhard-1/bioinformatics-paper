"""Derive equal-donor cell-type means from verified diagnostic profiles only."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METRICS = ["theta_variance_trace", "theta_pairwise_tv", "S_cv"]
QUOTAS = [5, 10, 20, 25, 40, 50, 100, 250]
TRAJECTORIES = {"minor": [5, 10, 25], "balanced": [20, 40, 100], "dominant": [50, 100, 250]}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review_directory", type=Path)
    args = parser.parse_args()
    directory = args.review_directory.resolve()
    review_path = directory / "review.json"
    review = json.loads(review_path.read_bytes())
    assert review["status"] == "passed bounded diagnostic numerical review; manual figure review pending"
    source = directory / "profile_dispersion_unique.csv"
    indexed = {item["path"]: item["sha256"] for item in review["artifacts"]}
    assert sha(source) == indexed[source.relative_to(ROOT).as_posix()]
    frame = pd.read_csv(source)
    assert len(frame) == 672 and not frame.duplicated(["donor", "cell_type", "reference_cells"]).any()
    assert frame.donor.nunique() == 14 and frame.cell_type.nunique() == 6
    assert set(frame.reference_cells) == set(QUOTAS)
    assert np.isfinite(frame[METRICS].to_numpy(float)).all()
    sizes = frame.groupby(["cell_type", "reference_cells"]).size()
    assert len(sizes) == 48 and sizes.eq(14).all()
    means = frame.groupby(["cell_type", "reference_cells"], observed=True, sort=True)[METRICS].mean().reset_index()
    means["equally_weighted_donors"] = 14
    rows = []
    for name, quotas in TRAJECTORIES.items():
        for metric in METRICS:
            table = means.pivot(index="cell_type", columns="reference_cells", values=metric)
            for cell_type, row in table.iterrows():
                low, middle, high = (float(row[q]) for q in quotas)
                rows.append({"cell_type": cell_type, "trajectory": name, "metric": metric,
                    "low_quota": quotas[0], "middle_quota": quotas[1], "high_quota": quotas[2],
                    "low_equal_donor_mean": low, "middle_equal_donor_mean": middle,
                    "high_equal_donor_mean": high, "endpoint_change_of_means": high-low,
                    "endpoint_ratio_of_means": high/low if low > 0 else np.nan,
                    "ratio_defined": low > 0, "endpoint_direction_of_means": "lower" if high < low else "higher" if high > low else "equal",
                    "first_adjacent_change_of_means": middle-low, "second_adjacent_change_of_means": high-middle,
                    "both_adjacent_means_nonincreasing": middle <= low and high <= middle,
                    "both_adjacent_means_strictly_decreasing": middle < low and high < middle,
                    "any_adjacent_mean_increase": middle > low or high > middle,
                    "equally_weighted_donors": 14})
    trajectories = pd.DataFrame(rows)
    assert len(means) == 48 and len(trajectories) == 54
    out = directory / "type_summaries"
    out.mkdir(exist_ok=False)
    shutil.copyfile(__file__, out / Path(__file__).name)
    means.to_csv(out/"type_equal_donor_means.csv", index=False, lineterminator="\n", float_format="%.17g")
    trajectories.to_csv(out/"type_mean_trajectories.csv", index=False, lineterminator="\n", float_format="%.17g")
    record = {"status": "passed equal-donor cell-type summary checks",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input": {"path": source.relative_to(ROOT).as_posix(), "sha256": sha(source), "rows": len(frame)},
        "frozen_review": {"path": review_path.relative_to(ROOT).as_posix(), "sha256": sha(review_path)},
        "producer": {"path": Path(__file__).resolve().relative_to(ROOT).as_posix(), "sha256": sha(__file__)},
        "groups": len(sizes), "donors_per_group_minimum": int(sizes.min()), "donors_per_group_maximum": int(sizes.max()),
        "cell_types": sorted(frame.cell_type.unique().tolist()), "quotas": QUOTAS,
        "type_mean_rows": len(means), "type_mean_trajectory_rows": len(trajectories),
        "estimand": "Within each cell type and quota, mean of the 14 donor-specific conditional diagnostic metrics. Trajectory ratios divide these equal-donor means; they are not averages of individual donor ratios.",
        "scope": "Descriptive summaries of three fixed draws; no population inference or causal attribution. No new fitting or figures; prior results and manifests unchanged.",
        "artifacts": [{"file": p.name, "sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(out.iterdir()) if p.is_file()]}
    (out/"verification.json").write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    assert sha(source) == record["input"]["sha256"] and sha(review_path) == record["frozen_review"]["sha256"]
    print(out.relative_to(ROOT).as_posix())
    print("48 type/quota means; 54 type trajectories; exactly 14 donors per mean; no fitting or figure changes")


if __name__ == "__main__":
    main()
