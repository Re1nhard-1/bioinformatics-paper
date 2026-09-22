"""Run the frozen D041 endpoint extension, with paired batching controls first."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import os
import shutil
import subprocess
import threading
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/MUSIC_MC_EXTENSION_PROTOCOL_V1.md"
ADAPTER = ROOT / "analysis/run_music_mc_extension_v1.R"
BASELINE = ROOT / "analysis/run_music_reference_composition.R"
RSCRIPT = Path("C:/Program Files/R/R-4.5.2/bin/Rscript.exe")
LIMIT_SECONDS = 7200
BLOCKS = set(range(3, 15))
ENV = dict(os.environ, LC_ALL="C", LANG="C", OPENBLAS_NUM_THREADS="1",
           OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n",
                          encoding="utf-8", newline="\n")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_adapter_derivation():
    old = BASELINE.read_text(encoding="utf-8")
    new = ADAPTER.read_text(encoding="utf-8")
    replacements = {
        "# Restricted-grid adapter for MUSIC_REFERENCE_COMPOSITION_PROTOCOL.md.\n"
        "# Derived from frozen run_music_external.R; official fitting code unchanged.\n"
        "# Full external adapter and batching controls. The official package is unchanged.":
        "# Endpoint Monte Carlo adapter for MUSIC_MC_EXTENSION_PROTOCOL_V1.md.\n"
        "# Derived from frozen run_music_reference_composition.R.\n"
        "# Official fitting calls, tracing and paired target batching are unchanged.",
        "refs<-refs[refs$triple_key==triple,,drop=FALSE];stopifnot(nrow(refs)==24L)":
        "refs<-refs[refs$triple_key==triple,,drop=FALSE];stopifnot(nrow(refs)==96L)",
        "cases<-cases[cases$triple_key==triple,,drop=FALSE];stopifnot(nrow(cases)==96L)":
        "cases<-cases[cases$triple_key==triple,,drop=FALSE];stopifnot(nrow(cases)==384L)",
        'cases$block=="0"': 'cases$block=="3"',
        "stopifnot(total==11520L)": "stopifnot(total==46080L)",
    }
    expected = old
    for before, after in replacements.items():
        require(expected.count(before) == 1, f"Unexpected frozen baseline: {before}")
        expected = expected.replace(before, after)
    require(new == expected, "R adapter differs beyond five allowed replacement sites")
    return "".join(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                      fromfile=rel(BASELINE), tofile=rel(ADAPTER)))


def validate_input(input_path):
    inp = json.loads(input_path.read_bytes())
    selection_path = ROOT / inp["selection_manifest"]
    sel = json.loads(selection_path.read_bytes())
    expected = {}

    def register(path, digest):
        key = rel(ROOT / path)
        if key in expected:
            require(expected[key] == digest, f"Conflicting manifest hashes: {key}")
        expected[key] = digest

    register(inp["selection_manifest"], inp["selection_manifest_sha256"])
    register(inp["source_input_manifest"], inp["source_input_sha256"])
    for item in inp["artifacts"] + sel["sources"] + sel["artifacts"] + inp["donor_matrices"]:
        register(item["path"], item["sha256"])
    for path, digest in expected.items():
        require(sha(ROOT / path) == digest, f"Input hash mismatch: {path}")
    required = [inp["reference_file"], inp["case_file"]]
    required += [str(Path(inp["data_directory"]) / name)
                 for name in ("cells.tsv", "genes.tsv", "targets.tsv", "targets.mtx.gz")]
    for path in required:
        require(rel(ROOT / path) in expected, f"Input missing from frozen hashes: {path}")
    refs = pd.read_csv(ROOT / inp["reference_file"], sep="\t", keep_default_na=False)
    cases = pd.read_csv(ROOT / inp["case_file"], sep="\t", keep_default_na=False)
    targets = pd.read_csv(ROOT / inp["data_directory"] / "targets.tsv", sep="\t",
                          keep_default_na=False)
    require(len(sel["donors"]) == len(set(sel["donors"])) == 14, "Expected 14 donors")
    require(len(inp["cell_types"]) == 6, "Expected six cell types")
    require(set(targets.held_out) == set(sel["donors"]), "Target donor set changed")
    require(len(targets) == 840 and targets.target_name.is_unique
            and targets.groupby("held_out").size().eq(60).all(), "Expected 60 targets/donor")
    require(set(refs.block) == set(cases.block) == BLOCKS, "Expected new blocks 3..14")
    require(set(refs.budget) == set(cases.budget) == {60, 300}, "Expected endpoint budgets")
    require(set(refs.level) == set(cases.level) == {"balanced", "ratio10"}, "Unexpected level")
    require(len(refs) == 1344 and refs.reference_id.is_unique, "Expected 1344 references")
    require(len(cases) == 5376 and cases.case_id.is_unique, "Expected 5376 logical cases")
    require(refs.triple_key.nunique() == 14 and refs.groupby("triple_key").size().eq(96).all(),
            "Expected 96 references/triple")
    require(cases.groupby("triple_key").size().eq(384).all(), "Expected 384 cases/triple")
    require(cases.reference_id.value_counts().eq(4).all()
            and set(cases.reference_id) == set(refs.reference_id), "Expected four cases/reference")
    require(cases.groupby("held_out").triple_id.nunique().eq(4).all(), "Four triples/fold required")
    by_ref = refs.set_index("reference_id")
    columns = ["triple_key", "block", "budget", "level", "dominant_donor"]
    for name in columns:
        require(np.array_equal(cases[name].to_numpy(),
                               by_ref.loc[cases.reference_id, name].to_numpy()),
                f"Case/reference disagreement: {name}")
    for _, group in refs.groupby(["triple_key", "block", "budget"]):
        donors = group.reference_donors.iloc[0].split("|")
        require(len(donors) == len(set(donors)) == 3
                and group.reference_donors.nunique() == 1, "Invalid reference donor set")
        unequal = group[group.level == "ratio10"]
        equal = group[group.level == "balanced"]
        require(len(group) == 4 and len(equal) == 1 and len(unequal) == 3,
                "Expected balanced and all three dominant rotations")
        require(equal.dominant_donor.iloc[0] == "" and set(unequal.dominant_donor) == set(donors),
                "Dominant rotations incomplete")
        require(equal.major_count.eq(equal.budget // 3).all()
                and equal.minor_count.eq(equal.budget // 3).all()
                and unequal.major_count.eq(5 * unequal.budget // 6).all()
                and unequal.minor_count.eq(unequal.budget // 12).all(), "Quota mismatch")
    for triple, group in cases.groupby("triple_key"):
        donors = set(refs.loc[refs.triple_key == triple, "reference_donors"].iloc[0].split("|"))
        require(group.held_out.nunique() == 4 and not donors.intersection(group.held_out),
                "Target/reference donor isolation failed")
    original = json.loads((ROOT / inp["source_input_manifest"]).read_bytes())
    old_cases = pd.read_csv(ROOT / original["case_file"], sep="\t", keep_default_na=False)
    identity = ["held_out", "triple_id", "triple_key"]
    require(cases[identity].drop_duplicates().sort_values(identity).reset_index(drop=True).equals(
        old_cases[identity].drop_duplicates().sort_values(identity).reset_index(drop=True)),
        "Reference grouping differs from the original external design")
    for name in ("targets.tsv", "targets.mtx.gz", "genes.tsv"):
        require(sha(ROOT / inp["data_directory"] / name)
                == sha(ROOT / original["data_directory"] / name), f"Fixed input changed: {name}")
    return inp, sel, selection_path, refs, cases, expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--workers", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate manifests/design/adapter without running R or creating a run")
    args = parser.parse_args()
    input_path = args.input.resolve()
    source_diff = check_adapter_derivation()
    inp, sel, selection_path, refs, cases, input_hashes = validate_input(input_path)
    official_manifest = ROOT / "environment/music/source_manifest.json"
    official = json.loads(official_manifest.read_bytes())
    require(official["commit"] == "f21fe67f5670d5e9fca0ad7550abaae3423eb59c",
            "Official MuSiC source commit changed")
    require(sha(official["archive"]) == official["sha256"], "Official archive hash mismatch")
    code = [Path(__file__), ADAPTER, PROTOCOL, BASELINE, ROOT / "environment/music/runtime.R"]
    code += sorted((Path(official["source_root"]) / "R").glob("*.R"))
    code += sorted(path for path in (ROOT / ".tools/R-library/MuSiC").rglob("*") if path.is_file())
    sources = {rel(path): sha(path) for path in code + [input_path, selection_path, official_manifest]}
    if args.validate_only:
        print(json.dumps(dict(status="inputs and R derivation validated; no fits executed",
                              references=1344, logical_cases=5376, prediction_rows=645120,
                              workers=args.workers, fitting_limit_seconds=LIMIT_SECONDS,
                              sources=sources, input_artifact_count=len(input_hashes)), indent=2))
        return
    out = ROOT / "results/music_mc_extension" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True, exist_ok=False)
    snapshot = out / "source"
    snapshot.mkdir()
    for path in code[:5]:
        shutil.copyfile(path, snapshot / path.name)
    (snapshot / "adapter_diff.patch").write_text(source_diff, encoding="utf-8", newline="\n")
    record = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
                  input_manifest=rel(input_path), input_manifest_sha256=sha(input_path),
                  sources=[dict(path=path, sha256=digest) for path, digest in sources.items()],
                  verified_input_artifacts=[dict(path=path, sha256=digest)
                                            for path, digest in sorted(input_hashes.items())],
                  planned_reference_matrices=1344, planned_logical_cases=5376,
                  planned_prediction_rows=645120, planned_weighted_fits=322560,
                  new_blocks=sorted(BLOCKS), workers=args.workers,
                  fitting_limit_seconds=LIMIT_SECONDS, rscript_sha256=sha(RSCRIPT),
                  control_fits_excluded_from_scientific_grid=True)
    save(out / "started.json", record)
    print("OUTPUT", rel(out), flush=True)
    start = time.monotonic()
    stop = threading.Event()
    folds, failures = [], []

    def run_r(folder, triple, mode, batch="pooled"):
        if stop.is_set():
            raise RuntimeError("Another worker failed; no further dispatch")
        folder.mkdir(parents=True, exist_ok=False)
        command = [str(RSCRIPT), "--vanilla", str(ADAPTER), str(input_path), str(folder),
                   triple, mode, batch]
        save(folder / "command.json", dict(command=command, cwd=str(ROOT),
                                           remaining_seconds=LIMIT_SECONDS - (time.monotonic() - start)))
        with (folder / "console.log").open("w", encoding="utf-8", newline="\n") as log:
            process = subprocess.Popen(command, cwd=ROOT, env=ENV, stdout=log, stderr=subprocess.STDOUT)
            try:
                while process.poll() is None:
                    if stop.is_set():
                        raise RuntimeError("Another worker failed; active process stopped")
                    if time.monotonic() - start >= LIMIT_SECONDS:
                        raise TimeoutError("7200-second fitting deadline exceeded")
                    time.sleep(0.25)
                if process.returncode:
                    raise subprocess.CalledProcessError(process.returncode, command)
            finally:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()

    try:
        first = cases.loc[(cases.held_out == sel["donors"][0]) & (cases.triple_id == 0),
                          "triple_key"].unique()
        require(len(first) == 1, "Control triple not unique")
        run_r(out / "controls", first[0], "controls")
        controls = json.loads((out / "controls/controls.json").read_bytes())
        require(controls["status"] == "controls completed" and len(controls["checks"]) == 2
                and controls["prediction_rows_per_path"] == 960,
                "Incomplete or malformed control metadata")
        keys = ["case_id", "method", "mixture_id"]
        numeric = ["mae", "rmse"] + [f"{prefix}_{i}" for prefix in ("true", "pred") for i in range(6)]
        pooled = pd.read_csv(out / "controls/pooled_predictions.csv", keep_default_na=False).sort_values(keys).reset_index(drop=True)
        separate = pd.read_csv(out / "controls/separate_predictions.csv", keep_default_na=False).sort_values(keys).reset_index(drop=True)
        require(len(pooled) == len(separate) == 960 and not pooled.duplicated(keys).any()
                and not separate.duplicated(keys).any() and pooled[keys].equals(separate[keys]),
                "Control export keys or counts differ")
        delta = float(np.max(abs(pooled[numeric].to_numpy() - separate[numeric].to_numpy())))
        require(np.isfinite(delta), "Nonfinite batching export comparison")
        equivalent = bool(controls["batch_equivalence_passed"] and delta <= 1e-10)
        batch = "pooled" if equivalent else "separate"
        record.update(batch_mode=batch, batch_controls=controls, batch_export_maximum_difference=delta,
                      control_fallback_used=not equivalent)
        save(out / "controls_verified.json", dict(record, passed=equivalent))
        print("CONTROLS", batch, "max difference", delta, flush=True)

        def worker(triple):
            try:
                run_r(out / triple, triple, "full", batch)
                fold = json.loads((out / triple / "fold.json").read_bytes())
                require(fold["status"] == "completed" and fold["prediction_rows"] == 46080
                        and fold["reference_matrices"] == 96 and fold["logical_cases"] == 384
                        and fold["weighted_fits"] == 23040, f"Incomplete triple: {triple}")
                return fold
            except Exception:
                stop.set()
                raise

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs = {pool.submit(worker, triple): triple for triple in sorted(refs.triple_key.unique())}
            for job in as_completed(jobs):
                try:
                    folds.append(job.result())
                    print("COMPLETED", len(folds), "/14", jobs[job], "elapsed",
                          round(time.monotonic() - start, 1), flush=True)
                except Exception as error:
                    failures.append(dict(triple_key=jobs[job], error=repr(error)))
                    stop.set()
        require(not failures and len(folds) == 14, "Incomplete extension; partial outputs retained")
        require(sum(fold["prediction_rows"] for fold in folds) == 645120, "Grid count mismatch")
        for path, digest in sources.items():
            require(sha(ROOT / path) == digest, f"Source changed during execution: {path}")
        record["status"] = "completed endpoint extension; descriptive review pending"
    except BaseException as error:
        stop.set()
        record.update(status="failed endpoint execution; partial outputs retained", error=repr(error))
        raise
    finally:
        record.update(folds=sorted(folds, key=lambda fold: fold["triple_key"]), failures=failures,
                      prediction_rows=sum(fold["prediction_rows"] for fold in folds),
                      elapsed_seconds=time.monotonic() - start,
                      finished_utc=datetime.now(timezone.utc).isoformat())
        record["artifacts"] = [dict(path=rel(path), sha256=sha(path), bytes=path.stat().st_size)
                               for path in sorted(out.rglob("*"))
                               if path.is_file() and path.name != "run.json"]
        save(out / "run.json", record)
        print(record["status"], "elapsed", round(record["elapsed_seconds"], 1), flush=True)


if __name__ == "__main__":
    main()
