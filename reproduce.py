"""Recalculate saved predictions; no fitting or network access."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
def execute(args):
    subprocess.run([sys.executable, "-X", "utf8"] + args, cwd=ROOT, check=True)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["original", "additional", "all"], default="all", nargs="?")
    a = p.parse_args()
    if a.stage in ("original", "all"):
        out = "recomputed/original_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        args = ["reproduction_integrated/code/summarize.py"]
        for stage in ["endpoints", "gradient", "thinning", "external", "alternative"]:
            args.extend(["--" + stage, "work/predictions/" + stage])
        execute(args + ["--output", out])
    if a.stage in ("additional", "all"):
        folder = ROOT / "results/music_mc_extension_review"
        before = set(folder.glob("*/review.json"))
        execute(["analysis/review_music_mc_archive.py", "results/music_mc_extension/20260918T155055175930Z/run.json", "results/music_external_review/20260917T133518483857Z/review.json"])
        new = set(folder.glob("*/review.json")) - before
        if len(new) != 1:
            raise ValueError("Expected exactly one new complete summary")
        actual = new.pop()
        execute(["analysis/check_music_mc_summary_v1.py", str(actual)])
        execute(["verify_recomputed_tables.py", str(actual.parent), "results/music_mc_extension_review/20260918T163322469416Z"])
    print("Saved-prediction recalculation passed; no MuSiC refitting was performed.")
if __name__ == "__main__":
    main()
