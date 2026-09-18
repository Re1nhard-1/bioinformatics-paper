"""Check every regenerated summary CSV against the archived full-precision table."""
import csv
import json
import math
from pathlib import Path
import sys

actual, expected = map(Path, sys.argv[1:3])
names = sorted(p.name for p in expected.glob("*.csv"))
if {p.name for p in actual.glob("*.csv")} != set(names):
    raise ValueError("Summary file set differs")
numeric = text = rows = 0
maximum = 0.0
for name in names:
    with (actual / name).open(encoding="utf-8", newline="") as f:
        a = list(csv.reader(f))
    with (expected / name).open(encoding="utf-8", newline="") as f:
        e = list(csv.reader(f))
    if len(a) != len(e) or a[0] != e[0]:
        raise ValueError("Summary shape differs: " + name)
    rows += len(a) - 1
    for i, (ar, er) in enumerate(zip(a[1:], e[1:]), 1):
        if len(ar) != len(er):
            raise ValueError("Summary row width differs")
        for x, y in zip(ar, er):
            try:
                xx, yy = float(x), float(y)
            except ValueError:
                if x != y:
                    raise ValueError("Summary text differs: " + name + ":" + str(i))
                text += 1
            else:
                if math.isnan(xx) and math.isnan(yy):
                    continue
                delta = abs(xx - yy)
                if not math.isfinite(delta) or delta > 1e-10:
                    raise ValueError("Summary number differs: " + name + ":" + str(i))
                maximum = max(maximum, delta)
                numeric += 1
result = dict(status="passed", csv_files=len(names), rows=rows, numeric_cells=numeric,
              matching_text_cells=text, maximum_absolute_difference=maximum)
(actual / "archive_table_comparison.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result))
