#!/usr/bin/env python3
"""First-difference exactly 14 rate columns in a freshly built model CSV."""
import argparse
import csv
import json
import math
from pathlib import Path

COUNTRIES = "AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()
RATE_COLUMNS = [f"{c}_r" for c in COUNTRIES]

def transform(source, output_dir):
    source, output_dir = Path(source), Path(output_dir)
    # Never overwrite inputs: rerunning from the same level CSV is deterministic.
    if source.resolve().parent == output_dir.resolve():
        raise ValueError("Output directory must differ from input directory")
    if source.name == "model_input_rate_diff.csv":
        raise ValueError("Already-differenced input is forbidden; use the level input")
    with source.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields):
            raise ValueError("Missing/duplicate headers")
        rows = list(reader)
    if not set(["Quarter"] + RATE_COLUMNS).issubset(fields):
        raise ValueError("Missing Quarter or one of the 14 rate columns")
    if len(rows) < 2:
        raise ValueError("At least two quarters are required")
    import re
    quarters = [r["Quarter"] for r in rows]
    if any(not re.fullmatch(r"[12][0-9]{3}Q[1-4]", q) for q in quarters):
        raise ValueError("Invalid quarters")
    index = [int(q[:4]) * 4 + int(q[-1]) for q in quarters]
    if any(b-a != 1 for a, b in zip(index, index[1:])):
        raise ValueError("Quarters must be unique, ordered and consecutive")
    for row in rows:
        if set(row) != set(fields):
            raise ValueError("Malformed CSV row")
        for field in fields[1:]:
            if not math.isfinite(float(row[field])):
                raise ValueError(f"Missing/non-finite {field} in {row['Quarter']}")
    if fields[0] != "Quarter":
        raise ValueError("Quarter must be the first column")
    result, audit = [], []
    for previous, current in zip(rows, rows[1:]):
        new = current.copy()
        for field in RATE_COLUMNS:
            before, now = float(previous[field]), float(current[field])
            difference = now - before
            new[field] = repr(difference)
            audit.append(dict(Quarter=current["Quarter"], country=field[:2],
                              previous_rate=before, current_rate=now,
                              rate_difference=difference))
        result.append(new)
    output_dir.mkdir(parents=True, exist_ok=True)
    def write(path, columns, data):
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(data)
    write(output_dir / "model_input_rate_diff.csv", fields, result)
    write(output_dir / "model_input_rate_level_aligned.csv", fields, rows[1:])
    write(output_dir / "rate_difference_audit.csv", list(audit[0]), audit)
    # Check serialized output, including exact preservation of every other field.
    with (output_dir / "model_input_rate_diff.csv").open(newline="") as f:
        saved = list(csv.DictReader(f))
    for i, row in enumerate(saved, 1):
        for field in fields:
            if field in RATE_COLUMNS:
                assert float(row[field]) == float(rows[i][field])-float(rows[i-1][field])
            else:
                assert row[field] == rows[i][field]
    note = {
        "rate_transform": "r_t - r_(t-1); no log, percent change or scaling",
        "rate_internal_name": "*_r is now a first difference",
        "source": str(source), "countries": COUNTRIES,
        "dropped_quarter": quarters[0],
        "first_quarter": quarters[1], "last_quarter": quarters[-1],
        "observations_before": len(rows), "observations_after": len(saved),
        "verification": "PASS: 14 rate differences and all retained non-rate values",
        "irf": "r IRFs describe rate changes. For rate-level IRFs, cumulatively sum each posterior draw from horizon 0 before computing quantiles.",
        "comparison": "Use model_input_rate_level_aligned.csv for a same-sample level baseline; refit both models."
    }
    (output_dir / "RATE_TRANSFORM.json").write_text(
        json.dumps(note, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(note, ensure_ascii=False, indent=2))
    return saved

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/model_input.csv")
    parser.add_argument("--output-dir", default="data/rate_diff")
    args = parser.parse_args()
    transform(args.input, args.output_dir)
