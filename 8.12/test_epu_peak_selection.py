#!/usr/bin/env python3
import csv,json
from pathlib import Path
from select_epu_peaks import select_peaks
s,_=select_peaks(["a","b","c","d"],[2,2,1,3],6,1)
assert [r["index"] for r in s]==[0,3]
s,_=select_peaks(["a","b","c","d"],[2,1,2,1],6,3)
assert [r["index"] for r in s]==[0]
folder=Path("data/rate_diff")
rules=json.loads((folder/"epu_peak_rules.json").read_text())
with (folder/"model_input_rate_diff.csv").open() as f:
    qq=[r["Quarter"] for r in csv.DictReader(f)][rules["p"]:]
with open("8.12/epu_wti_quarterly.csv") as f:
    risk={r["Quarter"]:r for r in csv.DictReader(f)}
with (folder/"epu_peak_dates.csv").open() as f: peaks=list(csv.DictReader(f))
for ec in ("US","CN"):
    selected=[r for r in peaks if r["epu_country"]==ec]
    assert len(selected)==len(set(int(q[:4])//5 for q in qq))
    for r in selected:
        lo,hi=map(int,r["period"].split("-"))
        eligible=[q for q in qq if lo<=int(q[:4])<=hi]
        expected=max(eligible,key=lambda q:(float(risk[q][ec+"_EPU"]),-qq.index(q)))
        assert r["Quarter"]==expected
        assert float(r["EPU"])==float(risk[expected][ec+"_EPU"])
print("PASS country-specific period maxima, earliest tie, local plateau, boundary and gap checks")
