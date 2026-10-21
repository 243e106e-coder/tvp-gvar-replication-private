"""Strict quarterly levels -> transformed inputs with a 2000Q1 effective start.

raw_macro.csv columns: Quarter and CC_GDP, CC_CPI, CC_RATE, CC_REER, CC_EQ.
All series must have observation-level source lineage reviewed separately.
No index reconstruction from growth rates, imputation or source splicing.
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

COUNTRIES = "AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()

def transform(raw, risk, p=1, gdp="growth"):
    if p not in (1,2) or gdp not in ("growth","log_level"):
        raise ValueError("Invalid specification")
    raw=raw.set_index("Quarter").copy()
    if raw.index.duplicated().any(): raise ValueError("Duplicate raw quarters")
    raw.index=pd.PeriodIndex(raw.index,freq="Q")
    raw=raw.sort_index()
    if not raw.index.equals(pd.period_range(raw.index.min(),raw.index.max(),freq="Q")):
        raise ValueError("Raw quarters are not contiguous")
    if raw.index.min()>pd.Period("1999Q1",freq="Q"):
        raise ValueError("Raw inputs must include 1999Q1 onward")
    out=pd.DataFrame(index=raw.index)
    for cc in COUNTRIES:
        for source,target in [("GDP","y"),("CPI","dp"),("RATE","r"),("REER","de"),("EQ","deq")]:
            s=pd.to_numeric(raw[f"{cc}_{source}"],errors="raise")
            if source!="RATE":
                if (s.dropna()<=0).any(): raise ValueError(f"Nonpositive {cc}_{source}")
                s=np.log(s)
                s=s if source=="GDP" and gdp=="log_level" else 100*s.diff()
            else: s=s.diff() # percentage points; never log interest rates
            out[f"{cc}_{target}"]=s
    # Keep p transformed presample rows so the estimator loses no 2000Q1 row.
    first=pd.Period("2000Q1",freq="Q")-p
    out=out.loc[first:]
    risk=risk.set_index("Quarter").copy()
    if risk.index.duplicated().any(): raise ValueError("Duplicate risk quarters")
    risk.index=pd.PeriodIndex(risk.index,freq="Q")
    needed=[f"{k}_months" for k in ("US_EPU","CN_EPU","WTI")]+[f"LN_{k}" for k in ("US_EPU","CN_EPU","WTI")]
    r=risk.reindex(out.index)
    ok=out.notna().all(axis=1) & np.isfinite(out).all(axis=1)
    ok &= r[needed].notna().all(axis=1)
    ok &= np.isfinite(r[needed]).all(axis=1)
    ok &= (r[["US_EPU_months","CN_EPU_months","WTI_months"]]==3).all(axis=1)
    if not ok.iloc[0]:
        missing=list(out.columns[out.iloc[0].isna()])+list(r.columns[r.iloc[0].isna()])
        raise ValueError(f"Missing presample {first}: {missing}. Never shift the start silently.")
    # Internal gaps terminate the common sample. Later isolated valid rows are
    # disclosed and are never bridged by dropping the gap.
    bad=np.flatnonzero(~ok.to_numpy());end=bad[0] if len(bad) else len(out)
    if end<=p+30: raise ValueError("Insufficient continuous common sample from 2000Q1")
    sample=out.iloc[:end].copy(); sample.index=sample.index.astype(str)
    sample.index.name="Quarter"
    audit={"effective_start":"2000Q1","input_start":str(first),"common_end":sample.index[-1],
           "p":p,"gdp":gdp,"transformations":{"y":"100*dlog(real GDP)" if gdp=="growth" else "log(real GDP)",
             "dp":"100*dlog(CPI)","r":"diff(rate), percentage points","de":"100*dlog(REER)","deq":"100*dlog(equity)"},
           "first_gap":str(out.index[end]) if end<len(out) else None,
           "valid_quarters_after_gap":int(ok.iloc[end:].sum())}
    return sample.reset_index(),audit

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--p",type=int,default=1)
    a.add_argument("--gdp",choices=["growth","log_level"],default="growth")
    a.add_argument("--root",default="data/country_epu_oil");v=a.parse_args();root=Path(v.root)
    d,meta=transform(pd.read_csv(root/"raw_macro.csv"),pd.read_csv(root/"risk_quarterly.csv"),v.p,v.gdp)
    d.to_csv(root/"model_input.csv",index=False)
    (root/"sample_audit.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps(meta,indent=2))
