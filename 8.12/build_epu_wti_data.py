#!/usr/bin/env python3
"""Reproduce frozen official EPU/WTI quarterly inputs; never fill missing months."""
import hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "epu_wti_sources"
SOURCES = {
 "us.xlsx": ("https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx",
             "US news-based EPU", "Main News Index"),
 "cn.xlsx": ("https://www.policyuncertainty.com/media/SCMP_China_Policy_Uncertainty_Data.xlsx",
             "China SCMP EPU; Factiva only, no ProQuest splice", "China EPU Index_Factiva"),
 "wti.xls": ("https://www.eia.gov/dnav/pet/hist_xls/RWTCm.xls",
             "EIA RWTC monthly Cushing OK WTI spot FOB USD/barrel", "Data 1")
}
def aggregate(frame, label):
    if frame.duplicated(["year","month"]).any():
        raise ValueError("Duplicate source months: "+label)
    if not frame.month.between(1,12).all():
        raise ValueError("Invalid month: "+label)
    if not np.isfinite(frame.value).all() or (frame.value<0).any():
        raise ValueError("Negative or missing observation: "+label)
    frame = frame.sort_values(["year","month"]).copy()
    frame["Quarter"] = frame.year.astype(str)+"Q"+((frame.month-1)//3+1).astype(str)
    q=frame.groupby("Quarter").agg(value=("value","mean"),months=("month","nunique"))
    # Only full quarters are eligible. A missing quarter inside model window is fatal downstream.
    q=q[q.months==3]
    if (q.value<=0).any(): raise ValueError("Cannot log a non-positive quarterly mean: "+label)
    q.columns=[label, label+"_months"]
    return frame,q

def build():
    tables=[]; monthly=[]; manifest=[]
    for name,(url,label,sheet) in SOURCES.items():
        p=RAW/name
        if not p.exists(): raise FileNotFoundError(p)
        if name!="wti.xls":
            d=pd.read_excel(p,sheet_name=sheet)
            vals=[pd.to_numeric(d.iloc[:,i],errors="coerce") for i in range(3)]
            # Only explicit non-data footer rows are excluded. Missing values within dated rows fail.
            dated=vals[0].notna() & vals[1].notna()
            f=pd.DataFrame({"year":vals[0][dated].astype(int),
                            "month":vals[1][dated].astype(int),"value":vals[2][dated]})
            key="US_EPU" if name=="us.xlsx" else "CN_EPU"
        else:
            d=pd.read_excel(p,sheet_name=sheet,header=2)
            dates=pd.to_datetime(d.iloc[:,0],errors="coerce")
            dated=dates.notna()
            f=pd.DataFrame({"year":dates[dated].dt.year,"month":dates[dated].dt.month,
                           "value":pd.to_numeric(d.iloc[:,1][dated],errors="coerce")})
            key="WTI"
        f,q=aggregate(f,key); tables.append(q)
        f["series"]=key; monthly.append(f)
        manifest.append(dict(file=name,url=url,description=label,sheet=sheet,
                             sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                             first_month=f"{f.iloc[0].year:04d}-{f.iloc[0].month:02d}",
                             last_month=f"{f.iloc[-1].year:04d}-{f.iloc[-1].month:02d}"))
    q=pd.concat(tables,axis=1,join="inner").sort_index()
    for key in ["US_EPU","CN_EPU","WTI"]: q["LN_"+key]=np.log(q[key])
    q.to_csv(ROOT/"epu_wti_quarterly.csv",index_label="Quarter",float_format="%.15g")
    pd.concat(monthly).to_csv(ROOT/"epu_wti_monthly_audit.csv",index=False,float_format="%.15g")
    payload=dict(downloaded_utc="2026-10-19",sources=manifest,
                 aggregation="Equal-weight arithmetic mean of exactly 3 monthly values; then natural log. EIA monthly averages, not daily-weighted quarterly averages.",
                 missing="No zero replacement, interpolation, source splice or extrapolation.",
                 china="SCMP Factiva version; China index is not US EPU even when placed in the US block.",
                 wti="Nominal spot price, not futures or Brent; not an identified oil-supply shock.",
                 first_quarter=q.index[0],last_quarter=q.index[-1],observations=len(q))
    (ROOT/"epu_wti_source_manifest.json").write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))
if __name__=="__main__": build()
