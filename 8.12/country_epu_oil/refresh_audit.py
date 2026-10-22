"""Refresh source-identifiable series; disclose unresolved lineage and stop MCMC.

This is an audit, not a claim that 70 macro series have been updated.
Downloads never substitute an unverified series or splice onto old observations.
"""
import hashlib,io,json,re,argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
import numpy as np
import pandas as pd

COUNTRIES="AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()
ROOT=Path("data/country_epu_oil")
WORKBOOK=Path("8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx")

def get(url,path):
    data=urlopen(Request(url,headers={"User-Agent":"Mozilla/5.0 research-data-audit"}),timeout=35).read()
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    return data,hashlib.sha256(data).hexdigest()

def quarterly(d,freq):
    d=d.copy();d.columns=["date","value"]
    d["date"]=pd.to_datetime(d.date,errors="raise")
    d["value"]=pd.to_numeric(d.value,errors="coerce")
    if d.date.duplicated().any(): raise ValueError("Duplicate source dates")
    d=d[d.date>="1999-01-01"].sort_values("date")
    d["Quarter"]=d.date.dt.to_period("Q").astype(str)
    expected=3 if freq=="M" else 1
    q=d.groupby("Quarter").value.agg(["mean","count"])
    q.loc[q["count"]!=expected,"mean"]=np.nan
    return d,q

def refresh_one(item):
    cc,var,sid,freq=item;url=f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={sid}"
    if var=="REER":
        area={"EA":"XM","UK":"GB"}.get(cc,cc)
        url=f"https://stats.bis.org/api/v1/data/WS_EER/M.R.B.{area}?startPeriod=1999-01&format=csv"
    elif var=="WTI":
        url="https://www.eia.gov/dnav/pet/hist_xls/RWTCm.xls"
    meta=dict(country=cc,variable=var,series_id=sid,url=url,frequency=freq)
    try:
        ext="xls" if var=="WTI" else "csv"
        data,sha=get(url,ROOT/"sources"/f"{sid}.{ext}")
        if data.lstrip().startswith(b"<"): raise ValueError("Source returned HTML instead of observations (maintenance/access page)")
        if var=="REER":
            parsed=pd.read_csv(io.BytesIO(data))
            if not ((parsed.FREQ=="M")&(parsed.EER_TYPE=="R")&(parsed.EER_BASKET=="B")&(parsed.REF_AREA==area)).all():
                raise ValueError("BIS series metadata mismatch")
            parsed=parsed[["TIME_PERIOD","OBS_VALUE"]]
        elif var=="WTI":
            parsed=pd.read_excel(io.BytesIO(data),sheet_name="Data 1",header=2).iloc[:,:2]
            parsed=parsed[pd.to_datetime(parsed.iloc[:,0],errors="coerce").notna()]
        else: parsed=pd.read_csv(io.BytesIO(data))
        d,q=quarterly(parsed,freq)
        q.to_csv(ROOT/"sources"/f"{sid}_quarterly.csv")
        v=q["mean"].dropna()
        if v.empty: raise ValueError("No complete quarters")
        meta.update(status="downloaded",sha256=sha,first=v.index[0],last=v.index[-1],
                    missing_quarters=int(q["mean"].isna().sum()),
                    meaning="Full-source refresh; no backfill from old workbook")
    except Exception as e: meta.update(status="failed",error=f"{type(e).__name__}: {e}")
    return meta

def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    manifest=[]
    # These identifiers are explicitly documented in the historical workbook.
    tasks=[(cc,"REER","RB"+{"EA":"XM","UK":"GB"}.get(cc,cc)+"BIS","M") for cc in COUNTRIES]
    tasks += [(cc,"RATE",sid,f) for cc,sid,f in [
        ("BR","IRSTCI01BRQ156N","Q"),("NO","IR3TIB01NOQ156N","Q"),
        ("ZA","IR3TIB01ZAQ156N","Q"),("JP","IRSTCI01JPM156N","M"),
        ("TR","IRSTCI01TRM156N","M")]]
    tasks += [(cc,"EQ",sid,f) for cc,sid,f in [
        ("BR","SPASTT01BRQ661N","Q"),("NO","SPASTT01NOQ661N","Q"),
        ("TR","SPASTT01TRQ661N","Q"),("ZA","SPASTT01ZAM661N","M")]]
    with ThreadPoolExecutor(max_workers=6) as pool: manifest=list(pool.map(refresh_one,tasks))
    # Risk inputs are refreshed as independent sources, preserving each range.
    risk_tables=[]
    for key,url,sheet in [
        ("US_EPU","https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx","Main News Index"),
        ("CN_EPU","https://www.policyuncertainty.com/media/China_Mainland_Paper_EPU.xlsx","EPU 2000 onwards")]:
        meta=dict(variable=key,url=url,sheet=sheet)
        try:
            data,sha=get(url,ROOT/"sources"/f"{key}.xlsx")
            d=pd.read_excel(io.BytesIO(data),sheet_name=sheet).iloc[:,:3]
            d.columns=["year","month","value"]
            yy=pd.to_numeric(d.year,errors="coerce");mm=pd.to_numeric(d.month,errors="coerce")
            d=d[yy.notna()&mm.notna()].copy()
            d["date"]=pd.to_datetime(dict(year=d.year,month=d.month,day=1))
            monthly,q=quarterly(d[["date","value"]],"M")
            # Genuine zero monthly EPU is retained; only positive quarterly means can be logged.
            if (monthly.value.dropna()<0).any(): raise ValueError("Negative EPU")
            q=q.rename(columns={"mean":key,"count":key+"_months"})
            q["LN_"+key]=np.log(q[key].where(q[key]>0))
            risk_tables.append(q)
            v=q[key].dropna();meta.update(status="downloaded",sha256=sha,first=v.index[0],last=v.index[-1])
        except Exception as e: meta.update(status="failed",error=str(e))
        manifest.append(meta)
    oil=refresh_one(("GL","WTI","MCOILWTICO","M"));manifest.append(oil)
    if oil["status"]=="downloaded":
        q=pd.read_csv(ROOT/"sources/MCOILWTICO_quarterly.csv",index_col="Quarter")
        q.columns=["WTI","WTI_months"];q["LN_WTI"]=np.log(q.WTI.where(q.WTI>0));risk_tables.append(q)
    if len(risk_tables)==3:
        # Outer union, so missing Chinese pre-2000 data stays visible.
        pd.concat(risk_tables,axis=1).sort_index().to_csv(ROOT/"risk_quarterly.csv",index_label="Quarter")
    pd.DataFrame(manifest).to_csv(ROOT/"refresh_manifest.csv",index=False)
    bykey={(m.get("country"),m["variable"]):m for m in manifest}
    inventory=[];old=pd.ExcelFile(WORKBOOK)
    for var in ["GDP","CPI","RATE","REER","EQ"]:
        d=pd.read_excel(old,"Raw_"+var).set_index("Quarter")
        for cc in COUNTRIES:
            v=pd.to_numeric(d[cc],errors="coerce").dropna();m=bykey.get((cc,var),{})
            inventory.append(dict(country=cc,variable=var,old_first=v.index[0],old_last=v.index[-1],
                old_missing=int(d[cc].isna().sum()),refresh_status=m.get("status","unresolved_source"),
                source_id=m.get("series_id",""),new_first=m.get("first",""),new_last=m.get("last",""),
                note="Do not use CBRT borrowing rate as an undocumented substitute for OECD interbank rate" if cc=="TR" and var=="RATE" else
                     "Need exact source identifier, unit, adjustment and extraction before updating" if not m else "Source download is not a full-panel validation"))
    pd.DataFrame(inventory).to_csv(ROOT/"macro_inventory.csv",index=False)
    blockers=["Historical workbook raw tables start at 2000Q1, not 1999Q1.",
       "Exact reproducible GDP/CPI and remaining equity/rate source mappings are unresolved; full 70-series refresh is not complete.",
       "New macro raw file must be rebuilt from verified sources before estimation; old processed workbook is not a substitute.",
       "National Cholesky ordering and inherited autonomous Global block require sensitivity checks."]
    report=dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                macro_series=70,macro_downloads=sum(m.get("status")=="downloaded" for m in manifest[:23]),
                all_downloads=sum(m.get("status")=="downloaded" for m in manifest),
                estimation_ready=False,blockers=blockers)
    (ROOT/"READINESS.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()
