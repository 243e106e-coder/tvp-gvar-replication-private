"""Representative quarterly level peaks; fixed rules, independent of IRFs."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

RULE = {"height_quantile":0.75,"minimum_prominence_sd":0.5,
        "minimum_distance_quarters":4,"neighbourhood_quarters":4,
        "endpoints":"excluded: cannot confirm a two-sided peak",
        "ties":"first quarter of a flat peak; earliest date wins equal-height selection",
        "interpretation":"statistical representative peaks; event labels require separate verification"}

def select(risk, effective_quarters):
    r=risk.set_index("Quarter").loc[effective_quarters]
    rows=[]
    for series in ("CN_EPU","US_EPU","WTI"):
        z=pd.to_numeric(r[series],errors="raise")
        if z.isna().any() or not np.isfinite(z).all() or not (z>0).all():
            raise ValueError("Invalid peak series")
        x=z.to_numpy(); scale=float(z.std(ddof=0)); threshold=float(z.quantile(.75))
        candidates=[]
        for i in range(1,len(x)-1):
            j=i
            while j+1<len(x) and x[j+1]==x[i]: j+=1
            if x[i]<=x[i-1] or j==len(x)-1 or x[j+1]>=x[i]: continue
            left=x[max(0,i-4):i]; right=x[j+1:min(len(x),j+5)]
            prominence=float(x[i]-max(left.min(),right.min()))
            if x[i]>=threshold and prominence>=.5*scale:
                candidates.append((i,prominence))
        accepted=[]
        for i,prom in sorted(candidates,key=lambda t:(-x[t[0]],t[0])):
            if all(abs(i-j)>=4 for j,_ in accepted): accepted.append((i,prom))
        for i,prom in sorted(accepted):
            rows.append(dict(series=series,shock="WTI" if series=="WTI" else "EPU",
                epu_country=series[:2] if series!="WTI" else "BOTH",Quarter=z.index[i],
                value=float(x[i]),prominence=prom,height_threshold=threshold,
                prominence_threshold=.5*scale,method="representative_local_peak",
                event_status="not yet verified"))
    return pd.DataFrame(rows,columns=["series","shock","epu_country","Quarter","value",
        "prominence","height_threshold","prominence_threshold","method","event_status"])

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--p",type=int,default=1);v=a.parse_args()
    root=Path("data/country_epu_oil")
    q=pd.read_csv(root/"model_input.csv").Quarter.tolist()[v.p:]
    if q[0]!="2000Q2": raise ValueError("Effective sample must start at 2000Q2")
    risk=pd.read_csv(root/"risk_quarterly.csv")
    peaks=select(risk,q)
    peaks.to_csv(root/"peak_dates.csv",index=False)
    risk.set_index("Quarter").loc[q,["CN_EPU","US_EPU","WTI"]].to_csv(root/"peak_selection_input.csv")
    (root/"peak_rules.json").write_text(json.dumps(dict(rule=RULE,
        WTI="Same dates in both models; separate model responses",
        selection="Quarterly level peaks, never IRF significance",first=q[0],last=q[-1]),indent=2))
    if set(peaks.series)!={"CN_EPU","US_EPU","WTI"}:
        raise ValueError("A series has no qualifying peak; inspect inputs and rules, do not silently force a date")
