"""Select separate EPU/WTI quarterly level peaks, without consulting any IRF."""
import argparse,json
from pathlib import Path
import pandas as pd

def select(risk, effective_quarters):
    r=risk.set_index("Quarter").loc[effective_quarters]
    result=[]
    for series in ("CN_EPU","US_EPU","WTI"):
        z=pd.to_numeric(r[series],errors="raise")
        if z.isna().any() or not (z>0).all(): raise ValueError("Invalid peak series")
        groups=z.groupby([int(q[:4])//5*5 for q in z.index])
        for year,values in groups:
            q=values.idxmax() # first quarter wins ties (input ordered)
            result.append(dict(series=series,shock="WTI" if series=="WTI" else "EPU",
                epu_country=series[:2] if series!="WTI" else "BOTH",Quarter=q,
                value=float(values.loc[q]),period=f"{year}-{min(year+4,int(z.index[-1][:4]))}",
                method="calendar5y_max",observed_quarters=len(values)))
    return pd.DataFrame(result)

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--p",type=int,default=1);v=a.parse_args()
    root=Path("data/country_epu_oil")
    q=pd.read_csv(root/"model_input.csv").Quarter.tolist()[v.p:]
    if q[0]!="2000Q1": raise ValueError("Effective sample must start at 2000Q1")
    peaks=select(pd.read_csv(root/"risk_quarterly.csv"),q)
    peaks.to_csv(root/"peak_dates.csv",index=False)
    (root/"peak_rules.json").write_text(json.dumps(dict(rule="Maximum quarterly mean within each calendar five-year bin; earliest tie. Partial bins disclosed.",
        WTI="WTI dates apply in both separately estimated systems",selection="Ex-post level peaks, not innovations, returns or IRF significance",first=q[0],last=q[-1]),indent=2))
