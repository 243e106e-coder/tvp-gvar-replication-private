#!/usr/bin/env python3
"""Country/date audit of historic stable-only IRFs; no model re-estimation."""
import argparse,csv,math,pathlib,statistics
SHOCKS=("GPR","VIX","EPU","OIL")
COUNTRIES="AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()
DATES=("2003Q1","2008Q3","2014Q3","2020Q1","2022Q1","2023Q4")
def read(p):
 if not p.is_file():raise FileNotFoundError(p)
 with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def write(p,rows):
 with p.open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--input",default="results/audit_artifacts");ap.add_argument("--out",default="results/irf_stage2")
 args=ap.parse_args();inp=pathlib.Path(args.input);out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
 stability=[];countries=[];summary=[]
 for shock in SHOCKS:
  stab=read(inp/shock/"stability.csv");irf=read(inp/shock/"irf_stable_only.csv")
  shares={};retained={}
  for row in stab:
   d=row["date"];tot=int(row["total"]);ok=int(row["stable"]);share=float(row["stable_share"])
   if d in shares or tot<=0 or not 0<=ok<=tot or abs(ok/tot-share)>1e-6:raise ValueError("Invalid stability record "+shock+" "+d)
   shares[d]=share;retained[d]=ok
   stability.append(dict(shock=shock,date=d,draws=tot,stable_draws=ok,excluded_draws=tot-ok,stable_share=share,below_50pct=share<.5))
  if set(shares)!=set(DATES):raise ValueError("Missing dates for "+shock)
  grouped={}
  for r in irf:
   parts=r["variable"].split("_",1)
   if len(parts)!=2 or parts[0] not in COUNTRIES or parts[1] not in ("y","deq"):continue
   d=r["date"];h=int(r["horizon"]);m,lo,hi=map(float,(r["median"],r["low90"],r["high90"]))
   if d not in shares or not 0<=h<=12 or not all(map(math.isfinite,(m,lo,hi))) or not lo<=m<=hi:raise ValueError("Invalid IRF row")
   if int(r["stable_draws"])!=retained[d]:raise ValueError("IRF stable draw count mismatch")
   grouped.setdefault((d,parts[0],parts[1]),[]).append((h,m,lo,hi))
  if len(grouped)!=6*14*2:raise ValueError("Incomplete GDP/stock data "+shock)
  for (d,cc,v),rows in sorted(grouped.items()):
   if sorted(r[0] for r in rows)!=list(range(13)):raise ValueError("Missing or repeated horizon")
   countries.append(dict(shock=shock,date=d,country=cc,variable=v,stable_share=shares[d],
    positive_horizons=sum(r[1]>0 for r in rows),negative_horizons=sum(r[1]<0 for r in rows),
    significant_positive_90=sum(r[2]>0 for r in rows),significant_negative_90=sum(r[3]<0 for r in rows),
    median_h0=next(r[1] for r in rows if r[0]==0),median_h4=next(r[1] for r in rows if r[0]==4),median_h12=next(r[1] for r in rows if r[0]==12)))
  for v in ("y","deq"):
   q=[r for r in countries if r["shock"]==shock and r["variable"]==v]
   summary.append(dict(shock=shock,variable=v,mean_stable_share=statistics.mean(shares.values()),
    min_stable_share=min(shares.values()),max_stable_share=max(shares.values()),
    positive_point_estimates=sum(r["positive_horizons"] for r in q),negative_point_estimates=sum(r["negative_horizons"] for r in q),
    significant_positive_90=sum(r["significant_positive_90"] for r in q),significant_negative_90=sum(r["significant_negative_90"] for r in q),
    total_observations=6*14*13))
 write(out/"stability_by_date.csv",stability);write(out/"gdp_stock_by_country_date.csv",countries);write(out/"gdp_stock_summary.csv",summary)
 (out/"LIMITATIONS.md").write_text("# Interpretation limits\nThese IRFs are conditional on draws with spectral radius < 1. Historical CSVs cannot establish full spectral-radius distributions, G condition numbers, MCMC convergence or ESS. These are 90% Bayesian intervals, not independent tests; counts are dependent and unadjusted for multiple comparisons. They do not validate historical trade source data or structural identification.\n",encoding="utf-8")
 print("PASS",len(stability),"stability records",len(countries),"country/date summaries")
if __name__=="__main__":main()
