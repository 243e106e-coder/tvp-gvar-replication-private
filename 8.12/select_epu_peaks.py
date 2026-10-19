#!/usr/bin/env python3
"""Select EPU peaks ex post from observable quarterly EPU, never from IRFs."""
import argparse,csv,json,math
from pathlib import Path

def select_peaks(quarters,values,count=6,gap=4):
    if len(quarters)!=len(values) or len(values)<2:
        raise ValueError("Need at least two aligned quarters")
    if count<1 or gap<1 or any(not math.isfinite(v) or v<0 for v in values):
        raise ValueError("Invalid peak parameters/data")
    # Plateau: keep its earliest quarter, if above adjacent plateau values.
    candidates=[];i=0
    while i<len(values):
        j=i
        while j+1<len(values) and values[j+1]==values[i]: j+=1
        left=values[i-1] if i else -math.inf
        right=values[j+1] if j+1<len(values) else -math.inf
        if values[i]>left and values[i]>right:
            candidates.append(dict(index=i,Quarter=quarters[i],EPU=values[i],
                                   boundary=i==0 or j==len(values)-1))
        i=j+1
    ranked=sorted(candidates,key=lambda r:(-r["EPU"],r["index"]))
    chosen=[]
    for rank,r in enumerate(ranked,1):
        r["height_rank"]=rank
        r["selected"]=len(chosen)<count and all(abs(r["index"]-v["index"])>=gap for v in chosen)
        if r["selected"]: chosen.append(r)
    return sorted(chosen,key=lambda r:r["index"]),sorted(ranked,key=lambda r:r["index"])

def build(input_file,quarter_file,output_dir,p=1,count=6,gap=4,mode="period5y"):
    with open(input_file,encoding="utf-8-sig",newline="") as f: macro=list(csv.DictReader(f))
    qq=[r["Quarter"] for r in macro][p:]
    if p not in (1,2) or len(set(qq))!=len(qq): raise ValueError("Bad p/quarters")
    with open(quarter_file,encoding="utf-8-sig",newline="") as f:
        risk={r["Quarter"]:r for r in csv.DictReader(f)}
    chosen_all=[];candidates_all=[]
    for ec in ("US","CN"):
        vals=[float(risk[q][ec+"_EPU"]) for q in qq]
        local_selected,candidates=select_peaks(qq,vals,count,gap)
        if mode=="period5y":
            selected=[]
            periods=sorted(set(int(q[:4])//5*5 for q in qq))
            for start in periods:
                eligible=[i for i,q in enumerate(qq) if start<=int(q[:4])<=start+4]
                best=max(eligible,key=lambda i:(vals[i],-i))
                selected.append(dict(index=best,Quarter=qq[best],EPU=vals[best],
                    boundary=best==0 or best==len(qq)-1,height_rank=0,selected=True,
                    period=f"{start}-{min(start+4,int(qq[-1][:4]))}",method="period5y"))
        else:
            selected=[dict(r,period="full_sample",method="local_top") for r in local_selected]
        candidates=[dict(r,period="full_sample",method="local_top") for r in candidates]
        for dest,src in ((chosen_all,selected),(candidates_all,candidates)):
            for row in src: dest.append(dict(epu_country=ec,**row))
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    cols=["epu_country","index","Quarter","EPU","boundary","height_rank","selected","period","method"]
    for name,rows in (("epu_peak_dates.csv",chosen_all),("epu_peak_candidates.csv",candidates_all)):
        with (out/name).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
    rules=dict(frequency="quarterly; peak of arithmetic mean of 3 monthly EPU values",
       sample_first=qq[0],sample_last=qq[-1],p=p,target_peaks_per_country=count,
       minimum_quarter_distance_for_local_mode=gap,selected_mode=mode,
       method="Default period5y: highest quarterly EPU in each calendar five-year bin (2000-04,2005-09,...); first/last bins may be incomplete. Ties earliest quarter. Alternative local_top: top local maxima separated >=4 quarters, up to six. Candidates table preserves local peaks.",
       basis="US dates from US EPU; CN dates from Davis-Liu-Sheng mainland EPU. Shared across GLOBAL/US placement and both shocks.",
       interpretation="Ex-post descriptive peak dates, not real-time signal or dates chosen from IRF significance.")
    (out/"epu_peak_rules.json").write_text(json.dumps(rules,indent=2)+"\n")
    print(json.dumps(dict(rules=rules,selected=chosen_all),indent=2))
if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--input",default="data/rate_diff/model_input_rate_diff.csv")
    a.add_argument("--quarter-file",default="8.12/epu_wti_quarterly.csv")
    a.add_argument("--output-dir",default="data/rate_diff")
    a.add_argument("--p",type=int,default=1)
    a.add_argument("--count",type=int,default=6)
    a.add_argument("--gap",type=int,default=4)
    a.add_argument("--mode",choices=["period5y","local_top"],default="period5y")
    v=a.parse_args();build(v.input,v.quarter_file,v.output_dir,v.p,v.count,v.gap,v.mode)
