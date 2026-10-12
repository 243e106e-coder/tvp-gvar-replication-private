#!/usr/bin/env python3
"""Read-only audit of historical trade W and independent-shock artifacts.

Usage:
 python3 8.12/audit_historical_weights_irfs.py \
   --artifacts results/audit_artifacts/ --out results/historical_audit

Supply four downloaded/unzipped GitHub Actions 37872702188 artifact folders
under --artifacts, e.g. GPR/stability.csv and GPR/irf_stable_only.csv.
Never fabricate missing records or interpret absence as a zero response.
"""
import argparse,csv,hashlib,json,math,pathlib,re,sys,zipfile
COUNTRIES="AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()
SHOCKS="GPR VIX EPU OIL".split()
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(65536),b""):h.update(b)
    return h.hexdigest()
def csvrows(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None
def write(p,rows,cols):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--artifacts",default="results/audit_artifacts")
    pa.add_argument("--out",default="results/historical_audit")
    pa.add_argument("--weight",default="8.12/Trade_Weights_14_Economies_2000_2014.csv")
    a=pa.parse_args();o=pathlib.Path(a.out);o.mkdir(parents=True,exist_ok=True)
    weight=pathlib.Path(a.weight)
    if not weight.is_file():sys.exit("ERROR: weight CSV missing")
    rr=csvrows(weight);issues=[]
    if len(rr)!=14 or set(r.get("Reporter","") for r in rr)!=set(COUNTRIES):
        issues.append("INVALID_14_COUNTRY_REPORTER_SET")
    if not all(set(COUNTRIES).issubset(r) for r in rr):
        issues.append("MISSING_PARTNER_COLUMNS")
    weights=[]
    for r in rr:
        country=r.get("Reporter","")
        values=[num(r.get(cc)) for cc in COUNTRIES]
        if any(x is None or x<0 for x in values):issues.append("INVALID_WEIGHT_ROW_"+country)
        if country in COUNTRIES and values[COUNTRIES.index(country)] not in (None,0.0):
            issues.append("NONZERO_WEIGHT_DIAGONAL_"+country)
        sm=sum(z for z in values if z is not None)
        if abs(sm-1)>1e-6:issues.append("ROW_SUM_NOT_ONE_"+country)
        weights.append(dict(reporter=country,rowsum=sm,min_value=min((z for z in values if z is not None),default=None),
            max_value=max((z for z in values if z is not None),default=None)))
    write(o/"weight_matrix_checks.csv",weights,["reporter","rowsum","min_value","max_value"])
    # Upstream mirror-completion is explicitly documented in code, but cannot
    # be detected or attributed to individual cells from the normalized W alone.
    issues.append("UPSTREAM_2013_2014_MIRROR_COMPLETION_DOCUMENTED_CELL_LINEAGE_UNKNOWN")
    records=[];signs=[];missing=[]
    root=pathlib.Path(a.artifacts)
    for shock in SHOCKS:
        files=list(root.rglob("*.csv")) if root.exists() else []
        def pick(name):
            eligible=[p for p in files if p.name==name and shock.lower() in "/".join(p.parts).lower()]
            if len(eligible)!=1:
                missing.append(dict(shock=shock,file=name,matches=len(eligible)))
                return None
            return eligible[0]
        sf=pick("stability.csv");ir=pick("irf_stable_only.csv")
        if sf:
            for r in csvrows(sf):
                n=num(r.get("total"));stable=num(r.get("stable"));share=num(r.get("stable_share"))
                if n is None or stable is None or n<=0 or stable<0 or stable>n:issues.append("INVALID_STABILITY_"+shock)
                if n and stable is not None and share is not None and abs(stable/n-share)>1e-6:
                    issues.append("STABLE_SHARE_MISMATCH_"+shock)
                records.append(dict(shock=shock,date=r.get("date",""),total=n,stable=stable,stable_share=share,
                     low_stable_share=(share is not None and share<0.5),source_sha256=digest(sf)))
        if ir:
            for r in csvrows(ir):
                var=r.get("variable",""); date=r.get("date","");h=r.get("horizon","")
                if not re.match(r"^(AU|BR|CA|CH|CN|EA|UK|JP|KR|NO|SG|TR|US|ZA)_(y|dp|r|de|deq)$",var):continue
                med=num(r.get("median"));lo=num(r.get("low90"));hi=num(r.get("high90"))
                if None in (med,lo,hi):issues.append("INVALID_IRF_ROW_"+shock);continue
                if lo>hi or not lo<=med<=hi:issues.append("BAD_IRF_QUANTILE_ORDER_"+shock)
                # Positivity alone is not 'wrong'; record significance separately.
                signs.append(dict(shock=shock,date=date,variable=var,horizon=h,median=med,low90=lo,high90=hi,
                    direction=("positive" if med>0 else "negative" if med<0 else "zero"),
                    excludes_zero_90=(lo>0 or hi<0),stable_draws=r.get("stable_draws","")))
    write(o/"stability_review.csv",records,
          ["shock","date","total","stable","stable_share","low_stable_share","source_sha256"])
    write(o/"irf_direction_significance.csv",signs,
          ["shock","date","variable","horizon","median","low90","high90","direction","excludes_zero_90","stable_draws"])
    write(o/"missing_or_duplicate_artifacts.csv",missing,["shock","file","matches"])
    (o/"audit_summary.json").write_text(json.dumps(dict(weight_file=str(weight),weight_sha256=digest(weight),
        weight_cells=len(rr)*14,weight_checks=weights,issues=sorted(set(issues)),
        stability_records=len(records),irf_records=len(signs),missing_artifacts=missing,
        important="No raw bilateral source data or fill flags in final W; cannot validate individual 2013-14 observations.",
        note="Only reports file integrity and stable-only IRF; no MCMC convergence proof and no causal identification proof."),
        ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if missing:print("INCOMPLETE: provide exact independent-shock artifacts for full result review",file=sys.stderr)
    print(f"PASS read-only audit: {len(records)} stability rows, {len(signs)} IRF rows; issues={len(set(issues))}")
if __name__=="__main__":main()
