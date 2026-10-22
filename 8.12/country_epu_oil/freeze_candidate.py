"""Freeze selected observations, transformations, source receipts and coverage."""
import gzip, hashlib, json, shutil
from pathlib import Path
import numpy as np
import pandas as pd
from build_input import transform,COUNTRIES

FILES=["raw_macro.csv","risk_quarterly.csv","verified_macro_inventory.csv","source_receipts.json",
       "rate_proxy_changes.csv","raw_missing_observations.csv","model_input.csv","sample_audit.json",
       "peak_dates.csv","peak_rules.json"]
def verify(root):
    root=Path(root);m=json.loads((root/"snapshot_hashes.json").read_text())
    for name,sha in m.items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha:raise ValueError("Snapshot hash mismatch: "+name)
    inv=pd.read_csv(root/"verified_macro_inventory.csv")
    if len(inv)!=70 or inv.duplicated(["country","variable"]).any():raise ValueError("Not 70 unique macro series")
    raw=pd.read_csv(root/"raw_macro.csv");risk=pd.read_csv(root/"risk_quarterly.csv")
    built,audit=transform(raw,risk,p=1,gdp="growth")
    saved=pd.read_csv(root/"model_input.csv")
    if not saved.Quarter.equals(built.Quarter) or not np.allclose(saved.iloc[:,1:],built.iloc[:,1:],rtol=1e-10,atol=1e-10):
        raise ValueError("Transformed snapshot is not reproducible from raw levels")
    if audit["effective_start"]!="2000Q2":raise ValueError("Incorrect effective start")
    print(json.dumps(dict(status="PASS frozen candidate",macro_series=70,rate_proxy_changes=4,
         actual_estimation_start="2000Q2",common_end=audit["common_end"],estimation_quarters=len(saved)-1,
         formal_original_rate_specification_ready=False,exploratory_candidate_ready=True),indent=2))

def freeze(source="data/country_epu_oil",target="8.12/country_epu_oil/data_20261022"):
    source=Path(source);target=Path(target);target.mkdir(parents=True,exist_ok=True)
    for f in FILES:shutil.copy2(source/f,target/f)
    # Persist all selected raw daily/monthly/quarterly observations in a compact,
    # deterministic archive. This does not republish the author's unrelated data.
    with (source/"source_observations.csv").open("rb") as src, (target/"source_observations.csv.gz").open("wb") as out:
        with gzip.GzipFile(filename="",mode="wb",fileobj=out,mtime=0) as gz:shutil.copyfileobj(src,gz)
    raw=pd.read_csv(source/"raw_macro.csv").set_index("Quarter")
    transformed,audit=transform(raw.reset_index(),pd.read_csv(source/"risk_quarterly.csv"),1,"growth")
    # Exact first-gap causes are disclosed; no claim that all latest observations
    # can enter a continuous complete-case TVP-GVAR sample.
    gap=audit["first_gap"];columns=raw.columns[raw.loc[gap].isna()].tolist() if gap else []
    readiness=dict(exploratory_candidate_ready=True,formal_original_rate_specification_ready=False,
      data_series=70,effective_start="2000Q2",common_end=audit["common_end"],
      estimation_quarters=len(transformed)-1,first_raw_gap=gap,first_raw_gap_variables=columns,
      rate_proxy_changes=["CN: 1-day repo","US: overnight call-money rate","SG: 1-year Treasury bill yield","TR: discount rate"],
      no_imputation=True,no_cross_series_splicing=True,
      limitations=["Four rate proxies differ from inherited definitions; exploratory candidate only.",
        "China GDP and repo are publisher-constructed quarterly series, not untouched NBS raw releases; see Atlanta Fed readme.",
        "Native quarterly Australian CPI is used as a published full-quarter statistic; monthly series require exactly 3 observed months.",
        "Euro-area HICP uses official changing composition; no fixed-membership aggregate was constructed.",
        "Source snapshots include observations after common_end; these are not bridged across gaps.",
        "Trade-weight lineage remains unresolved; unchanged weight matrix is an additional model limitation."])
    (target/"READINESS.json").write_text(json.dumps(readiness,indent=2)+"\n")
    names=FILES+["source_observations.csv.gz","READINESS.json"]
    hashes={f:hashlib.sha256((target/f).read_bytes()).hexdigest() for f in names}
    (target/"snapshot_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    verify(target)

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser();p.add_argument("--verify",action="store_true");v=p.parse_args()
    if v.verify:verify("8.12/country_epu_oil/data_20261022")
    else:freeze()
