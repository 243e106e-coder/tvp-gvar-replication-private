#!/usr/bin/env python3
"""Retrieve one historical posterior.rda, without modifying repository data."""
import argparse,io,json,os,pathlib,sys,zipfile
import importlib.util
spec=importlib.util.spec_from_file_location("historical_fetch","8.12/fetch_historical_artifacts.py")
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--shock",required=True,choices=("GPR","VIX","EPU","OIL"))
 p.add_argument("--run-id",type=int,default=37872702188)
 p.add_argument("--out",default="results/raw_posterior_input")
 a=p.parse_args()
 token=os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
 if not token:sys.exit("GITHUB_TOKEN missing")
 repo=os.environ.get("GITHUB_REPOSITORY","243e106e-coder/tvp-gvar-replication-private")
 base=f"https://api.github.com/repos/{repo}"
 historical=module.get(f"{base}/actions/runs/{a.run_id}",token)
 if historical["conclusion"]!="success":sys.exit("Historical run not successful")
 data=module.get(f"{base}/actions/runs/{a.run_id}/artifacts?per_page=100",token)
 name=f"independent-{a.shock}-{a.run_id}"
 found=[x for x in data["artifacts"] if x["name"]==name and not x.get("expired")]
 if len(found)!=1:sys.exit(f"Expected one live artifact {name}; found {len(found)}")
 raw=module.get(f'{base}/actions/artifacts/{found[0]["id"]}/zip',token,binary=True)
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  matches=[n for n in z.namelist() if pathlib.PurePosixPath(n).name=="posterior.rda"]
  if len(matches)!=1:sys.exit(f"Expected one posterior.rda, got {len(matches)}")
  data=z.read(matches[0])
  if not data:sys.exit("Empty posterior.rda")
 out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
 (out/"posterior.rda").write_bytes(data)
 meta={"shock":a.shock,"historical_run_id":a.run_id,"historical_sha":historical["head_sha"],
       "artifact_id":found[0]["id"],"posterior_bytes":len(data),"zip_member":matches[0]}
 (out/"posterior_source.json").write_text(json.dumps(meta,indent=2)+"\n")
 print("Retrieved historical posterior for",a.shock,"bytes",len(data))
if __name__=="__main__":main()
