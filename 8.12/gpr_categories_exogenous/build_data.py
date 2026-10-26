"""Frozen official GPR event categories, each used separately as exogenous z.
Monthly zeros are genuine observations. Average THREE months before logging.
No interpolation, epsilon, rebasing, source splicing, or post-macro-sample peaks.
"""
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
DATA=HERE/'data_20261026'
BASE=HERE.parent/'country_epu_oil/data_20261022'
COUNTRIES='AU BR CA CH CN EA UK JP KR NO SG TR US ZA'.split()
GLOBALS={'GPR':'Recent overall GPR','GPRT':'Recent geopolitical threats','GPRA':'Recent geopolitical acts',
         'GPRH':'Historical overall GPR','GPRHT':'Historical geopolitical threats','GPRHA':'Historical geopolitical acts'}
CATEGORIES={f'SHAREH_CAT_{i}':name for i,name in enumerate(['War threats','Peace threats','Military buildups','Nuclear threats','Terror threats','Beginning of war','Escalation of war','Terror acts'],1)}
SERIES=list(GLOBALS)+list(CATEGORIES)
PEAK_RULE={'height_quantile':.75,'minimum_prominence_sd':.5,'minimum_distance_quarters':4,'neighbourhood_quarters':4,
 'endpoints':'excluded; two-sided confirmation required','ties':'first quarter of plateau; higher peaks first, earlier wins ties',
 'selection_sample':'actual common estimation window ONLY; not full available GPR history','selection_target':'each raw quarterly category mean, never its IRF or significance'}

def quarterly(monthly):
 x=monthly.copy();idx=pd.DatetimeIndex(x.index)
 if idx.duplicated().any() or not (idx.day==1).all():raise ValueError('Invalid monthly dates')
 x=x.sort_index().reindex(pd.date_range(idx.min(),idx.max(),freq='MS'))
 x=x.apply(pd.to_numeric,errors='raise')
 if (x<0).any().any():raise ValueError('Negative risk index')
 # Native zero remains zero; neither monthly logs nor two-month means.
 return x.resample('QE').mean().where(x.resample('QE').count()==3)

def select_peaks(series):
 z=pd.to_numeric(series,errors='raise')
 if z.isna().any() or not np.isfinite(z).all() or not (z>0).all():raise ValueError('Invalid quarterly peak series')
 x=z.to_numpy();sd=float(z.std(ddof=0));threshold=float(z.quantile(.75));candidates=[]
 for i in range(1,len(x)-1):
  j=i
  while j+1<len(x) and x[j+1]==x[i]:j+=1
  if x[i]<=x[i-1] or j==len(x)-1 or x[j+1]>=x[i]:continue
  prominence=float(x[i]-max(x[max(0,i-4):i].min(),x[j+1:min(len(x),j+5)].min()))
  if x[i]>=threshold and prominence>=.5*sd:candidates.append((i,prominence))
 accepted=[]
 for i,prom in sorted(candidates,key=lambda t:(-x[t[0]],t[0])):
  if all(abs(i-j)>=4 for j,_ in accepted):accepted.append((i,prom))
 return [dict(Quarter=str(z.index[i]),value=float(x[i]),prominence=prom,height_threshold=threshold,prominence_threshold=.5*sd) for i,prom in sorted(accepted)]

def build():
 src=DATA/'source_data_gpr_export.xls';d=pd.read_excel(src)
 if d.month.duplicated().any():raise ValueError('Duplicate official month')
 labels=d[['var_name','var_label']].dropna().set_index('var_name').var_label.to_dict()
 numeric=d.drop(columns=['var_name','var_label']).set_index('month');numeric.index=pd.DatetimeIndex(numeric.index)
 # Store a complete inventory, not 112 falsely distinct event-risk categories.
 inventory=[]
 for col in numeric:
  s=pd.to_numeric(numeric[col],errors='raise');valid=s.dropna()
  role='event_category' if col in CATEGORIES else 'overall_or_threats_acts' if col in GLOBALS else 'country_specific' if col.startswith(('GPRC_','GPRHC_')) else 'article_count_not_risk' if col in ['N10','N3H'] else 'same_index_alternative_scale' if col in ['SHARE_GPR','SHARE_GPRH'] else 'search_definition_robustness_variant'
  inventory.append(dict(series=col,publisher_label=labels.get(col,''),role=role,estimated_this_batch=col in SERIES,first_month=str(valid.index.min().date()),last_month=str(valid.index.max().date()),observations=len(valid),zero_months=int((valid==0).sum())))
 pd.DataFrame(inventory).to_csv(DATA/'official_series_inventory.csv',index=False)
 monthly=numeric.loc['1999-01-01':,SERIES]
 if monthly.isna().any().any():raise ValueError('Main category monthly history contains gaps')
 monthly.to_csv(DATA/'gpr_monthly_selected.csv',index_label='month',float_format='%.16g')
 q=quarterly(monthly);q.index=q.index.to_period('Q');q.index.name='Quarter'
 q.to_csv(DATA/'gpr_quarterly_levels.csv',float_format='%.16g')
 macro=pd.read_csv(BASE/'model_input.csv');cols=['Quarter']+[f'{c}_{v}' for c in COUNTRIES for v in ['y','dp','r','de','deq']]
 macro=macro[cols]
 if macro.Quarter.iloc[0]!='2000Q1' or macro.Quarter.iloc[-1]!='2025Q3' or len(macro)!=103:raise ValueError('Unexpected frozen macro sample')
 macro.to_csv(DATA/'macro_input.csv',index=False,float_format='%.16g')
 used=q.reindex(pd.PeriodIndex(macro.Quarter,freq='Q'))
 if used.isna().any().any() or (used<=0).any().any():raise ValueError('Cannot use natural log of a missing/nonpositive quarterly category; no epsilon or shifted start')
 np.log(used).to_csv(DATA/'gpr_quarterly_logs.csv',float_format='%.16g')
 peaks=[];manifest=[];design=[]
 for col in SERIES:
  pp=select_peaks(used[col].iloc[1:])
  if not pp:raise ValueError('No qualifying peak: '+col+'; no arbitrary date substitution')
  for row in pp:peaks.append(dict(series=col,category_name=GLOBALS.get(col,CATEGORIES.get(col)),**row,context_status='statistical peak; historical event attribution not yet verified'))
  manifest.append(dict(series=col,label=GLOBALS.get(col,CATEGORIES.get(col)),publisher_label=labels[col],source_kind='historical 3-newspaper article-share category' if col in CATEGORIES else 'recent 10 newspapers' if col in ['GPR','GPRT','GPRA'] else 'historical 3 newspapers',input_start='2000Q1',effective_start='2000Q2',common_end='2025Q3',effective_n=102,peaks=len(pp),transform='ln(quarterly mean); zero months retained; quarterly values all positive',shock='one external GPR innovation, normalized to +10% level at h=0'))
  # Independent Python preflight of actual country conditional designs.
  z=np.log(used[col]).to_numpy();base=macro.set_index('Quarter');w=pd.read_csv(HERE.parent/'Trade_Weights_14_Economies_2000_2014.csv',index_col=0).loc[COUNTRIES,COUNTRIES].to_numpy()
  if not np.allclose(w.sum(axis=1),1,atol=1e-6) or np.any(w<0) or not np.allclose(np.diag(w),0):raise ValueError('Invalid trade weights')
  w=w/w.sum(axis=1)[:,None]
  for ci,c in enumerate(COUNTRIES):
   own=base[[f'{c}_{v}' for v in ['y','dp','r','de','deq']]].to_numpy();foreign=np.column_stack([base[[f'{partner}_{v}' for partner in COUNTRIES]].to_numpy()@w[ci] for v in ['y','dp','r','de','deq']]);ex=np.column_stack([foreign,z]);X=np.column_stack([np.ones(102),ex[1:],ex[:-1],own[:-1]])
   rank=np.linalg.matrix_rank(X)
   if rank!=X.shape[1]:raise ValueError('Rank-deficient design: '+col+'/'+c)
   design.append(dict(series=col,unit=c,equations=5,n=102,regressors=X.shape[1],rank=rank,condition_number=np.linalg.cond(X)))
  X=np.column_stack([np.ones(102),z[:-1]]);design.append(dict(series=col,unit='GL_external_path',equations=1,n=102,regressors=2,rank=np.linalg.matrix_rank(X),condition_number=np.linalg.cond(X)))
 pd.DataFrame(peaks).to_csv(DATA/'peak_dates.csv',index=False,float_format='%.16g');pd.DataFrame(manifest).to_csv(DATA/'run_manifest.csv',index=False);pd.DataFrame(design).to_csv(DATA/'design_preflight.csv',index=False)
 (DATA/'peak_rules.json').write_text(json.dumps(PEAK_RULE,indent=2)+'\n')
 meta=dict(source_url='https://www.matteoiacoviello.com/gpr_files/data_gpr_export.xls',source_page='https://www.matteoiacoviello.com/gpr.htm',download_date='2026-10-26',source_page_last_monthly_update='2026-10-01',latest_downloaded_month=str(monthly.index.max().date()),latest_risk_quarter=str(q.index.max()),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),macro_snapshot='country_epu_oil/data_20261022',macro_baseline_run=39647769567,macro_input_sha256=hashlib.sha256((BASE/'model_input.csv').read_bytes()).hexdigest(),categories=14,national_endogenous_variables=70,external_path_variables=1,estimation_quarters=102,common_end='2025Q3',rate_proxy_limitations_retained=True,trade_weight_lineage_unresolved=True,no_imputation=True,scope='6 aggregate/threats/acts + 8 historical event categories; country indices and search variants inventoried separately, not silently called event categories',exogeneity='GPR is external in all 14 national VARX blocks; auxiliary TVP AR(1) receives no macro or foreign feedback; not an external-instrument causal identification')
 (DATA/'source_manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
 hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.is_file() and p.name!='snapshot_hashes.json'}
 (DATA/'snapshot_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
 print(json.dumps(dict(series=SERIES,peaks_by_series={m['series']:m['peaks'] for m in manifest},effective_sample='2000Q2--2025Q3',all_designs_full_rank=True),indent=2))

def verify():
 m=json.loads((DATA/'snapshot_hashes.json').read_text())
 for name,sha in m.items():
  if hashlib.sha256((DATA/name).read_bytes()).hexdigest()!=sha:raise ValueError('Frozen file changed: '+name)
 logs=pd.read_csv(DATA/'gpr_quarterly_logs.csv',index_col=0);levels=pd.read_csv(DATA/'gpr_quarterly_levels.csv',index_col=0).reindex(logs.index)
 if not np.allclose(logs,np.log(levels),rtol=1e-12,atol=1e-12):raise ValueError('Logged levels mismatch')
 pp=pd.read_csv(DATA/'peak_dates.csv');
 for col in SERIES:
  rebuilt=select_peaks(levels[col].iloc[1:]);actual=pp[pp.series==col]
  if [r['Quarter'] for r in rebuilt]!=actual.Quarter.tolist():raise ValueError('Peak mismatch: '+col)
 print('PASS frozen GPR data, transformation and independent own-series peak dates')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');v=a.parse_args();verify() if v.verify else build()
