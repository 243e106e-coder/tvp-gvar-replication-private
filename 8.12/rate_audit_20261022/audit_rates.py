"""Audit the frozen 14-rate candidate; never fill or replace an observation.

ADF/KPSS diagnose transforms, not measurement comparability or causality.
Run against the SAME 2000Q2--2025Q3 sample, including real negative rates.
"""
import argparse,hashlib,io,json,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.tsa.stattools import adfuller,kpss

HERE=Path(__file__).resolve().parent
SNAPSHOT=HERE.parent/'country_epu_oil'/'data_20261022'
COUNTRIES='AU BR CA CH CN EA UK JP KR NO SG TR US ZA'.split()
DEFINITIONS={c:('3-month money-market benchmark','3 months','monthly published rate') for c in COUNTRIES}
DEFINITIONS.update({
 'BR':('call-money / interbank rate','overnight','monthly published rate'),
 'JP':('call-money / interbank rate','overnight','monthly published rate'),
 'US':('call-money / interbank rate','overnight','monthly published rate'),
 'CN':('publisher quarterly repo rate','1 day','publisher quarterly series'),
 'SG':('government Treasury-bill yield','1 year','quarterly mean of three monthly end-period yields'),
 'TR':('central-bank discount rate','administered; not a traded maturity','monthly published rate')})
EXCEPTIONS={'CN':'Different maturity and source construction from the 3-month group.',
'US':'Different maturity from the 3-month group; this is NOT the incomplete US 3-month series.',
'SG':'Government yield, 1-year maturity and monthly endpoints. Not SORA or SIBOR.',
'TR':'Administered discount rate. Not the overnight market/funding rate; high-priority unresolved definition.'}

def run_tests(raw,start='2000Q2',end='2025Q3'):
 rows=[]
 for c in COUNTRIES:
  for form in ['level','difference']:
   z=pd.to_numeric(raw[c+'_RATE'],errors='raise');z=z if form=='level' else z.diff()
   x=z.loc[start:end]
   if len(x)!=102 or x.isna().any() or not np.isfinite(x).all():raise ValueError('Common test sample is incomplete: '+c)
   for regression,maxlag in [('c',4),('c',8),('ct',4)]:
    with warnings.catch_warnings(record=True) as caught:
     warnings.simplefilter('always')
     ad=adfuller(x,maxlag=maxlag,regression=regression,autolag='AIC')
     kp=kpss(x,regression=regression,nlags='auto')
    rows.append(dict(country=c,form=form,sample_start=start,sample_end=end,sample_n=len(x),
     deterministic=regression,ADF_maxlag=maxlag,ADF_lag=ad[2],ADF_nobs=ad[3],ADF_stat=ad[0],ADF_p=ad[1],
     ADF_reject_unit_root_5pct=bool(ad[1]<.05),KPSS_stat=kp[0],KPSS_p=kp[1],KPSS_bandwidth=kp[2],
     KPSS_reject_stationarity_5pct=bool(kp[1]<.05),KPSS_p_is_bound=bool(kp[1] in [.01,.1]),
     warnings='; '.join(str(w.message) for w in caught)))
 return pd.DataFrame(rows)

def fred_table(path):
 a=pd.read_csv(path)
 if list(a.columns)!=['DATE','VALUE']:raise ValueError('Invalid selected observation schema')
 a['period']=pd.PeriodIndex(pd.to_datetime(a.DATE),freq='M');a['value']=pd.to_numeric(a.VALUE,errors='coerce')
 if a.period.duplicated().any():raise ValueError('Duplicate FRED observation')
 return a.set_index('period').value

def verify_sources(out):
 monthly=[];receipts=[]
 for c,filename,sid in [('CN','fred_CN_3m.csv','IR3TIB01CNM156N'),('US','fred_US_3m.csv','IR3TIB01USM156N'),('TR','fred_TR_overnight.csv','IRSTCI01TRM156N')]:
  p=HERE/'sources'/filename;s=fred_table(p);s=s.reindex(pd.period_range('1999-01','2026-08' if c!='US' else '2026-09',freq='M'))
  for date,v in s.items():monthly.append(dict(country=c,period=str(date),value=v,source_id=sid,source_url='https://fred.stlouisfed.org/data/'+sid))
  receipts.append(dict(file=filename,url='https://fred.stlouisfed.org/data/'+sid,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 observed=pd.DataFrame(monthly);observed.to_csv(out/'independent_source_checks.csv',index=False)
 observed[observed.value.isna()].to_csv(out/'confirmed_missing_months.csv',index=False)
 # The CBRT table is an effective-date schedule. Carrying a legally effective
 # policy rate is meaningful here ONLY; it is not a rule for missing market data.
 p=HERE/'sources/cbrt_overnight.csv';a=pd.read_csv(p)
 a['date']=pd.to_datetime(a.effective_date,format='%d.%m.%y');a['rate']=pd.to_numeric(a.borrowing);a=a[a.date<=pd.Timestamp('2026-10-22')].sort_values('date')
 a[['date','rate']].to_csv(out/'cbrt_effective_schedule.csv',index=False)
 months=pd.date_range('2002-02-28','2026-08-31',freq='ME');b=pd.merge_asof(pd.DataFrame({'date':months}),a[['date','rate']],on='date')
 b['period']=b.date.dt.to_period('M');original=fred_table(HERE/'sources/fred_TR_overnight.csv')
 b['OECD_published']=b.period.map(original);b['difference']=b.OECD_published-b.rate
 b.to_csv(out/'rejected_CBRT_month_end_substitution.csv',index=False)
 common=b.dropna(subset=['difference']);mismatch=int((common.difference.abs()>1e-6).sum())
 if mismatch==0:raise ValueError('Unexpected validation change; manually reassess CBRT substitute')
 receipts.append(dict(file=p.name,url='https://www.tcmb.gov.tr/wps/wcm/connect/EN/TCMB%2BEN/Main%2BMenu/Core%2BFunctions/Monetary%2BPolicy/Central%2BBank%2BInterest%2BRates/CBRT%2BINTEREST%2BRATES',sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 (out/'independent_source_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
 return mismatch

def main(out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 raw=pd.read_csv(SNAPSHOT/'raw_macro.csv').set_index('Quarter');raw.index=pd.PeriodIndex(raw.index,freq='Q')
 input=pd.read_csv(SNAPSHOT/'model_input.csv').set_index('Quarter');input.index=pd.PeriodIndex(input.index,freq='Q')
 if len(input)!=103 or str(input.index[0])!='2000Q1' or str(input.index[-1])!='2025Q3':raise ValueError('Unexpected frozen sample')
 inv=pd.read_csv(SNAPSHOT/'verified_macro_inventory.csv');obs=pd.read_csv(SNAPSHOT/'source_observations.csv.gz');regs=[];des=[]
 for c in COUNTRIES:
  s=raw[c+'_RATE'];expected=s.diff().reindex(input.index)
  if not np.allclose(input[c+'_r'],expected,rtol=1e-10,atol=1e-10):raise ValueError('Rate transformation mismatch: '+c)
  iv=inv[(inv.country==c)&(inv.variable=='RATE')].iloc[0];definition,tenor,aggregation=DEFINITIONS[c]
  regs.append(dict(country=c,source_id=iv.source_id,definition=definition,tenor=tenor,source_aggregation=aggregation,
   input_form='first difference in annual-percentage-rate points; never log',three_month_reference_group=c not in set(EXCEPTIONS)|{'BR','JP'},proxy_changed_from_inherited_spec=c in EXCEPTIONS,
   unresolved_definition=c in EXCEPTIONS,issue=EXCEPTIONS.get(c,''),first=iv['first'],last=iv['last']))
  for form,x in [('level',s),('difference',s.diff())]:
   z=x.loc['2000Q2':'2025Q3'];imax=z.abs().idxmax();des.append(dict(country=c,form=form,n=len(z),negative=int((z<0).sum()),mean=z.mean(),std=z.std(),minimum=z.min(),maximum=z.max(),largest_absolute_quarter=str(imax),largest_absolute_value=z.loc[imax]))
 pd.DataFrame(regs).to_csv(out/'rate_definition_registry.csv',index=False);pd.DataFrame(des).to_csv(out/'rate_descriptives.csv',index=False)
 tests=run_tests(raw);tests.to_csv(out/'rate_level_difference_tests.csv',index=False)
 bad=verify_sources(out)
 primary=tests[(tests.deterministic=='c')&(tests.ADF_maxlag==4)]
 info=dict(source_run=39647769567,sample_start='2000Q2',sample_end='2025Q3',sample_n=102,all_14_frozen_rate_differences_verified=True,
  chosen_default='retain first difference; no automatic variable or source replacement',homogeneous_14_country_rate_ready=False,
  unresolved_countries=list(EXCEPTIONS),CBRT_month_end_substitution_rejected_overlap_mismatches=bad,
  tests={f:dict(ADF_rejects=int(g.ADF_reject_unit_root_5pct.sum()),KPSS_rejects=int(g.KPSS_reject_stationarity_5pct.sum())) for f,g in primary.groupby('form')},
  important='Source comparability remains unresolved. Stationarity screens do not validate the proxies, identify policy shocks, or prove convergence.')
 (out/'DECISION.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',default=str(HERE/'results'));v=a.parse_args();main(v.out)
