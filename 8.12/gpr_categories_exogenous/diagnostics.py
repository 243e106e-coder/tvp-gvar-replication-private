"""Descriptive screens on the exact frozen 102-quarter estimation sample."""
import warnings
import pandas as pd
from statsmodels.tsa.stattools import adfuller,kpss
from build_data import DATA,SERIES
x=pd.read_csv(DATA/'macro_input.csv').set_index('Quarter').iloc[1:]
z=pd.read_csv(DATA/'gpr_quarterly_logs.csv').set_index('Quarter').iloc[1:]
x=pd.concat([x,z[SERIES]],axis=1);rows=[]
for col in x:
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');a=adfuller(x[col],maxlag=4,regression='c',autolag='AIC');k=kpss(x[col],regression='c',nlags='auto')
 rows.append(dict(variable=col,n=len(x),ADF_lag=a[2],ADF_stat=a[0],ADF_p=a[1],KPSS_stat=k[0],KPSS_p=k[1],KPSS_bandwidth=k[2],KPSS_bound=k[1] in [.01,.1],note='Descriptive screen, not automatic differencing or proof of exogeneity'))
out=DATA.parents[2]/'results/gpr_categories_exogenous/preflight'
out.mkdir(parents=True,exist_ok=True)
pd.DataFrame(rows).to_csv(out/'unit_root_screens.csv',index=False)
print('84 variable screens saved; no automatic transformation changes')
