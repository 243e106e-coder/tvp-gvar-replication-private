"""Plot already-computed draw-wise cumulative rate responses. No MCMC or
cumulative sums of marginal quantiles/medians are performed by this module.
"""
import argparse,json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
HERE=Path(__file__).resolve().parent
COUNTRIES='AU BR CA CH CN EA UK JP KR NO SG TR US ZA'.split()
KEYS=['epu_country','shock','date','filter','variable','horizon']
BANDS=['median','low68','high68','low90','high90']
def validate(a):
 if a.duplicated(KEYS+['response_kind']).any():raise ValueError('Duplicate rate summary')
 if set(a.response_kind)!= {'raw','rate_level_cumulative'}:raise ValueError('Missing raw or cumulative summaries')
 if not np.isfinite(a[BANDS]).all().all():raise ValueError('Nonfinite intervals')
 if not ((a.low90<=a.low68)&(a.low68<=a['median'])&(a['median']<=a.high68)&(a.high68<=a.high90)).all():raise ValueError('Invalid interval ordering')
 for key,g in a.groupby(['epu_country','shock','date','filter','response_kind']):
  if set(g.variable)!={c+'_r' for c in COUNTRIES} or not g.groupby('variable').horizon.apply(lambda x:set(x)==set(range(13))).all():raise ValueError('Incomplete countries or horizons: '+str(key))
 zero=a[a.horizon==0];raw=zero[zero.response_kind=='raw'].set_index(KEYS).sort_index();cum=zero[zero.response_kind=='rate_level_cumulative'].set_index(KEYS).sort_index()
 if not raw.index.equals(cum.index) or not np.allclose(raw[BANDS],cum[BANDS],rtol=1e-10,atol=1e-12):raise ValueError('Cumulative h=0 differs from original rate-change impact')
 return a

def summary(a):
 rows=[]
 for key,g in a.groupby(['epu_country','shock','date','filter','response_kind']):
  d=dict(zip(['model','shock','date','filter','response_kind'],key))
  for band in [68,90]:
   pos=g[f'low{band}']>0;neg=g[f'high{band}']<0
   d.update({f'positive_cells_{band}':int(pos.sum()),f'negative_cells_{band}':int(neg.sum()),f'zero_including_cells_{band}':int((~(pos|neg)).sum()),f'countries_any_{band}':len(set(g.loc[pos|neg,'variable'])),f'countries_none_{band}':14-len(set(g.loc[pos|neg,'variable']))})
  rows.append(d)
 return pd.DataFrame(rows)

def make_fig(g,title):
 fig,axs=plt.subplots(4,4,figsize=(14,8));fig.subplots_adjust(left=.06,right=.985,bottom=.09,top=.88,hspace=.72,wspace=.4)
 kind=g.response_kind.iloc[0];level=kind=='rate_level_cumulative'
 fig.suptitle(title,fontsize=17,color='#15334E');fig.text(.06,.925,('Cumulative response of interest-rate LEVEL' if level else 'Response of interest-rate CHANGE')+'; annual-rate percentage points; dark 68%, light 90%',fontsize=11,color='#245782')
 for i,c in enumerate(COUNTRIES):
  ax=axs.flat[i];z=g[g.variable==c+'_r'].sort_values('horizon');h=z.horizon.to_numpy();ax.fill_between(h,z.low90.to_numpy(),z.high90.to_numpy(),color='#D9E8FA');ax.fill_between(h,z.low68.to_numpy(),z.high68.to_numpy(),color='#6CBBD5');ax.plot(h,z['median'],color='#234F83',lw=1.2);ax.axhline(0,color='#667789',lw=.7,ls='--');ax.set_title(c,loc='left',fontsize=11);ax.set_xticks([0,4,8,12]);ax.tick_params(labelsize=8);ax.yaxis.get_offset_text().set_fontsize(7)
  for spine in ax.spines.values():spine.set_color('#CED8E4')
 for i in [14,15]:axs.flat[i].axis('off')
 fig.text(.6,.17,'Mixed national rate definitions remain unresolved.',fontsize=10,color='#923C32');fig.text(.6,.14,'SG: 1-year T-bill; TR: administered discount rate.',fontsize=10,color='#923C32');fig.text(.5,.025,'Quarters after a +10% own-impact shock (h=0 to 12)',ha='center',fontsize=11)
 return fig

def main(out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True);a=validate(pd.read_csv(HERE/'rate_irf_original_summaries.csv.gz'));summary(a).to_csv(out/'rate_response_significance.csv',index=False)
 with PdfPages(out/'Rate_Level_and_Change_68_90.pdf') as pdf:
  for key,g in a.groupby(['epu_country','shock','date','filter','response_kind']):
   title=f'{key[0]} EPU + WTI / {key[1]} shock / {key[2]} / {key[3]}'
   fig=make_fig(g,title);pdf.savefig(fig);name='_'.join(key)+'.png';fig.savefig(out/name,dpi=125);plt.close(fig)
 (out/'source.json').write_text(json.dumps(dict(run_id=39647769567,input='rate_irf_original_summaries.csv.gz',MCMC_rerun=False,cumulation='performed within each original posterior draw before marginal quantiles, in country_epu_oil/run.R at source commit ac8fcda87a3c34bd557ae401590cff7f622109a7',main_filter='valid_all; stable_only is a conditional-posterior sensitivity, not a convergence correction'),indent=2)+'\n')
 print('Validated and plotted 19 shock dates x 2 filters x 2 response definitions = 76 pages')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',default=str(HERE/'plots'));v=a.parse_args();main(v.out)
