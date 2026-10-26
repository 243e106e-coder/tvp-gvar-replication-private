"""One category, one own-peak date and one response variable per figure.
Nested 68/90% posterior intervals. Counts and country/horizon lists kept apart.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from build_data import SERIES
COUNTRIES='AU BR CA CH CN EA UK JP KR NO SG TR US ZA'.split()
LABELS={'y':'Real GDP growth (percentage points)','dp':'Quarterly inflation (percentage points)',
'r':'Interest-rate change (percentage points)','de':'Real effective appreciation (approx. percent)','deq':'Equity price return (approx. percent)'}
def ranges(vals):
 out=[]
 for v in sorted(set(map(int,vals))):
  if out and v==out[-1][-1]+1:out[-1].append(v)
  else:out.append([v])
 return ', '.join(str(v[0]) if len(v)==1 else f'{v[0]}-{v[-1]}' for v in out) or 'none'

def validate(z):
 if len(z)!=182 or set(z.country)!=set(COUNTRIES):raise ValueError('Missing country responses')
 if not z.groupby('country').horizon.apply(lambda x:set(x)==set(range(13))).all():raise ValueError('Missing horizons')
 if not ((z.low90<=z.low68)&(z.low68<=z['median'])&(z['median']<=z.high68)&(z.high68<=z.high90)).all():raise ValueError('Invalid nested intervals')

def main(series,root):
 out=Path(root)/series;d=pd.read_csv(out/'irf_summary.csv');d=d[d.variable.str.fullmatch(r'[A-Z]{2}_(y|dp|r|de|deq)')].copy();d['country']=d.variable.str[:2];d['macro']=d.variable.str.split('_').str[1]
 if set(d.series)!={series}:raise ValueError('Mixed category IRFs')
 pp=pd.read_csv(out/'selected_peak_dates.csv');own=set(pp.Quarter)
 if set(d.date)!=own:raise ValueError('IRFs do not use all and only own-category peak dates')
 rows=[];details=[];figdir=out/'figures';figdir.mkdir(exist_ok=True)
 with PdfPages(out/'IRF_68_90.pdf') as pdf:
  for key,z in d.groupby(['date','filter','response_kind','macro']):
   date,filt,kind,macro=key;validate(z);row=dict(series=series,date=date,filter=filt,response_kind=kind,macro=macro,total_cells=182,total_countries=14)
   for b in [68,90]:
    pos=z[f'low{b}']>0;neg=z[f'high{b}']<0;sig=pos|neg
    row.update({f'positive_cells_{b}':int(pos.sum()),f'negative_cells_{b}':int(neg.sum()),f'zero_including_cells_{b}':int((~sig).sum()),f'countries_any_{b}':len(set(z.loc[sig,'country'])),f'countries_none_{b}':14-len(set(z.loc[sig,'country']))})
   rows.append(row)
   for cc in COUNTRIES:
    a=z[z.country==cc];r=dict(series=series,date=date,filter=filt,response_kind=kind,macro=macro,country=cc)
    for b in [68,90]:
     pos=a[f'low{b}']>0;neg=a[f'high{b}']<0;r.update({f'positive_horizons_{b}':ranges(a.loc[pos,'horizon']),f'negative_horizons_{b}':ranges(a.loc[neg,'horizon']),f'zero_including_horizons_{b}':ranges(a.loc[~(pos|neg),'horizon']),f'zero_including_count_{b}':int((~(pos|neg)).sum())})
    details.append(r)
   fig,axs=plt.subplots(4,4,figsize=(14.4,8.2));fig.subplots_adjust(left=.06,right=.985,bottom=.1,top=.86,hspace=.72,wspace=.38)
   fig.suptitle(f'{series}: own peak {date} / {macro} / {filt}',fontsize=16,color='#15334E')
   label='Interest-rate LEVEL response (percentage points)' if kind=='rate_level_cumulative' else LABELS[macro]
   fig.text(.06,.913,label+'; +10% external GPR impact; dark 68%, light 90%',fontsize=10,color='#245782')
   for i,cc in enumerate(COUNTRIES):
    ax=axs.flat[i];a=z[z.country==cc].sort_values('horizon');h=a.horizon.to_numpy();ax.fill_between(h,a.low90.to_numpy(),a.high90.to_numpy(),color='#D9E8FA');ax.fill_between(h,a.low68.to_numpy(),a.high68.to_numpy(),color='#6CBBD5');ax.plot(h,a['median'],color='#234F83',lw=1.15);ax.axhline(0,color='#667789',ls='--',lw=.65);ax.set_title(cc,loc='left',fontsize=11);ax.set_xticks([0,4,8,12]);ax.tick_params(labelsize=8);ax.yaxis.get_offset_text().set_fontsize(7)
    for sp in ax.spines.values():sp.set_color('#CED8E4')
   for i in [14,15]:axs.flat[i].axis('off')
   fig.text(.57,.16,f"90%: {row['countries_any_90']}/14 countries have at least one zero-excluding horizon",fontsize=9,color='#245782')
   fig.text(.57,.13,'Exploratory rate proxies; external-path restrictions are assumed.',fontsize=9,color='#923C32')
   fig.text(.5,.035,'Quarters after shock (h=0 to 12)',ha='center',fontsize=10)
   pdf.savefig(fig);fig.savefig(figdir/f'{series}_{date}_{filt}_{kind}_{macro}.png',dpi=120);plt.close(fig)
 pd.DataFrame(rows).to_csv(out/'significance_counts.csv',index=False);pd.DataFrame(details).to_csv(out/'country_horizon_significance.csv',index=False)
 print('Plotted',series,len(rows),'pages, all own peak dates; countries and horizons reported separately')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--series',choices=SERIES,required=True);a.add_argument('--root',default='results/gpr_categories_exogenous');v=a.parse_args();main(v.series,v.root)
