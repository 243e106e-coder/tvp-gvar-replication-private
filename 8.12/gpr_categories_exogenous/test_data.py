import unittest
import numpy as np
import pandas as pd
from build_data import quarterly,select_peaks,SERIES,DATA,verify
class CategoryTests(unittest.TestCase):
 def test_zero_month_is_observation_not_missing(self):
  a=pd.DataFrame({'x':[0.,2.,4.]},index=pd.date_range('2000-01-01',periods=3,freq='MS'));self.assertEqual(quarterly(a).iloc[0,0],2.)
  a.iloc[0,0]=np.nan;self.assertTrue(pd.isna(quarterly(a).iloc[0,0]))
 def test_peaks_follow_own_series_and_exclude_endpoint(self):
  idx=pd.period_range('2000Q2',periods=24,freq='Q');x=pd.Series(1.,index=idx);y=x.copy();x.iloc[5]=10;y.iloc[13]=12;x.iloc[-1]=50
  self.assertEqual([p['Quarter'] for p in select_peaks(x)],[str(idx[5])]);self.assertEqual([p['Quarter'] for p in select_peaks(y)],[str(idx[13])])
 def test_plateau_has_one_first_date(self):
  x=pd.Series(1.,index=pd.period_range('2000Q2',periods=30,freq='Q'));x.iloc[10:13]=10
  self.assertEqual([p['Quarter'] for p in select_peaks(x)],[str(x.index[10])])
 def test_frozen_data_and_scope(self):
  verify();m=pd.read_csv(DATA/'run_manifest.csv');self.assertEqual(set(m.series),set(SERIES));self.assertEqual(len(m),14)
  q=pd.read_csv(DATA/'gpr_quarterly_logs.csv');self.assertEqual(len(q),103);self.assertEqual(q.Quarter.iloc[-1],'2025Q3')
if __name__=='__main__':unittest.main()
