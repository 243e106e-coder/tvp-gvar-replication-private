import unittest
import pandas as pd
from audit_rates import HERE,fred_table
from plot_rate_irfs import validate,summary
class RateAuditTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.a=pd.read_csv(HERE/'rate_irf_original_summaries.csv.gz')
 def test_true_source_gaps_are_not_treated_as_zero(self):
  for f,q in [('fred_CN_3m.csv','2006-01'),('fred_US_3m.csv','2020-04'),('fred_TR_overnight.csv','2023-04')]:self.assertTrue(pd.isna(fred_table(HERE/'sources'/f).loc[q]))
 def test_rejects_missing_country(self):
  with self.assertRaisesRegex(ValueError,'Incomplete'):validate(self.a[self.a.variable!='SG_r'])
 def test_rejects_fake_cumulative_impact(self):
  a=self.a.copy();i=a.index[(a.horizon==0)&(a.response_kind=='rate_level_cumulative')][0];a.loc[i,'median']+=1e-10
  # A material false impact must be rejected even when the band ordering holds.
  a.loc[i,'median']=self.a.loc[i,'median']+1e-3
  with self.assertRaises(ValueError):validate(a)
 def test_complete_original_results_and_counts(self):
  z=summary(validate(self.a))
  self.assertEqual(len(z),76)
  for b in [68,90]:
   self.assertTrue((z[f'positive_cells_{b}']+z[f'negative_cells_{b}']+z[f'zero_including_cells_{b}']==182).all())
   self.assertTrue((z[f'countries_any_{b}']+z[f'countries_none_{b}']==14).all())
if __name__=='__main__':unittest.main()
