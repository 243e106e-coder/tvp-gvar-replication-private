import tempfile,unittest
import pandas as pd
from refresh_macro_panel import Builder

class AggregationTests(unittest.TestCase):
    def test_missing_month_is_not_a_quarter(self):
        with tempfile.TemporaryDirectory() as t:
            b=Builder(t)
            b.add("US","CPI",pd.DataFrame({"period":["2025-10","2025-12"],"value":[100.,102.]}),"test","index","M")
            self.assertTrue(pd.isna(b.series["US_CPI"].loc[pd.Period("2025Q4")]))
    def test_native_quarter_does_not_require_synthetic_months(self):
        with tempfile.TemporaryDirectory() as t:
            b=Builder(t)
            b.add("AU","CPI",pd.DataFrame({"period":["2000-Q1","2000-Q2"],"value":[100.,101.]}),"test","index","Q")
            self.assertEqual(b.series["AU_CPI"].loc[pd.Period("2000Q1")],100.)
    def test_negative_rate_is_an_observation(self):
        with tempfile.TemporaryDirectory() as t:
            b=Builder(t)
            b.add("CH","RATE",pd.DataFrame({"period":["2020-01","2020-02","2020-03"],"value":[-.5,-.4,-.3]}),"test","percent","M")
            self.assertAlmostEqual(b.series["CH_RATE"].loc[pd.Period("2020Q1")],-.4)

if __name__=="__main__":unittest.main()
