import unittest
import numpy as np
import pandas as pd
from build_input import transform,COUNTRIES
from select_peaks import select

class PipelineTests(unittest.TestCase):
    def fixture(self):
        qq=pd.period_range("1999Q1","2012Q4",freq="Q").astype(str)
        raw=pd.DataFrame({"Quarter":qq});risk=pd.DataFrame({"Quarter":qq})
        for cc in COUNTRIES:
            for v in ("GDP","CPI","REER","EQ"):
                raw[f"{cc}_{v}"]=100*np.exp(.01*np.arange(len(qq)))
            raw[f"{cc}_RATE"]=-1+np.arange(len(qq))*.1
        for v in ("CN_EPU","US_EPU","WTI"):
            risk[v]=100.;risk["LN_"+v]=np.log(100.);risk[v+"_months"]=3
        return raw,risk

    def test_presample_and_units(self):
        raw,risk=self.fixture()
        for p in (1,2):
            out,a=transform(raw,risk,p)
            self.assertEqual(out.Quarter.iloc[p],"2000Q1")
            self.assertTrue(np.allclose(out.AU_y,1))
            self.assertTrue(np.allclose(out.AU_r,.1)) # negative rates are valid
            self.assertEqual(a["effective_start"],"2000Q1")

    def test_internal_gap_is_not_bridged(self):
        raw,risk=self.fixture();raw.loc[raw.Quarter=="2010Q2","US_CPI"]=np.nan
        out,a=transform(raw,risk)
        self.assertEqual(out.Quarter.iloc[-1],"2010Q1")
        self.assertEqual(a["first_gap"],"2010Q2")
        self.assertGreater(a["valid_quarters_after_gap"],0)

    def test_cannot_invent_chinese_presample(self):
        raw,risk=self.fixture();risk.loc[risk.Quarter<"2000Q1","LN_CN_EPU"]=np.nan
        with self.assertRaisesRegex(ValueError,"Missing presample"): transform(raw,risk)

    def test_shocks_have_own_dates(self):
        _,risk=self.fixture()
        risk.loc[risk.Quarter=="2001Q1","CN_EPU"]=400
        risk.loc[risk.Quarter=="2002Q2","US_EPU"]=500
        risk.loc[risk.Quarter=="2003Q3","WTI"]=600
        peaks=select(risk,risk.Quarter[risk.Quarter>="2000Q1"].tolist())
        first=peaks[peaks.period=="2000-2004"].set_index("series")
        self.assertEqual(first.loc["CN_EPU","Quarter"],"2001Q1")
        self.assertEqual(first.loc["US_EPU","Quarter"],"2002Q2")
        self.assertEqual(first.loc["WTI","Quarter"],"2003Q3")

if __name__=="__main__":unittest.main()
