"""Build a reproducible 70-series candidate from full single-source histories.

No interpolation, zero replacement, cross-series gap filling or old/new splicing.
The four explicitly disclosed alternative interest-rate proxies are a new
exploratory specification, NOT a claim to have completed the old rate measures.
"""
import argparse, hashlib, io, json, re, zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import quote
import numpy as np
import pandas as pd

COUNTRIES="AU BR CA CH CN EA UK JP KR NO SG TR US ZA".split()
AREAS=dict(zip(COUNTRIES,"AUS BRA CAN CHE CHN EA20 GBR JPN KOR NOR SGP TUR USA ZAF".split()))
OECD="https://sdmx.oecd.org/public/rest/v1/data/"
GDP_URL=OECD+"OECD.SDD.NAD,DSD_NAMAIN1@DF_QNA_EXPENDITURE_INDICES,1.1/Q.Y.AUS+BRA+CAN+CHE+EA+GBR+JPN+KOR+NOR+TUR+USA+ZAF.S1.S1.B1GQ._Z._Z._Z.IX.LR.N.T0102?startPeriod=1999-Q1&format=csvfile"
FIN_URL=OECD+"OECD.SDD.STES,DSD_STES@DF_FINMARK,4.0/AUS+BRA+CAN+CHE+CHN+EA20+GBR+JPN+KOR+NOR+TUR+USA+ZAF.M+Q.SHARE+IRSTCI+IR3TIB......?startPeriod=1999-01&format=csvfile"
OLD_CPI_URL=OECD+"OECD.SDD.TPS,DSD_PRICES@DF_PRICES_ALL,1.0/CHN+USA+JPN+KOR+SGP+GBR+EA19+CHE+CAN+AUS+NOR+BRA+TUR+ZAF.M.N.CPI.IX._T.N._Z?startPeriod=1999-01&format=csvfile"
NEW_CPI_URL=OECD+"OECD.SDD.TPS,DSD_PRICES_COICOP2018@DF_PRICES_C2018_ALL,1.0/AUS+BRA+CAN+CHE+CHN+EA20+GBR+JPN+KOR+NOR+TUR+USA+ZAF.M+Q.N.CPI.IX._T.N._Z?startPeriod=1999-01&format=csvfile"
CHINA_URL="https://www.atlantafed.org/-/media/Project/Atlanta/FRBA/Documents/research-and-data/data/china-macroeconomy/data/2026-09-data.zip"
IMF_URL="https://api.imf.org/external/sdmx/2.1/data/IMF.STA,MFS_IR,9.0.0/TUR+SGP+CHN+USA.MMRT_RT_PT_A_PT+IBRT_S3M_RT_PT_A_PT+DISR_RT_PT_A_PT.M?startPeriod=1999-01"
EA_CPI_URL="https://data-api.ecb.europa.eu/service/data/ICP/M.U2.N.000000.4.INX?startPeriod=1999-01&format=csvdata"
EA_EQ_URL="https://data-api.ecb.europa.eu/service/data/FM/M.U2.EUR.DS.EI.DJES50I.HSTA?startPeriod=1999-01&format=csvdata"
SG_IDS={"sg_gdp":"d_8ab3c73416b81c8a20ac8a93a1c80a40","sg_cpi":"d_bdaff844e3ef89d39fceb962ff8f0791","sg_rate":"d_5fe5a4bb4a1ecc4d8a56a095832e2b24"}
YAHOO={"AU":"^AXJO","CA":"^GSPTSE","CH":"^SSMI","CN":"000001.SS","UK":"^FTSE","JP":"^N225","KR":"^KS11","SG":"^STI","US":"^GSPC"}
RATE_CHANGES={
 "CN":("3-month interbank rate (missing 2006-01)","1-day repo, Atlanta Fed/Higgins-Zha; fraction multiplied by 100"),
 "US":("Inherited short-rate series; 3-month OECD comparator missing 2020-04","OECD immediate/call-money overnight rate; full single history"),
 "SG":("Inherited rate lacks a single-source SORA history before 2005","MAS/SingStat 1-year Treasury-bill yield, monthly end-period; full single history"),
 "TR":("OECD overnight interbank rate; nine missing months in 2023","IMF discount rate, percent per annum; full single history, not a fill for interbank rates")}

class Builder:
    def __init__(self,root,reuse=False):
        self.root=Path(root);self.raw=self.root/"sources_full";self.raw.mkdir(parents=True,exist_ok=True)
        self.reuse=reuse;self.receipts=[];self.series={};self.lineage=[];self.inventory=[]
    def fetch(self,name,url):
        p=self.raw/name
        if self.reuse and p.exists(): b=p.read_bytes();when=datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat()
        else:
            b=urlopen(Request(url,headers={"User-Agent":"Mozilla/5.0","Accept":"text/csv" if 'api.imf' in url else "*/*","Cache-Control":"no-cache"}),timeout=60).read()
            if b.lstrip().startswith((b"<!doctype",b"<!DOCTYPE")):raise ValueError("HTML response instead of data: "+url)
            p.write_bytes(b);when=datetime.now(timezone.utc).isoformat()
        sha=hashlib.sha256(b).hexdigest()
        self.receipts.append(dict(file=name,url=url,sha256=sha,bytes=len(b),retrieved_utc=when))
        return b
    def gov(self,name,ds):
        url=f"https://api-open.data.gov.sg/v1/public/api/datasets/{ds}/poll-download"
        p=self.raw/(name+".csv")
        if self.reuse and p.exists():return pd.read_csv(io.BytesIO(self.fetch(name+".csv",url)))
        m=json.loads(urlopen(Request(url,headers={"Cache-Control":"no-cache","User-Agent":"Mozilla/5.0"}),timeout=45).read())
        if m.get("data",{}).get("status")!="DOWNLOAD_SUCCESS":raise ValueError("Government download not ready: "+ds)
        b=urlopen(m["data"]["url"],timeout=45).read();p.write_bytes(b)
        # Persist the stable API URL, never the expiring signed URL/token.
        return pd.read_csv(io.BytesIO(self.fetch(name+".csv",url))) if self.reuse else self.gov_record(name,url,b)
    def gov_record(self,name,url,b):
        self.receipts.append(dict(file=name+".csv",url=url,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),retrieved_utc=datetime.now(timezone.utc).isoformat()))
        return pd.read_csv(io.BytesIO(b))
    def add(self,cc,var,d,source_id,unit,frequency="M",note=""):
        """Keep source observations and quarterly counts; missing stays missing."""
        d=d[["period","value"]].copy();d["value"]=pd.to_numeric(d.value,errors="coerce")
        if frequency=="Q":
            d["period"]=d.period.astype(str).str.replace("-Q","Q",regex=False)
            idx=pd.PeriodIndex(d.period,freq="Q")
        else:
            d["period"]=pd.to_datetime(d.period,errors="raise")
            idx=pd.PeriodIndex(d.period,freq="M" if frequency=="M" else "D")
        d.index=idx;d=d.sort_index();d=d[d.index>=pd.Period("1999-01-01",freq=idx.freq)]
        # Do not include the current, incomplete quarter even for daily markets.
        last_complete=pd.Timestamp.now(tz="UTC").tz_localize(None).to_period("Q")-1
        d=d[d.index.asfreq("Q")<=last_complete]
        if d.index.duplicated().any():raise ValueError("Duplicate observations inside completed quarters: "+source_id)
        key=cc+"_"+var
        if frequency=="Q":
            q=d.value.copy()
        else:
            # Reindex monthly sources before averaging; all three source months
            # must have real observations. Daily markets use all supplied daily
            # closes and require observations in each of the three months.
            d["quarter"]=d.index.asfreq("Q")
            q=d.groupby("quarter").value.mean()
            if frequency=="M":
                count=d.groupby("quarter").value.count()
            else:
                count=d[d.value.notna()].groupby("quarter").apply(lambda z:len(z.index.asfreq("M").unique()),include_groups=False)
            q=q.where(count==3)
        if var!="RATE" and (q.dropna()<=0).any():raise ValueError("Nonpositive level: "+key)
        q.index=pd.PeriodIndex(q.index,freq="Q");self.series[key]=q
        valid=q.dropna();self.inventory.append(dict(country=cc,variable=var,source_id=source_id,unit=unit,native_frequency=frequency,
          first=str(valid.index.min()),last=str(valid.index.max()),valid_quarters=len(valid),
          rate_proxy_change=cc in RATE_CHANGES and var=="RATE",note=note))
        for i,v in d.value.items():self.lineage.append(dict(country=cc,variable=var,period=str(i),value=None if pd.isna(v) else float(v),source_id=source_id,unit=unit))
    def build(self):
        gdp=pd.read_csv(io.BytesIO(self.fetch("oecd_gdp.csv",GDP_URL)))
        assert (gdp.ADJUSTMENT=="Y").all() and (gdp.PRICE_BASE=="LR").all() and (gdp.TRANSACTION=="B1GQ").all()
        for cc in [c for c in COUNTRIES if c not in ("CN","SG")]:
            a="EA" if cc=="EA" else AREAS[cc];d=gdp[gdp.REF_AREA==a]
            self.add(cc,"GDP",d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),"OECD_QNA_LR_"+a,"2020=100; calendar and seasonally adjusted","Q")
        china=zipfile.ZipFile(io.BytesIO(self.fetch("china_macro_202609.zip",CHINA_URL)))
        ch=pd.read_csv(china.open("outdata2609_hz_quarterly.csv"))
        qq=[str(int(x))+"Q"+str(int(round((x-int(x))*4))+1) for x in ch.q_dates_data]
        for var,values,unit in [("GDP",np.exp(ch.logrealGDP_va),"real GDP by value added; exp(publisher logrealGDP_va), RMB billion"),
                                ("RATE",100*ch.R1dRepo,"percent per annum; source fractions multiplied by 100")]:
            self.add("CN",var,pd.DataFrame(dict(period=qq,value=values)),"AtlantaFed_202609_"+("logrealGDP_va" if var=="GDP" else "R1dRepo"),unit,"Q",
              "Publisher-constructed/seasonally adjusted quarterly series; no additional interpolation or growth-rate extension. See publisher readme.")
        for name,ds in SG_IDS.items():
            d=self.gov(name,ds);d.iloc[:,0]=d.iloc[:,0].astype(str).str.strip()
            wanted={"sg_gdp":"GDP In Chained (2015) Dollars","sg_cpi":"All Items","sg_rate":"Government Securities - 1-Year Treasury Bills Yield"}[name]
            row=d[d.iloc[:,0]==wanted]
            if len(row)!=1:raise ValueError("Ambiguous SingStat row: "+wanted)
            row=row.iloc[0,1:]
            if name=="sg_gdp":
                periods=[re.sub(r"^(\d{4})([1-4])Q$",r"\1Q\2",str(x)) for x in row.index];var="GDP";freq="Q";unit="million SGD, chained 2015 dollars, seasonally adjusted"
            else:
                periods=pd.to_datetime(row.index,format="%Y%b");var="CPI" if name=="sg_cpi" else "RATE";freq="M";unit="2024=100, all items, NSA" if var=="CPI" else "percent per annum; monthly end-period observations"
            self.add("SG",var,pd.DataFrame(dict(period=periods,value=row.values)),"SingStat_"+ds,unit,freq)
        old=pd.read_csv(io.BytesIO(self.fetch("oecd_cpi_1999.csv",OLD_CPI_URL)))
        new=pd.read_csv(io.BytesIO(self.fetch("oecd_cpi_2018.csv",NEW_CPI_URL)))
        for cc in [c for c in COUNTRIES if c not in ("SG","EA")]:
            use_new=cc in ("AU","CA","CH","JP","NO","TR","ZA")
            d=(new if use_new else old);freq="Q" if cc=="AU" else "M"
            d=d[(d.REF_AREA==AREAS[cc])&(d.FREQ==freq)]
            self.add(cc,"CPI",d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),
                "OECD_CPI_"+("COICOP2018_" if use_new else "COICOP1999_")+AREAS[cc],
                "All-items price index, NSA; base period "+str(d.BASE_PER.unique().tolist()),freq,
                "Entire revised history from one classification; no old/new classification splice.")
        for var,name,url in [("CPI","ecb_hicp.csv",EA_CPI_URL),("EQ","ecb_eurostoxx50.csv",EA_EQ_URL)]:
            d=pd.read_csv(io.BytesIO(self.fetch(name,url)))
            self.add("EA",var,d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),
              "ECB_"+("ICP.M.U2.N.000000.4.INX" if var=="CPI" else "FM.M.U2.EUR.DS.EI.DJES50I.HSTA"),
              "Euro-area HICP overall index, changing composition, NSA" if var=="CPI" else "EURO STOXX 50 price index; published monthly averages","M")
        fin=pd.read_csv(io.BytesIO(self.fetch("oecd_financial.csv",FIN_URL)))
        rate_measure={cc:("IRSTCI" if cc in ("BR","JP","US") else "IR3TIB") for cc in COUNTRIES if cc not in ("CN","SG","TR")}
        for cc,m in rate_measure.items():
            d=fin[(fin.REF_AREA==AREAS[cc])&(fin.FREQ=="M")&(fin.MEASURE==m)]
            self.add(cc,"RATE",d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),"OECD_FINMARK_"+AREAS[cc]+"_"+m,"percent per annum","M")
        imf=pd.read_csv(io.BytesIO(self.fetch("imf_interest.csv",IMF_URL)))
        d=imf[(imf.COUNTRY=="TUR")&(imf.INDICATOR=="DISR_RT_PT_A_PT")].copy()
        d["period"]=d.TIME_PERIOD.str.replace("-M","-",regex=False);d["value"]=d.OBS_VALUE
        self.add("TR","RATE",d,"IMF_MFS_IR_TUR_DISR_RT_PT_A_PT_M","discount rate, percent per annum","M","Entire discount-rate series; never used to fill missing interbank observations.")
        for cc in ("BR","NO","TR","ZA"):
            d=fin[(fin.REF_AREA==AREAS[cc])&(fin.FREQ=="M")&(fin.MEASURE=="SHARE")]
            self.add(cc,"EQ",d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),"OECD_FINMARK_"+AREAS[cc]+"_SHARE","price index; source base "+str(d.BASE_PER.unique().tolist()),"M")
        for cc,ticker in YAHOO.items():
            url=f"https://query1.finance.yahoo.com/v8/finance/chart/{quote(ticker,safe='')}?period1=915148800&period2=1792627200&interval=1d"
            d=json.loads(self.fetch("yahoo_"+cc+".json",url))["chart"]["result"][0]
            dates=pd.to_datetime(d["timestamp"],unit="s",utc=True).tz_convert(d["meta"]["exchangeTimezoneName"]).tz_localize(None).normalize()
            self.add(cc,"EQ",pd.DataFrame(dict(period=dates,value=d["indicators"]["quote"][0]["close"])),"Yahoo_"+ticker,"price-index daily unadjusted close; mean of supplied daily closes per complete quarter","D")
        for cc in COUNTRIES:
            area={"EA":"XM","UK":"GB"}.get(cc,cc)
            url=f"https://stats.bis.org/api/v1/data/WS_EER/M.R.B.{area}?startPeriod=1999-01&format=csv"
            d=pd.read_csv(io.BytesIO(self.fetch("bis_reer_"+cc+".csv",url)))
            assert (d.FREQ=="M").all() and (d.EER_TYPE=="R").all() and (d.EER_BASKET=="B").all()
            self.add(cc,"REER",d.rename(columns={"TIME_PERIOD":"period","OBS_VALUE":"value"}),"BIS_WS_EER_M.R.B."+area,"real broad effective exchange rate, NSA; 2020=100","M")
        assert len(self.series)==70 and len(self.inventory)==70
        periods=pd.period_range("1999Q1",max(s.index.max() for s in self.series.values()),freq="Q")
        raw=pd.DataFrame({k:s.reindex(periods) for k,s in self.series.items()},index=periods)
        raw.index=raw.index.astype(str);raw.index.name="Quarter"
        raw.to_csv(self.root/"raw_macro.csv",float_format="%.15g")
        pd.DataFrame(self.inventory).to_csv(self.root/"verified_macro_inventory.csv",index=False)
        pd.DataFrame(self.lineage).to_csv(self.root/"source_observations.csv",index=False,float_format="%.15g")
        (self.root/"source_receipts.json").write_text(json.dumps(self.receipts,indent=2)+"\n")
        pd.DataFrame([dict(country=c,old_measure=o,candidate_measure=n,selection_reason="Use a complete single observed history; no cross-measure gap filling; exploratory proxy change") for c,(o,n) in RATE_CHANGES.items()]).to_csv(self.root/"rate_proxy_changes.csv",index=False)
        missing=[]
        for col in raw:
            for q in raw.index[raw[col].isna()]:missing.append(dict(Quarter=q,variable=col))
        pd.DataFrame(missing).to_csv(self.root/"raw_missing_observations.csv",index=False)
        print("Built 70 full-source candidate series. Four interest-rate proxy changes disclosed.")

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--root",default="data/country_epu_oil");a.add_argument("--reuse-local-cache",action="store_true")
    v=a.parse_args();Builder(v.root,v.reuse_local_cache).build()
