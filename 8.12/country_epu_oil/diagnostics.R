# Descriptive unit-root screens, no automatic model selection by p-value.
suppressPackageStartupMessages(library(tseries))
d<-read.csv("data/country_epu_oil/model_input.csv",check.names=FALSE)
risk<-read.csv("data/country_epu_oil/risk_quarterly.csv",check.names=FALSE)
ix<-match(d$Quarter,risk$Quarter)
for(v in c("LN_CN_EPU","LN_US_EPU","LN_WTI")) d[[v]]<-risk[[v]][ix]
# Screens use the actual estimation window, not the initial-lag observation.
d<-d[-1,,drop=FALSE]
rows<-lapply(names(d)[-1],function(v) {
  x<-d[[v]];stopifnot(all(is.finite(x)))
  ad<-tryCatch(suppressWarnings(adf.test(x)),error=function(e)e)
  kp<-tryCatch(suppressWarnings(kpss.test(x,null="Level")),error=function(e)e)
  data.frame(variable=v,n=length(x),sd=sd(x),min=min(x),max=max(x),
    adf_p=if(inherits(ad,"htest")) ad$p.value else NA_real_,
    kpss_p=if(inherits(kp,"htest")) kp$p.value else NA_real_,
    adf_status=if(inherits(ad,"htest")) "computed; tseries p-values may be bounded" else conditionMessage(ad),
    kpss_status=if(inherits(kp,"htest")) "computed; tseries p-values may be bounded" else conditionMessage(kp))
})
write.csv(do.call(rbind,rows),"data/country_epu_oil/unit_root_screen.csv",row.names=FALSE)
