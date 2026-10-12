#!/usr/bin/env Rscript
# Source: uploaded original Raw_* sheets transformed in model_log_levels_balanced.csv.
# Diagnostic only; does NOT estimate GVAR, infer cointegration, or change source data.
args <- commandArgs(trailingOnly = TRUE)
infile <- if(length(args)>=1) args[1] else 'model_log_levels_balanced.csv'
outdir <- if(length(args)>=2) args[2] else 'preflight_results'
if(!file.exists(infile)) stop('Input missing: ',infile)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
x <- read.csv(infile,check.names=FALSE,stringsAsFactors=FALSE)
countries <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
vars <- c('y','p','r','q','eq')
needed <- c('Quarter',as.vector(t(outer(countries, vars, paste, sep='_'))))
if(!identical(names(x),needed)) stop('Unexpected column names/order; stop rather than silently align.')
if(anyNA(x)) stop('Balanced input contains NA; stop.')
if(any(!is.finite(as.matrix(x[,-1])))) stop('Non-finite numeric values found.')
qi <- as.integer(substr(x$Quarter,1,4))*4L+as.integer(substr(x$Quarter,6,6))
if(any(diff(qi)!=1)) stop('Non-consecutive quarterly sample.')
summ <- do.call(rbind,lapply(names(x)[-1],function(v){z=x[[v]];data.frame(variable=v,mean=mean(z),sd=sd(z),min=min(z),max=max(z),lag1_acf=as.numeric(acf(z,plot=FALSE,lag.max=1)$acf[2]),stringsAsFactors=FALSE)}))
write.csv(summ,file.path(outdir,'descriptive_statistics.csv'),row.names=FALSE)
# ADF unit-root evidence; urca package is optional. Never manufacture p-values.
if(requireNamespace('urca',quietly=TRUE)){
  rows <- lapply(names(x)[-1],function(v){
    s <- x[[v]]
    res <- tryCatch(urca::ur.df(s,type='drift',lags=4,selectlags='AIC'),error=function(e)e)
    if(inherits(res,'error'))return(data.frame(variable=v,adf_tau2=NA_real_,critical_5pct=NA_real_,reject_unit_root_5pct=NA,error=conditionMessage(res)))
    st <- as.numeric(res@teststat[1,'tau2'])
    cv <- as.numeric(res@cval['tau2','5pct'])
    data.frame(variable=v,adf_tau2=st,critical_5pct=cv,reject_unit_root_5pct=(st<cv),error='')
  })
  write.csv(do.call(rbind,rows),file.path(outdir,'adf_unit_root_screen.csv'),row.names=FALSE)
} else {
  writeLines('urca unavailable: skipped ADF; install.packages("urca") then rerun. Do not conclude stationarity.',file.path(outdir,'adf_not_run.txt'))
}
writeLines(c('Preflight only: no TVP-GVAR was fit.',
 'Log-level I(1) variables can violate stationary-VAR stability assumptions.',
 'Evaluate integration order for each series and formal cointegration/model specification.',
 'Do NOT condition final claims on IRF sign selection or silent stable-draw rejection.',
 'No imputation, no zero-filling, no mirrored data.'),file.path(outdir,'READ_ME_FIRST.txt'))
cat('PASS preflight:',nrow(x),'quarters and',length(countries),'countries. Results at',outdir,'\n')
