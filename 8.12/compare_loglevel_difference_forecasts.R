#!/usr/bin/env Rscript
# Comparable one-quarter-ahead forecasting diagnostics for level AR(1),
# difference AR(1), and random walk. NOT a GVAR, TVP, or ECM estimator.
suppressPackageStartupMessages(library(readxl))
file <- "8.12/TVPGVAR_14C_5VAR_LOG_LEVELS.xlsx"
out <- "results/loglevels_model_comparison"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(file)) stop("Missing workbook")
x <- as.data.frame(readxl::read_excel(file,sheet="MODEL_BALANCED",.name_repair="minimal"),
                   check.names=FALSE)
countries<-c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
vars<-c("y","p","r","q","eq")
req<-c("Quarter",unlist(lapply(countries,function(i)paste0(i,"_",vars))))
if(!identical(names(x),req))stop("Workbook columns differ from documented schema")
x$Quarter<-as.character(x$Quarter)
if(anyDuplicated(x$Quarter)||any(!grepl("^[12][0-9]{3}Q[1-4]$",x$Quarter)))stop("Bad quarters")
qi<-as.integer(substr(x$Quarter,1,4))*4L+as.integer(substr(x$Quarter,6,6))
if(any(diff(qi)!=1L))stop("Non-consecutive sample")
for(v in req[-1]){
 x[[v]]<-suppressWarnings(as.numeric(x[[v]]))
 if(any(!is.finite(x[[v]])))stop("No imputation: missing/nonfinite ",v)
}
# Same quarters and same level-scale forecast target for all 3 specifications.
# Historical fit ends t-1; forecast observed level at t. A minimum of 60
# historical quarters is required for each origin.
min_train<-60L
if(nrow(x)<min_train+10L)stop("Too few observations")
pred_ar<-function(y,next_lag){
 n<-length(y)
 if(n<3L || any(!is.finite(y)))stop("Invalid AR sample")
 fit<-lm(y[-1L]~y[-n])
 co<-coef(fit)
 if(any(!is.finite(co)))return(NA_real_)
 unname(co[1L]+co[2L]*next_lag)
}
rows<-list()
for(country in countries)for(v in vars){
 name<-paste0(country,"_",v); y<-x[[name]]
 for(t in seq.int(min_train+1L,length(y))){
   hist<-y[seq_len(t-1L)]
   prev<-tail(hist,1L)
   p_level<-pred_ar(hist,prev)
   changes<-diff(hist)
   p_diff<-prev+pred_ar(changes,tail(changes,1L))
   p_rw<-prev
   if(!all(is.finite(c(p_level,p_diff,p_rw))))stop("Unusable fit: ",name," ",x$Quarter[t])
   actual<-y[t]
   for(model in c("level_AR1","difference_AR1","random_walk")){
     predicted<-switch(model,level_AR1=p_level,difference_AR1=p_diff,random_walk=p_rw)
     rows[[length(rows)+1L]]<-data.frame(country=country,variable=v,series=name,
       Quarter=x$Quarter[t],model=model,actual_level=actual,predicted_level=predicted,
       error=actual-predicted,sq_error=(actual-predicted)^2,abs_error=abs(actual-predicted))
   }
 }
}
pred<-do.call(rbind,rows)
write.csv(pred,file.path(out,"rolling_forecasts.csv"),row.names=FALSE)
group<-split(pred,paste(pred$series,pred$model,sep="|"))
metrics<-do.call(rbind,lapply(group,function(z)data.frame(country=z$country[1L],
 variable=z$variable[1L],series=z$series[1L],model=z$model[1L],
 n=nrow(z),test_start=z$Quarter[1L],test_end=tail(z$Quarter,1L),
 rmse=sqrt(mean(z$sq_error)),mae=mean(z$abs_error),bias=mean(z$error))))
rownames(metrics)<-NULL
write.csv(metrics,file.path(out,"forecast_metrics_same_level_target.csv"),row.names=FALSE)
# Residual dependence: in-sample diagnostics at final origin only; no structural claims.
res<-list()
for(country in countries)for(v in vars){
 nm<-paste0(country,"_",v); y<-x[[nm]]
 for(model in c("level_AR1","difference_AR1")){
   z<-if(model=="level_AR1")y else diff(y)
   fit<-lm(z[-1L]~z[-length(z)])
   e<-residuals(fit)
   rho1<-if(length(e)>2L)cor(e[-1L],e[-length(e)]) else NA_real_
   # AR1 persistence below 1 is not sufficient evidence for GVAR system stability.
   res[[length(res)+1L]]<-data.frame(series=nm,model=model,ar1=unname(coef(fit)[2L]),
      residual_sd=sd(e),residual_acf1=rho1)
 }
}
write.csv(do.call(rbind,res),file.path(out,"ar1_residual_diagnostics.csv"),row.names=FALSE)
writeLines(c(
 "This is a three-way *univariate benchmark* comparing forecasts on the SAME level target and SAME evaluation quarters.",
 "level_AR1 predicts y[t] from y[t-1]; difference_AR1 models delta y then adds prediction to y[t-1]; random_walk predicts y[t-1].",
 "For interest rates y is rate level in annual percentage points, not a log index.",
 "Do not compare raw AIC values from level versus differenced regressions: dependent variables differ.",
 "RMSE/MAE are common-target errors, not TVP-GVAR evidence and not structural IRFs.",
 "No cointegrated VECM or TVP-GVAR is estimated; this deliberately avoids false claims of model equivalence.",
 "Sample and timing are held constant, with expanding-window one-step-ahead forecasts and no future leakage.",
 "No zero filling, interpolation, mirroring or changes to original input."
 ),file.path(out,"READ_ME_FIRST.txt"))
cat("PASS: comparison complete; ",nrow(pred)," forecasts.\n",sep="")
