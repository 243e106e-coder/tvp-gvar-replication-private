#!/usr/bin/env Rscript
# Robust exploratory preflight for 14 economies x 5 log-level variables.
# No model estimation, source mutation, imputation, or sign-targeted selection.
suppressPackageStartupMessages(library(urca))
infile <- Sys.getenv("LEVEL_INPUT", "data/loglevels/model_semantic.csv")
outdir <- Sys.getenv("LEVEL_TEST_OUTPUT", "results/loglevels_robust_preflight")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
if(!file.exists(infile)) stop("Missing semantic input: ", infile)
d <- read.csv(infile, check.names=FALSE, stringsAsFactors=FALSE)
cc <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
vv <- c("y","p","r","q","eq")
cols <- c("Quarter",unlist(lapply(cc,function(i) paste0(i,"_",vv))))
if(!identical(names(d),cols)) stop("Unexpected input schema or order")
if(anyDuplicated(d$Quarter) || any(!grepl("^[12][0-9]{3}Q[1-4]$",d$Quarter))) stop("Invalid quarter labels")
qi <- as.integer(substr(d$Quarter,1,4))*4L+as.integer(substr(d$Quarter,6,6))
if(any(diff(qi)!=1)) stop("Non-consecutive quarters")
for(nm in cols[-1]){
  d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  if(any(!is.finite(d[[nm]]))) stop("Missing or nonfinite data: ",nm)
}
adf_one <- function(z,spec,cap) tryCatch({
  fit <- urca::ur.df(z,type=spec,lags=cap,selectlags="AIC")
  key <- if(spec=="trend") "tau3" else "tau2"
  stat <- as.numeric(fit@teststat[1,key]); cv <- as.numeric(fit@cval[key,"5pct"])
  data.frame(stat=stat,cv5=cv,reject5=stat<cv,error="")
},error=function(e)data.frame(stat=NA_real_,cv5=NA_real_,reject5=NA,error=conditionMessage(e)))
kpss_one <- function(z,spec) tryCatch({
  fit <- urca::ur.kpss(z,type=spec,lags="short")
  stat <- as.numeric(fit@teststat)[1]
  # urca::ur.kpss stores critical values under names such as "5pct";
  # match carefully instead of silently exporting NA for all tests.
  critical <- fit@cval
  if(is.matrix(critical)) {
    indices <- which(grepl("^5\\s*(pct|%)?$",trimws(colnames(critical)),ignore.case=TRUE))
    if(length(indices)==1L) critical <- critical[,indices] else {
      indices <- which(grepl("^5\\s*(pct|%)?$",trimws(rownames(critical)),ignore.case=TRUE))
      if(length(indices)!=1L) stop("KPSS critical value label not found")
      critical <- critical[indices,]
    }
  }
  critical_names <- names(critical)
  idx <- which(grepl("^5\\s*(pct|%)?$",trimws(critical_names),ignore.case=TRUE))
  if(length(idx)!=1L) stop("Cannot find unambiguous 5% KPSS critical value; available: ",paste(critical_names,collapse=","))
  cv <- as.numeric(critical[idx])
  if(length(stat)!=1L || length(cv)!=1L || !is.finite(stat) || !is.finite(cv)) stop("Invalid KPSS statistic/critical value")
  data.frame(stat=stat,cv5=cv,reject_stationarity5=stat>cv,error="")
},error=function(e)data.frame(stat=NA_real_,cv5=NA_real_,reject_stationarity5=NA,error=conditionMessage(e)))
adf <- list()
for(name in cols[-1]) for(transform in c("level","difference")){
  z <- if(transform=="level") d[[name]] else diff(d[[name]])
  for(spec in c("drift","trend")) for(cap in c(4L,8L)){
    adf[[length(adf)+1L]] <- cbind(data.frame(variable=name,transform=transform,deterministic=spec,max_lags=cap),adf_one(z,spec,cap))
  }
}
adf <- do.call(rbind,adf)
write.csv(adf,file.path(outdir,"adf_robustness.csv"),row.names=FALSE)
kpss <- list()
for(name in cols[-1]) for(transform in c("level","difference")){
  z <- if(transform=="level") d[[name]] else diff(d[[name]])
  for(spec in c("mu","tau")){
    kpss[[length(kpss)+1L]] <- cbind(data.frame(variable=name,transform=transform,kpss_null=if(spec=="mu") "level_stationary" else "trend_stationary"),kpss_one(z,spec))
  }
}
kpss <- do.call(rbind,kpss)
write.csv(kpss,file.path(outdir,"kpss_robustness.csv"),row.names=FALSE)
if(anyNA(kpss$cv5) || anyNA(kpss$reject_stationarity5)) stop("KPSS critical values not resolved: check implementation, do not publish incomplete output.")
joh <- list()
for(country in cc){
  Z <- as.matrix(d[paste0(country,"_",vv)])
  for(ecdet in c("const","trend")) for(K in c(2L,3L)){
    item <- tryCatch({
      fit <- urca::ca.jo(Z,type="trace",ecdet=ecdet,K=K,spec="transitory")
      n <- length(fit@teststat)
      stat <- as.numeric(fit@teststat); cv <- as.numeric(fit@cval[,"5pct"])
      data.frame(country=country,ecdet=ecdet,K=K,rank_hypothesis=paste0("r <= ",rev(seq_len(n)-1L)),
                 stat=stat,critical_5pct=cv,reject5=stat>cv,error="")
    },error=function(e)data.frame(country=country,ecdet=ecdet,K=K,rank_hypothesis=NA_character_,
      stat=NA_real_,critical_5pct=NA_real_,reject5=NA,error=conditionMessage(e)))
    joh[[length(joh)+1L]] <- item
  }
}
write.csv(do.call(rbind,joh),file.path(outdir,"johansen_corrected_exploratory.csv"),row.names=FALSE)
write.csv(subset(adf,variable %in% c("KR_p","TR_p")),file.path(outdir,"KR_TR_CPI_adf.csv"),row.names=FALSE)
write.csv(subset(kpss,variable %in% c("KR_p","TR_p")),file.path(outdir,"KR_TR_CPI_kpss.csv"),row.names=FALSE)
writeLines(c("EXPLORATORY EVIDENCE ONLY; not an estimated cointegrated TVP-GVAR.",
 "ADF null is a unit root; KPSS null is stationarity.",
 "Johansen hypotheses are labeled in descending order returned by urca.",
 "Johansen including potentially I(0) short rates is not an established cointegration inference.",
 "Assess structural breaks, lag sensitivity, integration orders, and deterministic assumptions.",
 "No missing-value replacement, no source modification, no IRF sign selection."),
 file.path(outdir,"METHODOLOGY_LIMITATIONS.txt"))
cat("PASS robust log-level diagnostics in ",outdir,"\n",sep="")
