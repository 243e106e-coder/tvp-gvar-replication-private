#!/usr/bin/env Rscript
# Exploratory integration-order and conditional cointegration screening.
# Never impute, replace raw data, or estimate a formal TVP-GVAR here.
suppressPackageStartupMessages(library(urca))
input <- "data/loglevels/model_semantic.csv"
outdir <- "results/loglevels_integration_screen"
if(!file.exists(input)) stop("Run build_five_loglevels_input.R first")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(input,check.names=FALSE)
countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
variables <- c("y","p","r","q","eq")
expected <- c("Quarter",unlist(lapply(countries,function(cc)paste0(cc,"_",variables))))
if(!identical(names(d),expected)) stop("Unexpected schema")
if(anyDuplicated(d$Quarter) || any(!grepl("^[12][0-9]{3}Q[1-4]$",d$Quarter))) stop("Bad quarters")
quarter_index <- as.integer(substr(d$Quarter,1,4))*4L+as.integer(substr(d$Quarter,6,6))
if(any(diff(quarter_index)!=1L)) stop("Nonconsecutive quarters")
for(nm in expected[-1]){
 d[[nm]]<-suppressWarnings(as.numeric(d[[nm]]))
 if(any(!is.finite(d[[nm]]))) stop("Missing or invalid observations: ",nm)
}
adf <- function(x,det){
  tryCatch({
    a<-urca::ur.df(x,type=det,lags=4,selectlags="AIC")
    key<-if(det=="trend")"tau3" else "tau2"
    stat<-as.numeric(a@teststat[1,key]); crit<-as.numeric(a@cval[key,"5pct"])
    c(reject=as.integer(stat<crit),stat=stat,critical=crit)
  },error=function(e)c(reject=NA_real_,stat=NA_real_,critical=NA_real_))
}
kpss <- function(x,det){
  tryCatch({
    k<-urca::ur.kpss(x,type=det,lags="short")
    z<-k@cval
    if(is.matrix(z)){
      if("5pct"%in%colnames(z)&&nrow(z)==1L) crit<-as.numeric(z[1,"5pct"])
      else if("5pct"%in%rownames(z)&&ncol(z)==1L) crit<-as.numeric(z["5pct",1])
      else stop("Ambiguous KPSS critical-value matrix")
    } else if(!is.null(names(z))&&"5pct"%in%names(z)) crit<-as.numeric(z["5pct"])
    else if(length(z)==4L) crit<-as.numeric(z[2L])
    else stop("Unknown KPSS critical values")
    statistic<-as.numeric(k@teststat)[1]
    if(!is.finite(crit)||!is.finite(statistic)) stop("Nonfinite KPSS")
    c(reject=as.integer(statistic>crit),stat=statistic,critical=crit)
  },error=function(e)c(reject=NA_real_,stat=NA_real_,critical=NA_real_))
}
screen <- list()
for(cc in countries) for(v in variables){
 nm<-paste0(cc,"_",v); x<-d[[nm]]
 for(level in 0:2){
   z<-if(level==0L)x else diff(x,differences=level)
   for(det in c("drift","trend")){
     a<-adf(z,det)
     k<-kpss(z,if(det=="drift")"mu" else "tau")
     screen[[length(screen)+1L]]<-data.frame(country=cc,variable=v,series=nm,difference_order=level,
        deterministic=det,adf_stat=a["stat"],adf_cv5=a["critical"],adf_reject_unit_root5=a["reject"],
        kpss_stat=k["stat"],kpss_cv5=k["critical"],kpss_reject_stationarity5=k["reject"])
   }
 }
}
results<-do.call(rbind,screen)
write.csv(results,file.path(outdir,"integration_tests_all.csv"),row.names=FALSE)
# Two tests give contrasting nulls; mark only well-aligned evidence, never automatic I(2).
classify <- function(name){
  z<-subset(results,series==name)
  ok <- function(order) {
    a<-subset(z,difference_order==order)
    if(anyNA(a[c("adf_reject_unit_root5","kpss_reject_stationarity5")]))return(FALSE)
    all(a$adf_reject_unit_root5==1L & a$kpss_reject_stationarity5==0L)
  }
  nonstat <- function(order){
    a<-subset(z,difference_order==order)
    if(anyNA(a[c("adf_reject_unit_root5","kpss_reject_stationarity5")]))return(FALSE)
    all(a$adf_reject_unit_root5==0L & a$kpss_reject_stationarity5==1L)
  }
  if(ok(0L)) return("I0_candidate")
  if(nonstat(0L)&&ok(1L)) return("I1_candidate")
  if(nonstat(0L)&&nonstat(1L)&&ok(2L)) return("I2_candidate_REQUIRES_VERIFICATION")
  "ambiguous_manual_review"
}
classification<-data.frame(series=expected[-1],stringsAsFactors=FALSE)
classification$country<-sub("_.*","",classification$series)
classification$variable<-sub(".*_","",classification$series)
classification$screen_label<-vapply(classification$series,classify,character(1))
write.csv(classification,file.path(outdir,"integration_classification_conservative.csv"),row.names=FALSE)
# Conditioned Johansen: Only test subsets where every member has conservative I1 evidence.
# At least two variables; short rate is NOT automatically included.
joh <- list()
for(cc in countries){
  subset1<-subset(classification,country==cc & screen_label=="I1_candidate")$variable
  if(length(subset1)<2L){
    joh[[length(joh)+1L]]<-data.frame(country=cc,variables=paste(subset1,collapse="|"),K=NA_integer_,ecdet=NA_character_,
      hypothesis=NA_character_,stat=NA_real_,critical5=NA_real_,reject5=NA,error="Fewer than two conservative I1 candidates")
    next
  }
  block<-as.matrix(d[paste0(cc,"_",subset1)])
  for(K in c(2L,3L))for(ecdet in c("const","trend")){
    val<-tryCatch({
      obj<-urca::ca.jo(block,type="trace",ecdet=ecdet,K=K,spec="transitory")
      n<-length(obj@teststat)
      data.frame(country=cc,variables=paste(subset1,collapse="|"),K=K,ecdet=ecdet,
        hypothesis=paste0("r <= ",rev(seq_len(n)-1L)),stat=as.numeric(obj@teststat),
        critical5=as.numeric(obj@cval[,"5pct"]),reject5=as.numeric(obj@teststat)>as.numeric(obj@cval[,"5pct"]),error="")
    },error=function(e)data.frame(country=cc,variables=paste(subset1,collapse="|"),K=K,ecdet=ecdet,
        hypothesis=NA_character_,stat=NA_real_,critical5=NA_real_,reject5=NA,error=conditionMessage(e)))
    joh[[length(joh)+1L]]<-val
  }
}
write.csv(do.call(rbind,joh),file.path(outdir,"johansen_on_screened_i1_subsets.csv"),row.names=FALSE)
write.csv(subset(results,country%in%c("KR","TR") & variable=="p"),file.path(outdir,"KR_TR_CPI_I0_I1_I2_screens.csv"),row.names=FALSE)
write.csv(subset(results,variable=="r"),file.path(outdir,"interest_rates_integration_screens.csv"),row.names=FALSE)
writeLines(c("Diagnostic screening only. ADF null=unit root; KPSS null=stationarity.",
 "I0 and I1 candidate labels demand agreement under TWO deterministic specifications; ambiguous does not mean nonstationary.",
 "I2 candidate is a SCREEN FLAG and MUST NOT be treated as a confirmed I(2) classification.",
 "Johansen subsets can vary across countries and hence cannot directly form a consistent 14-country error-correction GVAR.",
 "Johansen is not valid if unrecognized I(0), I(2), breaks or cross-country interdependence invalidate assumptions.",
 "Do not use these tests to force IRF signs; preserve original growth-rate baseline until an appropriate model is developed.",
 "No imputation, no zero filling, no mirroring, no changes to original inputs."),
 file.path(outdir,"READ_ME_FIRST.txt"))
cat("Completed conservative integration screening on ",nrow(d),"quarters.\n")
