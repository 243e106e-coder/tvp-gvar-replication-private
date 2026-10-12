#!/usr/bin/env Rscript
# Break-aware integration sensitivity; exploratory, not formal post-selection inference.
# Uses predetermined split dates and never alters or imputes source data.
suppressPackageStartupMessages(library(urca))
input <- "data/loglevels/model_semantic.csv"
out <- "results/break_aware_integration"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(input)) stop("Missing data: run build_five_loglevels_input.R")
x <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
vars <- c("y","p","r","q","eq")
req <- c("Quarter",unlist(lapply(countries,function(cc)paste0(cc,"_",vars))))
if(!identical(names(x),req)) stop("Unexpected column order")
if(anyDuplicated(x$Quarter) || any(!grepl("^[12][0-9]{3}Q[1-4]$",x$Quarter))) stop("Invalid quarter index")
idx <- as.integer(substr(x$Quarter,1,4))*4L+as.integer(substr(x$Quarter,6,6))
if(any(diff(idx)!=1L)) stop("Quarter sequence has gaps")
for(v in req[-1]){
  x[[v]] <- suppressWarnings(as.numeric(x[[v]]))
  if(any(!is.finite(x[[v]]))) stop("Nonfinite series: ",v)
}
crit5 <- function(obj) {
  z <- obj@cval
  if(is.matrix(z)) {
    if("5pct" %in% colnames(z) && nrow(z)==1L) return(as.numeric(z[1L,"5pct"]))
    if("5pct" %in% rownames(z) && ncol(z)==1L) return(as.numeric(z["5pct",1L]))
  }
  if(!is.null(names(z)) && "5pct" %in% names(z)) return(as.numeric(z["5pct"]))
  if(length(z)==4L) return(as.numeric(z[2L])) # standard ur.kpss order: 10%,5%,2.5%,1%
  stop("Unrecognized KPSS critical-value schema")
}
one <- function(z,det,lag) tryCatch({
  if(length(z)<32L) stop("Insufficient observations (<32)")
  fit <- urca::ur.df(z,type=det,lags=lag,selectlags="AIC")
  key <- if(det=="trend") "tau3" else "tau2"
  av <- as.numeric(fit@teststat[1L,key]); ac <- as.numeric(fit@cval[key,"5pct"])
  k <- urca::ur.kpss(z,type=if(det=="trend")"tau" else "mu",lags="short")
  kv <- as.numeric(k@teststat)[1L]; kc <- crit5(k)
  if(any(!is.finite(c(av,ac,kv,kc)))) stop("Nonfinite test output")
  data.frame(adf_stat=av,adf_critical5=ac,adf_reject_unit_root5=av<ac,
             kpss_stat=kv,kpss_critical5=kc,kpss_reject_stationarity5=kv>kc,error="")
},error=function(e)data.frame(adf_stat=NA_real_,adf_critical5=NA_real_,
    adf_reject_unit_root5=NA,kpss_stat=NA_real_,kpss_critical5=NA_real_,
    kpss_reject_stationarity5=NA,error=conditionMessage(e)))
# Prespecified alternative periods; not selected to maximize statistical significance.
periods <- list(
 full=rep(TRUE,nrow(x)),
 pre_2020=idx < (2020L*4L+1L),
 from_2020=idx >= (2020L*4L+1L),
 pre_2022=idx < (2022L*4L+1L),
 from_2022=idx >= (2022L*4L+1L)
)
focus <- c("KR_p","TR_p",paste0(countries,"_r"))
rows<-list()
for(v in focus) for(p in names(periods)){
  mask <- periods[[p]]
  series <- x[[v]][mask]
  quarters <- x$Quarter[mask]
  if(length(series)<32L) next
  for(order in 0:2){
    z <- if(order==0L) series else diff(series,differences=order)
    for(det in c("drift","trend")) for(lag in c(2L,4L)){
      ans <- one(z,det,lag)
      rows[[length(rows)+1L]] <- cbind(
        data.frame(series=v,period=p,first_quarter=quarters[1L],
          last_quarter=tail(quarters,1L),transform_order=order,
          deterministic=det,max_adf_lags=lag,n=length(z)),
        ans)
    }
  }
}
res <- do.call(rbind,rows)
write.csv(res,file.path(out,"cpi_and_interest_rate_subsample_tests.csv"),row.names=FALSE)
# Deterministic break regressions are descriptive, NOT valid post-selection stationarity tests.
br <- list()
for(v in c("KR_p","TR_p")){
  z <- x[[v]]
  for(cut in c("2020Q1","2022Q1")){
    j <- match(cut,x$Quarter)
    if(is.na(j)) next
    t <- seq_along(z)
    post <- as.integer(t>=j)
    slope <- pmax(0L,t-j+1L)
    fit <- lm(z~t+post+slope)
    cf <- summary(fit)$coefficients
    br[[length(br)+1L]] <- data.frame(series=v,cut=cut,n=length(z),
      level_shift=unname(coef(fit)["post"]),slope_shift=unname(coef(fit)["slope"]),
      adj_r2=summary(fit)$adj.r.squared,
      note="Prespecified-break OLS coefficients are descriptive, not causal nor unit-root evidence")
  }
}
write.csv(do.call(rbind,br),file.path(out,"cpi_prespecified_break_regressions.csv"),row.names=FALSE)
# Compare model candidates without automatically declaring integration order.
writeLines(c(
 "SCREENING REPORT: do not infer I(2) solely from ADF/KPSS; unit-root tests have limited power and breaks matter.",
 "Break dates 2020Q1 and 2022Q1 are prespecified sensitivity splits, not estimated optimal break dates.",
 "Post-2022 subsamples are short and skipped when fewer than 32 observations.",
 "Pre/post subsample ADF/KPSS are not equivalent to a formal structural-break-aware unit-root test.",
 "Full-country lag and trend choices should be justified on economic grounds, not IRF signs.",
 "For a future 14-country GVAR: evaluate comparability of integrated variables, exogeneity, long-run links, and global stability.",
 "No imputations, no zero filling, no mirrored data, and no modification to prior model files."
),file.path(out,"INTERPRETATION_README.txt"))
cat("PASS generated ",nrow(res),"subsample test rows in ",out,"\n",sep="")
