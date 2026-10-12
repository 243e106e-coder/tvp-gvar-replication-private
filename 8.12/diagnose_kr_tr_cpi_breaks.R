#!/usr/bin/env Rscript
# Targeted Korean / Turkish CPI diagnostics with structural breaks.
# Zivot-Andrews: endogenous ONE break, unit-root null; not a proof of I(0)/I(2).
# strucchange: breakpoints in deterministic trend regression, DESCRIPTIVE only.
suppressPackageStartupMessages({library(urca);library(strucchange)})
path <- "data/loglevels/model_semantic.csv"
out <- "results/loglevels_robust_preflight"
if(!file.exists(path)) stop("Run 8.12/build_five_loglevels_input.R first")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
dat <- read.csv(path,check.names=FALSE)
if(!all(c("Quarter","KR_p","TR_p") %in% names(dat))) stop("Missing required CPI series")
if(anyNA(dat[c("KR_p","TR_p")])) stop("CPI contains missing observations; no filling")
q <- as.character(dat$Quarter)
if(anyDuplicated(q)) stop("Duplicate quarter labels")
n <- nrow(dat)
if(n<45L) stop("Insufficient sample length for break diagnostic")
rows <- list()
breakrows <- list()
for(country in c("KR","TR")){
  z <- as.numeric(dat[[paste0(country,"_p")]])
  if(any(!is.finite(z))) stop("Nonfinite CPI: ",country)
  for(scale in c("log_level","first_difference")){
    series <- if(scale=="log_level") z else diff(z)
    quarters <- if(scale=="log_level") q else q[-1]
    for(model in c("intercept","trend","both")){
      one <- tryCatch({
        fit <- urca::ur.za(series,model=model,lag=4)
        val <- as.numeric(fit@teststat)[1]
        cv <- as.numeric(fit@cval)
        nms <- names(fit@cval)
        # ur.za critical values can be a named numeric vector or matrix; do not invent comparisons.
        if(is.matrix(fit@cval)){
          zcv <- fit@cval
          if("5pct" %in% colnames(zcv)) crit <- as.numeric(zcv[1,"5pct"])
          else if("5pct" %in% rownames(zcv)) crit <- as.numeric(zcv["5pct",1])
          else stop("Unknown ZA critical-value matrix layout")
        } else if(!is.null(nms) && "5pct" %in% nms) crit <- as.numeric(fit@cval["5pct"])
        else if(length(cv)==3L) crit <- cv[2] # ur.za documentation reports 1%, 5%, 10% in order
        else stop("Unknown ZA critical-value layout")
        bp <- as.integer(fit@bpoint)[1]
        if(!is.finite(val) || !is.finite(crit) || !is.finite(bp) || bp<1 || bp>length(quarters)) stop("Invalid ZA result")
        data.frame(country=country,transform=scale,model=model,statistic=val,critical_5pct=crit,
                   reject_unit_root_5pct=val<crit,break_index=bp,break_quarter=quarters[bp],error="")
      },error=function(e)data.frame(country=country,transform=scale,model=model,statistic=NA_real_,
        critical_5pct=NA_real_,reject_unit_root_5pct=NA,break_index=NA_integer_,break_quarter=NA_character_,
        error=conditionMessage(e)))
      rows[[length(rows)+1L]] <- one
    }
    t <- seq_along(series)
    b <- tryCatch({
      f <- strucchange::breakpoints(series~t,h=0.15,breaks=1)
      idx <- as.integer(f$breakpoints)[1]
      data.frame(country=country,transform=scale,break_index=idx,
                 break_quarter=if(is.finite(idx)) quarters[idx] else NA_character_,
                 note="Single candidate break in intercept/trend regression, not an endogenous unit-root test",error="")
    },error=function(e)data.frame(country=country,transform=scale,break_index=NA_integer_,
      break_quarter=NA_character_,note="Diagnostic failed",error=conditionMessage(e)))
    breakrows[[length(breakrows)+1L]] <- b
  }
}
za <- do.call(rbind,rows)
write.csv(za,file.path(out,"KR_TR_CPI_zivot_andrews.csv"),row.names=FALSE)
write.csv(do.call(rbind,breakrows),file.path(out,"KR_TR_CPI_trend_break_candidates.csv"),row.names=FALSE)
writeLines(c("ZA accommodates one endogenous break only, not every real-world regime change.",
 "The ZA null is a unit root, and critical values depend on the model and implementation.",
 "Short-run inflation persistence does not prove I(2).",
 "Trend-break candidates are descriptive, not independent evidence for cointegration.",
 "Never interpolate or substitute missing CPI observations."),
 file.path(out,"CPI_BREAK_INTERPRETATION.txt"))
if(any(nzchar(za$error))) warning("Some ZA specifications failed; consult error column, do not ignore.")
cat("CPI break diagnostics generated at ",out,"\n",sep="")
