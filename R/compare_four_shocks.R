#!/usr/bin/env Rscript
# Independent one-factor dominant-unit TVP-GVAR comparisons:
# GPR / VIX / EPU / OIL. Run from repository root.
# Prerequisites: installed repo dependencies (incl. threshtvp); 
# run 8.12/patch_BVAR_dominant_unit.R once beforehand.
# Use GDP log-level preprocessed macro input (data/model_input.csv).
# VIX: either US_vix in macro input, or processed quarterly CSV.

suppressPackageStartupMessages({
  library(compiler); library(snowfall); library(Matrix)
  library(mvtnorm); library(threshtvp); library(ggplot2)
})
source("R/BVAR_ttvp_dominant.r")
source("R/Datahandling.r")
source("R/auxilliary_functions_tvp.r")
source("8.12/dominant_gpr_vix_irf.R")

risk <- toupper(Sys.getenv("RISK", "GPR"))
stopifnot(risk %in% c("GPR","VIX","EPU","OIL"))
root <- Sys.getenv("COMPARE_OUTPUT_ROOT", "results/independent_shocks")
out <- file.path(root, tolower(risk))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
p <- as.integer(Sys.getenv("TVPGVAR_P","1"))
saves <- as.integer(Sys.getenv("TVPGVAR_SAVES","500"))
burns <- as.integer(Sys.getenv("TVPGVAR_BURNS","500"))
thin <- as.numeric(Sys.getenv("TVPGVAR_THIN","0.5"))
nhor <- as.integer(Sys.getenv("TVPGVAR_HORIZON","12"))
seed <- as.integer(Sys.getenv("TVPGVAR_SEED","20260816"))
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_SHOCK_PCT","10"))
stopifnot(p %in% 1:2, saves>=10, burns>=10, thin>0, thin<=1, nhor>=1, shock_pct>0)
RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
macrovars <- c("y","dp","r","de","deq")
macrocols <- unlist(lapply(countries, function(c) paste0(c,"_",macrovars)))
readq <- function(file, qcol, vcol) {
  if (!file.exists(file)) stop("Missing input: ",file)
  d <- read.csv(file, check.names=FALSE, fileEncoding="UTF-8-BOM")
  if (!all(c(qcol,vcol) %in% names(d))) stop("Missing columns ",qcol,", ",vcol," in ",file)
  if(anyDuplicated(d[[qcol]])) stop("Duplicate quarters in ",file)
  data.frame(Quarter=as.character(d[[qcol]]), value=as.numeric(d[[vcol]]))
}
align <- function(src, qq, label) {
  i <- match(qq,src$Quarter)
  if (anyNA(i)) stop(label,": missing quarters: ",paste(qq[is.na(i)], collapse=","))
  z <- src$value[i]
  if(any(!is.finite(z))) stop(label,": non-finite data in model window")
  z
}
macro_file <- Sys.getenv("TVPGVAR_MACRO_INPUT","data/model_input.csv")
if(!file.exists(macro_file)) stop("Run existing GDP log-level input builder first: ", macro_file)
base <- read.csv(macro_file,check.names=FALSE)
stopifnot("Quarter" %in% names(base), all(macrocols %in% names(base)))
if(anyDuplicated(base$Quarter)) stop("Duplicate macro quarters")
qq <- as.character(base$Quarter)
if(!all(grepl("^[0-9]{4}Q[1-4]$",qq))) stop("Invalid quarter labels")
if (risk=="GPR") {
  src <- readq("8.12/gpr_quarterly_processed.csv","Quarter","LN_GPR_QMEAN")
  note <- "log(quarterly mean GPR)"
} else if (risk=="VIX") {
  if ("US_vix" %in% names(base)) {
    src <- data.frame(Quarter=qq,value=as.numeric(base$US_vix))
  } else {
    src <- readq("data/vix_quarterly_processed.csv","Quarter","US_vix")
  }
  note <- "log(quarterly arithmetic mean of positive daily VIX closes)"
} else if (risk=="EPU") {
  src <- readq("8.12/Global_EPU_Quarterly_1997Q1_2026Q3.csv","Quarter","GEPU_current")
  metadata <- read.csv("8.12/Global_EPU_Quarterly_1997Q1_2026Q3.csv")
  if ("Months_available" %in% names(metadata)) {
    incomplete <- metadata$Quarter[is.na(metadata$Months_available) | metadata$Months_available<3]
    if(any(qq %in% incomplete)) stop("EPU contains incomplete quarterly observations within model sample")
  }
  if(any(src$value<=0,na.rm=TRUE)) stop("Cannot log EPU<=0")
  src$value <- log(src$value)
  note <- "log(GEPU_current, quarterly arithmetic mean of 3 months)"
} else {
  src <- readq("8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv","quarter","brent_usd_per_barrel")
  if(any(src$value<=0,na.rm=TRUE)) stop("Cannot log Brent<=0")
  src$value <- log(src$value)
  note <- "log(quarterly Brent USD per barrel)"
}
x <- base[,c("Quarter",macrocols)]
x$GL_shock <- align(src,qq,risk)
write.csv(x,file.path(out,"aligned_input.csv"),row.names=FALSE)
xglobal <- as.matrix(x[,-1]); storage.mode(xglobal)<-"double"
if(any(!is.finite(xglobal))) stop("Non-finite model input")
stopifnot(ncol(xglobal)==71L, tail(colnames(xglobal),1)=="GL_shock")

# Exactly 14x5 macro variables + one independent global variable;
# original country trade weights retained. No other risk is included.
w <- read.csv(Sys.getenv("TVPGVAR_WEIGHT_FILE", "data/trade_weights.csv"),
              check.names = FALSE, fileEncoding = "UTF-8-BOM",
              stringsAsFactors = FALSE)
names(w) <- sub("^\ufeff", "", names(w))
# The repository's official 14-economy trade-weight CSV uses 'Reporter',
# not 'Country', as the reporter-economy key.
if (!"Reporter" %in% names(w)) {
  stop("Trade-weight file must have a 'Reporter' column; found: ",
       paste(names(w), collapse = ", "))
}
if (anyNA(w$Reporter) || anyDuplicated(w$Reporter)) {
  stop("Missing or duplicate Reporter economies in trade weights")
}
rownames(w) <- as.character(w$Reporter)
stopifnot(all(countries %in% rownames(w)), all(countries %in% names(w)))
Wt <- as.matrix(w[countries, countries, drop = FALSE]); storage.mode(Wt) <- "double"
if (any(!is.finite(Wt)) || any(Wt < 0)) stop("Invalid trade weights")
if (any(abs(diag(Wt)) > 1e-10)) stop("Trade-weight diagonal must be zero")
if (any(rowSums(Wt) <= 0)) stop("Zero trade-weight row")
Wt <- Wt / rowSums(Wt)
cat("Trade weights loaded and validated (Reporter x partner): 14 x 14\n")
units<-c(countries,"GL"); K<-ncol(xglobal); gl_idx<-match("GL_shock",colnames(xglobal))
gW<-setNames(vector("list",length(units)),units)
for(cc in countries) {
  Wi<-matrix(0,11L,K)
  own<-match(paste0(cc,"_",macrovars),colnames(xglobal))
  Wi[cbind(1:5,own)]<-1
  for(k in 1:5) for(other in setdiff(countries,cc)) {
    j<-match(paste0(other,"_",macrovars[k]),colnames(xglobal))
    Wi[5+k,j]<-Wt[cc,other]
  }
  Wi[11,gl_idx]<-1
  rownames(Wi)<-c(paste0(cc,"_",macrovars),paste0("foreign_",macrovars),"global_shock")
  gW[[cc]]<-Wi
}
gW[["GL"]]<-matrix(0,1,K,dimnames=list("GL_shock",colnames(xglobal)))
gW[["GL"]][1,gl_idx]<-1
for(cc in units) colnames(gW[[cc]])<-colnames(xglobal)
Data.setup<-list(bigx=xglobal,gW=gW,new.data=xglobal,countries=units,
                 country_units=countries,quarters=qq)
Daten<-xglobal
cN<-units
# Override only the 2x2 Cholesky part: scalar dominant innovation.
# Rest of original GVAR transition reconstruction is reused.
dominant_gpr_struct_irf <- function(G,F,sig,x,units,horizon=12,
                                    shock_var="GL_shock",shock_pct=10) {
  j<-match("GL_shock",rownames(x))
  if(is.na(j)) stop("Missing GL_shock in global state")
  ui<-match("GL",units)
  v<-as.numeric(sig[[ui]])
  if(length(v)!=1L || !is.finite(v) || v<=0) stop("Invalid scalar shock variance")
  u<-rep(0,nrow(x)); u[j]<-sqrt(v)
  raw<-as.numeric(solve(G,u))
  if(!is.finite(raw[j]) || abs(raw[j])<1e-12) stop("Cannot normalize impact")
  impact<-raw*(log1p(shock_pct/100)/raw[j])
  ir<-matrix(NA_real_,nrow(x),horizon+1,dimnames=list(rownames(x),as.character(0:horizon)))
  ir[,1]<-impact
  for(h in seq_len(horizon)) {
    z<-rep(0,length(impact))
    for(k in seq_len(min(length(F),h))) z<-z+as.vector(F[[k]] %*% ir[,h-k+1])
    ir[,h+1]<-z
  }
  ir
}
# Original get_dominant_gpr_irf_t() dynamically calls the overridden scalar function.
# Its parameter names "shock_pct" remain the same.
ext.inst<-FALSE
shrink.parm<-list(B_1=as.numeric(Sys.getenv("TVPGVAR_B1","8")),
                  B_2=as.numeric(Sys.getenv("TVPGVAR_B2","0.01")),
                  kappa0=as.numeric(Sys.getenv("TVPGVAR_KAPPA0","-0.005")))
# Pre-flight current and first-lag foreign controls, domestic lags.
for(i in seq_along(cN)) {
  cc<-cN[i]; End<-xglobal[,substr(colnames(xglobal),1,2)==cc,drop=FALSE]
  Wex<-tvpgvar_extract_wex(gW[[cc]],xglobal,ncol(End))
  X<-cbind(1,Wex,tvpgvar_wex_lag(Wex,1L),mlag(End,p))
  X<-X[(p+1L):nrow(X),,drop=FALSE]
  if(any(!is.finite(X))) stop("Bad preflight ",cc)
}
# Scalar global (GL) block requires separate handling. The original BVAR
# assumes M >= 2 and drops dimensions for a one-variable endogenous block.
BVAR_scalar_GL <- function(i, gW, bigx, Daten, cN,
                           nsave, nburn, thin_chain, ext.inst, parms) {
  stopifnot(cN[[i]] == "GL")
  suppressPackageStartupMessages(library(threshtvp))
  Yraw <- matrix(as.numeric(bigx[, "GL_shock"]), ncol = 1L,
                 dimnames = list(NULL, "GL_shock"))
  lag_y <- mlag(Yraw, p)
  X <- cbind(constant = 1, lag_y)
  X <- X[(p+1L):nrow(X), , drop = FALSE]
  colnames(X) <- c("constant", paste0("Ylag", seq_len(p)))
  Y <- Yraw[(p+1L):nrow(Yraw), 1L]
  stopifnot(nrow(X) == length(Y), all(is.finite(X)), all(is.finite(Y)))
  # Same TVP, stochastic-volatility and shrinkage settings as country blocks.
  est <- threshtvp::estimate_tvp(
    Y, X, save = nsave, burn = nburn, p = p,
    sv_on = TRUE, thin = thin_chain, priorbtheta = parms,
    priormu = c(0, 10), h0prior = "stationary",
    grid.length = 150, thrsh.pct = 0.1, thrsh.pct.high = 1.5,
    TVS = TRUE, CPU = 1
  )
  Avec <- est$posterior$A
  Hd <- est$posterior$H
  if(length(dim(Avec)) != 3L) {
    stop("Scalar GL: unexpected posterior A shape: ",
         paste(dim(Avec), collapse = "x"))
  }
  Aa <- aperm(Avec, c(2, 3, 1))
  T <- nrow(X); K <- ncol(X)
  if(dim(Aa)[1L] != T || dim(Aa)[2L] != K) {
    stop("Scalar GL posterior A mismatch: ",
         paste(dim(Aa), collapse = "x"), " expected ", T, "x", K, "xDraws")
  }
  nd <- min(dim(Aa)[3L], as.integer(round(thin_chain * nsave)))
  if(nd < 1L) stop("Scalar GL has no retained draws")
  if(!is.matrix(Hd) || ncol(Hd) != T || nrow(Hd) < nd) {
    stop("Scalar GL unexpected posterior H shape: ",
         paste(dim(Hd), collapse = "x"))
  }
  alpha <- array(Aa[, , seq_len(nd), drop = FALSE],
                 dim = c(T, K, 1L, nd),
                 dimnames = list(NULL, colnames(X), "GL_shock", NULL))
  sig <- array(NA_real_, c(T, 1L, 1L, nd))
  for(j in seq_len(nd)) sig[, 1L, 1L, j] <- exp(Hd[j, ])
  mean_coef <- apply(alpha, c(1L, 2L), mean)
  if(!is.matrix(mean_coef)) mean_coef <- matrix(mean_coef, nrow = T)
  resid <- matrix(Y - rowSums(X * mean_coef), ncol = 1L)
  list(ALPHA = alpha, SIGMApost = sig,
       W = gW[[i]], cc.res = resid)
}
# Model's patched BVAR reads lag setting from TVPGVAR_P; held constant across runs.
BVAR<-cmpfun(BVAR)
CPU<-min(4L,max(1L,as.integer(Sys.getenv("TVPGVAR_CPU","2"))))
rng<-vector("list",length(cN));rng[[1]]<-.Random.seed
if(length(cN)>1) for(i in 2:length(cN)) rng[[i]]<-parallel::nextRNGStream(rng[[i-1]])
sfInit(parallel=TRUE,cpus=CPU)
sfExport(list=list("mlag","BVAR","datahandling","xglobal","gW","Daten","cN",
                   "bvartvpm","saves","burns","thin","ext.inst","shrink.parm","rng",
                   "tvpgvar_extract_wex","tvpgvar_wex_lag","tvpgvar_safe_inverse",
                   "BVAR_scalar_GL","p"))
predDens<-tryCatch(sfLapply(seq_along(cN),function(i) {
  assign(".Random.seed",rng[[i]],envir=.GlobalEnv)
  tryCatch({
    estimator <- if(cN[[i]] == "GL") BVAR_scalar_GL else BVAR
    estimator(i,gW=gW,bigx=xglobal,Daten=Daten,cN=cN,
              nsave=saves,nburn=burns,thin_chain=thin,
              ext.inst=ext.inst,parms=shrink.parm)
  }, error=function(e) {
    stop(paste0("UNIT=",cN[[i]],"; ERROR=",conditionMessage(e)),
         call.=FALSE)
  })
}),finally=sfStop())
save(predDens,Data.setup,file=file.path(out,"posterior.rda"))
A<-lapply(predDens,`[[`,"ALPHA");S<-lapply(predDens,`[[`,"SIGMApost")
globalG<-lapply(predDens,`[[`,"W")
nd<-dim(A[[1]])[4]
Tirf<-nrow(xglobal)-p
dates<-qq[-seq_len(p)]
stopifnot(all(vapply(A,function(z)dim(z)[1],integer(1))==Tirf))
chosen<-intersect(c("2003Q1","2008Q3","2014Q3","2020Q1","2022Q1","2023Q4"),dates)
rows<-list();stability<-list();n<-0L
for(d in chosen) {
 tt<-match(d,dates); IR<-array(NA_real_,c(K,nhor+1L,nd))
 rho<-rep(NA_real_,nd)
 for(dd in seq_len(nd)) {
   fit<-get_dominant_gpr_irf_t(tt,
      draw_i=lapply(A,function(z) {
        array(z[,,,dd,drop=FALSE], dim=dim(z)[1:3],
              dimnames=dimnames(z)[1:3])
      }),
      Sig_draw_i=lapply(S,function(z) {
        array(z[,,,dd,drop=FALSE], dim=dim(z)[1:3])
      }),
      x=t(xglobal),globalG=globalG,units=cN,
      horizon=nhor,shock_pct=shock_pct)
   IR[,,dd]<-fit$IRF_post;rho[dd]<-fit$max_eigen_modulus
 }
 stable<-is.finite(rho)&rho<1
 stability[[d]]<-data.frame(shock=risk,date=d,total=nd,stable=sum(stable),
                           stable_share=mean(stable))
 if(!any(stable)) next
 for(j in seq_len(K)) for(h in 0:nhor) {
   z<-IR[j,h+1L,stable]
   qs<-quantile(z,c(.05,.16,.5,.84,.95),names=FALSE)
   n<-n+1L
   rows[[n]]<-data.frame(shock=risk,date=d,variable=colnames(xglobal)[j],
                        horizon=h,median=qs[3],low68=qs[2],high68=qs[4],
                        low90=qs[1],high90=qs[5],stable_draws=sum(stable))
 }
}
ss<-do.call(rbind,stability)
write.csv(ss,file.path(out,"stability.csv"),row.names=FALSE)
if(length(rows)) {
 tab<-do.call(rbind,rows)
 write.csv(tab,file.path(out,"irf_stable_only.csv"),row.names=FALSE)
}
writeLines(c(paste("Shock:",risk),paste("Data:",note),
             paste("Macro input:",macro_file),paste("Sample:",head(qq,1),tail(qq,1)),
             paste("p:",p,"q:1, K:",K,"saves:",saves,"burns:",burns),
             paste("IRF normalization: +",shock_pct,"% log level impact",sep=""),
             "One independent scalar dominant shock; no other risk index in system.",
             "GDP is log level; DO NOT cumulate its IRF.",
             "Only posterior draws with global spectral radius < 1 retained.",
             "No missing-data imputation."),file.path(out,"provenance.txt"))
cat("COMPLETED:",risk,"output:",out,"\n")
