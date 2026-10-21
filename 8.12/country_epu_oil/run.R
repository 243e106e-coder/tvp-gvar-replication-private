#!/usr/bin/env Rscript
source("8.12/country_epu_oil/prepare.R")
epu_country<-Sys.getenv("EPU_COUNTRY","US")
placement<-"COUNTRY_OIL"
Data.setup<-prepare_epu_wti(epu_country,placement)
run_id<-paste(epu_country,placement,sep="_")
out<-file.path("results/country_epu_oil",run_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(Data.setup$aligned,file.path(out,"aligned_input.csv"),row.names=FALSE)
xglobal<-Data.setup$bigx;gW<-Data.setup$gW
Daten<-xglobal;cN<-Data.setup$countries
p<-as.integer(Sys.getenv("TVPGVAR_P","1"))
saves<-as.integer(Sys.getenv("TVPGVAR_SAVES","500"))
burns<-as.integer(Sys.getenv("TVPGVAR_BURNS","500"))
thin<-as.numeric(Sys.getenv("TVPGVAR_THIN","0.5"))
nhor<-as.integer(Sys.getenv("TVPGVAR_HORIZON","12"))
shock_pct<-as.numeric(Sys.getenv("TVPGVAR_SHOCK_PCT","10"))
stopifnot(p %in% 1:2,saves>=10,burns>=10,thin>0,thin<=1,nhor>=1,shock_pct>0)
source("8.12/country_epu_oil/patch_estimator.R")
source("R/BVAR_country_oil.r")
source("R/Datahandling.r")
source("R/auxilliary_functions_tvp.r")
source("8.12/country_epu_oil/irf.R")
design<-do.call(rbind,lapply(cN,function(cc) {
  End<-xglobal[,substr(colnames(xglobal),1,2)==cc,drop=FALSE]
  Wex<-tvpgvar_extract_wex(gW[[cc]],xglobal,ncol(End))
  X<-cbind(1,Wex,tvpgvar_wex_lag(Wex,1L),mlag(End,p))
  X<-X[-seq_len(p),,drop=FALSE]
  stopifnot(all(is.finite(X)),nrow(X)>ncol(X),qr(X)$rank==ncol(X))
  data.frame(unit=cc,equations=ncol(End),Wex=ncol(Wex),
             observations=nrow(X),regressors=ncol(X),rank=qr(X)$rank)
}))
write.csv(design,file.path(out,"design_audit.csv"),row.names=FALSE)
settings<-list(epu_country=epu_country,placement=placement,p=p,q=1L,
               saves=saves,burns=burns,thin=thin,horizon=nhor,shock_pct=shock_pct,
               seed=20260816,rate="first difference; level response is draw-wise cumulative",
               gdp="See data/country_epu_oil/sample_audit.json",
               order=paste0(epu_country, ":[EPU,y,dp,r,de,deq]; OO:[WTI]"),
               sample=range(Data.setup$quarters),
               interpretation="Recursive EPU and oil-price innovations, not isolated causal policy/supply shocks")
dput(settings,file=file.path(out,"settings.R"))
file.copy("data/country_epu_oil/sample_audit.json",out,overwrite=TRUE)
file.copy("data/country_epu_oil/refresh_manifest.csv",out,overwrite=TRUE)
writeLines(c("EPU is endogenous only in its economic country's national module.",
             "WTI is in a separate autonomous scalar OO module: no macro or EPU feedback into oil equations.",
             "Other countries receive trade-weighted foreign macro variables and common EPU/WTI controls.",
             "National Cholesky order EPU,y,dp,r,de,deq; blockwise innovations are uncorrelated by assumption.",
             "This restriction does not identify an oil-supply shock.",
             "GDP growth is the initial specification; log-level alternative is explicitly configured in build_input.py.",
             "Historical trade-weight provenance remains unresolved; no new mirror completion is performed.",
             "valid_all and stable_only summaries must be interpreted together."),file.path(out,"provenance.txt"))
peak_dir<-"data/country_epu_oil"
peaks<-read.csv(file.path(peak_dir,"peak_dates.csv"),stringsAsFactors=FALSE)
peaks<-peaks[peaks$series %in% c(paste0(epu_country,"_EPU"),"WTI"),]
chosen<-sort(unique(as.character(peaks$Quarter)))
available<-Data.setup$quarters[-seq_len(p)]
stopifnot(available[1]=="2000Q1",length(chosen)>0,all(chosen %in% available))
for(f in c("peak_dates.csv","peak_rules.json"))
  stopifnot(file.copy(file.path(peak_dir,f),out,overwrite=TRUE))
write.csv(peaks,file.path(out,"selected_peak_dates.csv"),row.names=FALSE)
if(Sys.getenv("PREFLIGHT_ONLY","false")=="true") {
  cat("PASS preflight",run_id,nrow(xglobal),"quarters",ncol(xglobal),"variables\n")
  quit(save="no",status=0)
}
suppressPackageStartupMessages({
  library(compiler);library(snowfall);library(Matrix);library(mvtnorm);library(threshtvp)
})
RNGkind("L'Ecuyer-CMRG");set.seed(settings$seed)
ext.inst<-FALSE
shrink.parm<-list(B_1=8,B_2=0.01,kappa0=-0.005)
BVAR<-cmpfun(BVAR)
CPU<-min(4L,max(1L,as.integer(Sys.getenv("TVPGVAR_CPU","2"))))
rng<-vector("list",length(cN));rng[[1]]<-.Random.seed
if(length(cN)>1) for(i in 2:length(cN)) rng[[i]]<-parallel::nextRNGStream(rng[[i-1]])
sfInit(parallel=TRUE,cpus=CPU)
sfExport(list=list("mlag","BVAR","datahandling","xglobal","gW","Daten","cN",
                   "bvartvpm","saves","burns","thin","ext.inst","shrink.parm","rng",
                   "tvpgvar_extract_wex","tvpgvar_wex_lag","tvpgvar_safe_inverse",
                   "p"))
predDens<-tryCatch(sfLapply(seq_along(cN),function(i) {
  assign(".Random.seed",rng[[i]],envir=.GlobalEnv)
  tryCatch({
    estimator <- BVAR
    estimator(i,gW=gW,bigx=xglobal,Daten=Daten,cN=cN,
              nsave=saves,nburn=burns,thin_chain=thin,
              ext.inst=ext.inst,parms=shrink.parm)
  }, error=function(e) {
    stop(paste0("UNIT=",cN[[i]],"; ERROR=",conditionMessage(e)),
         call.=FALSE)
  })
}),finally=sfStop())

save(predDens,Data.setup,settings,file=file.path(out,"posterior.rda"))
A<-lapply(predDens,`[[`,"ALPHA");S<-lapply(predDens,`[[`,"SIGMApost")
globalG<-lapply(predDens,`[[`,"W")
nd<-dim(A[[1]])[4];K<-ncol(xglobal)
stopifnot(all(vapply(A,function(z)dim(z)[4],integer(1))==nd),
          all(vapply(A,function(z)dim(z)[1],integer(1))==nrow(xglobal)-p))
dates<-Data.setup$quarters[-seq_len(p)]
stopifnot(all(chosen %in% dates))
shock_host<-Data.setup$host;active_risk_names<-Data.setup$risk_names
rows<-list();stability<-list();failures<-list();raw_by_date<-list()
add_summary<-function(IR,mask,date,shock,filter,response_kind="raw") {
  if(!any(mask)) return(NULL)
  do.call(rbind,lapply(seq_len(dim(IR)[1]),function(j)
    do.call(rbind,lapply(0:nhor,function(h) {
      z<-IR[j,h+1,mask]
      qs<-quantile(z,c(.05,.16,.5,.84,.95),names=FALSE)
      data.frame(epu_country=epu_country,placement=placement,date=date,shock=shock,
                 variable=dimnames(IR)[[1]][j],response_kind=response_kind,filter=filter,
                 horizon=h,median=qs[3],low68=qs[2],high68=qs[4],
                 low90=qs[1],high90=qs[5],draws=length(z))
    }))))
}
for(date in chosen) {
  tt<-match(date,dates);rho<-rep(NA_real_,nd);condition<-rho
  IRall<-array(NA_real_,c(K,nhor+1,nd,2),
              dimnames=list(colnames(xglobal),0:nhor,NULL,c("EPU","WTI")))
  for(dd in seq_len(nd)) {
    # If either injection fails, exclude this draw from BOTH shocks and log it.
    fit<-tryCatch({
      draw_i<-lapply(A,function(z) array(z[,,,dd,drop=FALSE],dim=dim(z)[1:3],dimnames=dimnames(z)[1:3]))
      sig_i<-lapply(S,function(z) array(z[,,,dd,drop=FALSE],dim=dim(z)[1:3]))
      answers<-lapply(1:2,function(sh) {
        active_shock<<-sh
        get_dominant_gpr_irf_t(tt,draw_i,sig_i,t(xglobal),globalG,cN,
                               horizon=nhor,shock_pct=shock_pct)
      })
      answers
    },error=function(e) {
      failures[[length(failures)+1L]]<<-data.frame(date=date,draw=dd,error=conditionMessage(e))
      NULL
    })
    if(is.null(fit)) next
    rho[dd]<-fit[[1]]$max_eigen_modulus
    condition[dd]<-fit[[1]]$G_condition_number
    for(sh in 1:2) IRall[,,dd,sh]<-fit[[sh]]$IRF_post
  }
  valid<-is.finite(rho)
  for(dd in seq_len(nd)) if(valid[dd]&&!all(is.finite(IRall[,,dd,]))) valid[dd]<-FALSE
  stable<-valid&rho<1
  stability[[date]]<-data.frame(epu_country=epu_country,placement=placement,date=date,
     total=nd,valid=sum(valid),failed=sum(!valid),stable=sum(stable),
     stable_share=sum(stable)/nd,
     stable_share_of_valid=if(any(valid))mean(rho[valid]<1)else NA_real_,
     median_rho=if(any(valid))median(rho[valid])else NA_real_,
     median_G_condition=if(any(valid))median(condition[valid])else NA_real_)
  for(sh in 1:2) {
    if(!any(peaks$Quarter==date & peaks$shock==c("EPU","WTI")[sh])) next
    IR<-array(IRall[,,,sh],dim=c(K,nhor+1,nd),dimnames=list(colnames(xglobal),0:nhor,NULL))
    for(filter in c("valid_all","stable_only")) {
      mask<-if(filter=="stable_only")stable else valid
      rows[[length(rows)+1L]]<-add_summary(IR,mask,date,c("EPU","WTI")[sh],filter)
      rr<-grep("_r$",colnames(xglobal))
      level<-IR[rr,,,drop=FALSE]
      for(j in seq_along(rr)) for(dd in which(mask)) level[j,,dd]<-cumsum(level[j,,dd])
      rows[[length(rows)+1L]]<-add_summary(level,mask,date,c("EPU","WTI")[sh],filter,"rate_level_cumulative")
    }
  }
  raw_by_date[[date]]<-list(IRF=IRall,rho=rho,G_condition=condition,valid=valid,stable=stable)
  cat("IRF complete",run_id,date,"stable",sum(stable),"/",nd,"\n")
}
save(raw_by_date,settings,file=file.path(out,"irf_draws.rda"))
ss<-do.call(rbind,stability);write.csv(ss,file.path(out,"stability.csv"),row.names=FALSE)
tab<-do.call(rbind,rows)
if(!is.null(tab)) {
  tab$significant68<-tab$low68>0|tab$high68<0
  tab$significant90<-tab$low90>0|tab$high90<0
}
if(!is.null(tab)) write.csv(tab,file.path(out,"irf_summary.csv"),row.names=FALSE)
if(length(failures)) write.csv(do.call(rbind,failures),file.path(out,"failed_draws.csv"),row.names=FALSE)
if(is.null(tab)) stop("No valid IRFs; posterior and failure diagnostics saved")
cat("COMPLETED",run_id,"\n")
