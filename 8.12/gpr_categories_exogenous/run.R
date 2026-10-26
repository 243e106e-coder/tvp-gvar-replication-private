#!/usr/bin/env Rscript
source("8.12/gpr_categories_exogenous/prepare.R")
series<-Sys.getenv("GPR_SERIES","GPR")
root<-"8.12/gpr_categories_exogenous/data_20261026"
manifest<-read.csv(file.path(root,"run_manifest.csv"),stringsAsFactors=FALSE)
stopifnot(series %in% manifest$series)
Data.setup<-prepare_external_gpr(series)
out<-file.path("results/gpr_categories_exogenous",series)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(Data.setup$aligned,file.path(out,"aligned_input.csv"),row.names=FALSE)
xglobal<-Data.setup$bigx;gW<-Data.setup$gW;Daten<-xglobal;cN<-Data.setup$countries
p<-as.integer(Sys.getenv("TVPGVAR_P","1"))
saves<-as.integer(Sys.getenv("TVPGVAR_SAVES","1000"));burns<-as.integer(Sys.getenv("TVPGVAR_BURNS","1000"))
thin<-as.numeric(Sys.getenv("TVPGVAR_THIN","0.5"));nhor<-as.integer(Sys.getenv("TVPGVAR_HORIZON","12"));shock_pct<-10
stopifnot(p==1L,saves>=10,burns>=10,thin>0,thin<=1,nhor==12L)
source("8.12/country_epu_oil/patch_estimator.R")
source("R/BVAR_country_oil.r");source("R/Datahandling.r");source("R/auxilliary_functions_tvp.r")
source("8.12/gpr_categories_exogenous/irf.R")
design<-do.call(rbind,lapply(cN,function(cc) {
  End<-xglobal[,startsWith(colnames(xglobal),paste0(cc,"_")),drop=FALSE]
  Wex<-tvpgvar_extract_wex(gW[[cc]],xglobal,ncol(End))
  X<-cbind(1,Wex,tvpgvar_wex_lag(Wex,1L),mlag(End,p))[-seq_len(p),,drop=FALSE]
  stopifnot(all(is.finite(X)),nrow(X)>ncol(X),qr(X)$rank==ncol(X))
  if(cc=="GL")stopifnot(ncol(Wex)==0L,ncol(X)==2L)
  data.frame(series=series,unit=cc,equations=ncol(End),external_controls=ncol(Wex),observations=nrow(X),regressors=ncol(X),rank=qr(X)$rank,condition_number=kappa(X,exact=TRUE))
}))
write.csv(design,file.path(out,"design_audit.csv"),row.names=FALSE)
bg<-do.call(rbind,lapply(cN,function(cc) {
  End<-xglobal[,startsWith(colnames(xglobal),paste0(cc,"_")),drop=FALSE]
  Wex<-tvpgvar_extract_wex(gW[[cc]],xglobal,ncol(End));X<-cbind(1,Wex,tvpgvar_wex_lag(Wex,1L),mlag(End,p))[-seq_len(p),,drop=FALSE];Y<-End[-seq_len(p),,drop=FALSE]
  U<-as.matrix(qr.resid(qr(X),Y))
  do.call(rbind,lapply(seq_len(ncol(Y)),function(j)do.call(rbind,lapply(c(4L,8L),function(lags) {
    u<-U[,j];ul<-mlag(matrix(u,ncol=1),lags);idx<-(lags+1L):length(u);Z<-cbind(X[idx,,drop=FALSE],ul[idx,,drop=FALSE]);qz<-qr(Z)
    if(qz$rank<ncol(Z))return(data.frame(unit=cc,variable=colnames(Y)[j],lags=lags,LM=NA_real_,pvalue=NA_real_,status="rank deficient auxiliary regression"))
    er<-qr.resid(qz,u[idx]);LM<-length(idx)*(1-sum(er^2)/sum((u[idx]-mean(u[idx]))^2))
    data.frame(unit=cc,variable=colnames(Y)[j],lags=lags,LM=LM,pvalue=pchisq(LM,lags,lower.tail=FALSE),status="conditional OLS screen; not exogeneity test or posterior residual check")
  }))))
}))
bg$pvalue_BH<-p.adjust(bg$pvalue,method="BH");write.csv(bg,file.path(out,"residual_BG_screen.csv"),row.names=FALSE)
peaks<-read.csv(file.path(root,"peak_dates.csv"),stringsAsFactors=FALSE);peaks<-peaks[peaks$series==series,,drop=FALSE];chosen<-as.character(peaks$Quarter)
stopifnot(length(chosen)>0,all(chosen %in% Data.setup$quarters[-1]),Data.setup$quarters[2]=="2000Q2")
write.csv(peaks,file.path(out,"selected_peak_dates.csv"),row.names=FALSE)
settings<-list(series=series,description=manifest$label[manifest$series==series],p=p,q=1L,saves=saves,burns=burns,thin=thin,horizon=nhor,shock_pct=shock_pct,
 seed=20261026L+match(series,manifest$series),sample=range(Data.setup$quarters),effective_sample=range(Data.setup$quarters[-1]),
 national_endogenous_variables=70L,external_path_variables=1L,rate="first difference; rate level IRF is draw-wise cumulative",risk="natural log of complete quarterly mean; +10% original level on impact",exogeneity="GPR and lag are external in every national VARX; separate autonomous TVP AR(1) has no macro feedback",interpretation="conditional external-process innovation, not an instrument-identified causal GPR shock",MCMC="one independent chain per category; ESS screens do not verify multi-chain convergence")
dput(settings,file=file.path(out,"settings.R"))
for(f in c("source_manifest.json","snapshot_hashes.json","peak_rules.json","run_manifest.csv"))stopifnot(file.copy(file.path(root,f),out,overwrite=TRUE))
for(f in c("rate_proxy_changes.csv","verified_macro_inventory.csv","READINESS.json"))stopifnot(file.copy(file.path("8.12/country_epu_oil/data_20261022",f),out,overwrite=TRUE))
writeLines(c("One category per fit, conditional on the same frozen macro sample; no EPU, WTI, VIX or other GPR category controls.","Country VARX coefficients condition on current and lagged external GPR; foreign macro controls retain inherited trade-weight assumptions.","Auxiliary external GPR path is an autonomous scalar TVP AR(1). Country innovations have zero effect on this path by imposed block exogeneity.","Path innovation exogeneity and zero cross-block innovation covariance are restrictions, not tested causal identification.","Each IRF freezes coefficients at the category's OWN quarterly level peak; future coefficients are not forecast.","Both valid_all and stable_only are saved. Stability uses the frozen stacked transition matrix; it is not convergence verification.","Four inherited rate-proxy changes and historical trade-weight lineage remain unresolved; see copied READINESS and rate registry."),file.path(out,"LIMITATIONS.txt"))
if(Sys.getenv("PREFLIGHT_ONLY","false")=="true"){cat("PASS external GPR preflight",series,nrow(xglobal)-1,"effective quarters\n");quit(save="no",status=0)}
suppressPackageStartupMessages({library(compiler);library(snowfall);library(Matrix);library(mvtnorm);library(threshtvp)})
RNGkind("L'Ecuyer-CMRG");set.seed(settings$seed);ext.inst<-FALSE;shrink.parm<-list(B_1=8,B_2=.01,kappa0=-.005);BVAR<-cmpfun(BVAR)
CPU<-min(4L,max(1L,as.integer(Sys.getenv("TVPGVAR_CPU","2"))));rng<-vector("list",length(cN));rng[[1]]<-.Random.seed
for(i in 2:length(cN))rng[[i]]<-parallel::nextRNGStream(rng[[i-1]])
sfInit(parallel=TRUE,cpus=CPU)
sfExport(list=list("mlag","BVAR","datahandling","xglobal","gW","Daten","cN","bvartvpm","saves","burns","thin","ext.inst","shrink.parm","rng","tvpgvar_extract_wex","tvpgvar_wex_lag","tvpgvar_safe_inverse","p"))
predDens<-tryCatch(sfLapply(seq_along(cN),function(i){assign(".Random.seed",rng[[i]],envir=.GlobalEnv);tryCatch(BVAR(i,gW=gW,bigx=xglobal,Daten=Daten,cN=cN,nsave=saves,nburn=burns,thin_chain=thin,ext.inst=ext.inst,parms=shrink.parm),error=function(e)stop(paste0("UNIT=",cN[[i]],"; ",conditionMessage(e)),call.=FALSE))}),finally=sfStop())
save(predDens,Data.setup,settings,file=file.path(out,"posterior.rda"))
A<-lapply(predDens,`[[`,"ALPHA");S<-lapply(predDens,`[[`,"SIGMApost")
ess<-do.call(rbind,lapply(seq_along(A),function(i) {
  a<-A[[i]];times<-unique(c(1L,ceiling(dim(a)[1]/2),dim(a)[1]))
  constant<-which(dimnames(a)[[2]]=="constant")[1]
  lag1<-which(dimnames(a)[[2]]=="Ylag1")
  do.call(rbind,lapply(times,function(ti) do.call(rbind,lapply(seq_len(dim(a)[3]),function(j) {
    do.call(rbind,lapply(c(constant,lag1[j]),function(k) {
      z<-as.numeric(a[ti,k,j,]);ok<-all(is.finite(z))&&sd(z)>0
      val<-if(ok)tryCatch(as.numeric(coda::effectiveSize(coda::mcmc(z))),error=function(e)NA_real_)else NA_real_
      data.frame(unit=cN[i],date=Data.setup$quarters[p+ti],equation=dimnames(a)[[3]][j],
        coefficient=if(k==constant)"constant"else"own_lag1",draws=length(z),ESS=val,
        status=if(ok)"single-chain ESS screen; not R-hat"else"constant or nonfinite trace")
    }))
  }))))
}))
write.csv(ess,file.path(out,"posterior_ESS_screen.csv"),row.names=FALSE)
globalG<-lapply(predDens,`[[`,"W");nd<-dim(A[[1]])[4];K<-ncol(xglobal);dates<-Data.setup$quarters[-seq_len(p)]
stopifnot(K==71L,all(vapply(A,function(z)dim(z)[4],integer(1))==nd),all(vapply(A,function(z)dim(z)[1],integer(1))==102L))
rows<-list();stability<-list();failures<-list();raw_by_date<-list()
add_summary<-function(IR,mask,date,response_kind="raw"){
 if(!any(mask))return(NULL)
 do.call(rbind,lapply(seq_len(dim(IR)[1]),function(j)do.call(rbind,lapply(0:nhor,function(h){qs<-quantile(IR[j,h+1,mask],c(.05,.16,.5,.84,.95),names=FALSE);data.frame(series=series,date=date,variable=dimnames(IR)[[1]][j],response_kind=response_kind,filter=filter,horizon=h,median=qs[3],low68=qs[2],high68=qs[4],low90=qs[1],high90=qs[5],draws=sum(mask))}))))
}
for(date in chosen){
 tt<-match(date,dates);rho<-rep(NA_real_,nd);condition<-rho;IRall<-array(NA_real_,c(K,nhor+1,nd),dimnames=list(colnames(xglobal),0:nhor,NULL))
 for(dd in seq_len(nd)){
  fit<-tryCatch({draw_i<-lapply(A,function(z)array(z[,,,dd,drop=FALSE],dim=dim(z)[1:3],dimnames=dimnames(z)[1:3]));sig_i<-lapply(S,function(z)array(z[,,,dd,drop=FALSE],dim=dim(z)[1:3]));get_dominant_gpr_irf_t(tt,draw_i,sig_i,t(xglobal),globalG,cN,horizon=nhor,shock_pct=shock_pct)},error=function(e){failures[[length(failures)+1L]]<<-data.frame(series=series,date=date,draw=dd,error=conditionMessage(e));NULL})
  if(is.null(fit))next
  rho[dd]<-fit$max_eigen_modulus;condition[dd]<-fit$G_condition_number;IRall[,,dd]<-fit$IRF_post
 }
 valid<-is.finite(rho);for(dd in seq_len(nd))if(valid[dd]&&!all(is.finite(IRall[,,dd])))valid[dd]<-FALSE
 stable<-valid&rho<1
 stability[[date]]<-data.frame(series=series,date=date,total=nd,valid=sum(valid),failed=sum(!valid),stable=sum(stable),stable_share=sum(stable)/nd,median_rho=if(any(valid))median(rho[valid])else NA_real_,median_G_condition=if(any(valid))median(condition[valid])else NA_real_)
 for(filter in c("valid_all","stable_only")){
  mask<-if(filter=="stable_only")stable else valid
  rows[[length(rows)+1L]]<-add_summary(IRall,mask,date)
  rr<-grep("_r$",colnames(xglobal));level<-IRall[rr,,,drop=FALSE]
  for(j in seq_along(rr))for(dd in which(mask))level[j,,dd]<-cumsum(level[j,,dd])
  rows[[length(rows)+1L]]<-add_summary(level,mask,date,"rate_level_cumulative")
 }
 raw_by_date[[date]]<-list(IRF=IRall,rho=rho,G_condition=condition,valid=valid,stable=stable)
 cat("IRF complete",series,date,"stable",sum(stable),"/",nd,"\n")
}
save(raw_by_date,settings,file=file.path(out,"irf_draws.rda"))
write.csv(do.call(rbind,stability),file.path(out,"stability.csv"),row.names=FALSE)
tab<-do.call(rbind,rows)
if(!is.null(tab)){tab$significant68<-tab$low68>0|tab$high68<0;tab$significant90<-tab$low90>0|tab$high90<0;write.csv(tab,file.path(out,"irf_summary.csv"),row.names=FALSE)}
if(length(failures))write.csv(do.call(rbind,failures),file.path(out,"failed_draws.csv"),row.names=FALSE)
if(is.null(tab))stop("No valid IRFs; saved posterior and diagnostics")
cat("COMPLETED external GPR",series,"\n")
