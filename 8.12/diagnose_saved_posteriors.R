#!/usr/bin/env Rscript
# Saved-posterior diagnostics; base R only; no estimation or data modification.
a <- commandArgs(TRUE)
if(length(a)<3) stop('Usage: Rscript diagnose_saved_posteriors.R RISK INPUT_DIR OUTPUT_DIR')
risk <- toupper(a[1]); input <- normalizePath(a[2], mustWork=TRUE); out <- a[3]
stopifnot(risk %in% c('GPR','VIX','EPU','OIL'))
parent <- normalizePath(dirname(out),mustWork=FALSE)
if(startsWith(paste0(parent,'/'),paste0(input,'/')) || normalizePath(out,mustWork=FALSE)==input) stop('Output cannot be inside input')
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop('Output directory must be empty')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
findfile <- function(n) {
 z <- list.files(input,recursive=TRUE,full.names=TRUE); z <- z[basename(z)==n]
 if(length(z)!=1) stop('Expected exactly one ',n,', found ',length(z)); z
}
paths <- vapply(c('posterior.rda','stability.csv','irf_stable_only.csv','provenance.txt','aligned_input.csv'),findfile,character(1))
hashes <- unname(tools::md5sum(paths))
write.csv(data.frame(file=paths,md5=hashes),file.path(out,'input_checksums.csv'),row.names=FALSE)
e <- new.env(parent=baseenv()); load(paths['posterior.rda'],envir=e)
if(!all(c('predDens','Data.setup') %in% ls(e))) stop('Saved objects missing')
pd <- e$predDens; ds <- e$Data.setup
x <- as.matrix(ds$bigx); qq <- as.character(ds$quarters); units <- as.character(ds$countries)
countries <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
vars <- c('y','dp','r','de','deq')
expected <- c(unlist(lapply(countries,function(c)paste0(c,'_',vars))),'GL_shock')
stopifnot(identical(units,c(countries,'GL')),identical(colnames(x),expected),ncol(x)==71,nrow(x)==length(qq),all(is.finite(x)))
A <- lapply(pd,`[[`,'ALPHA'); S <- lapply(pd,`[[`,'SIGMApost'); W <- lapply(pd,`[[`,'W')
stopifnot(length(A)==15,length(S)==15,length(W)==15)
TT <- dim(A[[1]])[1]; nd <- dim(A[[1]])[4]; p <- nrow(x)-TT
stopifnot(p %in% 1:2,nd>1)
dates <- qq[-seq_len(p)]
prov <- readLines(paths['provenance.txt'],warn=FALSE)
line <- grep('^IRF normalization:',prov,value=TRUE)
if(length(line)!=1 || !grepl('+10% log level impact',line,fixed=TRUE)) stop('Expected saved +10% normalization')
for(i in seq_along(A)) {
 stopifnot(length(dim(A[[i]]))==4,dim(A[[i]])[1]==TT,dim(A[[i]])[4]==nd,dim(S[[i]])[1]==TT,dim(S[[i]])[4]==nd)
 if(is.null(dimnames(A[[i]])[[2]]) || is.null(dimnames(A[[i]])[[3]])) stop('Missing ALPHA labels: ',units[i])
 k <- dim(A[[i]])[3]; wi <- as.matrix(W[[i]])
 stopifnot(ncol(wi)==71,nrow(wi)>=k,all(is.finite(wi)))
 if(!identical(colnames(wi),expected)) stop('Weight column order mismatch: ',units[i])
 if(i<=14) {
  stopifnot(k==5,nrow(wi)==11)
  if(max(abs(wi[11,]-diag(71)[71,]))>1e-12) stop('Global-control row mapping mismatch: ',units[i])
 } else stopifnot(k==1,nrow(wi)==1)
 own <- match(dimnames(A[[i]])[[3]],expected)
 if(anyNA(own) || max(abs(wi[seq_len(k),,drop=FALSE]-diag(71)[own,,drop=FALSE]))>1e-12) stop('Own-variable mapping mismatch: ',units[i])
}
aligned <- read.csv(paths['aligned_input.csv'],check.names=FALSE)
stopifnot(identical(as.character(aligned$Quarter),qq),identical(names(aligned)[-1],expected))
if(max(abs(as.matrix(aligned[,-1])-x))>1e-10) stop('Saved input differs from posterior input')
ss <- read.csv(paths['stability.csv'],stringsAsFactors=FALSE)
original <- read.csv(paths['irf_stable_only.csv'],stringsAsFactors=FALSE)
stopifnot(all(toupper(ss$shock)==risk),all(toupper(original$shock)==risk),all(ss$date %in% dates),all(ss$total==nd))
if(anyDuplicated(ss$date) || anyDuplicated(original[,c('date','variable','horizon')])) stop('Duplicate saved keys')
nh <- max(original$horizon); target <- log1p(.1); gl <- 71L
rho <- function(F) {
 k <- nrow(F[[1]]); pp <- length(F)
 C <- if(pp==1) F[[1]] else rbind(do.call(cbind,F),cbind(diag(k*(pp-1)),matrix(0,k*(pp-1),k)))
 max(Mod(eigen(C,only.values=TRUE)$values))
}
block <- function(i,tt,dd,tag,optional=FALSE) {
 z <- A[[i]]; ids <- which(dimnames(z)[[2]]==tag); k <- dim(z)[3]
 if(!length(ids)) {if(optional)return(matrix(0,k,0)); stop('Missing coefficient tag: ',units[i],'/',tag)}
 t(matrix(z[tt,ids,,dd,drop=FALSE],nrow=length(ids),ncol=k))
}
local <- function(i,tt,dd) {
 k <- dim(A[[i]])[3]; wi <- as.matrix(W[[i]]); ex <- nrow(wi)-k
 l0 <- block(i,tt,dd,'Wex',TRUE); stopifnot(ncol(l0)==ex)
 theta <- lapply(seq_len(p),function(j)block(i,tt,dd,paste0('Ylag',j)))
 lagx <- lapply(seq_len(p),function(j) {
  z <- block(i,tt,dd,paste0('Wexlag',j),TRUE)
  if(!ncol(z) && ex>0 && j>1)z <- matrix(0,k,ex)
  stopifnot(ncol(z)==ex,ncol(theta[[j]])==k); z
 })
 list(G=cbind(diag(k),-l0)%*%wi,H=lapply(seq_len(p),function(j)cbind(theta[[j]],lagx[[j]])%*%wi),ownG=cbind(diag(k),matrix(0,k,ex))%*%wi,ownH=lapply(theta,function(z)cbind(z,matrix(0,k,ex))%*%wi),theta=theta,l0=l0)
}
build <- function(L) list(G=do.call(rbind,lapply(L,`[[`,'G')),H=lapply(seq_len(p),function(j)do.call(rbind,lapply(L,function(z)z$H[[j]]))))
rows <- list(); locrows <- list(); abrows <- list(); impactrows <- list(); reconcile <- list(); traces <- list(); errors <- list()
append_error <- function(d,dd,stage,msg) {errors[[length(errors)+1L]] <<- data.frame(shock=risk,date=d,draw=dd,stage=stage,error=msg)}
for(d in ss$date) {
 cat('Diagnosing',risk,d,'all',nd,'saved draws\n'); flush.console()
 tt <- match(d,dates); IR <- array(NA_real_,c(71,nh+1,nd)); stable <- rep(FALSE,nd)
 for(dd in seq_len(nd)) {
  tryCatch({
   L <- lapply(seq_along(units),function(i)local(i,tt,dd)); M <- build(L)
   G <- M$G; H <- M$H; singular <- svd(G,nu=0,nv=0)$d
   cond <- max(singular)/min(singular); rc <- rcond(G)
   F <- lapply(H,function(h)solve(G,h)); rr <- rho(F); stable[dd] <- is.finite(rr)&&rr<1
   lr <- vapply(L,function(z)rho(z$theta),numeric(1))
   u <- numeric(71); variance <- as.numeric(S[[15]][tt,,,dd]); stopifnot(length(variance)==1,is.finite(variance),variance>0)
   u[gl] <- sqrt(variance); raw <- as.vector(solve(G,u)); stopifnot(is.finite(raw[gl]),raw[gl]>0)
   impact <- raw*(target/raw[gl]); injected <- u*(target/raw[gl])
   ir <- matrix(NA_real_,71,nh+1); ir[,1] <- impact
   for(h in seq_len(nh)) {z <- numeric(71); for(j in seq_len(min(p,h)))z <- z+as.vector(F[[j]]%*%ir[,h-j+1]);ir[,h+1] <- z}
   IR[,,dd] <- ir
   rows[[length(rows)+1]] <- data.frame(shock=risk,date=d,draw=dd,global_rho=rr,stable=stable[dd],G_condition=cond,G_rcond=rc,min_singular=min(singular),max_abs_F=max(vapply(F,function(z)max(abs(z)),numeric(1))),max_local_rho=max(lr),GL_rho=lr[15],solve_residual=max(abs(G%*%impact-injected)),GL_impact=impact[gl],GL_target_error=abs(impact[gl]-target))
   for(i in seq_along(units)) {
    locrows[[length(locrows)+1]] <- data.frame(shock=risk,date=d,draw=dd,unit=units[i],own_lag_rho=lr[i],global_stable=stable[dd])
    # Diagnostic replacement: remove this block's current and lagged external controls.
    alt <- L; alt[[i]]$G <- L[[i]]$ownG; alt[[i]]$H <- L[[i]]$ownH
    ar <- tryCatch({mm <- build(alt);rho(lapply(mm$H,function(h)solve(mm$G,h)))},error=function(e)NA_real_)
    abrows[[length(abrows)+1]] <- data.frame(shock=risk,date=d,draw=dd,unit=units[i],baseline_rho=rr,without_unit_external_rho=ar,rho_reduction=rr-ar,unstable_to_stable=is.finite(ar)&&rr>=1&&ar<1)
    # Across-draw traces of one own-lag coefficient and scalar global coefficient.
    traces[[length(traces)+1]] <- data.frame(shock=risk,date=d,draw=dd,unit=units[i],parameter='own_first_equation_lag1',value=L[[i]]$theta[[1]][1,1])
   }
   for(i in seq_along(countries)) {
    ix <- match(paste0(countries[i],'_deq'),expected)
    direct <- as.numeric(L[[i]]$l0[5,ncol(L[[i]]$l0)])*target
    impactrows[[length(impactrows)+1]] <- data.frame(shock=risk,date=d,draw=dd,country=countries[i],global_stable=stable[dd],equity_impact=impact[ix],direct_global_component=direct,foreign_network_component=impact[ix]-direct,equity_h1=ir[ix,2])
   }
  },error=function(e)append_error(d,dd,'reconstruction',conditionMessage(e)))
 }
 obs <- ss[ss$date==d,]; counted <- sum(stable)
 reconcile[[length(reconcile)+1]] <- data.frame(shock=risk,date=d,check='stable_draw_count',recomputed=counted,saved=obs$stable,abs_error=abs(counted-obs$stable),pass=counted==obs$stable)
 if(counted>0) {
  saved <- original[original$date==d,]
  for(j in seq_len(nrow(saved))) {
   ix <- match(saved$variable[j],expected); h <- saved$horizon[j]
   if(is.na(ix))stop('Unknown saved variable')
   qs <- quantile(IR[ix,h+1,stable],c(.05,.5,.95),names=FALSE)
   old <- as.numeric(unlist(saved[j,c('low90','median','high90')],use.names=FALSE)); err <- max(abs(qs-old))
   reconcile[[length(reconcile)+1]] <- data.frame(shock=risk,date=d,check=paste(saved$variable[j],h,sep='/'),recomputed=qs[2],saved=old[2],abs_error=err,pass=err<=1e-8*(1+max(abs(old))))
  }
 }
 rm(IR); gc(verbose=FALSE)
}
write <- function(z,n) {if(length(z))write.csv(do.call(rbind,z),file.path(out,n),row.names=FALSE)}
write(rows,'draw_stability_and_identification.csv');write(locrows,'local_block_stability.csv');write(abrows,'network_ablation.csv');write(impactrows,'equity_impact_decomposition.csv');write(reconcile,'original_csv_reconciliation.csv');write(traces,'selected_parameter_traces.csv');write(errors,'reconstruction_errors.csv')
if(!length(rows))stop('No reconstructed draws; see reconstruction_errors.csv')
r <- do.call(rbind,rows); la <- do.call(rbind,locrows); ab <- do.call(rbind,abrows)
summary <- do.call(rbind,lapply(split(la,interaction(la$date,la$unit,drop=TRUE)),function(z) {
 az <- ab[ab$date==z$date[1]&ab$unit==z$unit[1],]; bad <- is.finite(az$baseline_rho)&az$baseline_rho>=1
 data.frame(shock=risk,date=z$date[1],unit=z$unit[1],valid_draws=nrow(z),own_unstable_share=mean(z$own_lag_rho>=1),median_own_rho=median(z$own_lag_rho),baseline_unstable_draws=sum(bad),ablation_rescue_share=if(any(bad))mean(az$unstable_to_stable[bad])else NA_real_,median_rho_reduction=median(az$rho_reduction,na.rm=TRUE))
}))
write.csv(summary,file.path(out,'block_diagnostic_summary.csv'),row.names=FALSE)
# A rough single-chain autocorrelation diagnostic, not R-hat or convergence proof.
ess <- function(z) {
 n <- length(z); if(n<4 || sd(z)==0)return(c(lag1=NA,ESS=NA))
 ac <- as.numeric(acf(z,lag.max=min(n-1,100),plot=FALSE)$acf)[-1]
 np <- floor(length(ac)/2); pairs <- ac[2*seq_len(np)-1]+ac[2*seq_len(np)]
 first <- which(pairs<=0); if(length(first))pairs <- head(pairs,first[1]-1)
 c(lag1=ac[1],ESS=min(n,n/max(1,1+2*sum(pairs))))
}
tr <- do.call(rbind,traces)
et <- do.call(rbind,lapply(split(tr,interaction(tr$date,tr$unit,drop=TRUE)),function(z) {
 z <- z[order(z$draw),]; gap <- nrow(z)!=nd || !identical(z$draw,seq_len(nd)); v <- if(gap)c(lag1=NA,ESS=NA)else ess(z$value)
 data.frame(shock=risk,date=z$date[1],unit=z$unit[1],parameter=z$parameter[1],draws=nrow(z),chain_has_gaps=gap,lag1_acf=v[1],approx_single_chain_ESS=v[2])
}))
write.csv(et,file.path(out,'selected_parameter_autocorrelation.csv'),row.names=FALSE)
pdf(file.path(out,'posterior_diagnostic_plots.pdf'),width=11,height=8)
for(d in ss$date) {
 z <- r[r$date==d,]; par(mfrow=c(2,2))
 plot(z$draw,z$global_rho,type='l',main=paste(risk,d,'global spectral radius'),xlab='Saved draw',ylab='rho');abline(h=1,col='red',lty=2)
 plot(z$draw,log10(pmax(z$G_condition,1)),type='l',main='Contemporaneous G conditioning',xlab='Saved draw',ylab='log10 condition number')
 boxplot(own_lag_rho~unit,data=la[la$date==d,],las=2,main='Own-block dynamics (not full GVAR)',ylab='rho');abline(h=1,col='red',lty=2)
 boxplot(rho_reduction~unit,data=ab[ab$date==d,],las=2,main='Remove external controls: rho reduction',ylab='Baseline minus ablation');abline(h=0,lty=2)
}
dev.off()
writeLines(c('# Saved-posterior diagnostics',paste('Shock:',risk),'No MCMC, coefficient refit, missing-data imputation, or raw input changes.',
'Reconstruction follows dominant_gpr_vix_irf.R at source commit 1ff13aa298148620728452f5f97d33f5b9443a0b, with the scalar GL shock override used by compare_four_shocks.R.',
'First check original_csv_reconciliation.csv. Any failed check blocks interpretation; reconstruction_errors.csv records failed draws.',
'Own-block spectral radii omit external feedback; stable own blocks do not establish global stability.',
'Network ablation removes a unit equation block\x27s contemporaneous and lagged external controls while retaining all posterior coefficients. It is a diagnostic, not re-estimation, a causal contribution, or a remedy.',
'Equity impact = direct GL coefficient times normalized GL impact + contemporaneous foreign-network component. This is an algebraic decomposition; h1 includes lag propagation.',
'Selected parameter traces and approximate single-chain ESS cover one own-lag coefficient per unit/date only; they cannot establish whole-model convergence. R-hat requires multiple independent chains and is not reported.',
'Input order and equality are audited. This cannot independently certify the economic meaning or upstream transformation of the original stock-price series.',
'No claims that positive equity responses or low stability retention are resolved.'),file.path(out,'README.md'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
stopifnot(identical(hashes,unname(tools::md5sum(paths))))
rec <- do.call(rbind,reconcile)
if(length(errors) || any(!rec$pass))stop('Diagnostics saved, but reconstruction/reconciliation failed; inspect output before interpreting')
cat('PASS: existing IRFs and stable counts reconciled; input checksums unchanged.\n')
