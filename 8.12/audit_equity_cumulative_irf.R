#!/usr/bin/env Rscript
# Reconstruct the saved posterior IRF draws, then accumulate EQUITY RETURN IRFs
# DRAW BY DRAW (never sum posterior quantiles). Does not fit a model or alter inputs.
# Usage: Rscript 8.12/audit_equity_cumulative_irf.R GPR posterior-inputs output_dir
# Warning: cumulative sum equals log EQUITY PRICE IRF only when EQ_RETURN =
# Delta log(price). If return is an arithmetic return, this is not an exact price IRF.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop('Expected: RISK INPUT_DIR OUTPUT_DIR')
risk <- toupper(args[1]); input <- args[2]; out <- args[3]
stopifnot(risk %in% c('GPR','VIX','EPU','OIL'))
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Output must be empty')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
findfile <- function(n) {
  files <- list.files(input,recursive=TRUE,full.names=TRUE)
  match <- files[basename(files)==n]
  if(length(match)!=1L)stop('Expected one ',n,', found ',length(match))
  match
}
files <- vapply(c('posterior.rda','stability.csv','irf_stable_only.csv','aligned_input.csv'),findfile,character(1))
hashes <- tools::md5sum(files)
env <- new.env(parent=baseenv());load(files[['posterior.rda']],envir=env)
stopifnot(all(c('predDens','Data.setup') %in% ls(env)))
pd <- env$predDens; ds <- env$Data.setup
x <- as.matrix(ds$bigx); units <- as.character(ds$countries)
cn <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
vars <- c('y','dp','r','de','deq')
expected <- c(unlist(lapply(cn,function(cc)paste0(cc,'_',vars))), 'GL_shock')
stopifnot(identical(colnames(x),expected),identical(units,c(cn,'GL')))
A <- lapply(pd,`[[`,'ALPHA'); Sigma <- lapply(pd,`[[`,'SIGMApost'); W <- lapply(pd,`[[`,'W')
nd <- dim(A[[1]])[4]; TT <- dim(A[[1]])[1]; lag <- nrow(x)-TT
stopifnot(lag %in% c(1,2),nd>=1, length(A)==15, length(W)==15)
dates <- as.character(ds$quarters)[-(seq_len(lag))]
st <- read.csv(files[['stability.csv']]); original <- read.csv(files[['irf_stable_only.csv']])
if(anyDuplicated(st$date))stop('Duplicate dates in stability.csv')
nh <- max(original$horizon); equity_idx <- match(paste0(cn,'_deq'),expected)
# Every sampled VAR includes contemporaneous weighted foreign variables.
block <- function(i,t,draw,tag,optional=FALSE) {
  z <- A[[i]]; ids <- which(dimnames(z)[[2]]==tag); k <- dim(z)[3]
  if(!length(ids)) {
    if(optional)return(matrix(0,k,0))
    stop('Missing coefficient tag: ',units[i],'/',tag)
  }
  t(matrix(z[t,ids,,draw,drop=FALSE],nrow=length(ids),ncol=k))
}
local <- function(i,t,draw) {
  k <- dim(A[[i]])[3]; wi <- as.matrix(W[[i]]); ne <- nrow(wi)-k
  l0 <- block(i,t,draw,'Wex',TRUE)
  if(ncol(l0)!=ne)stop('Current foreign coefficient length mismatch: ',units[i])
  H <- lapply(seq_len(lag), function(j) {
    own <- block(i,t,draw,paste0('Ylag',j))
    ext <- block(i,t,draw,paste0('Wexlag',j),TRUE)
    if(ncol(ext)==0 && ne>0 && j>1)ext <- matrix(0,k,ne)
    if(ncol(ext)!=ne || ncol(own)!=k)stop('Lag coefficient length mismatch: ',units[i])
    cbind(own,ext) %*% wi
  })
  list(G=cbind(diag(k),-l0)%*%wi,H=H)
}
spectral_radius <- function(F) {
  k <- nrow(F[[1]]); p <- length(F)
  companion <- if(p==1) F[[1]] else rbind(
    do.call(cbind,F),
    cbind(diag(k*(p-1)),matrix(0,k*(p-1),k)))
  max(Mod(eigen(companion,only.values=TRUE)$values))
}
quant <- function(x)as.numeric(quantile(x,c(.05,.5,.95),na.rm=TRUE,names=FALSE))
rows <- list(); checks <- list(); target <- log1p(.1)
for(date in as.character(st$date)) {
  t <- match(date,dates)
  if(is.na(t))stop('Date not in posterior: ',date)
  cat('Processing ',risk,' ',date,' - ',nd,' saved posterior draws\n',sep='')
  ret <- array(NA_real_,c(14,nh+1,nd))
  stable <- logical(nd)
  for(draw in seq_len(nd)) {
    L <- lapply(seq_along(units),function(i)local(i,t,draw))
    G <- do.call(rbind,lapply(L,`[[`,'G'))
    H <- lapply(seq_len(lag),function(j)do.call(rbind,lapply(L,function(z)z$H[[j]])))
    F <- lapply(H,function(h)solve(G,h))
    rho <- spectral_radius(F)
    stable[draw] <- is.finite(rho) && rho<1
    if(!stable[draw])next
    variance <- as.numeric(Sigma[[15]][t,,,draw])
    if(!is.finite(variance)||variance<=0)stop('Bad GL variance')
    u <- numeric(71);u[71] <- sqrt(variance)
    z0 <- as.vector(solve(G,u))
    if(!is.finite(z0[71]) || z0[71]<=0)stop('Bad global shock normalization')
    impact <- z0*(target/z0[71])
    response <- matrix(NA_real_,71,nh+1)
    response[,1] <- impact
    for(h in seq_len(nh)) {
      z <- numeric(71)
      for(j in seq_len(min(lag,h)))z <- z + as.vector(F[[j]]%*%response[,h-j+1])
      response[,h+1] <- z
    }
    ret[,,draw] <- response[equity_idx,,drop=FALSE]
  }
  expected_stable <- st$stable[match(date,st$date)]
  if(sum(stable)!=expected_stable)stop('Stable draw mismatch ',date,': ',sum(stable),' vs ',expected_stable)
  # Reconcile return IRF posterior summaries with original published CSV.
  for(cc in seq_along(cn))for(h in 0:nh) {
    v <- ret[cc,h+1,stable]
    q <- quant(v)
    # Explicit subset, avoiding ambiguous non-standard evaluation of 'date'.
    saved <- original[original$date==date & original$variable==paste0(cn[cc],'_deq') & original$horizon==h,,drop=FALSE]
    if(nrow(saved)!=1L)stop('Missing/duplicate original equity IRF key')
    orig <- as.numeric(saved[1,c('low90','median','high90')])
    if(max(abs(orig-q))>1e-8*(1+max(abs(orig))))stop('Posterior equity IRF mismatch ',date,'/',cn[cc],'/',h)
    # This is exact for log returns; do not cumulate posterior medians instead.
    cum_draw <- colSums(matrix(ret[cc,seq_len(h+1),stable],nrow=h+1))
    cq <- quant(cum_draw)
    rows[[length(rows)+1L]] <- data.frame(
      shock=risk,date=date,country=cn[cc],horizon=h,stable_draws=sum(stable),
      return_low90=q[1],return_median=q[2],return_high90=q[3],
      cumulative_low90=cq[1],cumulative_median=cq[2],cumulative_high90=cq[3],
      return_median_positive=q[2]>0,cumulative_median_positive=cq[2]>0,
      cumulative_ci_excludes_zero=cq[1]>0||cq[3]<0)
  }
  checks[[length(checks)+1L]] <- data.frame(shock=risk,date=date,total=nd,stable=sum(stable),expected=expected_stable,pass=TRUE)
  rm(ret);gc(verbose=FALSE)
}
outdata <- do.call(rbind,rows)
write.csv(outdata,file.path(out,'equity_return_vs_cumulative_drawlevel.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,'posterior_reconciliation.csv'),row.names=FALSE)
agg <- aggregate(cbind(
  negative_return=as.integer(outdata$return_median<0),
  negative_cumulative=as.integer(outdata$cumulative_median<0),
  negative_cumulative_ci90=as.integer(outdata$cumulative_high90<0)),
  by=list(shock=outdata$shock,date=outdata$date,horizon=outdata$horizon),FUN=sum)
write.csv(agg,file.path(out,'equity_sign_counts.csv'),row.names=FALSE)
# PDF pages: each date, with fourteen countries and paired IRFs.
pdf(file.path(out,'equity_return_cumulative_comparison.pdf'),width=13,height=10)
for(date in unique(outdata$date)) {
  par(mfrow=c(4,4),mar=c(3,3,3,1),oma=c(0,0,3,0))
  for(cc in cn) {
    zz <- outdata[outdata$date==date & outdata$country==cc,]
    ylim <- range(c(zz$return_low90,zz$return_high90,zz$cumulative_low90,zz$cumulative_high90))
    plot(zz$horizon,zz$cumulative_median,type='l',col='blue',lwd=2,ylim=ylim,
         xlab='h (quarters)',ylab='IRF',main=cc)
    lines(zz$horizon,zz$cumulative_low90,col='blue',lty=3)
    lines(zz$horizon,zz$cumulative_high90,col='blue',lty=3)
    lines(zz$horizon,zz$return_median,col='black',lwd=1)
    abline(h=0,col='gray60',lty=2)
  }
  mtext(paste(risk,date,'Black = return; blue = cumulative conditional on log-return definition'),outer=TRUE,cex=1)
}
dev.off()
writeLines(c(
  'EQUITY IRF AUDIT — NO RE-ESTIMATION',
  'Source: saved TVP-GVAR posterior.rda for the specified shock; all stable draws only.',
  'All posterior return IRF quantiles and stable counts reconciled with saved original CSVs.',
  'Cumulative IRFs are calculated within each posterior draw, then quantiles are computed.',
  'IMPORTANT: the cumulated IRF is a log-stock-PRICE response ONLY IF EQ_RETURN equals Delta log(stock price).',
  'The upstream model builder maps deq to EQ_RETURN; it does not establish the EQ_RETURN calculation.',
  'If EQ_RETURN is 100*Delta log(price), divide cumulative output by 100 to get log-price units.',
  'If EQ_RETURN is arithmetic return, cumulation is not an exact transformed price-level IRF.',
  'Changing the estimated dependent variable to log price LEVEL requires a fresh model and stationarity/cointegration work.',
  'Positive cumulative median is not a failure; a 90% credible interval spanning zero is inconclusive.',
  'No inputs are changed; no imputation, sign restrictions, or artificial flips.'
),file.path(out,'README_METHOD.txt'))
if(!identical(unname(hashes),unname(tools::md5sum(files))))stop('Source file changed')
cat('PASS: completed draw-level cumulative equity IRFs in ',out,'\n',sep='')
