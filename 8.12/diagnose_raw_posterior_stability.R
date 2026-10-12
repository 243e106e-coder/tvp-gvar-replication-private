#!/usr/bin/env Rscript
# Read-only reconstruction of historical TVP-GVAR posterior diagnostics.
arg <- commandArgs(trailingOnly=TRUE)
if(length(arg)<3)stop("Usage: script posterior.rda SHOCK outdir")
path<-arg[1];risk<-arg[2];out<-arg[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source("8.12/dominant_gpr_vix_irf.R")
load(path)
stopifnot(exists("predDens"),exists("Data.setup"))
a<-lapply(predDens,function(z)z$ALPHA)
s<-lapply(predDens,function(z)z$SIGMApost)
w<-lapply(predDens,function(z)z$W)
units<-Data.setup$countries
x<-as.matrix(Data.setup$bigx)
p<-nrow(x)-dim(a[[1]])[1]
dates<-Data.setup$quarters[-seq_len(p)]
chosen<-intersect(c("2003Q1","2008Q3","2014Q3","2020Q1","2022Q1","2023Q4"),dates)
nd<-dim(a[[1]])[4]
nsel<-min(nd,as.integer(Sys.getenv("MAX_DIAGNOSTIC_DRAWS","40")))
draws<-unique(as.integer(round(seq(1,nd,length.out=nsel))))
rows<-list();i<-0L
for(date in chosen)for(k in draws){
 i<-i+1L; tt<-match(date,dates)
 rows[[i]]<-tryCatch({
  fit<-get_dominant_gpr_irf_t(tt,
    draw_i=lapply(a,function(z)array(z[,,,k,drop=FALSE],dim=dim(z)[1:3],dimnames=dimnames(z)[1:3])),
    Sig_draw_i=lapply(s,function(z)array(z[,,,k,drop=FALSE],dim=dim(z)[1:3])),
    x=t(x),globalG=w,units=units,horizon=12,shock_pct=10)
  data.frame(shock=risk,date=date,draw=k,rho=fit$max_eigen_modulus,
    stable=fit$max_eigen_modulus<1,
    G_condition=fit$G_condition_number,
    G_min_singular_value=fit$G_min_singular_value,
    max_abs_transition=fit$max_abs_transition,
    GDP_positive_h4=sum(fit$IRF_post[grep("_y$",rownames(fit$IRF_post)),5]>0),
    STOCK_positive_h4=sum(fit$IRF_post[grep("_deq$",rownames(fit$IRF_post)),5]>0),
    error="")
 },error=function(e)data.frame(shock=risk,date=date,draw=k,
     rho=NA_real_,stable=NA,G_condition=NA_real_,
     G_min_singular_value=NA_real_,max_abs_transition=NA_real_,
     GDP_positive_h4=NA_integer_,STOCK_positive_h4=NA_integer_,
     error=conditionMessage(e)))
}
z<-do.call(rbind,rows)
write.csv(z,file.path(out,"draw_level_diagnostics.csv"),row.names=FALSE)
writeLines(c(paste("Model historical shock:",risk),
  paste("Stored posterior draws:",nd),
  paste("Draw subset per date:",length(draws)),
  "Draw subset is systematic and diagnostic, not full posterior inference.",
  "Rho and G condition reconstructed with historical project function.",
  "This does not validate MCMC convergence or raw bilateral trade W provenance."),
  file.path(out,"README.txt"))
if(all(is.na(z$rho)))stop("All selected draws failed reconstruction")
cat("OK",nrow(z),"rows; valid",sum(is.finite(z$rho)),"\n")
