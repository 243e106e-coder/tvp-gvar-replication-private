# Functional smoke test of the actual sampler and scalar covariance assembly.
# Twenty draws are only a shape/numeric check, never convergence evidence.
Sys.setenv(TVPGVAR_P="1")
source("8.12/country_epu_oil/patch_estimator.R")
source("R/BVAR_country_oil.r")
source("R/auxilliary_functions_tvp.r")
suppressPackageStartupMessages({library(threshtvp);library(snowfall)})
set.seed(20261022)
x<-as.numeric(arima.sim(list(ar=.5),n=80,sd=.2))
big<-matrix(x,ncol=1,dimnames=list(NULL,"GL_wti"))
w<-matrix(1,1,1,dimnames=list("GL_wti","GL_wti"))
sfInit(parallel=FALSE,cpus=1)
z<-tryCatch(BVAR(1,gW=list(w),bigx=big,Daten=big,cN="GL",nsave=20,nburn=20,thin_chain=1,
            ext.inst=FALSE,parms=list(B_1=8,B_2=.01,kappa0=-.005)),finally=sfStop())
stopifnot(identical(dim(z$ALPHA),c(79L,2L,1L,20L)),
          identical(dim(z$SIGMApost),c(79L,1L,1L,20L)),
          all(is.finite(z$ALPHA)),all(is.finite(z$SIGMApost)),all(z$SIGMApost>0))
cat("PASS actual scalar GL sampler and covariance assembly; 20-draw functional smoke only\n")
