# Two separate systems: CN EPU in CN, or US EPU in US; WTI in OO.
# OO is an autonomous one-equation oil block. This weak-exogeneity assumption
# must be assessed before interpreting a WTI innovation as economic causality.
prepare_country_epu_oil <- function(epu_country, base, risk, wt) {
  stopifnot(epu_country %in% c("CN","US"))
  countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  vars <- c("y","dp","r","de","deq")
  units <- c(countries,"OO")
  risk_names <- c(paste0(epu_country,"_epu"),"OO_wti")
  stopifnot(!anyDuplicated(base$Quarter),!anyDuplicated(risk$Quarter))
  ix <- match(base$Quarter,risk$Quarter)
  if(anyNA(ix)) stop("Risk data do not cover all macro and presample quarters")
  stopifnot(all(risk[ix,c(paste0(epu_country,"_EPU_months"),"WTI_months")]==3))
  own <- setNames(lapply(units,function(cc) {
    if(cc=="OO") "OO_wti" else c(if(cc==epu_country) risk_names[1],paste0(cc,"_",vars))
  }),units)
  dat <- base[,c("Quarter",unlist(lapply(countries,function(cc)paste0(cc,"_",vars))))]
  dat[[risk_names[1]]] <- risk[[paste0("LN_",epu_country,"_EPU")]][ix]
  dat$OO_wti <- risk$LN_WTI[ix]
  dat <- dat[,c("Quarter",unlist(own,use.names=FALSE))]
  x <- as.matrix(dat[,-1]);storage.mode(x)<-"double"
  wt <- as.matrix(wt[countries,countries]);storage.mode(wt)<-"double"
  stopifnot(all(is.finite(x)),ncol(x)==72,all(is.finite(wt)),all(wt>=0),
            max(abs(diag(wt)))<1e-10,max(abs(rowSums(wt)-1))<1e-6)
  wt <- wt/rowSums(wt)
  gW <- setNames(lapply(units,function(cc) {
    k <- length(own[[cc]])
    controls <- if(cc=="OO") character() else c(paste0("foreign_",vars),
                    if(cc!=epu_country) risk_names[1],"OO_wti")
    Wi <- matrix(0,k+length(controls),ncol(x),dimnames=list(c(own[[cc]],controls),colnames(x)))
    Wi[cbind(seq_len(k),match(own[[cc]],colnames(x)))] <- 1
    if(cc!="OO") {
      for(v in seq_along(vars)) for(partner in setdiff(countries,cc))
        Wi[k+v,match(paste0(partner,"_",vars[v]),colnames(x))] <- wt[cc,partner]
      j <- k+5L
      if(cc!=epu_country) {j<-j+1L;Wi[j,match(risk_names[1],colnames(x))]<-1}
      Wi[nrow(Wi),match("OO_wti",colnames(x))] <- 1
    }
    Wi
  }),units)
  stopifnot(identical(unlist(own,use.names=FALSE),colnames(x)),nrow(gW$OO)==1L)
  list(bigx=x,gW=gW,countries=units,country_units=countries,quarters=dat$Quarter,
       aligned=dat,new.data=x,risk_names=risk_names,host=epu_country,
       epu_country=epu_country,placement="COUNTRY_OIL",Wtrade=wt)
}

prepare_epu_wti <- function(epu_country,placement="COUNTRY_OIL") {
  stopifnot(placement=="COUNTRY_OIL")
  base <- read.csv("data/country_epu_oil/model_input.csv",check.names=FALSE)
  risk <- read.csv("data/country_epu_oil/risk_quarterly.csv",check.names=FALSE)
  w <- read.csv("8.12/Trade_Weights_14_Economies_2000_2014.csv",check.names=FALSE)
  rownames(w)<-w$Reporter
  prepare_country_epu_oil(epu_country,base,risk,w[,-1])
}
