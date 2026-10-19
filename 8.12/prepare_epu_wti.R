# Shared data/mapping builder for the four EPU+WTI experiments.
prepare_epu_wti <- function(epu_country, placement) {
  stopifnot(epu_country %in% c("US","CN"), placement %in% c("GLOBAL","US"))
  countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  vars <- c("y","dp","r","de","deq")
  base <- read.csv("data/rate_diff/model_input_rate_diff.csv",check.names=FALSE)
  risk <- read.csv("8.12/epu_wti_quarterly.csv",check.names=FALSE)
  stopifnot(!anyDuplicated(base$Quarter),!anyDuplicated(risk$Quarter))
  qi <- function(q) as.integer(substr(q,1,4))*4L+as.integer(substr(q,6,6))
  stopifnot(all(grepl("^[12][0-9]{3}Q[1-4]$",base$Quarter)),all(diff(qi(base$Quarter))==1))
  ix <- match(base$Quarter,risk$Quarter)
  if(anyNA(ix)) stop("EPU/WTI does not cover the common macro sample")
  stopifnot(all(risk[ix,c("US_EPU_months","CN_EPU_months","WTI_months")]==3))
  host <- if(placement=="GLOBAL") "GL" else "US"
  risk_names <- paste0(host,c("_epu","_wti"))
  risk_values <- as.matrix(risk[ix,c(paste0("LN_",epu_country,"_EPU"),"LN_WTI")])
  units <- if(placement=="GLOBAL") c(countries,"GL") else countries
  own <- setNames(lapply(units,function(cc) {
    if(cc=="GL") risk_names else if(cc=="US"&&placement=="US")
      c(risk_names,paste0(cc,"_",vars)) else paste0(cc,"_",vars)
  }),units)
  macro_cols <- unlist(lapply(countries,function(cc) paste0(cc,"_",vars)))
  dat <- base[,c("Quarter",macro_cols)]
  for(j in 1:2) dat[[risk_names[j]]] <- risk_values[,j]
  ordered <- unlist(own,use.names=FALSE)
  dat <- dat[,c("Quarter",ordered)]
  x <- as.matrix(dat[,-1]); storage.mode(x)<-"double"
  stopifnot(ncol(x)==72,all(is.finite(x)))
  w <- read.csv("8.12/Trade_Weights_14_Economies_2000_2014.csv",
                fileEncoding="UTF-8-BOM",check.names=FALSE)
  stopifnot("Reporter"%in%names(w),!anyDuplicated(w$Reporter))
  rownames(w)<-w$Reporter
  wt<-as.matrix(w[countries,countries]);storage.mode(wt)<-"double"
  stopifnot(all(is.finite(wt)),all(wt>=0),all(abs(diag(wt))<1e-10),
            all(abs(rowSums(wt)-1)<1e-6))
  wt<-wt/rowSums(wt)
  K<-ncol(x);gW<-setNames(vector("list",length(units)),units)
  for(cc in units) {
    k<-length(own[[cc]])
    nex<-if(cc=="GL") 0L else if(cc=="US"&&placement=="US") 5L else 7L
    Wi<-matrix(0,k+nex,K,dimnames=list(NULL,colnames(x)))
    Wi[cbind(seq_len(k),match(own[[cc]],colnames(x)))]<-1
    labels<-own[[cc]]
    if(cc!="GL") {
      for(v in seq_along(vars)) for(partner in setdiff(countries,cc)) {
        Wi[k+v,match(paste0(partner,"_",vars[v]),colnames(x))]<-wt[cc,partner]
      }
      labels<-c(labels,paste0("foreign_",vars))
      if(nex==7L) {
        Wi[cbind(k+6:7,match(risk_names,colnames(x)))]<-1
        labels<-c(labels,"common_epu","common_wti")
      }
    }
    rownames(Wi)<-labels;gW[[cc]]<-Wi
    stopifnot(identical(colnames(x)[substr(colnames(x),1,2)==cc],own[[cc]]))
  }
  # Own rows stacked in exactly the order of xglobal; no duplication.
  stopifnot(identical(unlist(own,use.names=FALSE),colnames(x)),
            sum(vapply(own,length,integer(1)))==K)
  if(placement=="GLOBAL") stopifnot(nrow(gW$GL)==2,all(gW$GL[,!startsWith(colnames(x),"GL_")]==0))
  list(bigx=x,gW=gW,countries=units,country_units=countries,quarters=dat$Quarter,
       aligned=dat,new.data=x,host=host,risk_names=risk_names,
       epu_country=epu_country,placement=placement,Wtrade=wt)
}
