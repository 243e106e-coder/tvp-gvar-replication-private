# Each national equation is a conditional TVP-VARX: z and z(-1) are
# external controls. The auxiliary GL row only generates an autonomous TVP
# AR(1) shock path, with NO national variables, trade weights or feedback.
prepare_external_gpr <- function(series, base=NULL, logged=NULL, weights=NULL) {
  root<-"8.12/gpr_categories_exogenous/data_20261026"
  if(is.null(base)) base<-read.csv(file.path(root,"macro_input.csv"),check.names=FALSE)
  if(is.null(logged)) logged<-read.csv(file.path(root,"gpr_quarterly_logs.csv"),check.names=FALSE)
  if(is.null(weights)) {
    w<-read.csv("8.12/Trade_Weights_14_Economies_2000_2014.csv",check.names=FALSE)
    rownames(w)<-w$Reporter;weights<-w[,-1]
  }
  countries<-c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  vars<-c("y","dp","r","de","deq");units<-c(countries,"GL")
  stopifnot(length(series)==1L,series %in% names(logged),!anyDuplicated(base$Quarter),!anyDuplicated(logged$Quarter))
  ix<-match(base$Quarter,logged$Quarter)
  if(anyNA(ix)) stop("External GPR does not cover macro sample and initial lag")
  own<-setNames(lapply(units,function(cc)if(cc=="GL") "GL_gpr" else paste0(cc,"_",vars)),units)
  dat<-base[,c("Quarter",unlist(own[countries],use.names=FALSE)),drop=FALSE]
  dat$GL_gpr<-logged[[series]][ix]
  x<-as.matrix(dat[,-1]);storage.mode(x)<-"double"
  wt<-as.matrix(weights[countries,countries]);storage.mode(wt)<-"double"
  stopifnot(ncol(x)==71L,all(is.finite(x)),all(is.finite(wt)),all(wt>=0),
            max(abs(diag(wt)))<1e-10,max(abs(rowSums(wt)-1))<1e-6)
  wt<-wt/rowSums(wt)
  gW<-setNames(lapply(units,function(cc) {
    k<-length(own[[cc]]);controls<-if(cc=="GL")character()else c(paste0("foreign_",vars),"external_gpr")
    W<-matrix(0,k+length(controls),ncol(x),dimnames=list(c(own[[cc]],controls),colnames(x)))
    W[cbind(seq_len(k),match(own[[cc]],colnames(x)))]<-1
    if(cc!="GL") {
      for(v in seq_along(vars))for(partner in setdiff(countries,cc))W[k+v,match(paste0(partner,"_",vars[v]),colnames(x))]<-wt[cc,partner]
      W[nrow(W),match("GL_gpr",colnames(x))]<-1
    }
    W
  }),units)
  stopifnot(nrow(gW$GL)==1L,all(gW$GL[1,seq_len(70)]==0),gW$GL[1,71]==1,
            identical(colnames(x),unlist(own,use.names=FALSE)))
  list(bigx=x,gW=gW,countries=units,country_units=countries,quarters=dat$Quarter,aligned=dat,
       new.data=x,Wtrade=wt,series=series,risk_names="GL_gpr",placement="EXOGENOUS_GPR")
}
