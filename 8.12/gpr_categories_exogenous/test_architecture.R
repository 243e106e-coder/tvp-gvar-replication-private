source("8.12/gpr_categories_exogenous/prepare.R")
source("8.12/gpr_categories_exogenous/irf.R")
# Analytical VARX: every national own lag=.5, exogenous path AR=.6,
# current external coefficient=.2, lag external coefficient=-.1.
# Thus national response at h is .2*ln(1.1)*.6^h exactly, with no feedback.
for(series in c("GPR","SHAREH_CAT_2")) {
  d<-prepare_external_gpr(series);x<-d$bigx;K<-ncol(x);aa<-ss<-vector("list",15)
  for(i in seq_along(d$countries)) {
    cc<-d$countries[i];own<-colnames(x)[startsWith(colnames(x),paste0(cc,"_"))]
    k<-length(own);ex<-nrow(d$gW[[i]])-k;tags<-c("constant",rep("Wex",ex),rep("Wexlag1",ex),rep("Ylag1",k))
    a<-array(0,c(1,length(tags),k),dimnames=list(NULL,tags,own))
    a[1,which(tags=="Ylag1"),]<-diag(if(cc=="GL") .6 else .5,k)
    if(cc!="GL") {a[1,tail(which(tags=="Wex"),1),]<-.2;a[1,tail(which(tags=="Wexlag1"),1),]<--.1}
    aa[[i]]<-a;ss[[i]]<-array(diag(k),c(1,k,k))
  }
  z<-get_dominant_gpr_irf_t(1,aa,ss,t(x),d$gW,d$countries,horizon=4,shock_pct=10)
  stopifnot(max(abs(z$IRF_post[K,]-log1p(.1)*.6^(0:4)))<1e-12,
            max(abs(z$IRF_post[-K,]-matrix(rep(.2*log1p(.1)*.6^(0:4),each=70),70)))<1e-12)
  # Deliberately inject forbidden macro feedback into external G: must fail.
  G<-diag(K);G[K,1]<-.1;dimnames(G)<-list(colnames(x),colnames(x))
  bad<-try(dominant_gpr_struct_irf(G,list(diag(.5,K)),lapply(1:15,function(i)if(i==15)matrix(1,1)else diag(5)),t(x),d$countries),silent=TRUE)
  stopifnot(inherits(bad,"try-error"))
}
cat("PASS external VARX placement, autonomous AR path, analytical propagation, feedback rejection and +10% normalization\n")
