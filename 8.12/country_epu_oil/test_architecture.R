source("8.12/country_epu_oil/prepare.R")
source("8.12/country_epu_oil/irf.R")
# Analytical two-channel test: a WTI innovation reaches the national EPU
# equation contemporaneously via G, and then propagates with lag 0.5.
for(ec in c("CN","US")) for(p in 1:2) {
  cc<-c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  vars<-c("y","dp","r","de","deq")
  base<-as.data.frame(matrix(1,20,71));names(base)<-c("Quarter",unlist(lapply(cc,function(c)paste0(c,"_",vars))))
  base$Quarter<-paste0(2000+rep(0:4,each=4),"Q",rep(1:4,5))
  risk<-data.frame(Quarter=base$Quarter,US_EPU_months=3,CN_EPU_months=3,WTI_months=3,
                   LN_CN_EPU=1,LN_US_EPU=1,LN_WTI=1)
  w<-matrix(1/13,14,14,dimnames=list(cc,cc));diag(w)<-0
  d<-prepare_country_epu_oil(ec,base,risk,w);x<-d$bigx;K<-ncol(x)
  stopifnot(sum(startsWith(colnames(x),paste0(ec,"_")))==6,
            identical(colnames(x)[K],"OO_wti"),
            !any(grepl("_wti$",colnames(x)[1:(K-1)])))
  aa<-ss<-vector("list",length(d$countries))
  for(i in seq_along(d$countries)) {
    unit<-d$countries[i];own<-colnames(x)[startsWith(colnames(x),paste0(unit,"_"))]
    k<-length(own);ex<-nrow(d$gW[[i]])-k
    tags<-c("constant",rep("Wex",ex),rep("Wexlag1",ex),
            unlist(lapply(seq_len(p),function(j)rep(paste0("Ylag",j),k))))
    a<-array(0,c(1,length(tags),k),dimnames=list(NULL,tags,own))
    a[1,which(tags=="Ylag1"),]<-diag(.5,k)
    if(unit==ec) {
      oilrow<-which(rownames(d$gW[[i]])=="OO_wti")-k
      a[1,which(tags=="Wex")[oilrow],1]<-.2
      # Make reduced lag F exactly .5 I after contemporaneous solve.
      a[1,which(tags=="Wexlag1")[oilrow],1]<--.1
    }
    aa[[i]]<-a;ss[[i]]<-array(diag(k),c(1,k,k))
  }
  active_risk_names<-d$risk_names
  for(sh in 1:2) {
    active_shock<-sh
    ans<-get_dominant_gpr_irf_t(1,aa,ss,t(x),d$gW,d$countries,horizon=3,shock_pct=10)
    j<-match(active_risk_names[sh],colnames(x))
    stopifnot(max(abs(ans$IRF_post[j,]-log1p(.1)*.5^(0:3)))<1e-12)
    if(sh==2) stopifnot(max(abs(ans$IRF_post[match(paste0(ec,"_epu"),colnames(x)),]-.2*log1p(.1)*.5^(0:3)))<1e-12)
  }
}
cat("PASS: own country EPU, separate oil scalar, contemporaneous transmission, p=1/2 and shock normalization\n")
