# Deterministic stack/normalization tests for BOTH placements and p=1/2.
source("8.12/prepare_epu_wti.R")
source("8.12/epu_wti_irf.R")
for(ec in c("US","CN")) for(pl in c("GLOBAL","US")) for(lags in 1:2) {
  d<-prepare_epu_wti(ec,pl);x<-d$bigx;units<-d$countries;K<-ncol(x)
  aa<-ss<-vector("list",length(units))
  for(i in seq_along(units)) {
    cc<-units[i];own<-colnames(x)[startsWith(colnames(x),paste0(cc,"_"))]
    k<-length(own);ex<-nrow(d$gW[[cc]])-k
    tags<-c("constant",rep("Wex",ex),rep("Wexlag1",ex),
            unlist(lapply(seq_len(lags),function(j)rep(paste0("Ylag",j),k))))
    A<-array(0,c(1,length(tags),k),dimnames=list(NULL,tags,own))
    A[1,which(tags=="Ylag1"),]<-diag(.5,k)
    aa[[i]]<-A;ss[[i]]<-array(diag(k),c(1,k,k))
  }
  shock_host<-d$host;active_risk_names<-d$risk_names
  for(sh in 1:2) {
    active_shock<-sh
    ans<-get_dominant_gpr_irf_t(1,aa,ss,t(x),d$gW,units,horizon=3,shock_pct=10)
    j<-match(active_risk_names[sh],colnames(x))
    stopifnot(max(abs(ans$IRF_post[j,]-log1p(.1)*(.5)^(0:3)))<1e-12,
              max(abs(ans$IRF_post[-j,]))<1e-12,
              abs(ans$max_eigen_modulus-.5)<1e-10)
  }
}
cat("PASS: all four mappings, p=1/2 stacking, both shocks, +10% normalization and analytical AR(1) paths\n")
