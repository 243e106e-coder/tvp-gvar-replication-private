source("8.12/dominant_gpr_vix_irf.R")
# A national EPU innovation uses the first column of its six-equation Cholesky
# factor. WTI uses the separate GL scalar innovation. Cross-block innovation
# covariance is zero in this inherited blockwise identification.
dominant_gpr_struct_irf <- function(G,F,sig,x,units,horizon=12,shock_var=NULL,shock_pct=10) {
  target <- active_risk_names[active_shock]
  host <- substr(target,1,2);ui<-match(host,units)
  stopifnot(!is.na(ui),active_shock %in% 1:2)
  k<-sum(startsWith(rownames(x),paste0(host,"_")))
  S<-matrix(sig[[ui]],k,k);S<-(S+t(S))/2
  if(any(!is.finite(S))) stop("Nonfinite host covariance")
  L<-t(chol(S))
  own<-rownames(x)[startsWith(rownames(x),paste0(host,"_"))]
  jlocal<-match(target,own)
  u<-unlist(lapply(seq_along(units),function(i) {
    ki<-sum(startsWith(rownames(x),paste0(units[i],"_")))
    if(i==ui) L[,jlocal] else rep(0,ki)
  }),use.names=FALSE)
  raw<-as.vector(solve(G,u));j<-match(target,rownames(x))
  if(!is.finite(raw[j])||abs(raw[j])<1e-12) stop("Cannot normalize impact")
  ir<-matrix(0,nrow(x),horizon+1,dimnames=list(rownames(x),0:horizon))
  ir[,1]<-raw*log1p(shock_pct/100)/raw[j]
  for(h in seq_len(horizon)) for(lag in seq_len(min(length(F),h)))
    ir[,h+1]<-ir[,h+1]+as.vector(F[[lag]]%*%ir[,h-lag+1])
  stopifnot(all(is.finite(ir)))
  ir
}
