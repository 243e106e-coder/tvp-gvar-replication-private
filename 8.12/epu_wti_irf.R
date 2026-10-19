# Generalized recursive injection for GLOBAL (2 equations) and US (7).
# The inherited stack/transition reconstruction is unchanged.
source("8.12/dominant_gpr_vix_irf.R")
dominant_gpr_struct_irf <- function(G,F,sig,x,units,horizon=12,
                                    shock_var=NULL,shock_pct=10) {
  ui<-match(shock_host,units)
  S<-as.matrix(sig[[ui]]); S<-(S+t(S))/2
  if(any(!is.finite(S))) stop("Invalid host covariance")
  # No automatic covariance ridge or generalized inverse in IRF construction.
  L<-t(chol(S))
  k<-if(shock_host=="GL") 2L else 7L
  stopifnot(all(dim(S)==c(k,k)),active_shock %in% c(1L,2L))
  blocks<-lapply(seq_along(sig),function(i) {
    if(i==ui) L[,active_shock] else rep(0,nrow(sig[[i]]))
  })
  u<-unlist(blocks,use.names=FALSE);stopifnot(length(u)==nrow(x))
  raw<-as.numeric(solve(G,u))
  j<-match(active_risk_names[active_shock],rownames(x))
  if(is.na(j)||!is.finite(raw[j])||abs(raw[j])<1e-12) stop("Cannot normalize shock impact")
  impact<-raw*(log1p(shock_pct/100)/raw[j])
  ir<-matrix(0,nrow(x),horizon+1,dimnames=list(rownames(x),as.character(0:horizon)))
  ir[,1]<-impact
  for(h in seq_len(horizon)) for(lag in seq_len(min(length(F),h)))
    ir[,h+1]<-ir[,h+1]+as.vector(F[[lag]]%*%ir[,h-lag+1])
  stopifnot(all(is.finite(ir)))
  ir
}
