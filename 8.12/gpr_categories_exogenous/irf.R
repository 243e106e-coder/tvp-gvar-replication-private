source("8.12/dominant_gpr_vix_irf.R")
# Override only the old two-variable GPR/VIX shock injection. No VIX/EPU/WTI
# is present in these separate 14 fits. No ridge, pseudo-inverse or forced
# stability. Contemporaneous national innovations are set to zero.
dominant_gpr_struct_irf <- function(G,F,sig,x,units,horizon=12,shock_var="GL_gpr",shock_pct=10) {
  K<-nrow(x);j<-match("GL_gpr",rownames(x));ui<-match("GL",units)
  stopifnot(!is.na(j),!is.na(ui),j==K,all(dim(G)==c(K,K)))
  target<-numeric(K);target[j]<-1
  if(max(abs(G[j,]-target))>1e-10) stop("GPR path is not external: contemporaneous macro feedback")
  for(f in F) if(max(abs(f[j,-j]))>1e-10) stop("GPR path is not external: lagged macro feedback")
  variance<-as.numeric(sig[[ui]])
  if(length(variance)!=1L||!is.finite(variance)||variance<=0)stop("Invalid external-path innovation variance")
  u<-numeric(K);u[j]<-sqrt(variance)
  raw<-as.numeric(solve(G,u))
  if(!is.finite(raw[j])||abs(raw[j])<1e-12)stop("Cannot normalize external GPR shock")
  ans<-matrix(0,K,horizon+1,dimnames=list(rownames(x),0:horizon))
  ans[,1]<-raw*log1p(shock_pct/100)/raw[j]
  for(h in seq_len(horizon))for(lag in seq_len(min(length(F),h))) ans[,h+1]<-ans[,h+1]+as.numeric(F[[lag]]%*%ans[,h-lag+1])
  stopifnot(all(is.finite(ans)),abs(ans[j,1]-log1p(shock_pct/100))<1e-12)
  ans
}
