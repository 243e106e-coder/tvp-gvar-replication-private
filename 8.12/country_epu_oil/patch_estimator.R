source("8.12/patch_BVAR_dominant_unit.R")
s<-paste(readLines("R/BVAR_ttvp_dominant.r"),collapse="\n")
replace_once<-function(old,new) {
  if(length(gregexpr(old,s,fixed=TRUE)[[1]])!=1L||!grepl(old,s,fixed=TRUE))
    stop("Estimator patch does not match source: ",old)
  s<<-sub(old,new,s,fixed=TRUE)
}
# Keep the one-equation OO block and p=1 lag matrices two-dimensional.
replace_once('End <- xglobal[,substr(colnames(xglobal),1,2)==Names[[nr]]]',
             'End <- xglobal[,substr(colnames(xglobal),1,2)==Names[[nr]],drop=FALSE]')
replace_once('X1 <- X1[(p+1):nrow(X1),]', 'X1 <- X1[(p+1):nrow(X1),,drop=FALSE]')
replace_once('Y <- Yraw[(p+1):nrow(Yraw),]', 'Y <- Yraw[(p+1):nrow(Yraw),,drop=FALSE]')
replace_once('A0 <- A0_store[nn,,,ii]', 'A0 <- matrix(A0_store[nn,,,ii],M,M)')
# diag(x) interprets a length-one vector as a matrix SIZE. Always provide M
# so the scalar oil variance remains a 1x1 covariance, not an identity of an
# accidental size determined by exp(log-variance).
replace_once('diag(exp(H_store[nn,,ii]))', 'diag(as.numeric(exp(H_store[nn,,ii])),nrow=M,ncol=M)')
# Reject singular A0 rather than silently applying the inherited generalized inverse.
replace_once('ans <- MASS::ginv(A)', 'stop("Singular A0; no generalized inverse in this specification: ",context)')
writeLines(s,"R/BVAR_country_oil.r")
