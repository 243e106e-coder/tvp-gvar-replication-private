#!/usr/bin/env Rscript
x<-read.csv('8.12/VIXCLS_daily_fred.csv',check.names=FALSE)
dc<-intersect(c('DATE','observation_date','Date'),names(x)); vc<-intersect(c('VIXCLS','vixcls'),names(x))
if(length(dc)!=1||length(vc)!=1) stop('Unexpected VIX columns')
d<-as.Date(x[[dc]]);v<-suppressWarnings(as.numeric(x[[vc]]));ok<-!is.na(d)&is.finite(v)&v>0
if(!any(ok)) stop('No positive daily VIX values')
q<-paste0(format(d[ok],'%Y'),'Q',((as.integer(format(d[ok],'%m'))-1L)%/%3L+1L))
a<-tapply(v[ok],q,mean)
y<-data.frame(Quarter=names(a),US_vix=log(as.numeric(a)))
y<-y[order(as.integer(substr(y$Quarter,1,4))*4+as.integer(substr(y$Quarter,6,6))),]
dir.create('data/loglevels',recursive=TRUE,showWarnings=FALSE)
write.csv(y,'data/loglevels/vix_quarterly_processed.csv',row.names=FALSE)
cat('PASS VIX quarterly log means\n')
