#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
root<-"results/epu_wti"
files<-list.files(root,pattern="^irf_summary.csv$",recursive=TRUE,full.names=TRUE)
expected<-c("US_GLOBAL","US_US","CN_GLOBAL","CN_US")
status<-data.frame(model=expected,has_irf=vapply(expected,function(m)
                  file.exists(file.path(root,m,"irf_summary.csv")),logical(1)))
write.csv(status,file.path(root,"model_completion.csv"),row.names=FALSE)
if(!length(files)) stop("No model IRFs; inspect job logs and posterior artifacts")
all<-do.call(rbind,lapply(files,read.csv))
all$variable<-sub("^(GL|US)_(epu|wti)$","COMMON_\\2",all$variable)
write.csv(all,file.path(root,"all_irf_summaries.csv"),row.names=FALSE)
pf<-list.files(root,pattern="^selected_peak_dates.csv$",recursive=TRUE,full.names=TRUE)
if(length(pf)) {
  peaks<-unique(do.call(rbind,lapply(pf,read.csv)))
  write.csv(peaks,file.path(root,"selected_epu_peak_dates.csv"),row.names=FALSE)
  q<-read.csv("8.12/epu_wti_quarterly.csv")
  pdf(file.path(root,"EPU_series_with_selected_peaks.pdf"),width=11,height=5)
  for(ec in c("US","CN")) {
    z<-data.frame(Quarter=q$Quarter,EPU=q[[paste0(ec,"_EPU")]])
    z<-subset(z,Quarter>="2000Q3"&Quarter<="2025Q3")
    z$time<-as.integer(substr(z$Quarter,1,4))+(as.integer(substr(z$Quarter,6,6))-1)/4
    pp<-subset(peaks,epu_country==ec)
    pp$time<-as.integer(substr(pp$Quarter,1,4))+(as.integer(substr(pp$Quarter,6,6))-1)/4
    print(ggplot(z,aes(time,EPU))+geom_line()+
       geom_point(data=pp,color="red",size=2)+
       geom_text(data=pp,aes(label=Quarter),vjust=-.6,size=3)+
       theme_bw()+labs(title=paste(ec,"EPU: selected period peaks"),
          subtitle=if(ec=="CN") "Davis-Liu-Sheng mainland overall EPU; quarterly mean" else "US news-based EPU; quarterly mean",
          x="Year",y="EPU index"))
  }
  dev.off()
}
sf<-list.files(root,pattern="^stability.csv$",recursive=TRUE,full.names=TRUE)
if(length(sf)) write.csv(do.call(rbind,lapply(sf,read.csv)),file.path(root,"all_stability.csv"),row.names=FALSE)
global<-subset(all,placement=="GLOBAL")
us<-subset(all,placement=="US")
keys<-c("epu_country","date","shock","variable","response_kind","filter","horizon")
paired<-merge(global,us,by=keys,suffixes=c("_global","_us"))
if(nrow(paired)) {
  paired$median_us_minus_global<-paired$median_us-paired$median_global
  paired$sign_flip<-sign(paired$median_us)!=sign(paired$median_global)
  paired$significant_global<-paired$low90_global>0|paired$high90_global<0
  paired$significant_us<-paired$low90_us>0|paired$high90_us<0
  write.csv(paired,file.path(root,"placement_comparison.csv"),row.names=FALSE)
}
# PDFs: one page per date/shock/macro variable, 14 countries, both placements.
plot_data<-subset(all,filter=="stable_only"&grepl("^[A-Z]{2}_(y|dp|r|de|deq)$",variable))
plot_data$country<-substr(plot_data$variable,1,2)
plot_data$macro<-sub("^[A-Z]{2}_","",plot_data$variable)
labels<-c(y="log GDP response",dp="CPI inflation response",r="interest-rate change response",
          de="REER change response",deq="equity return response")
for(ec in c("US","CN")) {
  d<-subset(plot_data,epu_country==ec)
  if(!nrow(d)) next
  pdf(file.path(root,paste0(ec,"_EPU_WTI_placement_comparison.pdf")),width=12,height=9)
  for(kind in unique(d$response_kind)) for(date_i in unique(d$date))
    for(sh in c("EPU","WTI")) for(v in unique(d$macro)) {
      z<-subset(d,date==date_i&shock==sh&macro==v&response_kind==kind)
      if(!nrow(z)) next
      title<-if(kind=="rate_level_cumulative")"Cumulative rate-level response" else labels[[v]]
      print(ggplot(z,aes(horizon,median,color=placement,fill=placement))+
        geom_hline(yintercept=0,color="grey60",linewidth=.3)+
        geom_ribbon(aes(ymin=low90,ymax=high90),alpha=.12,color=NA)+
        geom_line(linewidth=.55)+facet_wrap(~country,ncol=4,scales="free_y")+
        theme_bw(base_size=10)+labs(title=paste(ec,"EPU + WTI |",date_i,"| +10%",sh,"shock"),
          subtitle=paste(title,"; stable posterior draws only; shaded 90% intervals"),
          x="Quarters after shock",y=title))
    }
  dev.off()
}
writeLines(c("Compare placements only within a given EPU series and common sample.",
 "IRF dates are country-specific EPU period maxima; GLOBAL and US placements share dates within each EPU country.",
 "China EPU is Davis-Liu-Sheng MAINLAND overall EPU (not SCMP, not TPU).",
 "Both EPU and WTI shocks are normalized to +10% own-variable impact.",
 "Placement differences are descriptive differences between posterior medians, not a posterior credible interval for their difference.",
 "Stable-only posterior conditioning differs by specification; inspect all_stability.csv and valid_all summaries.",
 "500 save/500 burn is an exploratory run; assess MCMC convergence before publication.",
 "The US-hosted China index remains CHINESE EPU. This is an architecture sensitivity test.",
 "WTI innovations mix supply/demand/news mechanisms; do not call them identified oil-supply shocks.",
 "If a model has no stable draws for a date, its stable-only page/row is absent; see completion and stability diagnostics."),
 file.path(root,"COMPARISON_NOTES.txt"))
if(!all(status$has_irf)) stop("Incomplete four-model comparison; partial outputs saved")
cat("PASS four-model comparison\n")
