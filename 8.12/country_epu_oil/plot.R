suppressPackageStartupMessages(library(ggplot2))
root<-"results/country_epu_oil"
for(ec in c("CN","US")) {
  out<-file.path(root,paste0(ec,"_COUNTRY_OIL"));f<-file.path(out,"irf_summary.csv")
  if(!file.exists(f)) stop("Missing completed model: ",f)
  d<-read.csv(f);d<-d[grepl("^[A-Z]{2}_(y|dp|r|de|deq)$",d$variable),]
  d$country<-substr(d$variable,1,2);d$macro<-sub("^[A-Z]{2}_","",d$variable)
  peaks<-read.csv(file.path(out,"selected_peak_dates.csv"))
  # Both filters are exported. Each page is ONE shock/date/response across 14
  # country panels. Nested intervals come from the same posterior draws.
  pdf(file.path(out,"IRF_68_90.pdf"),width=12,height=9)
  for(sh in unique(d$shock)) for(date_i in unique(d$date[d$shock==sh]))
    for(filter_i in unique(d$filter)) for(kind in unique(d$response_kind)) for(v in unique(d$macro)) {
      z<-subset(d,shock==sh&date==date_i&filter==filter_i&response_kind==kind&macro==v)
      if(!nrow(z)) next
      pp<-peaks[peaks$Quarter==date_i&peaks$shock==sh,]
      # r is estimated in first differences. Only the draw-wise cumulative
      # summary represents the response of the interest-rate LEVEL.
      ylabel <- if(v == "r") {
        if(kind == "rate_level_cumulative") "Interest-rate LEVEL response (percentage points)"
        else "Interest-rate CHANGE response (percentage points)"
      } else "Response in transformed variable units"
      print(ggplot(z,aes(horizon,median))+
        geom_hline(yintercept=0,color="grey50",linewidth=.3)+
        geom_ribbon(aes(ymin=low90,ymax=high90),fill="#66B2FF",alpha=.25)+
        geom_ribbon(aes(ymin=low68,ymax=high68),fill="#0875D1",alpha=.4)+
        geom_line(color="#003E85",linewidth=.65)+facet_wrap(~country,ncol=4,scales="free_y")+
        theme_bw()+labs(title=paste(ec,"EPU + WTI |",sh,"shock |",date_i,"|",v),
          subtitle=paste("Exploratory rate-proxy candidate | Peak series:",paste(pp$series,collapse=","),
             "|",filter_i,"|",kind,"| dark 68%, light 90% credible intervals"),
          x="Quarters after +10% impact",y=ylabel))
    }
  dev.off()
  # Country counts and horizon counts are separate quantities.
  counts<-aggregate(cbind(significant68,significant90)~date+shock+macro+response_kind+filter+horizon,d,sum)
  names(counts)[names(counts)=="significant68"]<-"countries_significant68"
  names(counts)[names(counts)=="significant90"]<-"countries_significant90"
  write.csv(counts,file.path(out,"significant_country_counts_by_horizon.csv"),row.names=FALSE)
}
