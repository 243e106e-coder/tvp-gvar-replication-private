#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
dirs <- file.path("results/independent_shocks",c("gpr","vix","epu","oil"))
paths <- file.path(dirs,"irf_stable_only.csv")
if(!all(file.exists(paths))) stop("Missing independent-model IRFs: ",paste(paths[!file.exists(paths)],collapse=","))
d <- do.call(rbind,lapply(paths,read.csv))
d$shock <- factor(d$shock,levels=c("GPR","VIX","EPU","OIL"))
d <- d[grepl("_(y|dp|r|de|deq)$",d$variable),]
dir.create("results/independent_shocks/plots",recursive=TRUE,showWarnings=FALSE)
for(v in unique(sub("^.*_","",d$variable))) {
 z<-d[endsWith(d$variable,paste0("_",v)),]
 g<-ggplot(z,aes(horizon,median,color=shock,fill=shock))+
   geom_hline(yintercept=0,linetype=2,color="grey55")+
   geom_ribbon(aes(ymin=low90,ymax=high90),alpha=.10,color=NA)+
   geom_line(linewidth=.6)+
   facet_grid(date~variable,scales="free_y")+
   labs(title=paste("Independent global shocks:",v),
        subtitle="Each shock: +10% in its own log index; stable posterior draws only",
        x="Quarter after shock",y="Log-level response",color="Shock",fill="Shock")+
   theme_bw(base_size=9)+theme(legend.position="bottom")
 ggsave(file.path("results/independent_shocks/plots",paste0("compare_",v,".pdf")),
        g,width=20,height=16,limitsize=FALSE)
}
