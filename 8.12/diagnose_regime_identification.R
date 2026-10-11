#!/usr/bin/env Rscript
# TVP-GVAR post-estimation comparison inspired by Lhuissier & Tripier (2021).
# This is NOT an MS-SVAR, MS-DSGE, or newly identified structural-shock estimate.
# It reads existing independent-shock IRFs and makes descriptive diagnostics.
# No data imputation, no MCMC, no modification to input artifacts.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
input_root <- Sys.getenv('IRF_INPUT_ROOT', 'results/independent_shocks')
output_root <- Sys.getenv('IRF_DIAG_OUT', 'results/regime_identification_diagnostics')
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
shocks <- c('GPR','VIX','EPU','OIL')
countries <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
macros <- c('y','dp','r','de','deq')
files <- setNames(file.path(input_root,tolower(shocks),'irf_stable_only.csv'),shocks)
missing <- files[!file.exists(files)]
if (length(missing)) stop('Missing IRF CSV(s): ',paste(missing,collapse=', '))
read_shock <- function(s) {
  z <- read.csv(files[[s]], check.names=FALSE)
  required <- c('shock','date','variable','horizon','median','low90','high90','stable_draws')
  if(!all(required %in% names(z))) stop(s,': missing ',paste(setdiff(required,names(z)),collapse=', '))
  z$shock <- s
  z
}
irf <- do.call(rbind,lapply(shocks,read_shock))
if(anyDuplicated(irf[c('shock','date','variable','horizon')])) stop('Duplicated IRF keys')
if(any(!is.finite(irf$median)) || any(!is.finite(irf$low90)) || any(!is.finite(irf$high90))) stop('Nonfinite IRF entries')
if(any(irf$low90>irf$median | irf$median>irf$high90)) stop('Posterior quantiles out of order')
irf$country <- sub('_.*$','',irf$variable)
irf$macro <- sub('^.*_','',irf$variable)
irf <- subset(irf,country %in% countries & macro %in% macros)
if(nrow(irf)==0) stop('No country IRF records matched expected variable names')
irf$positive <- irf$median>0
irf$negative <- irf$median<0
irf$above_zero90 <- irf$low90>0
irf$below_zero90 <- irf$high90<0
irf$excludes_zero90 <- irf$above_zero90|irf$below_zero90
irf$shock <- factor(irf$shock,levels=shocks)
write.csv(irf,file.path(output_root,'country_irfs_annotated.csv'),row.names=FALSE)

# All dates are named historical snapshots, NOT inferred latent Markov regimes.
# In particular, 2003Q1 and 2023Q4 must not be automatically called tranquil.
counts <- aggregate(cbind(
  positive=as.integer(irf$positive),negative=as.integer(irf$negative),
  positive_ci90=as.integer(irf$above_zero90),negative_ci90=as.integer(irf$below_zero90)),
  by=list(shock=irf$shock,date=irf$date,macro=irf$macro,horizon=irf$horizon),FUN=sum)
counts$n_economies <- as.integer(aggregate(irf$country,
  by=list(shock=irf$shock,date=irf$date,macro=irf$macro,horizon=irf$horizon),
  FUN=function(x) length(unique(x)))$x)
write.csv(counts,file.path(output_root,'sign_and_credible_interval_counts.csv'),row.names=FALSE)

# A *necessary* sign-pattern diagnostic, NOT a sign-restricted identification:
# Lhuissier & Tripier identify uncertainty with GDP down, VIX up, credit spread up,
# plus a >=50% VIX 1-step FEVD restriction. Current separate-shock IRFs have
# neither a jointly modeled credit spread nor the required joint FEVD/rotations.
# Accordingly this report NEVER marks identification as passed.
focus <- subset(irf,date %in% c('2008Q3','2020Q1','2022Q1') & horizon %in% c(0,1,4,8))
gdp <- subset(focus,macro=='y')
gdp$gdp_negative_median <- gdp$median<0
gdp$gdp_negative_credible90 <- gdp$high90<0
write.csv(gdp[,c('shock','date','horizon','country','median','low90','high90','gdp_negative_median','gdp_negative_credible90')],
  file.path(output_root,'gdp_recessionary_sign_audit.csv'),row.names=FALSE)

# Country / date / horizon equity puzzle (GPR, VIX, EPU)
equity <- subset(focus,macro=='deq' & shock %in% c('GPR','VIX','EPU'))
write.csv(equity[,c('shock','date','horizon','country','median','low90','high90','positive','above_zero90')],
  file.path(output_root,'equity_positive_response_audit.csv'),row.names=FALSE)

# Four-shock model stability: purely descriptive diagnostic from saved stability CSV.
stfiles <- setNames(file.path(input_root,tolower(shocks),'stability.csv'),shocks)
if(all(file.exists(stfiles))) {
  st <- do.call(rbind,lapply(shocks,function(s){z<-read.csv(stfiles[[s]]);z$shock<-s;z}))
  if(all(c('date','stable_share','total','stable') %in% names(st))) {
    st$stable_percent <- st$stable_share*100
    write.csv(st,file.path(output_root,'stable_posterior_shares.csv'),row.names=FALSE)
    ps<-ggplot(st,aes(date,stable_percent,group=shock,color=shock))+geom_line()+geom_point()+
      geom_hline(yintercept=50,linetype=2,color='grey50')+theme_bw()+
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      labs(title='Stable posterior share by independent shock',x='Snapshot',y='Stable draws (%)')
    ggsave(file.path(output_root,'posterior_stability.pdf'),ps,width=10,height=5)
  }
} else {
  warning('Some stability CSVs unavailable; omitting stability chart')
}

# Compare IRFs across actual snapshots; no regime inference is made.
for (s in shocks) for (v in macros) {
  dat<-subset(irf,shock==s & macro==v)
  dat$country<-factor(dat$country,levels=countries)
  plot<-ggplot(dat,aes(horizon,median,color=date,fill=date))+
    geom_hline(yintercept=0,linetype=2,color='grey50')+
    geom_ribbon(aes(ymin=low90,ymax=high90),alpha=.05,color=NA)+
    geom_line(linewidth=.42)+facet_wrap(~country,ncol=4,scales='free_y')+
    theme_bw(base_size=9)+theme(legend.position='bottom')+
    labs(title=paste('Historical snapshot comparison:',s,'shock /',v),
         subtitle='Existing TVP-GVAR conditional IRFs; not Markov regimes',
         x='Horizon (quarters)',y='Posterior response')
  ggsave(file.path(output_root,paste0('snapshot_',tolower(s),'_',v,'.pdf')),plot,width=13,height=10)
}

writeLines(c(
 'DIAGNOSTICS INSPIRED BY LHUISSIER & TRIPIER (2021)',
 'These reports reuse previously saved TVP-GVAR independent-shock CSVs.',
 'Snapshots such as 2008Q3, 2020Q1, 2022Q1 are not estimated latent Markov regimes.',
 'The paper uses MS-SVAR in GDP log-levels, VIX and BAA-AAA credit spread.',
 'Its uncertainty shock uses impact signs (GDP negative, VIX positive, spread positive)',
 'and a VIX one-step forecast-error variance contribution >=50%.',
 'Current independent-shock outputs cannot test the joint spread sign or FEVD restriction.',
 'The GDP sign audit is descriptive and must NOT be labeled structural identification.',
 'The current baseline shock is a 10% log-index impact; the paper uses a VIX 10-point rise.',
 'No new Bayesian posterior is calculated. Bands are from the stable draws already saved.',
 'Results conditioned on stable draws may differ from the full posterior.',
 'Do not interpret a positive equity median as statistically supported when CI crosses zero.',
 'Data are neither imputed nor overwritten.'),file.path(output_root,'METHODOLOGY_LIMITATIONS.txt'))
cat('Diagnostics generated in ',output_root,'\n',sep='')
