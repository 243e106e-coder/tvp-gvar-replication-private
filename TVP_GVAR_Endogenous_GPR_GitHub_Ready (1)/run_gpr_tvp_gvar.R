rm(list = ls())
set.seed(20260714)

suppressPackageStartupMessages({
  library(compiler); library(snowfall); library(Matrix); library(mvtnorm)
  library(threshtvp); library(ggplot2)
})
source("R/BVAR_ttvp.r")
source("R/Datahandling.r")
source("R/stacking_irf_efficient.r")
source("R/auxilliary_functions_tvp.r")
source("R/gpr_girf.R")
source("R/prepare_data.R")

dir.create("results", showWarnings = FALSE)
Data.setup <- prepare_gvar_data()
xglobal <- Data.setup$bigx
gW <- Data.setup$gW
Daten <- Data.setup$new.data
cN <- Data.setup$countries

CPU <- max(1L, min(4L, parallel::detectCores() - 1L))
saves <- as.integer(Sys.getenv("TVPGVAR_SAVES", "50"))
burns <- as.integer(Sys.getenv("TVPGVAR_BURNS", "50"))
thin <- as.numeric(Sys.getenv("TVPGVAR_THIN", "0.1"))
nhor <- as.integer(Sys.getenv("TVPGVAR_HORIZON", "12"))
ext.inst <- FALSE
shrink.parm <- list(B_1 = 3, B_2 = 0.03, kappa0 = -0.1/20)
BVAR <- cmpfun(BVAR)

sfInit(parallel = TRUE, cpus = CPU)
sfExport(list = list("mlag","BVAR","datahandling","xglobal","gW","Daten",
                     "cN","bvartvpm","saves","burns","thin","ext.inst",
                     "shrink.parm"))
predDens <- sfLapply(seq_along(cN), function(i) {
  BVAR(i, gW = gW, bigx = xglobal, Daten = Daten, cN = cN,
       nsave = saves, nburn = burns, thin_chain = thin,
       ext.inst = ext.inst, parms = shrink.parm)
})
sfStop()
save(predDens, Data.setup, file = "results/predDens_gpr_endogenous.rda")

Sigma.posterior <- A.list <- globalG <- vector("list", length(predDens))
for (i in seq_along(predDens)) {
  globalG[[i]] <- predDens[[i]]$W
  Sigma.posterior[[i]] <- predDens[[i]]$SIGMApost
  A.list[[i]] <- predDens[[i]]$ALPHA
}

thin.fac <- round(thin * saves)
Tirf <- nrow(xglobal) - 1L
IRF_post <- array(NA_real_, c(Tirf, ncol(xglobal), nhor + 1L, thin.fac),
                  dimnames = list(Data.setup$quarters[-1], colnames(xglobal),
                                  0:nhor, NULL))
for (irep in seq_len(thin.fac)) {
  A.i <- rapply(A.list, classes = "array", how = "list", f = function(z) z[,,,irep])
  S.i <- rapply(Sigma.posterior, classes = "array", how = "list", f = function(z) z[,,,irep])
  for (tt in seq_len(Tirf)) {
    IRF_post[tt,,,irep] <- get.gpr.irfa.t(tt, A.i, S.i, t(xglobal), globalG,
                                          horz = nhor)$IRF_post
  }
}
save(IRF_post, file = "results/irf_gpr_endogenous.rda")

med <- apply(IRF_post, c(1,2,3), median, na.rm = TRUE)
lo  <- apply(IRF_post, c(1,2,3), quantile, probs = 0.16, na.rm = TRUE)
hi  <- apply(IRF_post, c(1,2,3), quantile, probs = 0.84, na.rm = TRUE)

plot_dates <- c("2007Q3","2008Q4","2020Q2","2022Q1")
plot_dates <- plot_dates[plot_dates %in% dimnames(IRF_post)[[1]]]
vars <- setdiff(colnames(xglobal), "US_gpr")
dd <- do.call(rbind, lapply(plot_dates, function(d) {
  ti <- match(d, dimnames(IRF_post)[[1]])
  do.call(rbind, lapply(vars, function(v) data.frame(
    date=d, variable=v, horizon=0:nhor,
    median=med[ti,v,], low=lo[ti,v,], high=hi[ti,v,]
  )))
}))
write.csv(dd, "results/gpr_irf_summary.csv", row.names = FALSE)

for (vv in c("y","dp","de","r","deq")) {
  z <- dd[grepl(paste0("_",vv,"$"), dd$variable), ]
  p <- ggplot(z, aes(horizon, median, color=date, fill=date)) +
    geom_hline(yintercept=0, color="grey55", linewidth=.35) +
    geom_ribbon(aes(ymin=low, ymax=high), alpha=.10, color=NA) +
    geom_line(linewidth=.7) + facet_wrap(~variable, scales="free_y", ncol=2) +
    theme_minimal(base_size=10) +
    labs(title=paste0("Responses to a one-standard-deviation global GPR shock: ",vv),
         x="Quarters after shock", y="Response", color="Shock date", fill="Shock date")
  ggsave(paste0("results/GPR_IRF_",vv,".png"), p, width=10, height=13, dpi=220)
}

