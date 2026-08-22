#!/usr/bin/env Rscript

rm(list = ls())
seed <- as.integer(Sys.getenv("TVPGVAR_SEED", "20260816"))
RNGkind("L'Ecuyer-CMRG"); set.seed(seed)

suppressPackageStartupMessages({
  library(compiler)
  library(snowfall)
  library(Matrix)
  library(mvtnorm)
  library(threshtvp)
  library(ggplot2)
})

source("R/BVAR_ttvp_dominant.r")
source("R/Datahandling.r")
source("R/auxilliary_functions_tvp.r")
source("R/prepare_data_structural_dominant_gpr_vix.R")
source("R/dominant_gpr_vix_irf.R")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

Data.setup <- prepare_gvar_data()
xglobal <- Data.setup$bigx
gW <- Data.setup$gW
Daten <- Data.setup$new.data
cN <- Data.setup$countries
country_units <- Data.setup$country_units

# -----------------------------------------------------------------------------
# Experiment parameters
# -----------------------------------------------------------------------------
CPU <- max(1L, min(4L, parallel::detectCores() - 1L))
saves <- as.integer(Sys.getenv("TVPGVAR_SAVES", "500"))
burns <- as.integer(Sys.getenv("TVPGVAR_BURNS", "500"))
thin <- as.numeric(Sys.getenv("TVPGVAR_THIN", "0.5"))
nhor <- as.integer(Sys.getenv("TVPGVAR_HORIZON", "12"))
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10"))
lag_order <- as.integer(Sys.getenv("TVPGVAR_P", "2"))
near_unit_threshold <- as.numeric(Sys.getenv("TVPGVAR_NEAR_UNIT_THRESHOLD", "0.98"))
ext.inst <- FALSE
shrink.parm <- list(
  B_1 = as.numeric(Sys.getenv("TVPGVAR_B1", "8")),
  B_2 = as.numeric(Sys.getenv("TVPGVAR_B2", "0.01")),
  kappa0 = as.numeric(Sys.getenv("TVPGVAR_KAPPA0", "-0.005"))
)

if (!lag_order %in% c(1L, 2L)) stop("TVPGVAR_P must be 1 or 2.")
if (saves < 10L || burns < 10L) stop("Too few MCMC draws.")
if (!is.finite(thin) || thin <= 0 || thin > 1) stop("Bad thin fraction.")
if (!is.finite(shock_pct) || shock_pct <= 0) stop("Bad GPR shock percent.")

# -----------------------------------------------------------------------------
# Hard architecture audit
# -----------------------------------------------------------------------------
expected_us <- paste0("US_", c("y","dp","r","de","deq"))
us_cols <- colnames(xglobal)[startsWith(colnames(xglobal), "US_")]
if (!identical(us_cols, expected_us)) stop("US block is not macro-only.")
if (!identical(tail(colnames(xglobal), 2L), c("GL_gpr", "GL_vix"))) {
  stop("Dominant ordering must be GL_gpr -> GL_vix.")
}

mapping_audit <- do.call(rbind, lapply(seq_along(cN), function(i) {
  cc <- cN[[i]]; Wi <- as.matrix(gW[[i]])
  own_k <- if (cc == "GL") 2L else 5L
  data.frame(
    unit = cc,
    role = if (cc == "GL") "dominant" else "country",
    own_variables = own_k,
    mapping_rows = nrow(Wi),
    foreign_global_rows = nrow(Wi) - own_k,
    loads_on_GL_gpr = FALSE,
    loads_on_GL_vix = FALSE,
    stringsAsFactors = FALSE
  )
}))
# R has no inline-if syntax; overwrite the two logical columns safely below.
for (ii in seq_len(nrow(mapping_audit))) {
  cc <- mapping_audit$unit[ii]; Wi <- as.matrix(gW[[cc]])
  own_k <- if (cc == "GL") 2L else 5L
  if (nrow(Wi) > own_k) {
    fr <- Wi[(own_k + 1L):nrow(Wi), , drop = FALSE]
    mapping_audit$loads_on_GL_gpr[ii] <- any(abs(fr[, match("GL_gpr", colnames(xglobal))]) > 1e-14)
    mapping_audit$loads_on_GL_vix[ii] <- any(abs(fr[, match("GL_vix", colnames(xglobal))]) > 1e-14)
  } else {
    mapping_audit$loads_on_GL_gpr[ii] <- FALSE
    mapping_audit$loads_on_GL_vix[ii] <- FALSE
  }
}

if (mapping_audit$foreign_global_rows[mapping_audit$unit == "GL"] != 0L) {
  stop("Dominant unit must have no country/foreign regressors.")
}
if (any(!mapping_audit$loads_on_GL_gpr[mapping_audit$unit != "GL"]) ||
    any(!mapping_audit$loads_on_GL_vix[mapping_audit$unit != "GL"])) {
  stop("Every country must receive both dominant risk variables.")
}
write.csv(mapping_audit, "results/dominant_unit_mapping_audit.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Fast design preflight
# -----------------------------------------------------------------------------
design_preflight <- do.call(rbind, lapply(seq_along(cN), function(i) {
  cc <- cN[[i]]
  End <- as.matrix(xglobal[, substr(colnames(xglobal), 1, 2) == cc, drop = FALSE])
  Wi <- as.matrix(gW[[i]])
  Wex <- tvpgvar_extract_wex(Wi, xglobal, ncol(End))
  Ylag <- mlag(End, lag_order)
  Wexlag <- tvpgvar_wex_lag(Wex, 1L)
  X <- cbind(Wex, Wexlag, Ylag)
  X <- X[(lag_order + 1L):nrow(X), , drop = FALSE]
  X <- cbind(constant = 1, X)
  Y <- End[(lag_order + 1L):nrow(End), , drop = FALSE]
  storage.mode(X) <- "double"; storage.mode(Y) <- "double"
  if (nrow(X) != nrow(Y) || any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("Preflight failed for ", cc)
  }
  data.frame(
    unit = cc,
    role = if (cc == "GL") "dominant" else "country",
    p = lag_order,
    q = if (ncol(Wex) > 0L) 1L else 0L,
    Wex = ncol(Wex),
    Wexlag1 = ncol(Wexlag),
    equations = ncol(Y),
    regressors = ncol(X),
    observations = nrow(X),
    design_rank = qr(X)$rank,
    full_column_rank = qr(X)$rank == ncol(X),
    stringsAsFactors = FALSE
  )
}))
write.csv(design_preflight, "results/dominant_design_preflight.csv", row.names = FALSE)
cat("Dominant-unit design preflight passed for all ", length(cN), " units.\n", sep = "")

# -----------------------------------------------------------------------------
# MCMC estimation, unit by unit
# -----------------------------------------------------------------------------
rng_streams <- vector("list", length(cN)); rng_streams[[1L]] <- .Random.seed
if (length(cN) > 1L) for (ii in 2:length(cN)) {
  rng_streams[[ii]] <- parallel::nextRNGStream(rng_streams[[ii - 1L]])
}

BVAR <- cmpfun(BVAR)
sfInit(parallel = TRUE, cpus = CPU)
on.exit(try(sfStop(), silent = TRUE), add = TRUE)
sfExport(list = list(
  "mlag","BVAR","datahandling","xglobal","gW","Daten","cN","bvartvpm",
  "saves","burns","thin","ext.inst","shrink.parm","rng_streams",
  "tvpgvar_extract_wex","tvpgvar_wex_lag","tvpgvar_safe_inverse"
))

predDens <- sfLapply(seq_along(cN), function(i) {
  assign(".Random.seed", rng_streams[[i]], envir = .GlobalEnv)
  tryCatch(
    BVAR(i, gW = gW, bigx = xglobal, Daten = Daten, cN = cN,
         nsave = saves, nburn = burns, thin_chain = thin,
         ext.inst = ext.inst, parms = shrink.parm),
    error = function(e) stop("Unit ", cN[[i]], " failed: ", conditionMessage(e), call. = FALSE)
  )
})
sfStop()

save(predDens, Data.setup, file = "results/predDens_dominant_gpr_vix.rda")

# -----------------------------------------------------------------------------
# Regressor structure audit
# -----------------------------------------------------------------------------
reg_diag <- do.call(rbind, lapply(seq_along(predDens), function(i) {
  rn <- dimnames(predDens[[i]]$ALPHA)[[2]]
  data.frame(
    unit = cN[[i]],
    role = if (cN[[i]] == "GL") "dominant" else "country",
    Wex = sum(rn == "Wex"),
    Wexlag1 = sum(rn == "Wexlag1"),
    Wexlag2 = sum(rn == "Wexlag2"),
    Ylag1 = sum(rn == "Ylag1"),
    Ylag2 = sum(rn == "Ylag2"),
    stringsAsFactors = FALSE
  )
}))
expected_wex <- ifelse(reg_diag$unit == "GL", 0L, 7L)
expected_y <- ifelse(reg_diag$unit == "GL", 2L, 5L)
if (any(reg_diag$Wex != expected_wex) || any(reg_diag$Wexlag1 != expected_wex) ||
    any(reg_diag$Wexlag2 != 0L) || any(reg_diag$Ylag1 != expected_y) ||
    any(reg_diag$Ylag2 != if (lag_order == 2L) expected_y else 0L)) {
  print(reg_diag); stop("Estimated regressor structure mismatch.")
}
reg_diag$domestic_lag_order <- lag_order
reg_diag$foreign_lag_order <- ifelse(reg_diag$unit == "GL", 0L, 1L)
write.csv(reg_diag, "results/dominant_lag_spec_diagnostic.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Residual diagnostics
# -----------------------------------------------------------------------------
serial_diag <- do.call(rbind, lapply(seq_along(predDens), function(i) {
  cc <- cN[[i]]; res <- as.matrix(predDens[[i]]$cc.res)
  own <- colnames(xglobal)[substr(colnames(xglobal), 1, 2) == cc]
  do.call(rbind, lapply(seq_along(own), function(j) {
    z <- as.numeric(res[, j]); z <- z[is.finite(z)]
    lag_use <- min(4L, max(1L, floor(length(z)/5L)))
    bt <- tryCatch(Box.test(z, lag=lag_use, type="Ljung-Box", fitdf=0), error=function(e) NULL)
    data.frame(unit=cc, variable=own[j], n=length(z), test_lag=lag_use,
               p_value=if (is.null(bt)) NA_real_ else bt$p.value,
               reject_5pct=if (is.null(bt)) NA else bt$p.value < 0.05,
               stringsAsFactors=FALSE)
  }))
}))
write.csv(serial_diag, "results/residual_serial_correlation_diagnostic.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# IRF + stability posterior reconstruction
# -----------------------------------------------------------------------------
A.list <- lapply(predDens, `[[`, "ALPHA")
S.list <- lapply(predDens, `[[`, "SIGMApost")
globalG <- lapply(predDens, `[[`, "W")

n_irf_draws <- dim(A.list[[1]])[4]
if (any(vapply(A.list, function(z) dim(z)[4], integer(1)) != n_irf_draws)) {
  stop("Inconsistent retained draw counts.")
}
Tirf <- nrow(xglobal) - lag_order
irf_dates <- Data.setup$quarters[-seq_len(lag_order)]
if (any(vapply(A.list, function(z) dim(z)[1], integer(1)) != Tirf)) {
  stop("Coefficient paths are not aligned across units.")
}

IRF_post <- array(NA_real_, c(Tirf, ncol(xglobal), nhor+1L, n_irf_draws),
                  dimnames=list(irf_dates, colnames(xglobal), 0:nhor, NULL))
stability_rho <- matrix(NA_real_, Tirf, n_irf_draws, dimnames=list(irf_dates, NULL))
G_condition <- G_min_singular <- G_max_singular <- max_abs_transition <- stability_rho
local_rho_array <- array(NA_real_, c(Tirf, length(cN), n_irf_draws),
                         dimnames=list(irf_dates, cN, NULL))

for (dd in seq_len(n_irf_draws)) {
  A.i <- lapply(A.list, function(z) z[, , , dd])
  S.i <- lapply(S.list, function(z) z[, , , dd])
  for (tt in seq_len(Tirf)) {
    ans <- get_dominant_gpr_irf_t(
      tt=tt, draw_i=A.i, Sig_draw_i=S.i, x=t(xglobal), globalG=globalG,
      units=cN, horizon=nhor, shock_pct=shock_pct
    )
    IRF_post[tt, , , dd] <- ans$IRF_post
    stability_rho[tt, dd] <- ans$max_eigen_modulus
    G_condition[tt, dd] <- ans$G_condition_number
    G_min_singular[tt, dd] <- ans$G_min_singular_value
    G_max_singular[tt, dd] <- ans$G_max_singular_value
    max_abs_transition[tt, dd] <- ans$max_abs_transition
    local_rho_array[tt, , dd] <- ans$local_rho[cN]
  }
  cat("IRF draw ", dd, "/", n_irf_draws, " complete\n", sep = "")
}

if (any(!is.finite(stability_rho)) || any(!is.finite(local_rho_array))) {
  stop("Non-finite stability diagnostics.")
}
stable_mask <- stability_rho < 1
save(IRF_post, stability_rho, stable_mask, G_condition, G_min_singular,
     G_max_singular, max_abs_transition, local_rho_array,
     file="results/irf_dominant_gpr_vix.rda")

# -----------------------------------------------------------------------------
# Stability summaries
# -----------------------------------------------------------------------------
stability_summary <- do.call(rbind, lapply(seq_len(Tirf), function(tt) {
  r <- stability_rho[tt, ]
  data.frame(date=irf_dates[tt], draws=length(r), stable_draws=sum(r<1),
             stable_share=mean(r<1), median_rho=median(r),
             p95_rho=as.numeric(quantile(r,.95)), max_rho=max(r),
             stringsAsFactors=FALSE)
}))
write.csv(stability_summary, "results/stability_summary_by_date.csv", row.names=FALSE)

local_summary <- do.call(rbind, lapply(seq_along(cN), function(i) {
  r <- as.numeric(local_rho_array[,i,])
  data.frame(unit=cN[i], role=if(cN[i]=="GL") "dominant" else "country",
             observations=length(r), stable_share=mean(r<1), median_rho=median(r),
             p95_rho=as.numeric(quantile(r,.95)), max_rho=max(r), stringsAsFactors=FALSE)
}))
write.csv(local_summary, "results/local_stability_summary_by_country.csv", row.names=FALSE)

G_warn <- as.numeric(Sys.getenv("TVPGVAR_G_CONDITION_WARN", "10000"))
global_source <- do.call(rbind, lapply(seq_len(n_irf_draws), function(dd) {
  do.call(rbind, lapply(seq_len(Tirf), function(tt) {
    lr <- local_rho_array[tt,,dd]
    all_local <- all(lr < 1)
    gill <- !is.finite(G_condition[tt,dd]) || G_condition[tt,dd] > G_warn
    gs <- stability_rho[tt,dd] < 1
    src <- if (gs) "global_stable" else if (!all_local && gill) {
      "local_dynamics_and_G_conditioning"
    } else if (!all_local) "local_dynamics" else if (gill) "near_singular_G" else "cross_country_feedback"
    data.frame(date=irf_dates[tt], draw=dd, global_rho=stability_rho[tt,dd],
               global_stable=gs, all_local_stable=all_local,
               worst_local_country=names(which.max(lr)), worst_local_rho=max(lr),
               G_condition_number=G_condition[tt,dd], G_ill_conditioned=gill,
               likely_source=src, stringsAsFactors=FALSE)
  }))
}))
write.csv(global_source, "results/global_stability_source_diagnostic.csv", row.names=FALSE)

Gdiag <- do.call(rbind, lapply(seq_len(n_irf_draws), function(dd) data.frame(
  date=irf_dates, draw=dd, condition_number=G_condition[,dd],
  min_singular_value=G_min_singular[,dd], max_singular_value=G_max_singular[,dd],
  stringsAsFactors=FALSE)))
write.csv(Gdiag, "results/G_condition_diagnostic.csv", row.names=FALSE)

# -----------------------------------------------------------------------------
# IRF summaries
# -----------------------------------------------------------------------------
selected_dates <- intersect(c("2003Q1","2008Q3","2014Q3","2020Q1","2022Q1","2023Q4"), irf_dates)

summarize_irf <- function(stable_only=FALSE) {
  rows <- list()
  for (date in selected_dates) {
    tt <- match(date, irf_dates)
    draws <- seq_len(n_irf_draws)
    if (stable_only) draws <- draws[stable_mask[tt,draws]]
    if (!length(draws)) next
    for (v in colnames(xglobal)) {
      vi <- match(v, colnames(xglobal))
      for (h in 0:nhor) {
        z <- IRF_post[tt,vi,h+1L,draws]
        rows[[length(rows)+1L]] <- data.frame(
          date=date, variable=v, horizon=h, median=median(z),
          low68=as.numeric(quantile(z,.16)), high68=as.numeric(quantile(z,.84)),
          low90=as.numeric(quantile(z,.05)), high90=as.numeric(quantile(z,.95)),
          stable_only=stable_only, draws_used=length(draws), stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

irf_all <- summarize_irf(FALSE)
irf_stable <- summarize_irf(TRUE)
write.csv(irf_all, "results/gpr_structural_irf_summary.csv", row.names=FALSE)
write.csv(irf_stable, "results/gpr_structural_irf_summary_stable_only.csv", row.names=FALSE)

# GDP sign diagnostics at selected horizons.
gdp_vars <- paste0(country_units, "_y")
sel_h <- intersect(c(0L,1L,4L,8L,12L), 0:nhor)
make_gdp_sign <- function(stable_only=FALSE) {
  rows <- list()
  for (date in selected_dates) {
    tt <- match(date, irf_dates)
    draws <- seq_len(n_irf_draws)
    if (stable_only) draws <- draws[stable_mask[tt,draws]]
    if (!length(draws)) next
    for (v in gdp_vars) for (h in sel_h) {
      z <- IRF_post[tt,match(v,colnames(xglobal)),h+1L,draws]
      qs <- quantile(z,c(.05,.16,.5,.84,.95),names=FALSE)
      rows[[length(rows)+1L]] <- data.frame(
        date=date, variable=v, horizon=h, median=qs[3], low68=qs[2], high68=qs[4],
        low90=qs[1], high90=qs[5], positive_median=qs[3]>0,
        sig68=(qs[2]>0 || qs[4]<0), sig90=(qs[1]>0 || qs[5]<0),
        stable_only=stable_only, draws_used=length(draws), stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}
gdp_all <- make_gdp_sign(FALSE); gdp_stable <- make_gdp_sign(TRUE)
write.csv(gdp_all, "results/gdp_sign_diagnostic.csv", row.names=FALSE)
write.csv(gdp_stable, "results/gdp_sign_diagnostic_stable_only.csv", row.names=FALSE)

# Cumulative GDP response (appropriate for interpreting cumulative response of a
# growth/difference variable; retain raw IRF summary for primary interpretation).
make_cum_gdp <- function(stable_only=FALSE) {
  rows <- list()
  for (date in selected_dates) {
    tt <- match(date, irf_dates)
    draws <- seq_len(n_irf_draws)
    if (stable_only) draws <- draws[stable_mask[tt,draws]]
    if (!length(draws)) next
    for (v in gdp_vars) {
      vi <- match(v,colnames(xglobal))
      Z <- IRF_post[tt,vi,,draws,drop=FALSE]
      Z <- matrix(Z,nrow=nhor+1L,ncol=length(draws))
      CZ <- apply(Z,2,cumsum)
      for (h in 0:nhor) {
        z <- CZ[h+1L,]
        qs <- quantile(z,c(.05,.16,.5,.84,.95),names=FALSE)
        rows[[length(rows)+1L]] <- data.frame(
          date=date, variable=v, horizon=h, median=qs[3], low68=qs[2], high68=qs[4],
          low90=qs[1], high90=qs[5], stable_only=stable_only,
          draws_used=length(draws), stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind,rows)
}
write.csv(make_cum_gdp(FALSE), "results/gpr_structural_cumulative_gdp.csv", row.names=FALSE)
write.csv(make_cum_gdp(TRUE), "results/gpr_structural_cumulative_gdp_stable_only.csv", row.names=FALSE)

# Shock normalization validation.
target <- log1p(shock_pct/100)
shock_validation <- do.call(rbind,lapply(seq_len(Tirf),function(tt){
  z <- IRF_post[tt,match("GL_gpr",colnames(xglobal)),1,]
  data.frame(date=irf_dates[tt],target_log_jump=target,median_h0=median(z),
             min_h0=min(z),max_h0=max(z),max_abs_error=max(abs(z-target)),
             stringsAsFactors=FALSE)
}))
write.csv(shock_validation,"results/gpr_shock_validation.csv",row.names=FALSE)

# Selected-date stability convenience file.
write.csv(stability_summary[stability_summary$date %in% selected_dates,,drop=FALSE],
          "results/stability_summary_selected_dates.csv",row.names=FALSE)

# Architecture / provenance notes.
writeLines(c(
  "Separate dominant-unit GPR + VIX TVP-GVAR",
  paste0("Profile: ", Sys.getenv("TVPGVAR_PROFILE", "P2_strong_dominant_GPR_VIX")),
  paste0("Sample: ", Data.setup$quarters[1], " - ", tail(Data.setup$quarters,1)),
  paste0("Observations before lags: ", nrow(xglobal)),
  paste0("Country domestic lag p: ", lag_order),
  "Foreign macro/global lag q: 1 for country blocks",
  "Dominant unit: z_t = [GL_gpr, GL_vix]'",
  "Dominant unit receives no country macro variables.",
  "Country Wex = 5 trade-weighted foreign macro variables + current GPR + current VIX.",
  "Country Wexlag1 = the same 7 variables at t-1.",
  "GPR shock identification: recursive Cholesky inside dominant 2x2 block, GPR -> VIX.",
  "No global 72-variable Cholesky ordering is imposed.",
  "Paper-inspired architecture, not exact replication: the reference paper uses [GPR,HWWI] and BIC-selected one lag; this experiment uses [GPR,VIX] and the selected TVPGVAR_P setting."
), "results/dominant_model_architecture_summary.txt")

cat("\nDominant-unit TVP-GVAR run complete.\n")
cat("Global stable share: ", mean(stable_mask), "\n", sep="")
