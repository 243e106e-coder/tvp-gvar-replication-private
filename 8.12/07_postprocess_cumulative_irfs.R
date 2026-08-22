#!/usr/bin/env Rscript

# =============================================================================
# POST-PROCESS EXISTING TVP-GVAR IRFs (NO MCMC)
#
# Current model transformation:
#   y   = log real GDP level            -> keep IRF as-is
#   dp  = change in log CPI / inflation -> cumulative IRF across horizons
#   r   = interest-rate level           -> keep IRF as-is
#   de  = change in log REER            -> cumulative IRF across horizons
#   deq = change in log equity price    -> cumulative IRF across horizons
#
# IMPORTANT:
# Cumulation is performed DRAW BY DRAW first. Quantiles are computed only after
# cumulation. Do not sum posterior medians from the summary CSV.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

input_rda <- Sys.getenv(
  "TVPGVAR_IRF_RDA",
  "prior_artifact/results/irf_dominant_gpr_vix.rda"
)
out_dir <- Sys.getenv("TVPGVAR_POST_DIR", "results/cumulative_postprocess")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_rda)) stop("Missing draw-level IRF file: ", input_rda)
load(input_rda)

need <- c("IRF_post", "stable_mask")
missing_obj <- need[!vapply(need, exists, logical(1), inherits = FALSE)]
if (length(missing_obj)) stop("IRF RDA is missing: ", paste(missing_obj, collapse = ", "))

if (length(dim(IRF_post)) != 4L) stop("IRF_post must be [date, variable, horizon, draw].")
if (!all(dim(stable_mask) == c(dim(IRF_post)[1], dim(IRF_post)[4]))) {
  stop("stable_mask dimensions do not match IRF_post date/draw dimensions.")
}

irf_dates <- dimnames(IRF_post)[[1]]
var_names <- dimnames(IRF_post)[[2]]
horizons <- dimnames(IRF_post)[[3]]
if (is.null(irf_dates) || is.null(var_names)) stop("IRF_post must have date and variable dimnames.")
if (is.null(horizons)) horizons <- as.character(seq_len(dim(IRF_post)[3]) - 1L)
h_num <- suppressWarnings(as.integer(horizons))
if (anyNA(h_num)) h_num <- seq_len(dim(IRF_post)[3]) - 1L

# Selected diagnostic dates. If one is absent, silently skip it.
selected_dates <- c("2003Q1", "2008Q3", "2014Q3", "2020Q1", "2022Q1", "2023Q4")
selected_dates <- intersect(selected_dates, irf_dates)
if (!length(selected_dates)) selected_dates <- irf_dates

suffix_of <- function(v) {
  if (grepl("_deq$", v)) return("deq")
  if (grepl("_dp$",  v)) return("dp")
  if (grepl("_de$",  v)) return("de")
  if (grepl("_r$",   v)) return("r")
  if (grepl("_y$",   v)) return("y")
  return(NA_character_)
}

country_of <- function(v) sub("_(y|dp|r|de|deq)$", "", v)

country_vars <- var_names[vapply(var_names, function(v) !is.na(suffix_of(v)), logical(1))]
if (!length(country_vars)) stop("No country macro variables found in IRF_post.")

# -----------------------------------------------------------------------------
# Transform draw-level IRFs
# -----------------------------------------------------------------------------
transformed <- IRF_post
cum_suffix <- c("dp", "de", "deq")

for (v in country_vars) {
  vv <- match(v, var_names)
  suf <- suffix_of(v)
  if (suf %in% cum_suffix) {
    # IRF_post[, vv, , ] is date x horizon x draw.
    for (tt in seq_len(dim(IRF_post)[1])) {
      m <- IRF_post[tt, vv, , , drop = FALSE]
      m <- matrix(m, nrow = dim(IRF_post)[3], ncol = dim(IRF_post)[4])
      transformed[tt, vv, , ] <- apply(m, 2L, cumsum)
    }
  }
}

# -----------------------------------------------------------------------------
# Summary helper: median + 68% and 90% credible intervals
# -----------------------------------------------------------------------------
qfun <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(low90=NA, low68=NA, median=NA, high68=NA, high90=NA, n=0))
  qq <- quantile(x, probs = c(.05, .16, .50, .84, .95), na.rm = TRUE, names = FALSE)
  c(low90=qq[1], low68=qq[2], median=qq[3], high68=qq[4], high90=qq[5], n=length(x))
}

summarise_array <- function(stable_only = FALSE) {
  rows <- vector("list", length(selected_dates) * length(country_vars) * length(h_num))
  rr <- 0L
  for (dd in selected_dates) {
    tt <- match(dd, irf_dates)
    keep_draw <- if (stable_only) which(stable_mask[tt, ]) else seq_len(dim(IRF_post)[4])
    if (!length(keep_draw)) next
    for (v in country_vars) {
      vv <- match(v, var_names)
      suf <- suffix_of(v)
      for (hh in seq_along(h_num)) {
        z <- transformed[tt, vv, hh, keep_draw]
        qs <- qfun(as.numeric(z))
        rr <- rr + 1L
        rows[[rr]] <- data.frame(
          date = dd,
          country = country_of(v),
          variable = v,
          concept = suf,
          representation = if (suf %in% cum_suffix) "cumulative" else "direct",
          horizon = h_num[hh],
          median = qs[["median"]],
          low68 = qs[["low68"]],
          high68 = qs[["high68"]],
          low90 = qs[["low90"]],
          high90 = qs[["high90"]],
          draws = as.integer(qs[["n"]]),
          stable_only = stable_only,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows[seq_len(rr)])
}

sum_all <- summarise_array(FALSE)
sum_stable <- summarise_array(TRUE)

write.csv(sum_all,
          file.path(out_dir, "mixed_transform_irf_summary_all_draws.csv"),
          row.names = FALSE)
write.csv(sum_stable,
          file.path(out_dir, "mixed_transform_irf_summary_stable_only.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Sign diagnostics for the transformed economic effects
# -----------------------------------------------------------------------------
sign_diag <- do.call(rbind, lapply(split(sum_stable, list(sum_stable$date, sum_stable$concept, sum_stable$horizon), drop=TRUE), function(d) {
  data.frame(
    date = d$date[1], concept = d$concept[1], horizon = d$horizon[1],
    representation = d$representation[1],
    positive_median = sum(d$median > 0, na.rm = TRUE),
    negative_median = sum(d$median < 0, na.rm = TRUE),
    zero_median = sum(d$median == 0, na.rm = TRUE),
    significant_positive_68 = sum(d$low68 > 0, na.rm = TRUE),
    significant_negative_68 = sum(d$high68 < 0, na.rm = TRUE),
    significant_positive_90 = sum(d$low90 > 0, na.rm = TRUE),
    significant_negative_90 = sum(d$high90 < 0, na.rm = TRUE),
    economies = nrow(d), stringsAsFactors = FALSE
  )
}))
write.csv(sign_diag, file.path(out_dir, "mixed_transform_sign_diagnostic_stable_only.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Plot stable-only results: one variable x date, 14 country facets.
# 90% ribbon = lighter; 68% ribbon = darker.
# No fixed colors are required; base ggplot defaults are used.
# -----------------------------------------------------------------------------
concept_titles <- c(
  y   = "GDP log-level response (direct IRF)",
  dp  = "CPI price-level response (cumulative inflation IRF)",
  r   = "Interest-rate response (direct IRF)",
  de  = "REER level response (cumulative change IRF)",
  deq = "Equity-price level response (cumulative return IRF)"
)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (ccpt in names(concept_titles)) {
  for (dd in selected_dates) {
    d <- subset(sum_stable, concept == ccpt & date == dd)
    if (!nrow(d)) next
    p <- ggplot(d, aes(x = horizon, y = median)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
      geom_ribbon(aes(ymin = low90, ymax = high90), alpha = 0.15) +
      geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.30) +
      geom_line(linewidth = 0.55) +
      facet_wrap(~ country, scales = "free_y", ncol = 4) +
      labs(
        title = paste0(concept_titles[[ccpt]], " — ", dd),
        subtitle = "Stable posterior draws; median with 68% and 90% credible intervals",
        x = "Horizon (quarters)", y = "Response in model units"
      ) +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

    ggsave(
      filename = file.path(plot_dir, paste0(ccpt, "_", dd, "_stable.png")),
      plot = p, width = 13, height = 9, dpi = 220
    )
  }
}

# -----------------------------------------------------------------------------
# Small text audit so the output cannot be misread later.
# -----------------------------------------------------------------------------
audit <- c(
  "TVP-GVAR post-processing only: NO MCMC re-estimation was performed.",
  paste0("Source IRF RDA: ", input_rda),
  paste0("Dates summarized: ", paste(selected_dates, collapse = ", ")),
  "GDP y: direct IRF because current run uses log GDP level.",
  "CPI dp: draw-by-draw cumulative IRF because dp is a differenced/log-change variable.",
  "Interest rate r: direct IRF because r is in levels.",
  "REER de: draw-by-draw cumulative IRF because de is a differenced/log-change variable.",
  "Equity deq: draw-by-draw cumulative IRF because deq is a differenced/log-change variable.",
  "Credible intervals are computed AFTER transformation/cumulation.",
  "Do not interpret a single-period dp/de/deq response as a permanent level effect."
)
writeLines(audit, file.path(out_dir, "README_interpretation.txt"))
cat(paste(audit, collapse = "\n"), "\n")
cat("\nOutputs saved under: ", out_dir, "\n", sep = "")
