#!/usr/bin/env Rscript

# Preview IRFs from an EXISTING P2_strong GitHub Actions artifact.
# This script DOES NOT re-estimate the TVP-GVAR.
# It reads the stable-only posterior summaries already produced by
# run_gpr_tvp_gvar_structural.R and creates preview figures.

suppressPackageStartupMessages({
  library(ggplot2)
})

source_dir <- Sys.getenv("TVPGVAR_IRF_SOURCE_DIR", "source-artifact")
out_dir <- Sys.getenv("TVPGVAR_IRF_PREVIEW_OUT", "irf-preview")
source_run_id <- Sys.getenv("TVPGVAR_SOURCE_RUN_ID", "")
source_artifact <- Sys.getenv(
  "TVPGVAR_SOURCE_ARTIFACT",
  "tvp-gvar-stability-v2-P2_strong"
)
key_countries_text <- Sys.getenv(
  "TVPGVAR_PREVIEW_COUNTRIES",
  "US,JP,CN,KR,EA,SG"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "key_overview"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "by_variable"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "cumulative_gdp"), showWarnings = FALSE, recursive = TRUE)

find_one <- function(filename, required = TRUE) {
  hits <- list.files(
    source_dir,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(hits)) {
    if (required) stop("Required file not found in artifact: ", filename)
    return(NA_character_)
  }
  if (length(hits) > 1L) {
    stop("Multiple copies found for ", filename, ": ", paste(hits, collapse = ", "))
  }
  hits[[1L]]
}

read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

irf_file <- find_one("gpr_structural_irf_summary_stable_only.csv")
cum_file <- find_one("gpr_structural_cumulative_gdp_stable_only.csv")
stab_file <- find_one("stability_summary_selected_dates.csv")
param_file <- find_one("stability_profile_parameters.csv", required = FALSE)
shock_file <- find_one("gpr_shock_validation.csv", required = FALSE)

irf <- read_csv(irf_file)
cumgdp <- read_csv(cum_file)
stab <- read_csv(stab_file)

required_irf <- c(
  "date", "variable", "horizon", "median",
  "low68", "high68", "low90", "high90", "draws_used"
)
missing_irf <- setdiff(required_irf, names(irf))
if (length(missing_irf)) {
  stop("IRF summary is missing columns: ", paste(missing_irf, collapse = ", "))
}

if (!all(c("date", "draws", "stable_draws", "near_unit_draws") %in% names(stab))) {
  stop("stability_summary_selected_dates.csv does not contain expected draw-count columns.")
}

# Metadata from the actual completed grid cell, not hard-coded lag labels.
profile_label <- "P2_strong"
spec_label <- ""
param_lines <- character(0)

if (!is.na(param_file)) {
  params <- read_csv(param_file)
  if (nrow(params) != 1L) stop("Expected exactly one row in stability_profile_parameters.csv.")
  profile_label <- params$profile[[1L]]
  if (!identical(profile_label, "P2_strong")) {
    stop("Wrong artifact/profile supplied. Expected P2_strong, found: ", profile_label)
  }
  spec_label <- sprintf(
    "GVAR(%s,%s); B1=%s; B2=%s; kappa0=%s",
    params$p[[1L]], params$q[[1L]],
    params$B_1[[1L]], params$B_2[[1L]], params$kappa0[[1L]]
  )
  param_lines <- c(
    paste0("Profile: ", profile_label),
    paste0("Specification: GVAR(", params$p[[1L]], ",", params$q[[1L]], ")"),
    paste0("B_1: ", params$B_1[[1L]]),
    paste0("B_2: ", params$B_2[[1L]]),
    paste0("kappa0: ", params$kappa0[[1L]]),
    paste0("Recorded saves: ", params$saves[[1L]]),
    paste0("Recorded burns: ", params$burns[[1L]]),
    paste0("Recorded thin: ", params$thin[[1L]]),
    paste0("GPR lag mode: ", params$gpr_lag_mode[[1L]])
  )
}

shock_label <- "GPR shock"
if (!is.na(shock_file)) {
  shock <- read_csv(shock_file)
  if ("target_log_jump" %in% names(shock)) {
    z <- shock$target_log_jump[is.finite(shock$target_log_jump)]
    if (length(z)) {
      shock_pct <- 100 * expm1(median(z))
      shock_label <- sprintf("GPR shock = %.1f%%", shock_pct)
    }
  }
}

# Parse country / variable suffix from names such as JP_y, US_r, EA_deq.
irf$country <- sub("_.*$", "", irf$variable)
irf$suffix <- sub("^[^_]+_", "", irf$variable)

cumgdp$country <- sub("_.*$", "", cumgdp$variable)
cumgdp$suffix <- sub("^[^_]+_", "", cumgdp$variable)

variable_labels <- c(
  y = "Real GDP growth",
  dp = "Inflation",
  r = "Short-term interest rate",
  de = "REER change",
  deq = "Equity return"
)

macro_suffixes <- names(variable_labels)
irf <- irf[irf$suffix %in% macro_suffixes, , drop = FALSE]
if (!nrow(irf)) stop("No macro IRF rows were found after filtering y/dp/r/de/deq.")

irf$variable_label <- unname(variable_labels[irf$suffix])
key_countries <- trimws(strsplit(key_countries_text, ",", fixed = TRUE)[[1L]])
available_countries <- sort(unique(irf$country))
key_countries <- key_countries[key_countries %in% available_countries]
if (!length(key_countries)) {
  stop("None of TVPGVAR_PREVIEW_COUNTRIES are present in the IRF summary.")
}

plot_dates <- unique(as.character(irf$date))
plot_dates <- plot_dates[plot_dates %in% as.character(stab$date)]

date_meta <- function(d) {
  z <- stab[as.character(stab$date) == d, , drop = FALSE]
  if (!nrow(z)) return("")
  sprintf(
    "stable draws: %d/%d (%.1f%%); near-unit stable draws: %d",
    as.integer(z$stable_draws[[1L]]),
    as.integer(z$draws[[1L]]),
    100 * as.numeric(z$stable_draws[[1L]]) / as.numeric(z$draws[[1L]]),
    as.integer(z$near_unit_draws[[1L]])
  )
}

theme_preview <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

plot_irf <- function(z, title, subtitle, facet) {
  p <- ggplot(z, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = low90, ymax = high90),
                alpha = 0.12, na.rm = TRUE) +
    geom_ribbon(aes(ymin = low68, ymax = high68),
                alpha = 0.25, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_line(linewidth = 0.65, na.rm = TRUE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Horizon (quarters)",
      y = "Response",
      caption = "Line: posterior median; darker band: 68%; lighter band: 90%. Stable draws only."
    ) +
    theme_preview

  if (facet == "country") {
    p <- p + facet_wrap(~country, ncol = 4)
  } else if (facet == "key_grid") {
    p <- p + facet_grid(variable_label ~ country, scales = "free_y")
  }
  p
}

# ---------------------------------------------------------------------------
# 1) Quick-look overview: six key countries x five variables, one image/date.
# ---------------------------------------------------------------------------
for (d in plot_dates) {
  z <- irf[as.character(irf$date) == d & irf$country %in% key_countries, , drop = FALSE]
  z$country <- factor(z$country, levels = key_countries)
  z$variable_label <- factor(
    z$variable_label,
    levels = unname(variable_labels[macro_suffixes])
  )

  p <- plot_irf(
    z,
    title = paste0("P2 strong stable-only IRF preview — ", d),
    subtitle = paste(shock_label, date_meta(d), spec_label, sep = " | "),
    facet = "key_grid"
  )

  ggsave(
    filename = file.path(out_dir, "key_overview", paste0("IRF_key_", d, ".png")),
    plot = p,
    width = 18, height = 12, dpi = 170
  )
}

# ---------------------------------------------------------------------------
# 2) All 14 countries: one image per date x macro variable.
# ---------------------------------------------------------------------------
for (d in plot_dates) {
  for (s in macro_suffixes) {
    z <- irf[
      as.character(irf$date) == d & irf$suffix == s,
      , drop = FALSE
    ]
    if (!nrow(z)) next

    p <- plot_irf(
      z,
      title = paste0(variable_labels[[s]], " response — ", d),
      subtitle = paste(
        profile_label, shock_label, date_meta(d),
        "stable posterior draws only",
        sep = " | "
      ),
      facet = "country"
    )

    ggsave(
      filename = file.path(
        out_dir, "by_variable",
        paste0("IRF_", d, "_", s, "_all_countries.png")
      ),
      plot = p,
      width = 14, height = 10, dpi = 170
    )
  }
}

# ---------------------------------------------------------------------------
# 3) Cumulative real-GDP response, all countries.
#    Existing structural runner defines this as cumulative y response.
# ---------------------------------------------------------------------------
required_cum <- c(
  "date", "variable", "horizon", "median",
  "low68", "high68", "low90", "high90", "draws_used"
)
missing_cum <- setdiff(required_cum, names(cumgdp))
if (length(missing_cum)) {
  stop("Cumulative GDP summary is missing columns: ",
       paste(missing_cum, collapse = ", "))
}

for (d in intersect(plot_dates, unique(as.character(cumgdp$date)))) {
  z <- cumgdp[as.character(cumgdp$date) == d, , drop = FALSE]
  if (!nrow(z)) next

  p <- plot_irf(
    z,
    title = paste0("Cumulative real-GDP response — ", d),
    subtitle = paste(
      profile_label, shock_label, date_meta(d),
      "stable posterior draws only",
      sep = " | "
    ),
    facet = "country"
  )

  ggsave(
    filename = file.path(
      out_dir, "cumulative_gdp",
      paste0("Cumulative_GDP_", d, "_all_countries.png")
    ),
    plot = p,
    width = 14, height = 10, dpi = 170
  )
}

# ---------------------------------------------------------------------------
# 4) Compact CSVs for quick numerical inspection at standard horizons.
# ---------------------------------------------------------------------------
key_horizons <- intersect(c(0, 1, 4, 8, 12), sort(unique(irf$horizon)))
quick <- irf[
  irf$country %in% key_countries & irf$horizon %in% key_horizons,
  c("date", "country", "suffix", "variable_label", "horizon",
    "median", "low68", "high68", "low90", "high90", "draws_used"),
  drop = FALSE
]
write.csv(
  quick,
  file.path(out_dir, "IRF_key_countries_key_horizons.csv"),
  row.names = FALSE
)

cum_quick <- cumgdp[
  cumgdp$country %in% key_countries & cumgdp$horizon %in% key_horizons,
  c("date", "country", "horizon",
    "median", "low68", "high68", "low90", "high90", "draws_used"),
  drop = FALSE
]
write.csv(
  cum_quick,
  file.path(out_dir, "Cumulative_GDP_key_countries_key_horizons.csv"),
  row.names = FALSE
)

# Preserve the exact selected-date stability counts next to the figures.
write.csv(
  stab,
  file.path(out_dir, "stability_summary_selected_dates_SOURCE.csv"),
  row.names = FALSE
)

readme <- c(
  "P2_strong IRF preview",
  "=====================",
  "",
  "Purpose: quick visual sanity check using an EXISTING 500-save grid run.",
  "No TVP-GVAR re-estimation is performed by this plotting script.",
  "IRFs are taken from gpr_structural_irf_summary_stable_only.csv.",
  "Unstable date/draw slices were already set to NA by the structural runner.",
  "Near-unit but stable draws are retained, consistent with the source model code.",
  "",
  paste0("Source workflow run ID: ", source_run_id),
  paste0("Source artifact: ", source_artifact),
  param_lines,
  paste0("Shock normalization: ", shock_label),
  paste0("Key preview countries: ", paste(key_countries, collapse = ", ")),
  "",
  "Figure conventions:",
  "- solid line: stable-only posterior median",
  "- darker ribbon: 68% posterior interval",
  "- lighter ribbon: 90% posterior interval",
  "- dashed horizontal line: zero response",
  "",
  "This is a PREVIEW diagnostic. Do not treat the 500-save figures as final paper output.",
  "Use the final larger MCMC run for publication-quality posterior inference."
)
writeLines(readme, file.path(out_dir, "IRF_PREVIEW_README.txt"))

cat("IRF preview complete.\n")
cat("Output directory:", normalizePath(out_dir), "\n")
cat("Source IRF file:", irf_file, "\n")
cat("Dates:", paste(plot_dates, collapse = ", "), "\n")
cat("Key countries:", paste(key_countries, collapse = ", "), "\n")
