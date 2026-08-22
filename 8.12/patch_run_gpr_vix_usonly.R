#!/usr/bin/env Rscript

# =============================================================================
# Patch the generated GPR+VIX run into the "US-endogenous-GPR only" experiment.
#
# Input:
#   8.12/run_gpr_vix_tvp_gvar_structural.R
#
# Output:
#   8.12/run_gpr_vix_usonly_tvp_gvar_structural.R
#
# The existing successful GPR+VIX experiment is never overwritten.
# =============================================================================

src <- "8.12/run_gpr_vix_tvp_gvar_structural.R"
out <- "8.12/run_gpr_vix_usonly_tvp_gvar_structural.R"

if (!file.exists(src)) {
  stop("Missing generated GPR+VIX run script: ", src)
}

text <- paste(
  readLines(src, warn = FALSE),
  collapse = "\n"
)

replace_once_fixed <- function(x, old, new, label) {
  pos <- gregexpr(old, x, fixed = TRUE)[[1]]
  n <- if (identical(pos[[1]], -1L)) 0L else length(pos)

  if (n != 1L) {
    stop(
      label,
      ": expected exactly one match, found ",
      n
    )
  }

  sub(old, new, x, fixed = TRUE)
}


# -----------------------------------------------------------------------------
# 1. Use the new data mapping: no direct global_gpr in non-US local systems
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'source("R/prepare_data_structural_gpr_vix.R")',
  'source("R/prepare_data_structural_gpr_vix_usonly.R")',
  "prepare-data source replacement"
)


# -----------------------------------------------------------------------------
# 2. Keep existing 7x7 structural GPR+VIX shock, then override Wex lag alignment
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'source("R/gpr_vix_overrides.R")',
  paste(
    'source("R/gpr_vix_overrides.R")',
    'source("R/gpr_vix_usonly_overrides.R")',
    sep = "\n"
  ),
  "US-only override insertion"
)


# -----------------------------------------------------------------------------
# 3. Feedback/exogeneity diagnostic labels
# -----------------------------------------------------------------------------
old_feedback <- paste(
  '  expected_labels <- if (cc == "US") {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq")',
  '  } else {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",',
  '      "global_gpr", "global_vix")',
  '  }',
  sep = "\n"
)

new_feedback <- paste(
  '  expected_labels <- if (cc == "US") {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq")',
  '  } else {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",',
  '      "global_vix")',
  '  }',
  sep = "\n"
)

text <- replace_once_fixed(
  text,
  old_feedback,
  new_feedback,
  "feedback-label replacement"
)


# -----------------------------------------------------------------------------
# 4. Estimated regressor dimensions
# -----------------------------------------------------------------------------
old_dims <- paste(
  'expected_lag <- if (gpr_lag_mode == "current_only") {',
  '  ifelse(cN == "US", 5L, 6L)',
  '} else {',
  '  ifelse(cN == "US", 5L, 7L)',
  '}',
  'expected_wex <- ifelse(cN == "US", 5L, 7L)',
  'expected_ylag <- ifelse(cN == "US", 7L, 5L)',
  'expected_ylag2 <- if (lag_order == 2L) expected_ylag else rep(0L, length(cN))',
  sep = "\n"
)

new_dims <- paste(
  'expected_lag <- ifelse(cN == "US", 5L, 6L)',
  'expected_wex <- ifelse(cN == "US", 5L, 6L)',
  'expected_ylag <- ifelse(cN == "US", 7L, 5L)',
  'expected_ylag2 <- if (lag_order == 2L) expected_ylag else rep(0L, length(cN))',
  sep = "\n"
)

text <- replace_once_fixed(
  text,
  old_dims,
  new_dims,
  "regressor-dimension replacement"
)


# -----------------------------------------------------------------------------
# 5. Separate posterior filename
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'save(predDens, Data.setup, file = "results/predDens_gpr_vix_structural.rda")',
  'save(predDens, Data.setup, file = "results/predDens_gpr_vix_usonly_structural.rda")',
  "posterior-output replacement"
)


# -----------------------------------------------------------------------------
# 6. Clear experiment banner
# -----------------------------------------------------------------------------
old_banner <- paste(
  'cat("\\n============================================================\\n")',
  'cat(" P2 GPR + VIX STRUCTURAL DIAGNOSTIC\\n")',
  'cat(" Identification: GPR -> VIX -> y -> dp -> r -> de -> deq\\n")',
  'cat(" Non-US direct lagged GPR: controlled by TVPGVAR_GPR_LAG_MODE\\n")',
  'cat("============================================================\\n\\n")',
  sep = "\n"
)

new_banner <- paste(
  'cat("\\n============================================================\\n")',
  'cat(" P2 STRONG: GPR ENDOGENOUS IN US ONLY + GLOBAL VIX\\n")',
  'cat(" US recursive order: GPR -> VIX -> y -> dp -> r -> de -> deq\\n")',
  'cat(" Non-US direct current GPR: REMOVED\\n")',
  'cat(" Non-US direct lagged GPR: REMOVED\\n")',
  'cat(" Non-US global VIX: retained current + q=1 lag\\n")',
  'cat("============================================================\\n\\n")',
  sep = "\n"
)

text <- replace_once_fixed(
  text,
  old_banner,
  new_banner,
  "experiment-banner replacement"
)


# -----------------------------------------------------------------------------
# 7. Add explicit structural specification audit after Data.setup is created
# -----------------------------------------------------------------------------
anchor <- 'cN <- Data.setup$countries'

audit <- paste(
  anchor,
  '',
  '# Audit: GPR is endogenous in US only; non-US foreign/global rows must not',
  '# load on the US_gpr column.',
  'gpr_col <- match("US_gpr", colnames(xglobal))',
  'if (is.na(gpr_col)) stop("US_gpr missing from global vector.")',
  'spec_audit <- do.call(rbind, lapply(seq_along(cN), function(i) {',
  '  cc <- cN[[i]]',
  '  Wi <- as.matrix(gW[[i]])',
  '  k_i <- if (cc == "US") 7L else 5L',
  '  foreign_rows <- if (nrow(Wi) > k_i) Wi[(k_i + 1L):nrow(Wi), , drop = FALSE] else matrix(0, 0, ncol(Wi))',
  '  gpr_loading <- if (cc == "US") NA_real_ else max(abs(foreign_rows[, gpr_col]), na.rm = TRUE)',
  '  data.frame(',
  '    country = cc,',
  '    local_mapping_rows = nrow(Wi),',
  '    endogenous_variables = k_i,',
  '    foreign_global_rows = nrow(Wi) - k_i,',
  '    max_nonUS_foreign_loading_on_US_gpr = gpr_loading,',
  '    stringsAsFactors = FALSE',
  '  )',
  '}))',
  'if (any(spec_audit$country != "US" & spec_audit$max_nonUS_foreign_loading_on_US_gpr > 1e-14)) {',
  '  print(spec_audit)',
  '  stop("US-only-GPR specification audit failed.")',
  '}',
  'write.csv(spec_audit, "results/gpr_usonly_mapping_audit.csv", row.names = FALSE)',
  'cat("GPR US-only mapping audit passed.\\n")',
  sep = "\n"
)

text <- replace_once_fixed(
  text,
  anchor,
  audit,
  "US-only specification audit insertion"
)


# -----------------------------------------------------------------------------
# 8. Static validation
# -----------------------------------------------------------------------------
required_tokens <- c(
  'source("R/prepare_data_structural_gpr_vix_usonly.R")',
  'source("R/gpr_vix_usonly_overrides.R")',
  'predDens_gpr_vix_usonly_structural.rda',
  'expected_wex <- ifelse(cN == "US", 5L, 6L)',
  'Non-US direct current GPR: REMOVED',
  'gpr_usonly_mapping_audit.csv'
)

for (tok in required_tokens) {
  if (!grepl(tok, text, fixed = TRUE)) {
    stop(
      "Generated US-only-GPR run script failed validation. Missing token: ",
      tok
    )
  }
}

writeLines(
  strsplit(text, "\n", fixed = TRUE)[[1]],
  out
)

cat("Created: ", out, "\n", sep = "")
cat("Source left unchanged: ", src, "\n", sep = "")
