#!/usr/bin/env Rscript

# =============================================================================
# Create a GPR+VIX version of the CURRENT repository main structural run script.
#
# Source:
#   8.12/run_gpr_tvp_gvar_structural.R
#
# Output:
#   8.12/run_gpr_vix_tvp_gvar_structural.R
#
# The original script is never overwritten.
# =============================================================================

src <- "8.12/run_gpr_tvp_gvar_structural.R"
out <- "8.12/run_gpr_vix_tvp_gvar_structural.R"

if (!file.exists(src)) {
  stop("Missing source run script: ", src)
}

text <- paste(
  readLines(
    src,
    warn = FALSE
  ),
  collapse = "\n"
)

replace_once_fixed <- function(
    x,
    old,
    new,
    label) {

  pos <- gregexpr(
    old,
    x,
    fixed = TRUE
  )[[1]]

  n <- if (
    identical(
      pos[[1]],
      -1L
    )
  ) {
    0L
  } else {
    length(pos)
  }

  if (n != 1L) {
    stop(
      label,
      ": expected exactly one match, found ",
      n
    )
  }

  sub(
    old,
    new,
    x,
    fixed = TRUE
  )
}


# -----------------------------------------------------------------------------
# 1. Source the VIX-aware data preparation
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'source("R/prepare_data_structural.R")',
  'source("R/prepare_data_structural_gpr_vix.R")',
  "prepare-data source replacement"
)


# -----------------------------------------------------------------------------
# 2. Source runtime VIX overrides AFTER the original IRF helper
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'source("R/gpr_structural_irf.R")',
  paste(
    'source("R/gpr_structural_irf.R")',
    'source("R/gpr_vix_overrides.R")',
    sep = "\n"
  ),
  "GPR-VIX override insertion"
)


# -----------------------------------------------------------------------------
# 3. Enforce GPR -> VIX -> macro ordering
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'expected_us <- c("US_gpr", "US_y", "US_dp", "US_r", "US_de", "US_deq")',
  'expected_us <- c("US_gpr", "US_vix", "US_y", "US_dp", "US_r", "US_de", "US_deq")',
  "US recursive-order replacement"
)


# -----------------------------------------------------------------------------
# 4. Update stationary feedback diagnostic labels
# -----------------------------------------------------------------------------
old_feedback <- paste(
  '  expected_labels <- if (cc == "US") {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq")',
  '  } else {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",',
  '      "global_gpr")',
  '  }',
  sep = "\n"
)

new_feedback <- paste(
  '  expected_labels <- if (cc == "US") {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq")',
  '  } else {',
  '    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",',
  '      "global_gpr", "global_vix")',
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
# 5. Update regressor-block dimension assertions
#
# current_only:
#   US Wexlag1:      5 foreign macro
#   non-US Wexlag1:  5 foreign macro + VIX = 6
#
# current_and_lag:
#   non-US Wexlag1:  5 foreign macro + GPR + VIX = 7
#
# Wex:
#   US = 5
#   non-US = 7
#
# Domestic lag blocks:
#   US = 7
#   non-US = 5
# -----------------------------------------------------------------------------
old_dims <- paste(
  'expected_lag <- if (gpr_lag_mode == "current_only") {',
  '  rep(5L, length(cN))',
  '} else {',
  '  ifelse(cN == "US", 5L, 6L)',
  '}',
  'expected_wex <- ifelse(cN == "US", 5L, 6L)',
  'expected_ylag <- ifelse(cN == "US", 6L, 5L)',
  'expected_ylag2 <- if (lag_order == 2L) expected_ylag else rep(0L, length(cN))',
  sep = "\n"
)

new_dims <- paste(
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

text <- replace_once_fixed(
  text,
  old_dims,
  new_dims,
  "regressor-dimension replacement"
)


# -----------------------------------------------------------------------------
# 6. Use a separate posterior file name for the experiment
# -----------------------------------------------------------------------------
text <- replace_once_fixed(
  text,
  'save(predDens, Data.setup, file = "results/predDens_gpr_structural.rda")',
  'save(predDens, Data.setup, file = "results/predDens_gpr_vix_structural.rda")',
  "posterior-output replacement"
)


# -----------------------------------------------------------------------------
# 7. Write an explicit experiment banner at startup
# -----------------------------------------------------------------------------
banner_anchor <- 'dir.create("results", showWarnings = FALSE, recursive = TRUE)'

banner_text <- paste(
  banner_anchor,
  'cat("\\n============================================================\\n")',
  'cat(" P2 GPR + VIX STRUCTURAL DIAGNOSTIC\\n")',
  'cat(" Identification: GPR -> VIX -> y -> dp -> r -> de -> deq\\n")',
  'cat(" Non-US direct lagged GPR: controlled by TVPGVAR_GPR_LAG_MODE\\n")',
  'cat("============================================================\\n\\n")',
  sep = "\n"
)

text <- replace_once_fixed(
  text,
  banner_anchor,
  banner_text,
  "experiment-banner insertion"
)


# -----------------------------------------------------------------------------
# 8. Static sanity checks
# -----------------------------------------------------------------------------
required_tokens <- c(
  'source("R/prepare_data_structural_gpr_vix.R")',
  'source("R/gpr_vix_overrides.R")',
  '"US_vix"',
  '"global_vix"',
  'predDens_gpr_vix_structural.rda'
)

for (tok in required_tokens) {
  if (!grepl(tok, text, fixed = TRUE)) {
    stop(
      "Generated run script failed validation. Missing token: ",
      tok
    )
  }
}

if (
  grepl(
    'expected_us <- c("US_gpr", "US_y"',
    text,
    fixed = TRUE
  )
) {
  stop(
    "Old US ordering still remains in generated run script."
  )
}

writeLines(
  strsplit(
    text,
    "\n",
    fixed = TRUE
  )[[1]],
  out
)

cat("Created: ", out, "\n", sep = "")
cat("Original left unchanged: ", src, "\n", sep = "")
