#!/usr/bin/env Rscript

# Patch the existing R/BVAR_ttvp.r without overwriting it.
# Purpose:
#   current_only    : non-US equations include current global GPR but NOT a
#                     separate direct lagged-GPR regressor.
#   current_and_lag : reproduce the original unrestricted specification.
#
# The patched estimator is written to R/BVAR_ttvp_gprlag.r.

src <- "R/BVAR_ttvp.r"
out <- "R/BVAR_ttvp_gprlag.r"

if (!file.exists(src)) stop("Missing source file: ", src)

x <- readLines(src, warn = FALSE)
text <- paste(x, collapse = "\n")

# Current repository uses p=1. Replace the lag construction of Wex only;
# endogenous Y lags and all other lag operations remain untouched.
pat <- "mlag\\(\\s*Wex\\s*,\\s*(?:1|1L|p)\\s*\\)"
hits <- gregexpr(pat, text, perl = TRUE)[[1]]
nhit <- if (identical(hits[1], -1L)) 0L else length(hits)
if (nhit < 1L) {
  stop(
    "Could not locate mlag(Wex, 1/p) in R/BVAR_ttvp.r. ",
    "The upstream BVAR file may have changed; no patch was applied."
  )
}

text <- gsub(
  pat,
  "tvpgvar_wex_lag(Wex, 1L, Names[[nr]], gpr_lag_mode)",
  text,
  perl = TRUE
)

helper <- paste0(
'\n# -----------------------------------------------------------------------------\n',
'# Added by 8.12/patch_BVAR_gpr_lag.R\n',
'# Non-US Wex ordering is required to be:\n',
'# foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq, global_gpr.\n',
'# US Wex contains only the five foreign macro variables.\n',
'# -----------------------------------------------------------------------------\n',
'tvpgvar_wex_lag <- function(Wex, lag_order = 1L, country_name,\n',
'                             mode = c("current_only", "current_and_lag")) {\n',
'  mode <- match.arg(mode)\n',
'  if (as.integer(lag_order) != 1L) {\n',
'    stop("This TVP-GVAR GPR patch is designed for lag_order = 1.")\n',
'  }\n',
'  L <- mlag(Wex, lag_order)\n',
'  if (mode == "current_only" && !identical(as.character(country_name), "US")) {\n',
'    if (ncol(L) != 6L) {\n',
'      stop("Expected six non-US Wex variables (5 foreign macro + global GPR), found ", ncol(L))\n',
'    }\n',
'    # GPR is the sixth/last non-US Wex variable. Keep its contemporaneous term\n',
'    # in Wex, but remove only its direct lagged regressor from Wexlag1.\n',
'    L <- L[, 1:5, drop = FALSE]\n',
'  }\n',
'  L\n',
'}\n\n'
)

text <- paste0(helper, text)
writeLines(strsplit(text, "\n", fixed = TRUE)[[1]], out)

cat("Created patched estimator:", out, "\n")
cat("Replaced Wex lag construction occurrences:", nhit, "\n")
cat("Original estimator left unchanged:", src, "\n")
