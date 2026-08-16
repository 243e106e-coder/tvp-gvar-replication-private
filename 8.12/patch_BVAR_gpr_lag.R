#!/usr/bin/env Rscript

# Patch the existing R/BVAR_ttvp.r without overwriting it.
#
# Modes:
#   current_only    : for non-US equations, keep current global GPR in Wex,
#                     but remove the separate direct lagged-GPR regressor.
#   current_and_lag : reproduce the original unrestricted specification.
#
# Output:
#   R/BVAR_ttvp_gprlag.r
#
# This patch is intentionally narrow and is matched to the current repository
# implementation where:
#   p = 1
#   pwex = 1
#   non-US Wex = 5 foreign macro variables + global GPR (6 columns)
#   US Wex     = 5 foreign macro variables

src <- "R/BVAR_ttvp.r"
out <- "R/BVAR_ttvp_gprlag.r"

if (!file.exists(src)) stop("Missing source file: ", src)

x <- readLines(src, warn = FALSE)
text <- paste(x, collapse = "\n")

# -----------------------------------------------------------------------------
# Helper inserted ahead of BVAR().
# It computes the ordinary Wex lag and, only in current_only mode for non-US
# blocks, removes the sixth lagged column (global GPR). Current GPR remains in
# contemporaneous Wex, and the US block is untouched.
# -----------------------------------------------------------------------------
helper <- paste0(
  '# -----------------------------------------------------------------------------\n',
  '# Added by 8.12/patch_BVAR_gpr_lag.R\n',
  '# -----------------------------------------------------------------------------\n',
  'tvpgvar_wex_lag <- function(Wex, lag_order = 1L, country_name,\n',
  '                             mode = c("current_only", "current_and_lag")) {\n',
  '  mode <- match.arg(mode)\n',
  '  lag_order <- as.integer(lag_order)\n',
  '  if (length(lag_order) != 1L || is.na(lag_order) || lag_order != 1L) {\n',
  '    stop("This GPR-lag patch currently requires pwex = 1.")\n',
  '  }\n',
  '  L <- mlag(Wex, lag_order)\n',
  '  cc <- as.character(country_name)\n',
  '  if (mode == "current_only" && !identical(cc, "US")) {\n',
  '    if (ncol(L) != 6L) {\n',
  '      stop("Expected six non-US Wex lags (5 foreign macro + global GPR), found ", ncol(L),\n',
  '           " for ", cc, ".")\n',
  '    }\n',
  '    L <- L[, 1:5, drop = FALSE]\n',
  '  }\n',
  '  L\n',
  '}\n\n'
)

# -----------------------------------------------------------------------------
# Patch 1: the current repository uses exactly
#   Wexlag <- mlag(Wex,pwex)
# Replace only this Wex-lag construction. Y lags and all other lag operations
# remain unchanged.
# -----------------------------------------------------------------------------
pat_lag <- 'Wexlag\\s*<-\\s*mlag\\(\\s*Wex\\s*,\\s*pwex\\s*\\)'
hit_lag <- gregexpr(pat_lag, text, perl = TRUE)[[1]]
n_lag <- if (identical(hit_lag[1], -1L)) 0L else length(hit_lag)

if (n_lag != 1L) {
  stop(
    "Expected exactly one `Wexlag <- mlag(Wex,pwex)` occurrence in ", src,
    ", found ", n_lag, ". The upstream BVAR file may have changed; no patch was written."
  )
}

text <- sub(
  pat_lag,
  'Wexlag <- tvpgvar_wex_lag(Wex, pwex, Names[[nr]], gpr_lag_mode)',
  text,
  perl = TRUE
)

# -----------------------------------------------------------------------------
# Patch 2: original code creates Wexlag1 names using ncol(Wex). Under
# current_only, non-US Wex still has 6 columns but Wexlag has only 5, so the
# name count must use ncol(Wexlag).
# -----------------------------------------------------------------------------
pat_names <- 'rep\\(paste\\("Wexlag",ii,sep=""\\),ncol\\(Wex\\)\\)'
hit_names <- gregexpr(pat_names, text, perl = TRUE)[[1]]
n_names <- if (identical(hit_names[1], -1L)) 0L else length(hit_names)

if (n_names != 1L) {
  stop(
    "Expected exactly one Wexlag name-construction occurrence using ncol(Wex) in ", src,
    ", found ", n_names, ". The upstream BVAR file may have changed; no patch was written."
  )
}

text <- sub(
  pat_names,
  'rep(paste("Wexlag",ii,sep=""),ncol(Wexlag))',
  text,
  perl = TRUE
)

# Insert helper before BVAR().
text <- paste0(helper, text)

# -----------------------------------------------------------------------------
# Sanity checks on the generated source before writing it.
# -----------------------------------------------------------------------------
if (!grepl('tvpgvar_wex_lag\\(Wex, pwex, Names\\[\\[nr\\]\\], gpr_lag_mode\\)', text, perl = TRUE)) {
  stop("Internal patch validation failed: patched Wexlag call not found.")
}
if (!grepl('rep\\(paste\\("Wexlag",ii,sep=""\\),ncol\\(Wexlag\\)\\)', text, perl = TRUE)) {
  stop("Internal patch validation failed: Wexlag name count was not changed.")
}
if (grepl('Wexlag\\s*<-\\s*mlag\\(\\s*Wex\\s*,\\s*pwex\\s*\\)', text, perl = TRUE)) {
  stop("Internal patch validation failed: original Wexlag construction remains.")
}

writeLines(strsplit(text, "\n", fixed = TRUE)[[1]], out)

cat("Created patched estimator:", out, "\n")
cat("Original estimator left unchanged:", src, "\n")
cat("Patched Wex lag construction occurrences:", n_lag, "\n")
cat("Patched Wexlag naming occurrences:", n_names, "\n")
cat("Supported modes: current_only, current_and_lag\n")
