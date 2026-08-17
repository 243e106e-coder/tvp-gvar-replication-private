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
#   p is selected at runtime from TVPGVAR_P (supported: 1 or 2)
#   pwex = 1, so foreign variables remain first-order even when p = 2
#   non-US Wex = 5 foreign macro variables + global GPR (6 columns)
#   US Wex     = 5 foreign macro variables

src <- "R/BVAR_ttvp.r"
out <- "R/BVAR_ttvp_gprlag.r"

if (!file.exists(src)) stop("Missing source file: ", src)

x <- readLines(src, warn = FALSE)
text <- paste(x, collapse = "\n")

# -----------------------------------------------------------------------------
# Patch 0: make the domestic lag order runtime-configurable.
#
# The original estimator hard-codes p <<- 1.  The generated estimator accepts
# TVPGVAR_P=1 or TVPGVAR_P=2, while pwex remains fixed at one.  This gives the
# intended GVAR(p, q) comparison with p in {1,2} and q=1.
# -----------------------------------------------------------------------------
pat_p <- 'p\\s*<<-\\s*1\\s*#number of lags of the dependent variable'
hit_p <- gregexpr(pat_p, text, perl = TRUE)[[1]]
n_p <- if (identical(hit_p[1], -1L)) 0L else length(hit_p)

if (n_p != 1L) {
  stop(
    "Expected exactly one hard-coded domestic lag line in ", src,
    ", found ", n_p, ". The upstream BVAR file may have changed; no patch was written."
  )
}

p_runtime <- paste(
  'p <<- as.integer(Sys.getenv("TVPGVAR_P", "1"))',
  'if (length(p) != 1L || is.na(p) || !p %in% c(1L, 2L)) {',
  '  stop("TVPGVAR_P must be 1 or 2.")',
  '}',
  '# number of lags of the dependent variable',
  sep = "\n  "
)

text <- sub(pat_p, p_runtime, text, perl = TRUE)

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
# Patch 2: separate domestic-lag and foreign-lag naming loops.
#
# The upstream code creates both name vectors inside `for (ii in 1:p)`. That
# only works accidentally when p=q=1. With p=2 and q=1 it would generate twice
# as many Wex-lag names as Wex-lag columns. The replacement below uses p for
# domestic Y lags and pwex for foreign lags. Under current_only, the foreign
# name count is based on the already-trimmed Wexlag matrix.
# -----------------------------------------------------------------------------
pat_names <- paste0(
  'for \\(ii in 1:p\\)\\{\\s*',
  'nameslags <- c\\(nameslags,rep\\(paste\\("Ylag",ii,sep=""\\),ncol\\(Yraw\\)\\)\\)\\s*',
  'wexnameslags <- c\\(wexnameslags,rep\\(paste\\("Wexlag",ii,sep=""\\),ncol\\(Wex\\)\\)\\)\\s*',
  '\\}'
)
hit_names <- gregexpr(pat_names, text, perl = TRUE)[[1]]
n_names <- if (identical(hit_names[1], -1L)) 0L else length(hit_names)

if (n_names != 1L) {
  stop(
    "Expected exactly one combined Y/Wex lag-name loop in ", src,
    ", found ", n_names, ". The upstream BVAR file may have changed; no patch was written."
  )
}

names_runtime <- paste(
  'for (ii in seq_len(p)) {',
  '  nameslags <- c(nameslags, rep(paste("Ylag", ii, sep = ""), ncol(Yraw)))',
  '}',
  'for (ii in seq_len(pwex)) {',
  '  block_width <- ncol(Wexlag) / pwex',
  '  if (!isTRUE(all.equal(block_width, as.integer(block_width)))) {',
  '    stop("Wexlag columns are not divisible by pwex.")',
  '  }',
  '  wexnameslags <- c(',
  '    wexnameslags,',
  '    rep(paste("Wexlag", ii, sep = ""), as.integer(block_width))',
  '  )',
  '}',
  sep = "\n  "
)

text <- sub(
  pat_names,
  names_runtime,
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
if (!grepl('Sys.getenv\\("TVPGVAR_P", "1"\\)', text, perl = TRUE)) {
  stop("Internal patch validation failed: runtime domestic lag order was not inserted.")
}
if (!grepl('for \\(ii in seq_len\\(pwex\\)\\)', text, perl = TRUE) ||
    !grepl('block_width <- ncol\\(Wexlag\\) / pwex', text, perl = TRUE)) {
  stop("Internal patch validation failed: p/q lag-name loops were not separated.")
}
if (grepl('Wexlag\\s*<-\\s*mlag\\(\\s*Wex\\s*,\\s*pwex\\s*\\)', text, perl = TRUE)) {
  stop("Internal patch validation failed: original Wexlag construction remains.")
}

writeLines(strsplit(text, "\n", fixed = TRUE)[[1]], out)

cat("Created patched estimator:", out, "\n")
cat("Original estimator left unchanged:", src, "\n")
cat("Patched Wex lag construction occurrences:", n_lag, "\n")
cat("Separated domestic/foreign lag-name loops:", n_names, "\n")
cat("Patched domestic lag-order occurrences:", n_p, "\n")
cat("Domestic lag order is selected at runtime with TVPGVAR_P=1 or 2.\n")
cat("Supported modes: current_only, current_and_lag\n")
