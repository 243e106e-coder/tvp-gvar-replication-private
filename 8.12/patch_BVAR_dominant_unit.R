#!/usr/bin/env Rscript

# =============================================================================
# Create a TVP-BVAR estimator that supports a separate dominant unit with ZERO
# weakly-exogenous regressors.
#
# Starting point: R/BVAR_ttvp.r
# Output:         R/BVAR_ttvp_dominant.r
#
# Changes relative to upstream:
#   1) runtime domestic lag p = TVPGVAR_P (1 or 2), same for all 15 units;
#   2) foreign lag q remains 1;
#   3) a unit whose W matrix has only own rows (GL) receives a valid 0-column
#      Wex/Wexlag matrix instead of the upstream invalid (k+1):k slice;
#   4) p/q coefficient-name loops are separated;
#   5) singular contemporaneous matrices use a numeric generalized inverse;
#   6) design matrices are checked before expensive TVP estimation.
# =============================================================================

src <- "R/BVAR_ttvp.r"
out <- "R/BVAR_ttvp_dominant.r"
if (!file.exists(src)) stop("Missing source estimator: ", src)

text <- paste(readLines(src, warn = FALSE), collapse = "\n")

count_regex <- function(pattern, x) {
  h <- gregexpr(pattern, x, perl = TRUE)[[1]]
  if (identical(h[1], -1L)) 0L else length(h)
}

# Runtime p.
pat_p <- 'p\\s*<<-\\s*1\\s*#number of lags of the dependent variable'
if (count_regex(pat_p, text) != 1L) stop("Could not uniquely patch p.")
text <- sub(
  pat_p,
  paste(
    'p <<- as.integer(Sys.getenv("TVPGVAR_P", "2"))',
    'if (length(p) != 1L || is.na(p) || !p %in% c(1L, 2L)) {',
    '  stop("TVPGVAR_P must be 1 or 2.")',
    '}',
    '# number of lags of the dependent variable',
    sep = "\n  "
  ),
  text,
  perl = TRUE
)

helpers <- paste0(
  '# -----------------------------------------------------------------------------\n',
  '# Added by patch_BVAR_dominant_unit.R\n',
  '# -----------------------------------------------------------------------------\n',
  'tvpgvar_extract_wex <- function(W, xglobal, k_i) {\n',
  '  all <- as.matrix(W) %*% t(as.matrix(xglobal))\n',
  '  if (nrow(all) < k_i) stop("W has fewer rows than own variables.")\n',
  '  if (nrow(all) == k_i) {\n',
  '    return(matrix(0, nrow = nrow(xglobal), ncol = 0L))\n',
  '  }\n',
  '  z <- t(all[(k_i + 1L):nrow(all), , drop = FALSE])\n',
  '  storage.mode(z) <- "double"\n',
  '  z\n',
  '}\n\n',
  'tvpgvar_wex_lag <- function(Wex, lag_order = 1L) {\n',
  '  lag_order <- as.integer(lag_order)\n',
  '  if (lag_order != 1L) stop("Foreign lag q must equal 1.")\n',
  '  Wex <- as.matrix(Wex)\n',
  '  if (ncol(Wex) == 0L) return(matrix(0, nrow = nrow(Wex), ncol = 0L))\n',
  '  mlag(Wex, lag_order)\n',
  '}\n\n',
  'tvpgvar_safe_inverse <- function(A, context = "") {\n',
  '  A <- as.matrix(A); storage.mode(A) <- "double"\n',
  '  if (nrow(A) != ncol(A) || any(!is.finite(A))) stop("Invalid A0 ", context)\n',
  '  ans <- try(solve(A), silent = TRUE)\n',
  '  if (inherits(ans, "try-error") || !is.matrix(ans) || any(!is.finite(ans))) {\n',
  '    ans <- MASS::ginv(A)\n',
  '  }\n',
  '  ans <- as.matrix(ans); storage.mode(ans) <- "double"\n',
  '  if (any(!is.finite(ans))) stop("Could not invert A0 ", context)\n',
  '  ans\n',
  '}\n\n'
)

# Safe Wex extraction.
pat_wex <- paste0(
  'all\\s*<-\\s*W%\\*%t\\(xglobal\\)\\s*',
  'Wex\\s*<-\\s*all\\[\\(ncol\\(End\\)\\+1\\):nrow\\(all\\),\\]\\s*',
  'Wex\\s*<-\\s*t\\(Wex\\)\\s*',
  'class\\(Wex\\)\\s*<-\\s*"numeric"'
)
if (count_regex(pat_wex, text) != 1L) stop("Could not uniquely patch Wex extraction.")
text <- sub(
  pat_wex,
  'Wex <- tvpgvar_extract_wex(W, xglobal, ncol(End))',
  text,
  perl = TRUE
)

# Safe Wex lag.
pat_lag <- 'Wexlag\\s*<-\\s*mlag\\(Wex,pwex\\)'
if (count_regex(pat_lag, text) != 1L) stop("Could not uniquely patch Wexlag.")
text <- sub(pat_lag, 'Wexlag <- tvpgvar_wex_lag(Wex, pwex)', text, perl = TRUE)

# Separate lag naming loops.
pat_names <- paste0(
  'for \\(ii in 1:p\\)\\{\\s*',
  'nameslags <- c\\(nameslags,rep\\(paste\\("Ylag",ii,sep=""\\),ncol\\(Yraw\\)\\)\\)\\s*',
  'wexnameslags <- c\\(wexnameslags,rep\\(paste\\("Wexlag",ii,sep=""\\),ncol\\(Wex\\)\\)\\)\\s*',
  '\\}'
)
if (count_regex(pat_names, text) != 1L) stop("Could not uniquely patch lag names.")
text <- sub(
  pat_names,
  paste(
    'for (ii in seq_len(p)) {',
    '  nameslags <- c(nameslags, rep(paste0("Ylag", ii), ncol(Yraw)))',
    '}',
    'for (ii in seq_len(pwex)) {',
    '  block_width <- if (pwex > 0L) ncol(Wexlag) / pwex else 0L',
    '  if (!isTRUE(all.equal(block_width, as.integer(block_width)))) stop("Bad Wexlag width.")',
    '  wexnameslags <- c(wexnameslags, rep(paste0("Wexlag", ii), as.integer(block_width)))',
    '}',
    sep = "\n  "
  ),
  text,
  perl = TRUE
)

# Repair A0 inverse fallback.
pat_inv <- paste0(
  'A0inv\\s*<-\\s*try\\(solve\\(A0\\),silent=TRUE\\)\\s*',
  'if\\s*\\(is\\(A0inv,"try-error"\\)\\)\\s*A0inv\\s*<-\\s*ginv\\(A0inv\\)'
)
if (count_regex(pat_inv, text) != 1L) stop("Could not uniquely patch A0 inverse.")
text <- sub(
  pat_inv,
  paste(
    'A0inv <- tvpgvar_safe_inverse(',
    '  A0, context = paste0("unit=", Names[[nr]], ",draw=", ii, ",time=", nn)',
    ')',
    sep = "\n      "
  ),
  text,
  perl = TRUE
)

# Numeric/finite design guard.
pat_design <- 'X\\s*<-\\s*X1\\s*\\n\\s*Z\\s*<-\\s*kronecker\\(diag\\(N\\),X\\)'
if (count_regex(pat_design, text) != 1L) stop("Could not uniquely patch design guard.")
text <- sub(
  pat_design,
  paste(
    'X <- as.matrix(X1); storage.mode(X) <- "double"',
    'Y <- as.matrix(Y); storage.mode(Y) <- "double"',
    'if (nrow(X) != nrow(Y) || any(!is.finite(X)) || any(!is.finite(Y))) {',
    '  stop("Invalid design for ", Names[[nr]], ": X=", paste(dim(X), collapse="x"),',
    '       ", Y=", paste(dim(Y), collapse="x"))',
    '}',
    'Z <- kronecker(diag(N), X)',
    sep = "\n  "
  ),
  text,
  perl = TRUE
)

text <- paste0(helpers, text)

required <- c(
  'tvpgvar_extract_wex(W, xglobal, ncol(End))',
  'tvpgvar_wex_lag(Wex, pwex)',
  'Sys.getenv("TVPGVAR_P", "2")',
  'tvpgvar_safe_inverse('
)
for (tok in required) if (!grepl(tok, text, fixed = TRUE)) stop("Patch validation failed: ", tok)
if (grepl('ginv(A0inv)', text, fixed = TRUE)) stop("Broken inverse fallback remains.")

writeLines(strsplit(text, "\n", fixed = TRUE)[[1]], out)
cat("Created dominant-unit estimator: ", out, "\n", sep = "")
