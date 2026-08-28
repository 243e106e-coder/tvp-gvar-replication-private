#!/usr/bin/env Rscript

# ============================================================
# Exact reconciliation:
# 05_country_specific_lag_selection.R
# vs
# 06_cn_za_instability_diagnostic.R
#
# Cases:
#   CN p=2
#   ZA p=2
#   JP p=1
#
# This script does NOT change the model. It sources both existing
# diagnostics in isolated environments, intercepts the lag-selection
# companion-matrix input, and compares it with the A matrices produced
# by the instability diagnostic.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_05 <- "8.12/05_country_specific_lag_selection.R"
SCRIPT_06 <- "8.12/06_cn_za_instability_diagnostic.R"
OUT_DIR   <- "8.12/lag_reconciliation"

CASES <- data.frame(
  country = c("CN","ZA","JP"),
  p = c(2L,2L,1L),
  stringsAsFactors = FALSE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg <- function(...) cat(sprintf(...), "\n")

if(!file.exists(SCRIPT_05)) stopf("Missing %s", SCRIPT_05)
if(!file.exists(SCRIPT_06)) stopf("Missing %s", SCRIPT_06)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
companion_from_A <- function(A) {
  p <- length(A)
  if(p < 1L) stopf("No lag matrices supplied")
  k <- nrow(A[[1]])
  if(p == 1L) return(A[[1]])

  rbind(
    do.call(cbind, A),
    cbind(
      diag(k * (p - 1L)),
      matrix(0, k * (p - 1L), k)
    )
  )
}

eigen_table <- function(C, source, country, p) {
  ev <- eigen(C, only.values = TRUE)$values
  o <- order(Mod(ev), decreasing = TRUE)
  data.frame(
    source = source,
    country = country,
    p = p,
    rank = seq_along(o),
    real = Re(ev[o]),
    imag = Im(ev[o]),
    modulus = Mod(ev[o]),
    angle = Arg(ev[o]),
    stringsAsFactors = FALSE
  )
}

matrix_long <- function(M, source, country, p, object, lag = NA_integer_) {
  nr <- nrow(M); nc <- ncol(M)
  rn <- rownames(M); cn <- colnames(M)
  if(is.null(rn)) rn <- paste0("r", seq_len(nr))
  if(is.null(cn)) cn <- paste0("c", seq_len(nc))

  z <- expand.grid(
    row = seq_len(nr),
    col = seq_len(nc),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  z$row_name <- rn[z$row]
  z$col_name <- cn[z$col]
  z$value <- M[cbind(z$row, z$col)]
  z$source <- source
  z$country <- country
  z$p <- p
  z$object <- object
  z$lag <- lag

  z[, c("source","country","p","object","lag",
        "row","col","row_name","col_name","value")]
}

safe_maxdiff <- function(A, B) {
  if(is.null(A) || is.null(B)) return(NA_real_)
  if(!identical(dim(A), dim(B))) return(Inf)
  max(abs(as.numeric(A) - as.numeric(B)), na.rm = TRUE)
}

safe_rho <- function(C) {
  max(Mod(eigen(C, only.values = TRUE)$values))
}

safe_qid <- function(x) {
  if(is.null(x)) return(integer())
  if(is.data.frame(x) && "qid" %in% names(x)) return(as.integer(x$qid))
  integer()
}

# ------------------------------------------------------------
# Source both scripts in isolated environments.
# This intentionally runs their existing main diagnostics first.
# ------------------------------------------------------------
msg("Sourcing lag-selection script...")
e05 <- new.env(parent = globalenv())
sys.source(SCRIPT_05, envir = e05)

msg("Sourcing instability-diagnostic script...")
e06 <- new.env(parent = globalenv())
sys.source(SCRIPT_06, envir = e06)

required05 <- c("fit_country","companion_radius","panel","foreign","globals","VARS")
miss05 <- setdiff(required05, ls(e05))
if(length(miss05)) {
  writeLines(ls(e05), file.path(OUT_DIR, "00_objects_env05.txt"))
  stopf("05 environment missing: %s. See 00_objects_env05.txt",
        paste(miss05, collapse=", "))
}

required06 <- c("fit_local_varx","panel","VARS")
miss06 <- setdiff(required06, ls(e06))
if(length(miss06)) {
  writeLines(ls(e06), file.path(OUT_DIR, "00_objects_env06.txt"))
  stopf("06 environment missing: %s. See 00_objects_env06.txt",
        paste(miss06, collapse=", "))
}

# Save object inventories for reproducibility.
writeLines(ls(e05), file.path(OUT_DIR, "00_objects_env05.txt"))
writeLines(ls(e06), file.path(OUT_DIR, "00_objects_env06.txt"))

# ------------------------------------------------------------
# Intercept 05 companion_radius() to capture the exact Bdom list
# passed by fit_country(). The returned rho is kept identical.
# ------------------------------------------------------------
capture05 <- new.env(parent = emptyenv())
capture05$key <- NULL
capture05$Bdom <- list()

e05$companion_radius <- function(Bdom, k, p) {
  key <- capture05$key
  if(is.null(key)) stop("capture key not set", call. = FALSE)

  capture05$Bdom[[key]] <- lapply(Bdom, function(M) {
    M <- as.matrix(M)
    M
  })

  C <- companion_from_A(Bdom)
  safe_rho(C)
}

# ------------------------------------------------------------
# Compare cases
# ------------------------------------------------------------
summary_rows <- list()
eig_rows <- list()
matrix_rows <- list()
sample_rows <- list()
raw_rows <- list()
ret_rows <- list()

for(ii in seq_len(nrow(CASES))) {
  cc <- CASES$country[ii]
  p  <- CASES$p[ii]
  key <- paste(cc, p, sep="_p")

  msg("Reconciling %s p=%d ...", cc, p)

  # ---- 05: run fit_country with intercepted companion_radius
  capture05$key <- key
  ans05 <- e05$fit_country(
    cc = cc,
    p = p,
    panel = e05$panel,
    foreign = e05$foreign,
    globals = e05$globals
  )

  A05 <- capture05$Bdom[[key]]
  if(is.null(A05))
    stopf("Failed to capture Bdom from 05 for %s p=%d", cc, p)

  C05 <- companion_from_A(A05)

  # ---- 06: native fit
  ans06 <- e06$fit_local_varx(e06$panel, cc, p)
  A06 <- ans06$A
  C06 <- companion_from_A(A06)

  # ---- matrix dimensions and exact numerical comparison
  maxA <- NA_real_
  if(length(A05) == length(A06)) {
    diffs <- vapply(seq_along(A05),
                    function(L) safe_maxdiff(A05[[L]], A06[[L]]),
                    numeric(1))
    maxA <- max(diffs)
  } else {
    diffs <- rep(Inf, max(length(A05), length(A06)))
    maxA <- Inf
  }

  maxC <- safe_maxdiff(C05, C06)

  rho05 <- safe_rho(C05)
  rho06 <- safe_rho(C06)

  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    country = cc,
    p = p,
    rho_05 = rho05,
    rho_06 = rho06,
    rho_difference = rho05 - rho06,
    nlag_05 = length(A05),
    nlag_06 = length(A06),
    companion_nrow_05 = nrow(C05),
    companion_nrow_06 = nrow(C06),
    max_abs_A_difference = maxA,
    max_abs_companion_difference = maxC,
    A_equal_1e_10 = is.finite(maxA) && maxA < 1e-10,
    companion_equal_1e_10 = is.finite(maxC) && maxC < 1e-10,
    stringsAsFactors = FALSE
  )

  # Export each domestic lag matrix from both implementations.
  for(L in seq_along(A05)) {
    matrix_rows[[length(matrix_rows)+1L]] <-
      matrix_long(A05[[L]], "05_lag_selection", cc, p, "A", L)
  }
  for(L in seq_along(A06)) {
    matrix_rows[[length(matrix_rows)+1L]] <-
      matrix_long(A06[[L]], "06_instability", cc, p, "A", L)
  }

  matrix_rows[[length(matrix_rows)+1L]] <-
    matrix_long(C05, "05_lag_selection", cc, p, "companion")
  matrix_rows[[length(matrix_rows)+1L]] <-
    matrix_long(C06, "06_instability", cc, p, "companion")

  eig_rows[[length(eig_rows)+1L]] <-
    eigen_table(C05, "05_lag_selection", cc, p)
  eig_rows[[length(eig_rows)+1L]] <-
    eigen_table(C06, "06_instability", cc, p)

  # ---- compare base domestic data before regression
  q05 <- safe_qid(e05$panel[[cc]])
  q06 <- as.integer(e06$panel$qid)

  common_q <- intersect(q05, q06)
  sample_rows[[length(sample_rows)+1L]] <- data.frame(
    country = cc,
    p = p,
    source = c("05_lag_selection","06_instability"),
    n_base_quarters = c(length(q05), length(q06)),
    first_qid = c(if(length(q05)) min(q05) else NA,
                  if(length(q06)) min(q06) else NA),
    last_qid = c(if(length(q05)) max(q05) else NA,
                 if(length(q06)) max(q06) else NA),
    n_common_base_quarters = length(common_q),
    stringsAsFactors = FALSE
  )

  if(length(common_q)) {
    m05 <- match(common_q, q05)
    m06 <- match(common_q, q06)

    X05 <- as.matrix(e05$panel[[cc]][m05, e05$VARS, drop=FALSE])
    X06 <- as.matrix(e06$panel$X[m06, cc, e06$VARS, drop=FALSE][,1,])

    for(v in intersect(colnames(X05), colnames(X06))) {
      d <- X05[,v] - X06[,v]
      raw_rows[[length(raw_rows)+1L]] <- data.frame(
        country = cc,
        p = p,
        variable = v,
        n_common = sum(is.finite(X05[,v]) & is.finite(X06[,v])),
        max_abs_difference = if(any(is.finite(d))) max(abs(d), na.rm=TRUE) else NA,
        mean_abs_difference = if(any(is.finite(d))) mean(abs(d), na.rm=TRUE) else NA,
        stringsAsFactors = FALSE
      )
    }
  }

  # ---- record return-object structure, useful if estimation rows differ
  ret_rows[[length(ret_rows)+1L]] <- data.frame(
    country = cc,
    p = p,
    source = "05_lag_selection",
    return_names = paste(names(ans05), collapse=" | "),
    stringsAsFactors = FALSE
  )
  ret_rows[[length(ret_rows)+1L]] <- data.frame(
    country = cc,
    p = p,
    source = "06_instability",
    return_names = paste(names(ans06), collapse=" | "),
    stringsAsFactors = FALSE
  )

  # If 06 exposes regression row indices, save exact qids.
  if(!is.null(ans06$rows)) {
    rr <- as.integer(ans06$rows)
    qq <- e06$panel$qid[rr]
    write.csv(
      data.frame(country=cc, p=p, row_index=rr, qid=qq),
      file.path(OUT_DIR, sprintf("sample_06_%s_p%d.csv",cc,p)),
      row.names=FALSE
    )
  }

  # Save 05 returned diagnostics unchanged for inspection.
  if(!is.null(ans05$diag)) {
    write.csv(as.data.frame(ans05$diag),
              file.path(OUT_DIR, sprintf("returned_diag_05_%s_p%d.csv",cc,p)),
              row.names=FALSE)
  }
}

summary_df <- do.call(rbind, summary_rows)
eig_df <- do.call(rbind, eig_rows)
matrix_df <- do.call(rbind, matrix_rows)
sample_df <- do.call(rbind, sample_rows)
raw_df <- if(length(raw_rows)) do.call(rbind, raw_rows) else data.frame()
ret_df <- do.call(rbind, ret_rows)

write.csv(summary_df,
          file.path(OUT_DIR, "01_exact_reconciliation_summary.csv"),
          row.names=FALSE)

write.csv(sample_df,
          file.path(OUT_DIR, "02_base_sample_comparison.csv"),
          row.names=FALSE)

write.csv(raw_df,
          file.path(OUT_DIR, "03_domestic_data_comparison.csv"),
          row.names=FALSE)

write.csv(matrix_df,
          file.path(OUT_DIR, "04_A_and_companion_matrices_long.csv"),
          row.names=FALSE)

write.csv(eig_df,
          file.path(OUT_DIR, "05_eigenvalue_comparison.csv"),
          row.names=FALSE)

write.csv(ret_df,
          file.path(OUT_DIR, "06_return_object_structure.csv"),
          row.names=FALSE)

# Focused element-by-element A differences.
adiff <- list()
kk <- 1L
for(ii in seq_len(nrow(CASES))) {
  cc <- CASES$country[ii]
  p <- CASES$p[ii]
  for(L in seq_len(p)) {
    x <- subset(matrix_df,
                country==cc & p==p & object=="A" & lag==L &
                  source=="05_lag_selection")
    y <- subset(matrix_df,
                country==cc & p==p & object=="A" & lag==L &
                  source=="06_instability")

    keyx <- paste(x$row, x$col, sep="_")
    keyy <- paste(y$row, y$col, sep="_")
    common <- intersect(keyx, keyy)

    if(length(common)) {
      xx <- x[match(common,keyx),]
      yy <- y[match(common,keyy),]
      adiff[[kk]] <- data.frame(
        country=cc,
        p=p,
        lag=L,
        row=xx$row,
        col=xx$col,
        row_name_05=xx$row_name,
        col_name_05=xx$col_name,
        row_name_06=yy$row_name,
        col_name_06=yy$col_name,
        value_05=xx$value,
        value_06=yy$value,
        difference=xx$value-yy$value,
        abs_difference=abs(xx$value-yy$value),
        stringsAsFactors=FALSE
      )
      kk <- kk + 1L
    }
  }
}

if(length(adiff)) {
  adiff_df <- do.call(rbind, adiff)
  adiff_df <- adiff_df[order(adiff_df$country,
                             adiff_df$p,
                             -adiff_df$abs_difference),]
  write.csv(adiff_df,
            file.path(OUT_DIR, "07_A_elementwise_difference.csv"),
            row.names=FALSE)
}

# Automatic diagnosis.
diagnosis <- character()
diagnosis <- c(
  diagnosis,
  "Exact reconciliation of 05 lag-selection vs 06 instability diagnostic.",
  "",
  "Decision rules:"
)

for(i in seq_len(nrow(summary_df))) {
  z <- summary_df[i,]
  diagnosis <- c(
    diagnosis,
    sprintf(
      "%s p=%d: rho05=%.10f, rho06=%.10f, max|A05-A06|=%g, max|C05-C06|=%g",
      z$country, z$p, z$rho_05, z$rho_06,
      z$max_abs_A_difference, z$max_abs_companion_difference
    )
  )

  if(isTRUE(z$A_equal_1e_10) && isTRUE(z$companion_equal_1e_10)) {
    diagnosis <- c(
      diagnosis,
      "  -> A and companion matrices are numerically identical.",
      "     Any previously reported rho discrepancy comes from an older run/output,"
      ,"     not from the current implementations."
    )
  } else {
    diagnosis <- c(
      diagnosis,
      "  -> Current implementations construct different domestic/companion matrices.",
      "     Inspect 03_domestic_data_comparison.csv first.",
      "     If raw domestic data match, inspect 07_A_elementwise_difference.csv;",
      "     the discrepancy is then in regression sample/regressors/coefficient extraction."
    )
  }
}

diagnosis <- c(
  diagnosis,
  "",
  "Do not change lag order, weights, GPR, Brent, identification, or transformations",
  "until the two current implementations reconcile."
)

writeLines(diagnosis,
           file.path(OUT_DIR, "README_reconciliation.txt"))

msg("")
msg("=== Exact reconciliation summary ===")
print(summary_df, row.names=FALSE)
msg("")
msg("Outputs written to %s", OUT_DIR)
