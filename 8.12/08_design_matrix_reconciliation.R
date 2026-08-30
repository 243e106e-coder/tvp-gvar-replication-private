#!/usr/bin/env Rscript

# ============================================================
# Exact OLS design-matrix reconciliation
# 05_country_specific_lag_selection.R
# vs 06_cn_za_instability_diagnostic.R
#
# IMPORTANT:
# - Does NOT change the model.
# - Captures the ACTUAL X and y passed to stats::lm.fit().
# - Main focus: CN p=2; also ZA p=2 and JP p=1 as controls.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_05 <- "8.12/05_country_specific_lag_selection.R"
SCRIPT_06 <- "8.12/06_cn_za_instability_diagnostic.R"
OUT_DIR   <- "8.12/design_matrix_reconciliation"

CASES <- data.frame(
  country = c("CN", "ZA", "JP"),
  p       = c(2L, 2L, 1L),
  stringsAsFactors = FALSE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg   <- function(...) cat(sprintf(...), "\n")

if(!file.exists(SCRIPT_05)) stopf("Missing %s", SCRIPT_05)
if(!file.exists(SCRIPT_06)) stopf("Missing %s", SCRIPT_06)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

qid_label <- function(qid) {
  yr <- (qid - 1L) %/% 4L
  qq <- qid - 4L * yr
  sprintf("%dQ%d", yr, qq)
}

norm_col <- function(x) {
  x <- tolower(gsub("[^a-z0-9]+", "_", as.character(x)))
  gsub("^_+|_+$", "", x)
}

as_num_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

same_numeric <- function(a, b, tol = 1e-12) {
  if(length(a) != length(b)) return(FALSE)
  ok_na <- is.na(a) == is.na(b)
  if(!all(ok_na)) return(FALSE)
  ii <- which(!is.na(a) & !is.na(b))
  if(!length(ii)) return(TRUE)
  max(abs(a[ii] - b[ii])) <= tol
}

max_abs_diff <- function(a, b) {
  if(length(a) != length(b)) return(Inf)
  z <- abs(as.numeric(a) - as.numeric(b))
  if(!length(z) || all(is.na(z))) return(NA_real_)
  max(z, na.rm = TRUE)
}

# Find exact quarter rows by matching a captured dependent variable vector
# against the country domestic series. We use all 5 equations jointly where
# possible; this is robust to repeated values in one variable.
infer_rows_from_Y <- function(Ycap, Xdom, qid, vars) {
  Ycap <- as_num_matrix(Ycap)
  Xdom <- as_num_matrix(Xdom)

  if(ncol(Ycap) == 1L) {
    # One equation was captured. Try each domestic variable and every
    # contiguous window of matching length; retain all exact/near matches.
    n <- nrow(Ycap)
    candidates <- list()
    kk <- 1L
    for(j in seq_len(ncol(Xdom))) {
      if(nrow(Xdom) < n) next
      for(s in seq_len(nrow(Xdom) - n + 1L)) {
        idx <- s:(s+n-1L)
        a <- Xdom[idx, j]
        b <- Ycap[,1]
        good <- is.finite(a) & is.finite(b)
        if(sum(good) < max(5L, floor(0.8*n))) next
        d <- max(abs(a[good] - b[good]), na.rm=TRUE)
        if(is.finite(d) && d < 1e-10) {
          candidates[[kk]] <- list(var=colnames(Xdom)[j], idx=idx, diff=d)
          kk <- kk + 1L
        }
      }
    }
    if(length(candidates)) {
      z <- candidates[[1L]]
      return(list(rows=z$idx, qid=qid[z$idx], matched_by=z$var,
                  maxdiff=z$diff, ambiguous=length(candidates)>1L))
    }
    return(list(rows=integer(), qid=integer(), matched_by=NA_character_,
                maxdiff=NA_real_, ambiguous=TRUE))
  }

  list(rows=integer(), qid=integer(), matched_by=NA_character_,
       maxdiff=NA_real_, ambiguous=TRUE)
}

# ------------------------------------------------------------
# Capture actual stats::lm.fit(x, y) calls.
#
# trace() is used only during one targeted fit and removed immediately.
# ------------------------------------------------------------

capture_lmfit_calls <- function(expr, label) {
  cap <- new.env(parent = emptyenv())
  cap$calls <- list()

  assign(".DMR_CAPTURE_ENV", cap, envir = .GlobalEnv)

  tracer_expr <- quote({
    .ce <- get(".DMR_CAPTURE_ENV", envir = .GlobalEnv)
    .idx <- length(.ce$calls) + 1L

    .xx <- tryCatch(as.matrix(x), error=function(e) NULL)
    .yy <- tryCatch(as.matrix(y), error=function(e) {
      tryCatch(matrix(as.numeric(y), ncol=1L), error=function(e2) NULL)
    })

    .ce$calls[[.idx]] <- list(
      x = .xx,
      y = .yy,
      nrow_x = if(is.null(.xx)) NA_integer_ else nrow(.xx),
      ncol_x = if(is.null(.xx)) NA_integer_ else ncol(.xx),
      nrow_y = if(is.null(.yy)) NA_integer_ else nrow(.yy),
      ncol_y = if(is.null(.yy)) NA_integer_ else ncol(.yy),
      colnames_x = if(is.null(.xx)) character() else colnames(.xx)
    )
  })

  ns <- asNamespace("stats")
  trace("lm.fit", tracer = tracer_expr, print = FALSE, where = ns)
  on.exit({
    try(untrace("lm.fit", where = ns), silent = TRUE)
    if(exists(".DMR_CAPTURE_ENV", envir=.GlobalEnv, inherits=FALSE))
      rm(".DMR_CAPTURE_ENV", envir=.GlobalEnv)
  }, add = TRUE)

  result <- eval.parent(substitute(expr))

  untrace("lm.fit", where = ns)
  if(exists(".DMR_CAPTURE_ENV", envir=.GlobalEnv, inherits=FALSE))
    rm(".DMR_CAPTURE_ENV", envir=.GlobalEnv)

  list(result=result, calls=cap$calls, label=label)
}

select_main_design <- function(calls) {
  if(!length(calls)) return(NULL)

  good <- which(vapply(calls, function(z)
    !is.null(z$x) && !is.null(z$y) &&
      nrow(z$x) >= 20L && ncol(z$x) >= 5L &&
      nrow(z$x) == nrow(z$y),
    logical(1)))

  if(!length(good)) return(NULL)

  # Country VAR equations normally repeat the same X for 5 equations.
  # Prefer the largest-column design, then largest sample.
  score <- vapply(calls[good], function(z)
    100000 * ncol(z$x) + nrow(z$x), numeric(1))

  calls[[good[which.max(score)]]]
}

# ------------------------------------------------------------
# Source current scripts in isolated environments
# ------------------------------------------------------------

msg("Sourcing 05 lag-selection script...")
e05 <- new.env(parent = globalenv())
sys.source(SCRIPT_05, envir = e05)

msg("Sourcing 06 instability script...")
e06 <- new.env(parent = globalenv())
sys.source(SCRIPT_06, envir = e06)

need05 <- c("fit_country","panel","foreign","globals","VARS")
need06 <- c("fit_local_varx","panel","VARS")

m05 <- setdiff(need05, ls(e05))
m06 <- setdiff(need06, ls(e06))
if(length(m05)) stopf("05 missing objects: %s", paste(m05, collapse=", "))
if(length(m06)) stopf("06 missing objects: %s", paste(m06, collapse=", "))

# ------------------------------------------------------------
# Extract domestic base series and qids from each implementation
# ------------------------------------------------------------

get_dom05 <- function(cc) {
  z <- e05$panel[[cc]]
  if(is.null(z)) stopf("05 panel[[%s]] missing", cc)
  q <- as.integer(z$qid)
  X <- as.matrix(z[, e05$VARS, drop=FALSE])
  colnames(X) <- e05$VARS
  list(qid=q, X=X)
}

get_dom06 <- function(cc) {
  i <- match(cc, e06$COUNTRIES)
  if(is.na(i)) stopf("06 country %s missing", cc)
  q <- as.integer(e06$panel$qid)
  X <- e06$panel$X[, i, , drop=FALSE][,1,]
  X <- as.matrix(X)
  colnames(X) <- e06$VARS
  list(qid=q, X=X)
}

# ------------------------------------------------------------
# Reconstruct exact quarter labels for captured designs.
# 06 exposes rows directly. For 05, infer rows from captured y.
# ------------------------------------------------------------

align05_qid <- function(main05, cc) {
  dom <- get_dom05(cc)
  inf <- infer_rows_from_Y(main05$y, dom$X, dom$qid, e05$VARS)
  if(length(inf$qid) == nrow(main05$x)) return(inf)

  # Fallback: fit_country uses common quarterly sample and p/q lags.
  # If no NA trimming remains, the last n observations are the estimation sample.
  n <- nrow(main05$x)
  if(n <= length(dom$qid)) {
    idx <- (length(dom$qid)-n+1L):length(dom$qid)
    return(list(rows=idx, qid=dom$qid[idx],
                matched_by="fallback_last_n", maxdiff=NA_real_,
                ambiguous=TRUE))
  }

  list(rows=integer(), qid=integer(), matched_by="unresolved",
       maxdiff=NA_real_, ambiguous=TRUE)
}

# ------------------------------------------------------------
# Column matching
# ------------------------------------------------------------

canonical_name <- function(x) {
  z <- norm_col(x)

  # normalize common naming differences
  z <- gsub("intercept|const|constant", "const", z)
  z <- gsub("gpr_l0|gpr_0|gpr_current", "gpr_0", z)
  z <- gsub("oil_l0|oil_0|oil_current|brent_0|brent_l0", "oil_0", z)
  z <- gsub("brent", "oil", z)

  # foreign naming
  z <- gsub("star_?([a-z0-9]+)_l([0-9]+)", "\\1_star_l\\2", z)
  z <- gsub("([a-z0-9]+)_foreign_l([0-9]+)", "\\1_star_l\\2", z)
  z <- gsub("([a-z0-9]+)_ast_l([0-9]+)", "\\1_star_l\\2", z)

  z
}

compare_designs <- function(X05, q05, X06, q06, cc, p) {
  X05 <- as_num_matrix(X05)
  X06 <- as_num_matrix(X06)

  c05 <- colnames(X05)
  c06 <- colnames(X06)
  if(is.null(c05)) c05 <- paste0("V",seq_len(ncol(X05)))
  if(is.null(c06)) c06 <- paste0("V",seq_len(ncol(X06)))

  can05 <- canonical_name(c05)
  can06 <- canonical_name(c06)

  common_q <- intersect(q05, q06)
  out <- list()
  kk <- 1L

  # Compare matched columns by canonical names.
  common_c <- intersect(can05, can06)

  for(cn in common_c) {
    j05 <- which(can05 == cn)[1L]
    j06 <- which(can06 == cn)[1L]

    m05 <- match(common_q, q05)
    m06 <- match(common_q, q06)

    a <- X05[m05, j05]
    b <- X06[m06, j06]
    d <- a - b

    finite <- is.finite(a) & is.finite(b)

    first_bad <- NA_integer_
    if(any(finite & abs(d) > 1e-10))
      first_bad <- which(finite & abs(d) > 1e-10)[1L]

    out[[kk]] <- data.frame(
      country=cc,
      p=p,
      canonical_column=cn,
      column_05=c05[j05],
      column_06=c06[j06],
      n_common_quarters=length(common_q),
      n_finite_pairs=sum(finite),
      max_abs_difference=if(any(finite)) max(abs(d[finite]),na.rm=TRUE) else NA_real_,
      mean_abs_difference=if(any(finite)) mean(abs(d[finite]),na.rm=TRUE) else NA_real_,
      n_different_1e10=sum(finite & abs(d)>1e-10),
      first_different_qid=if(is.na(first_bad)) NA_integer_ else common_q[first_bad],
      first_different_quarter=if(is.na(first_bad)) NA_character_ else qid_label(common_q[first_bad]),
      stringsAsFactors=FALSE
    )
    kk <- kk + 1L
  }

  cmp <- if(length(out)) do.call(rbind,out) else data.frame()

  inventory <- rbind(
    data.frame(country=cc,p=p,source="05_lag_selection",
               position=seq_along(c05),column=c05,
               canonical=can05,stringsAsFactors=FALSE),
    data.frame(country=cc,p=p,source="06_instability",
               position=seq_along(c06),column=c06,
               canonical=can06,stringsAsFactors=FALSE)
  )

  missing <- rbind(
    data.frame(country=cc,p=p,
               present_in="05_only",
               canonical=setdiff(can05,can06),
               stringsAsFactors=FALSE),
    data.frame(country=cc,p=p,
               present_in="06_only",
               canonical=setdiff(can06,can05),
               stringsAsFactors=FALSE)
  )

  list(compare=cmp, inventory=inventory, missing=missing, common_q=common_q)
}

# ------------------------------------------------------------
# Run cases
# ------------------------------------------------------------

summary_rows <- list()
call_rows <- list()
colcmp_rows <- list()
inventory_rows <- list()
missing_rows <- list()
rowdetail_rows <- list()

for(ii in seq_len(nrow(CASES))) {
  cc <- CASES$country[ii]
  p  <- CASES$p[ii]

  msg("Capturing actual OLS design: %s p=%d", cc, p)

  cap05 <- capture_lmfit_calls(
    e05$fit_country(cc, p, e05$panel, e05$foreign, e05$globals),
    sprintf("05_%s_p%d",cc,p)
  )

  cap06 <- capture_lmfit_calls(
    e06$fit_local_varx(e06$panel, cc, p),
    sprintf("06_%s_p%d",cc,p)
  )

  main05 <- select_main_design(cap05$calls)
  main06 <- select_main_design(cap06$calls)

  if(is.null(main05)) stopf("Could not capture main lm.fit design from 05 for %s p=%d",cc,p)
  if(is.null(main06)) stopf("Could not capture main lm.fit design from 06 for %s p=%d",cc,p)

  a05 <- align05_qid(main05,cc)

  # 06 gives exact used rows in return object.
  ans06 <- cap06$result
  if(!is.null(ans06$rows)) {
    q06 <- e06$panel$qid[as.integer(ans06$rows)]
  } else {
    n06 <- nrow(main06$x)
    q06 <- tail(e06$panel$qid,n06)
  }

  q05 <- a05$qid

  # If capture selected one of 5 equation calls, X is still the exact design.
  cmp <- compare_designs(main05$x,q05,main06$x,q06,cc,p)

  if(nrow(cmp$compare)) colcmp_rows[[length(colcmp_rows)+1L]] <- cmp$compare
  inventory_rows[[length(inventory_rows)+1L]] <- cmp$inventory
  if(nrow(cmp$missing)) missing_rows[[length(missing_rows)+1L]] <- cmp$missing

  maxdiff <- if(nrow(cmp$compare))
    max(cmp$compare$max_abs_difference,na.rm=TRUE) else NA_real_
  if(!is.finite(maxdiff)) maxdiff <- NA_real_

  first_bad <- NA_character_
  if(nrow(cmp$compare)) {
    bad <- which(is.finite(cmp$compare$max_abs_difference) &
                   cmp$compare$max_abs_difference > 1e-10)
    if(length(bad)) {
      z <- cmp$compare[bad[order(cmp$compare$first_different_qid[bad],
                                na.last=TRUE)][1L],]
      first_bad <- sprintf("%s @ %s",
                           z$canonical_column,
                           z$first_different_quarter)
    }
  }

  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    country=cc,
    p=p,
    nrow_X05=nrow(main05$x),
    ncol_X05=ncol(main05$x),
    nrow_X06=nrow(main06$x),
    ncol_X06=ncol(main06$x),
    first_q05=if(length(q05)) qid_label(min(q05)) else NA_character_,
    last_q05=if(length(q05)) qid_label(max(q05)) else NA_character_,
    first_q06=if(length(q06)) qid_label(min(q06)) else NA_character_,
    last_q06=if(length(q06)) qid_label(max(q06)) else NA_character_,
    qid_alignment_method_05=a05$matched_by,
    qid_alignment_ambiguous_05=a05$ambiguous,
    n_common_columns=length(intersect(canonical_name(colnames(main05$x)),
                                      canonical_name(colnames(main06$x)))),
    n_05_only=length(setdiff(canonical_name(colnames(main05$x)),
                            canonical_name(colnames(main06$x)))),
    n_06_only=length(setdiff(canonical_name(colnames(main06$x)),
                            canonical_name(colnames(main05$x)))),
    max_abs_design_difference=maxdiff,
    first_detected_difference=first_bad,
    stringsAsFactors=FALSE
  )

  # Record all lm.fit calls to make sure we captured the intended one.
  for(src in c("05","06")) {
    calls <- if(src=="05") cap05$calls else cap06$calls
    if(length(calls)) for(k in seq_along(calls)) {
      z <- calls[[k]]
      call_rows[[length(call_rows)+1L]] <- data.frame(
        country=cc,p=p,source=src,call_index=k,
        nrow_x=z$nrow_x,ncol_x=z$ncol_x,
        nrow_y=z$nrow_y,ncol_y=z$ncol_y,
        x_columns=paste(z$colnames_x,collapse=" | "),
        stringsAsFactors=FALSE
      )
    }
  }

  # Full row-by-row matched design differences.
  common_q <- cmp$common_q
  common_c <- intersect(canonical_name(colnames(main05$x)),
                        canonical_name(colnames(main06$x)))

  if(length(common_q) && length(common_c)) {
    m05q <- match(common_q,q05)
    m06q <- match(common_q,q06)
    can05 <- canonical_name(colnames(main05$x))
    can06 <- canonical_name(colnames(main06$x))

    for(cn in common_c) {
      j05 <- which(can05==cn)[1L]
      j06 <- which(can06==cn)[1L]
      a <- main05$x[m05q,j05]
      b <- main06$x[m06q,j06]

      rowdetail_rows[[length(rowdetail_rows)+1L]] <- data.frame(
        country=cc,p=p,qid=common_q,quarter=qid_label(common_q),
        canonical_column=cn,
        value_05=a,value_06=b,
        difference=a-b,
        abs_difference=abs(a-b),
        stringsAsFactors=FALSE
      )
    }
  }

  # Save actual matrices separately for exact inspection.
  write.csv(
    data.frame(qid=q05,quarter=qid_label(q05),main05$x,
               check.names=FALSE),
    file.path(OUT_DIR,sprintf("X_actual_05_%s_p%d.csv",cc,p)),
    row.names=FALSE
  )

  write.csv(
    data.frame(qid=q06,quarter=qid_label(q06),main06$x,
               check.names=FALSE),
    file.path(OUT_DIR,sprintf("X_actual_06_%s_p%d.csv",cc,p)),
    row.names=FALSE
  )
}

summary_df <- do.call(rbind,summary_rows)
calls_df <- if(length(call_rows)) do.call(rbind,call_rows) else data.frame()
colcmp_df <- if(length(colcmp_rows)) do.call(rbind,colcmp_rows) else data.frame()
inventory_df <- do.call(rbind,inventory_rows)
missing_df <- if(length(missing_rows)) do.call(rbind,missing_rows) else data.frame()
rowdetail_df <- if(length(rowdetail_rows)) do.call(rbind,rowdetail_rows) else data.frame()

if(nrow(colcmp_df)) {
  colcmp_df <- colcmp_df[
    order(colcmp_df$country,colcmp_df$p,
          -ifelse(is.na(colcmp_df$max_abs_difference),-Inf,
                  colcmp_df$max_abs_difference)),
    ,drop=FALSE]
}

if(nrow(rowdetail_df)) {
  rowdetail_df <- rowdetail_df[
    order(rowdetail_df$country,rowdetail_df$p,
          -ifelse(is.na(rowdetail_df$abs_difference),-Inf,
                  rowdetail_df$abs_difference)),
    ,drop=FALSE]
}

write.csv(summary_df,
          file.path(OUT_DIR,"01_design_matrix_summary.csv"),
          row.names=FALSE)

write.csv(calls_df,
          file.path(OUT_DIR,"02_lmfit_call_inventory.csv"),
          row.names=FALSE)

write.csv(inventory_df,
          file.path(OUT_DIR,"03_design_column_inventory.csv"),
          row.names=FALSE)

write.csv(missing_df,
          file.path(OUT_DIR,"04_unmatched_columns.csv"),
          row.names=FALSE)

write.csv(colcmp_df,
          file.path(OUT_DIR,"05_columnwise_design_comparison.csv"),
          row.names=FALSE)

write.csv(rowdetail_df,
          file.path(OUT_DIR,"06_rowwise_design_difference.csv"),
          row.names=FALSE)

# Top differences only.
topdiff <- if(nrow(rowdetail_df)) {
  subset(rowdetail_df,is.finite(abs_difference) & abs_difference > 1e-10)
} else data.frame()

if(nrow(topdiff)) {
  topdiff <- head(topdiff,200L)
}
write.csv(topdiff,
          file.path(OUT_DIR,"07_top_200_design_differences.csv"),
          row.names=FALSE)

# ------------------------------------------------------------
# Automatic interpretation
# ------------------------------------------------------------

txt <- c(
  "DESIGN-MATRIX RECONCILIATION",
  "============================",
  "",
  "This diagnostic captures the ACTUAL X matrices passed to stats::lm.fit().",
  "It does not alter lags, weights, GPR, Brent, identification, or transformations.",
  ""
)

for(i in seq_len(nrow(summary_df))) {
  z <- summary_df[i,]
  txt <- c(
    txt,
    sprintf("%s p=%d",z$country,z$p),
    sprintf("  05 X: %d x %d; sample %s to %s",
            z$nrow_X05,z$ncol_X05,z$first_q05,z$last_q05),
    sprintf("  06 X: %d x %d; sample %s to %s",
            z$nrow_X06,z$ncol_X06,z$first_q06,z$last_q06),
    sprintf("  common columns: %d; 05-only: %d; 06-only: %d",
            z$n_common_columns,z$n_05_only,z$n_06_only),
    sprintf("  max absolute design difference: %s",
            format(z$max_abs_design_difference,digits=12)),
    sprintf("  first detected difference: %s",
            ifelse(is.na(z$first_detected_difference),
                   "NONE",z$first_detected_difference))
  )

  zz <- subset(colcmp_df,
               country==z$country & p==z$p &
                 is.finite(max_abs_difference) &
                 max_abs_difference > 1e-10)

  if(nrow(zz)) {
    zz <- zz[order(-zz$max_abs_difference),,drop=FALSE]
    take <- head(zz,8L)
    txt <- c(txt,"  largest differing columns:")
    for(j in seq_len(nrow(take))) {
      txt <- c(txt,
        sprintf("    %s: max diff=%g; first=%s",
                take$canonical_column[j],
                take$max_abs_difference[j],
                take$first_different_quarter[j]))
    }
  } else {
    txt <- c(txt,"  No matched design column differs above 1e-10.")
  }
  txt <- c(txt,"")
}

txt <- c(
  txt,
  "Interpretation order:",
  "1. If sample endpoints differ -> row trimming / complete-case alignment problem.",
  "2. If 04_unmatched_columns.csv is non-empty -> regressor specification/naming differs.",
  "3. If domestic lag columns match but star columns differ -> foreign-variable construction differs.",
  "4. If star columns match but GPR/Oil differs -> global-series alignment differs.",
  "5. If all X columns match -> inspect dependent-variable alignment / coefficient extraction.",
  "",
  "Do not choose the final country-specific lag vector until this reconciliation is resolved."
)

writeLines(txt,file.path(OUT_DIR,"README_design_matrix_reconciliation.txt"))

msg("")
msg("=== DESIGN MATRIX SUMMARY ===")
print(summary_df,row.names=FALSE)
msg("")
msg("Outputs: %s",OUT_DIR)
