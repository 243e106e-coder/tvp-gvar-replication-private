#!/usr/bin/env Rscript

# ============================================================
# 08_design_matrix_reconciliation.R  v2.0
#
# Direct design-matrix reconciliation:
#   05_country_specific_lag_selection.R
#   vs
#   06_cn_za_instability_diagnostic.R
#
# NO MODEL CHANGES.
# This script reproduces the design matrices exactly from the
# two currently implemented fit functions and compares them.
#
# Focus:
#   CN p=2
#   ZA p=2
#   JP p=1
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_05 <- "8.12/05_country_specific_lag_selection.R"
SCRIPT_06 <- "8.12/06_cn_za_instability_diagnostic.R"
OUT_DIR   <- "8.12/design_matrix_reconciliation"

CASES <- data.frame(
  country = c("CN","ZA","JP"),
  p = c(2L,2L,1L),
  stringsAsFactors = FALSE
)

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...), call.=FALSE)
msg <- function(...) cat(sprintf(...), "\n")

if(!file.exists(SCRIPT_05)) stopf("Missing %s", SCRIPT_05)
if(!file.exists(SCRIPT_06)) stopf("Missing %s", SCRIPT_06)

qid_label <- function(qid) {
  yr <- (qid-1L)%/%4L
  qq <- qid - 4L*yr
  sprintf("%dQ%d", yr, qq)
}

lag_matrix_local <- function(X,L) {
  X <- as.matrix(X)
  z <- matrix(NA_real_,nrow(X),ncol(X))
  if(L<nrow(X))
    z[(L+1L):nrow(X),] <- X[1L:(nrow(X)-L),,drop=FALSE]
  colnames(z) <- colnames(X)
  z
}

lag_vec_local <- function(x,L) {
  out <- rep(NA_real_,length(x))
  if(L<length(x))
    out[(L+1L):length(x)] <- x[1L:(length(x)-L)]
  out
}

canon <- function(x) {
  z <- tolower(gsub("[^a-z0-9]+","_",x))
  z <- gsub("^_+|_+$","",z)
  z <- gsub("^xstar_","",z)
  z
}

# ------------------------------------------------------------
# Source current scripts in isolated environments.
# Their main sections will execute and recreate their outputs;
# this diagnostic does not alter their model specifications.
# ------------------------------------------------------------

msg("Sourcing 05...")
e05 <- new.env(parent=globalenv())
sys.source(SCRIPT_05,envir=e05)

msg("Sourcing 06...")
e06 <- new.env(parent=globalenv())
sys.source(SCRIPT_06,envir=e06)

need05 <- c("panel","foreign","globals","VARS","P_CANDIDATES")
need06 <- c("panel","VARS","Q_FIXED")
miss05 <- setdiff(need05,ls(e05))
miss06 <- setdiff(need06,ls(e06))

if(length(miss05)) stopf("05 missing objects: %s",paste(miss05,collapse=", "))
if(length(miss06)) stopf("06 missing objects: %s",paste(miss06,collapse=", "))

# ------------------------------------------------------------
# EXACT reproduction of design used by 05 fit_country()
# ------------------------------------------------------------

build_05 <- function(cc,p) {
  Xi  <- as.matrix(e05$panel[[cc]][,e05$VARS,drop=FALSE])
  Xs  <- as.matrix(e05$foreign[[cc]][,e05$VARS,drop=FALSE])
  qid <- as.integer(e05$panel[[cc]]$qid)
  n   <- nrow(Xi)

  # IMPORTANT: this is exactly what current 05 uses:
  # idx starts after max(P_CANDIDATES), not after current p.
  idx <- seq.int(max(e05$P_CANDIDATES)+1L,n)

  Y <- Xi[idx,,drop=FALSE]
  colnames(Y) <- e05$VARS

  Z <- matrix(1,length(idx),1)
  colnames(Z) <- "const"

  # domestic lags
  for(L in seq_len(p)) {
    XL <- lag_matrix_local(Xi,L)[idx,,drop=FALSE]
    colnames(XL) <- paste0(e05$VARS,"_L",L)
    Z <- cbind(Z,XL)
  }

  # IMPORTANT:
  # current 05 includes contemporaneous foreign-star AND L1.
  Xs0 <- Xs[idx,,drop=FALSE]
  colnames(Xs0) <- paste0(e05$VARS,"_star_0")

  Xs1 <- lag_matrix_local(Xs,1L)[idx,,drop=FALSE]
  colnames(Xs1) <- paste0(e05$VARS,"_star_L1")

  Z <- cbind(Z,Xs0,Xs1)

  # globals current + L1
  G <- as.matrix(
    e05$globals[
      match(qid,e05$globals$qid),
      c("gpr","oil"),
      drop=FALSE
    ]
  )

  G0 <- G[idx,,drop=FALSE]
  colnames(G0) <- c("gpr_0","oil_0")

  G1 <- lag_matrix_local(G,1L)[idx,,drop=FALSE]
  colnames(G1) <- c("gpr_L1","oil_L1")

  Z <- cbind(Z,G0,G1)

  keep <- complete.cases(Y) & complete.cases(Z)

  list(
    X=Z[keep,,drop=FALSE],
    Y=Y[keep,,drop=FALSE],
    qid=qid[idx][keep],
    raw_idx=idx,
    keep=keep
  )
}

# ------------------------------------------------------------
# EXACT reproduction of design used by 06 fit_local_varx()
# ------------------------------------------------------------

build_06 <- function(cc,p) {
  i <- match(cc,e06$COUNTRIES)
  if(is.na(i)) stopf("06 country missing: %s",cc)

  Y <- e06$panel$X[,i,,drop=FALSE][,1,]
  Xs <- e06$panel$Xstar[,i,,drop=FALSE][,1,]

  Y <- as.matrix(Y)
  Xs <- as.matrix(Xs)

  colnames(Y) <- e06$VARS
  colnames(Xs) <- e06$VARS

  Tn <- nrow(Y)
  maxlag <- max(p,e06$Q_FIXED,1L)
  rows <- (maxlag+1L):Tn

  D <- data.frame(const=rep(1,length(rows)))

  # domestic lags
  for(L in seq_len(p)) {
    for(v in seq_along(e06$VARS)) {
      D[[paste0(e06$VARS[v],"_L",L)]] <-
        lag_vec_local(Y[,v],L)[rows]
    }
  }

  # IMPORTANT:
  # current 06 includes ONLY foreign-star L1.
  for(v in seq_along(e06$VARS)) {
    D[[paste0(e06$VARS[v],"_star_L1")]] <-
      lag_vec_local(Xs[,v],1L)[rows]
  }

  # globals current + L1
  D$gpr_0  <- e06$panel$gpr[rows]
  D$gpr_L1 <- lag_vec_local(e06$panel$gpr,1L)[rows]
  D$oil_0  <- e06$panel$oil[rows]
  D$oil_L1 <- lag_vec_local(e06$panel$oil,1L)[rows]

  YY <- Y[rows,,drop=FALSE]

  ok <- complete.cases(D) & complete.cases(YY)

  list(
    X=as.matrix(D[ok,,drop=FALSE]),
    Y=YY[ok,,drop=FALSE],
    qid=as.integer(e06$panel$qid[rows][ok]),
    raw_idx=rows,
    keep=ok
  )
}

# ------------------------------------------------------------
# Comparison helpers
# ------------------------------------------------------------

column_inventory <- function(a,b,cc,p) {
  A <- data.frame(
    country=cc,p=p,source="05_lag_selection",
    position=seq_len(ncol(a$X)),
    column=colnames(a$X),
    stringsAsFactors=FALSE
  )
  B <- data.frame(
    country=cc,p=p,source="06_instability",
    position=seq_len(ncol(b$X)),
    column=colnames(b$X),
    stringsAsFactors=FALSE
  )
  rbind(A,B)
}

compare_rows <- function(a,b,cc,p) {
  uq <- sort(unique(c(a$qid,b$qid)))
  data.frame(
    country=cc,p=p,
    qid=uq,
    quarter=qid_label(uq),
    in_05=uq %in% a$qid,
    in_06=uq %in% b$qid,
    stringsAsFactors=FALSE
  )
}

compare_common_columns <- function(a,b,cc,p) {
  ca <- canon(colnames(a$X))
  cb <- canon(colnames(b$X))

  common <- intersect(ca,cb)
  cq <- intersect(a$qid,b$qid)

  out <- list()
  k <- 1L

  for(nm in common) {
    ia <- which(ca==nm)[1L]
    ib <- which(cb==nm)[1L]
    ra <- match(cq,a$qid)
    rb <- match(cq,b$qid)

    va <- a$X[ra,ia]
    vb <- b$X[rb,ib]
    d <- va-vb
    good <- is.finite(va)&is.finite(vb)

    bad <- which(good & abs(d)>1e-10)

    out[[k]] <- data.frame(
      country=cc,p=p,
      canonical_column=nm,
      column_05=colnames(a$X)[ia],
      column_06=colnames(b$X)[ib],
      n_common_quarters=length(cq),
      max_abs_difference=if(any(good)) max(abs(d[good])) else NA_real_,
      mean_abs_difference=if(any(good)) mean(abs(d[good])) else NA_real_,
      n_different_1e10=length(bad),
      first_different_quarter=if(length(bad)) qid_label(cq[bad[1]]) else NA_character_,
      stringsAsFactors=FALSE
    )
    k <- k+1L
  }

  if(length(out)) do.call(rbind,out) else data.frame()
}

compare_Y <- function(a,b,cc,p) {
  cq <- intersect(a$qid,b$qid)
  if(!length(cq)) return(data.frame())

  ra <- match(cq,a$qid)
  rb <- match(cq,b$qid)

  out <- list()
  for(j in seq_along(e05$VARS)) {
    va <- a$Y[ra,j]
    vb <- b$Y[rb,j]
    d <- va-vb
    good <- is.finite(va)&is.finite(vb)

    out[[j]] <- data.frame(
      country=cc,p=p,
      equation=e05$VARS[j],
      n_common_quarters=length(cq),
      max_abs_difference=if(any(good)) max(abs(d[good])) else NA_real_,
      n_different_1e10=sum(good & abs(d)>1e-10),
      stringsAsFactors=FALSE
    )
  }
  do.call(rbind,out)
}

# ------------------------------------------------------------
# Main reconciliation
# ------------------------------------------------------------

summaries <- list()
inventories <- list()
rowcmp <- list()
colcmp <- list()
ycmp <- list()
unmatched <- list()

for(ii in seq_len(nrow(CASES))) {
  cc <- CASES$country[ii]
  pp <- CASES$p[ii]

  msg("Building direct designs: %s p=%d",cc,pp)

  a <- build_05(cc,pp)
  b <- build_06(cc,pp)

  write.csv(
    data.frame(qid=a$qid,quarter=qid_label(a$qid),a$X,
               check.names=FALSE),
    file.path(OUT_DIR,sprintf("X_05_%s_p%d.csv",cc,pp)),
    row.names=FALSE
  )

  write.csv(
    data.frame(qid=b$qid,quarter=qid_label(b$qid),b$X,
               check.names=FALSE),
    file.path(OUT_DIR,sprintf("X_06_%s_p%d.csv",cc,pp)),
    row.names=FALSE
  )

  ca <- canon(colnames(a$X))
  cb <- canon(colnames(b$X))
  only05 <- setdiff(ca,cb)
  only06 <- setdiff(cb,ca)

  inventories[[length(inventories)+1L]] <- column_inventory(a,b,cc,pp)
  rowcmp[[length(rowcmp)+1L]] <- compare_rows(a,b,cc,pp)
  ccx <- compare_common_columns(a,b,cc,pp)
  if(nrow(ccx)) colcmp[[length(colcmp)+1L]] <- ccx
  yy <- compare_Y(a,b,cc,pp)
  if(nrow(yy)) ycmp[[length(ycmp)+1L]] <- yy

  if(length(only05)) {
    unmatched[[length(unmatched)+1L]] <- data.frame(
      country=cc,p=pp,present_in="05_only",
      canonical_column=only05,stringsAsFactors=FALSE
    )
  }
  if(length(only06)) {
    unmatched[[length(unmatched)+1L]] <- data.frame(
      country=cc,p=pp,present_in="06_only",
      canonical_column=only06,stringsAsFactors=FALSE
    )
  }

  max_x_diff <- if(nrow(ccx) && any(is.finite(ccx$max_abs_difference)))
    max(ccx$max_abs_difference,na.rm=TRUE) else NA_real_

  summaries[[length(summaries)+1L]] <- data.frame(
    country=cc,p=pp,
    nobs_05=nrow(a$X),
    nobs_06=nrow(b$X),
    regressors_05=ncol(a$X),
    regressors_06=ncol(b$X),
    first_q_05=qid_label(min(a$qid)),
    last_q_05=qid_label(max(a$qid)),
    first_q_06=qid_label(min(b$qid)),
    last_q_06=qid_label(max(b$qid)),
    identical_qid=identical(as.integer(a$qid),as.integer(b$qid)),
    n_common_columns=length(intersect(ca,cb)),
    n_05_only=length(only05),
    n_06_only=length(only06),
    max_abs_difference_common_X=max_x_diff,
    stringsAsFactors=FALSE
  )
}

summary_df <- do.call(rbind,summaries)
inventory_df <- do.call(rbind,inventories)
rows_df <- do.call(rbind,rowcmp)
cols_df <- if(length(colcmp)) do.call(rbind,colcmp) else data.frame()
y_df <- if(length(ycmp)) do.call(rbind,ycmp) else data.frame()
unmatched_df <- if(length(unmatched)) do.call(rbind,unmatched) else data.frame()

write.csv(summary_df,
          file.path(OUT_DIR,"01_design_matrix_summary.csv"),
          row.names=FALSE)
write.csv(inventory_df,
          file.path(OUT_DIR,"02_design_column_inventory.csv"),
          row.names=FALSE)
write.csv(rows_df,
          file.path(OUT_DIR,"03_estimation_row_comparison.csv"),
          row.names=FALSE)
write.csv(unmatched_df,
          file.path(OUT_DIR,"04_unmatched_regressors.csv"),
          row.names=FALSE)
write.csv(cols_df,
          file.path(OUT_DIR,"05_common_column_value_comparison.csv"),
          row.names=FALSE)
write.csv(y_df,
          file.path(OUT_DIR,"06_dependent_variable_comparison.csv"),
          row.names=FALSE)

# ------------------------------------------------------------
# Explicit structural-specification audit
# ------------------------------------------------------------

spec <- data.frame(
  item=c(
    "domestic_lags",
    "foreign_star_contemporaneous",
    "foreign_star_lag1",
    "gpr_current",
    "gpr_lag1",
    "oil_current",
    "oil_lag1",
    "row_start_rule"
  ),
  implementation_05=c(
    "L1..Lp",
    "YES",
    "YES",
    "YES",
    "YES",
    "YES",
    "YES",
    "max(P_CANDIDATES)+1"
  ),
  implementation_06=c(
    "L1..Lp",
    "NO",
    "YES",
    "YES",
    "YES",
    "YES",
    "YES",
    "max(p,Q_FIXED,1)+1"
  ),
  same=c(
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),
  stringsAsFactors=FALSE
)

write.csv(spec,
          file.path(OUT_DIR,"07_structural_specification_difference.csv"),
          row.names=FALSE)

# ------------------------------------------------------------
# README / automatic conclusion
# ------------------------------------------------------------

readme <- c(
  "DESIGN-MATRIX RECONCILIATION v2",
  "================================",
  "",
  "This script directly reproduces the implemented design matrices of",
  "05_country_specific_lag_selection.R and",
  "06_cn_za_instability_diagnostic.R.",
  "",
  "No model setting is changed by this diagnostic.",
  "",
  "KNOWN IMPLEMENTATION DIFFERENCE FOUND BEFORE ESTIMATION:",
  "",
  "05 fit_country():",
  "  Z includes Xs[t] AND Xs[t-1].",
  "",
  "06 fit_local_varx():",
  "  D includes ONLY Xs[t-1].",
  "",
  "Therefore the two local systems are NOT the same VARX specification.",
  "Their domestic-lag coefficients and companion spectral radii are not",
  "expected to match exactly, even when the domestic raw data are identical.",
  "",
  "A second implementation difference exists for row initialization:",
  "  05 starts at max(P_CANDIDATES)+1 = 3 for BOTH p=1 and p=2.",
  "  06 starts at max(p,Q_FIXED,1)+1.",
  "Thus for p=1, 05 starts one quarter later than 06.",
  "For p=2, this row-start rule is the same.",
  "",
  "Interpretation:",
  "  - CN p=2 and ZA p=2: the main specification difference is",
  "    contemporaneous foreign-star variables included by 05 but omitted by 06.",
  "  - JP p=1: both the foreign-star specification and row-start rule differ.",
  "",
  "Do NOT use the earlier 05-vs-06 spectral-radius discrepancy as evidence",
  "of a data error. The estimators were not estimating the same local model.",
  "",
  "Next methodological decision should be explicit:",
  "choose the intended GVAR local specification for foreign variables, then",
  "make BOTH lag-selection and instability diagnostics use that same design.",
  ""
)

writeLines(readme,
           file.path(OUT_DIR,"README_design_matrix_reconciliation.txt"))

msg("")
msg("=== SUMMARY ===")
print(summary_df,row.names=FALSE)
msg("")
msg("=== STRUCTURAL DIFFERENCES ===")
print(spec,row.names=FALSE)
msg("")
msg("DONE. Outputs: %s",OUT_DIR)
