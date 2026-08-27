#!/usr/bin/env Rscript

# ============================================================
# Country-specific lag selection for GPR + Brent GVAR / TVP-GVAR
# Clean rebuild
#
# Compares p_i = 1 vs 2 country by country, holding q_i = 1.
# Outputs:
#   01_country_p1_p2_model_diagnostics.csv
#   02_equation_residual_serial_correlation.csv
#   03_country_lag_recommendation.csv
#   04_selected_country_lag_vector.R
#   README_country_specific_lag_selection.txt
#
# NOTE:
# This is a pre-estimation lag/specification diagnostic.
# It does not replace posterior TVP-GVAR stability diagnostics.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# --------------------------- CONFIG ---------------------------

MACRO_PATH  <- "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"
MACRO_SHEET <- 1

WEIGHT_PATH <- "8.12/Trade_Weights_14_Economies_2000_2014.csv"
GPR_PATH    <- "8.12/gpr_quarterly_processed.csv"
OIL_PATH    <- "8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"

OUT_DIR <- "8.12/country_specific_lag_selection"

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS      <- c("y","dp","r","de","deq")

SAMPLE_START <- "2000Q2"
SAMPLE_END   <- "2025Q3"

P_CANDIDATES <- c(1L, 2L)
Q_FIXED      <- 1L
LB_LAG       <- 4L
LB_ALPHA     <- 0.05
STABILITY_CUTOFF <- 1.0
BORDERLINE_RHO   <- 0.98

GPR_COL_EXACT <- "LN_GPR_QMEAN"
GDP_DLOG_DIVISOR <- 100

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------- UTILITIES --------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg   <- function(...) cat(sprintf(...), "\n")

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", ifelse(is.na(x), "", x)))
}

trim_chr <- function(x) trimws(as.character(x))

quarter_label <- function(qid) {
  yr <- (qid - 1L) %/% 4L
  qq <- qid - 4L * yr
  sprintf("%dQ%d", yr, qq)
}

safe_date_parse <- function(x) {
  sx <- trimws(as.character(x))
  out <- rep(as.Date(NA), length(sx))

  specs <- list(
    list(re="^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", fmt="%Y-%m-%d"),
    list(re="^[0-9]{4}/[0-9]{1,2}/[0-9]{1,2}$", fmt="%Y/%m/%d"),
    list(re="^[0-9]{4}\\.[0-9]{1,2}\\.[0-9]{1,2}$", fmt="%Y.%m.%d")
  )

  for (sp in specs) {
    ii <- which(is.na(out) & grepl(sp$re, sx))
    if (length(ii)) out[ii] <- suppressWarnings(as.Date(sx[ii], format=sp$fmt))
  }

  # YYYY-MM / YYYY/MM -> first day of month
  ii <- which(is.na(out) & grepl("^[0-9]{4}[-/][0-9]{1,2}$", sx))
  if (length(ii)) {
    z <- gsub("/", "-", sx[ii])
    z <- paste0(z, "-01")
    out[ii] <- suppressWarnings(as.Date(z, format="%Y-%m-%d"))
  }

  out
}

quarter_id <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    return(4L * (lt$year + 1900L) + (lt$mon %/% 3L) + 1L)
  }

  # Excel serial dates only if the vector is plausibly an Excel-date column.
  if (is.numeric(x)) {
    xx <- as.numeric(x)
    finite <- xx[is.finite(xx)]
    if (length(finite) && mean(finite >= 20000 & finite <= 70000) >= 0.80) {
      d <- as.Date(xx, origin="1899-12-30")
      lt <- as.POSIXlt(d)
      return(4L * (lt$year + 1900L) + (lt$mon %/% 3L) + 1L)
    }
  }

  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))

  # 2000Q1, 2000-Q1, 2000 Q1, 2000/Q1
  m <- regexec("^([0-9]{4})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L
  if (any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    qq <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L * yr + qq
  }

  # Never call as.Date() on arbitrary text.
  need <- which(is.na(out))
  if (length(need)) {
    d <- safe_date_parse(sx[need])
    ok2 <- !is.na(d)
    if (any(ok2)) {
      lt <- as.POSIXlt(d[ok2])
      out[need[ok2]] <- 4L * (lt$year + 1900L) + (lt$mon %/% 3L) + 1L
    }
  }

  out
}

quarter_parse_rate <- function(z) {
  q <- tryCatch(quarter_id(z), error=function(e) rep(NA_integer_, length(z)))
  if (!length(q)) return(0)
  mean(!is.na(q))
}

detect_date_col <- function(df, label="table") {
  if (!ncol(df)) stopf("%s has zero columns.", label)

  nn <- norm_name(names(df))
  preferred <- c("quarter","date","time","period","qtr","季度","日期","时间")
  hit <- match(norm_name(preferred), nn, nomatch=0L)
  hit <- hit[hit > 0L]

  if (length(hit)) {
    dcol <- names(df)[hit[1]]
    rr <- quarter_parse_rate(df[[dcol]])
    if (rr >= 0.50) {
      msg("[%s] date column: %s (parse %.1f%%)", label, dcol, 100*rr)
      return(dcol)
    }
  }

  rates <- vapply(df, quarter_parse_rate, numeric(1))
  rates[!is.finite(rates)] <- 0
  best <- which.max(rates)

  if (length(best) && rates[best] >= 0.80) {
    msg("[%s] date column by parse rate: %s (%.1f%%)",
        label, names(df)[best], 100*rates[best])
    return(names(df)[best])
  }

  ord <- order(rates, decreasing=TRUE)
  top <- head(ord, min(8L, length(ord)))
  diag <- paste(sprintf("%s=%.2f", names(df)[top], rates[top]), collapse=", ")
  stopf("Could not detect date/quarter column in %s. Best rates: %s", label, diag)
}

find_col_exact_norm <- function(df, candidates) {
  nn <- norm_name(names(df))
  cc <- norm_name(candidates)
  idx <- match(cc, nn, nomatch=0L)
  idx <- idx[idx > 0L]
  if (length(idx)) names(df)[idx[1]] else NA_character_
}

# ---------------------- MACRO EXCEL READER --------------------

ALIASES <- list(
  y      = c("y","gdploglevel","gdplog","reallogdp","gdp"),
  y_dlog = c("gdpdlog","ydlog","realgdpdlog","gdp_dlog"),
  dp     = c("dp","cpidlog","inflation","cpiinflation","cpi_dlog"),
  r      = c("r","rate","interestrate","shorttermrate","shortrate",
             "short_rate","interest_rate"),
  de     = c("de","reerdlog","exchangeratedlog","fxdlog","erdlog",
             "reer_dlog","exchange_rate_dlog"),
  deq    = c("deq","eqdlog","equitydlog","stockdlog","stockret",
             "equityret","eq_dlog","equity_dlog")
)

fill_right <- function(x) {
  x <- trim_chr(x)
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  last <- NA_character_
  for (i in seq_along(x)) {
    if (!is.na(x[i]) && nzchar(x[i])) last <- x[i]
    else if (!is.na(last)) x[i] <- last
  }
  x
}

country_var_candidates <- function(cc, aliases) {
  unique(c(
    paste0(cc, "_", aliases),
    paste0(aliases, "_", cc),
    paste0(cc, aliases),
    paste0(aliases, cc)
  ))
}

mapping_score <- function(nms) {
  nn <- norm_name(nms)
  score <- 0L
  for (cc in COUNTRIES) {
    for (vv in c("y","y_dlog","dp","r","de","deq")) {
      if (any(nn %in% norm_name(country_var_candidates(cc, ALIASES[[vv]])))) {
        score <- score + 1L
      }
    }
  }
  score
}

build_two_row_header <- function(raw, h1, h2) {
  if (h1 >= h2 || h2 >= nrow(raw)) return(NULL)

  a <- trim_chr(unlist(raw[h1, , drop=TRUE], use.names=FALSE))
  b <- trim_chr(unlist(raw[h2, , drop=TRUE], use.names=FALSE))
  af <- fill_right(a)

  nms <- character(ncol(raw))
  for (j in seq_len(ncol(raw))) {
    aj <- ifelse(is.na(af[j]), "", af[j])
    bj <- ifelse(is.na(b[j]), "", b[j])

    if (j == 1L && norm_name(bj) %in% norm_name(c("quarter","date","time","period","qtr","季度","日期","时间"))) {
      nms[j] <- bj
    } else if (j == 1L && norm_name(aj) %in% norm_name(c("quarter","date","time","period","qtr","季度","日期","时间"))) {
      nms[j] <- aj
    } else if (nzchar(aj) && nzchar(bj) && norm_name(aj) != norm_name(bj)) {
      nms[j] <- paste0(aj, "_", bj)
    } else if (nzchar(bj)) {
      nms[j] <- bj
    } else if (nzchar(aj)) {
      nms[j] <- aj
    } else {
      nms[j] <- paste0("V", j)
    }
  }

  dat <- raw[(h2+1L):nrow(raw), , drop=FALSE]
  names(dat) <- make.unique(nms, sep="_dup")
  as.data.frame(dat, check.names=FALSE)
}

read_macro_table <- function(path, sheet=1) {
  if (!file.exists(path)) stopf("Macro file not found: %s", path)
  if (!requireNamespace("readxl", quietly=TRUE)) stopf("Package readxl is required.")

  std <- as.data.frame(
    readxl::read_excel(path, sheet=sheet, .name_repair="unique"),
    check.names=FALSE
  )

  std_score <- mapping_score(names(std))
  std_date <- tryCatch(detect_date_col(std, "macro-standard"), error=function(e) NA_character_)
  msg("[macro] one-row header mapping score = %d", std_score)

  if (!is.na(std_date) && std_score >= 20L) {
    msg("[macro] using ordinary one-row header")
    return(std)
  }

  raw <- as.data.frame(
    readxl::read_excel(path, sheet=sheet, col_names=FALSE,
                       .name_repair="minimal", guess_max=500),
    check.names=FALSE
  )

  best <- NULL
  best_score <- -Inf
  max_header <- min(6L, nrow(raw)-1L)

  for (h1 in seq_len(max_header)) {
    for (h2 in seq_len(max_header)) {
      if (h1 >= h2) next
      cand <- build_two_row_header(raw, h1, h2)
      if (is.null(cand)) next

      sc <- mapping_score(names(cand))
      dc <- tryCatch(detect_date_col(cand, sprintf("macro-h%d-h%d",h1,h2)),
                     error=function(e) NA_character_)
      dr <- if (!is.na(dc)) quarter_parse_rate(cand[[dc]]) else 0
      total <- sc + 20*dr

      if (total > best_score) {
        best_score <- total
        best <- list(data=cand, h1=h1, h2=h2, score=sc, dcol=dc, drate=dr)
      }
    }
  }

  if (is.null(best) || is.na(best$dcol) || best$score < 20L) {
    stopf("Could not reconstruct the macro Excel header. Best mapping score=%s",
          ifelse(is.null(best), "none", best$score))
  }

  msg("[macro] using two-row header rows %d/%d; mapping score=%d; date=%s",
      best$h1, best$h2, best$score, best$dcol)
  best$data
}

extract_macro_long <- function(raw) {
  dcol <- detect_date_col(raw, "macro")
  qid <- quarter_id(raw[[dcol]])

  nn <- norm_name(names(raw))
  out <- list()

  for (cc in COUNTRIES) {
    pick <- function(vv) {
      cand <- country_var_candidates(cc, ALIASES[[vv]])
      idx <- match(norm_name(cand), nn, nomatch=0L)
      idx <- idx[idx > 0L]
      if (length(idx)) names(raw)[idx[1]] else NA_character_
    }

    cy  <- pick("y")
    cyd <- pick("y_dlog")
    cdp <- pick("dp")
    cr  <- pick("r")
    cde <- pick("de")
    ceq <- pick("deq")

    if (any(is.na(c(cdp,cr,cde,ceq))) || (is.na(cy) && is.na(cyd))) {
      stopf(
        "Macro mapping failed for %s: y=%s y_dlog=%s dp=%s r=%s de=%s deq=%s",
        cc, cy, cyd, cdp, cr, cde, ceq
      )
    }

    y <- if (!is.na(cy)) suppressWarnings(as.numeric(raw[[cy]])) else NULL
    if (is.null(y)) {
      gd <- suppressWarnings(as.numeric(raw[[cyd]])) / GDP_DLOG_DIVISOR
      y <- rep(NA_real_, length(gd))
      ok <- which(is.finite(gd) & !is.na(qid))
      if (length(ok)) {
        # Reconstruct an indexed log level separately over consecutive available observations.
        y[ok] <- log(100) + cumsum(replace(gd[ok], is.na(gd[ok]), 0))
      }
    }

    tmp <- data.frame(
      qid=qid,
      country=cc,
      y=y,
      dp=suppressWarnings(as.numeric(raw[[cdp]])),
      r=suppressWarnings(as.numeric(raw[[cr]])),
      de=suppressWarnings(as.numeric(raw[[cde]])),
      deq=suppressWarnings(as.numeric(raw[[ceq]]))
    )

    out[[cc]] <- tmp[!is.na(tmp$qid), ]
    msg("[macro] %s mapped: y=%s; dp=%s; r=%s; de=%s; deq=%s",
        cc, ifelse(!is.na(cy), cy, cyd), cdp, cr, cde, ceq)
  }

  do.call(rbind, out)
}

# --------------------- WEIGHTS / GLOBAL DATA ------------------

read_weights <- function(path) {
  if (!file.exists(path)) stopf("Weights file not found: %s", path)
  w0 <- read.csv(path, check.names=FALSE)

  # Allow first column to contain row country codes.
  first <- toupper(trimws(as.character(w0[[1]])))
  if (all(COUNTRIES %in% first)) {
    rownames(w0) <- first
    w0 <- w0[, -1, drop=FALSE]
  }

  # Columns must contain country codes after normalization.
  cn <- toupper(trimws(names(w0)))
  if (!all(COUNTRIES %in% cn)) {
    stopf("Weight columns do not contain all 14 economies. Columns: %s",
          paste(names(w0), collapse=", "))
  }
  w0 <- w0[, match(COUNTRIES, cn), drop=FALSE]

  if (is.null(rownames(w0)) || !all(COUNTRIES %in% toupper(rownames(w0)))) {
    if (nrow(w0) != length(COUNTRIES)) stopf("Weights must have 14 rows.")
    rownames(w0) <- COUNTRIES
  } else {
    w0 <- w0[match(COUNTRIES, toupper(rownames(w0))), , drop=FALSE]
  }

  W <- as.matrix(data.frame(lapply(w0, as.numeric), check.names=FALSE))
  rownames(W) <- COUNTRIES
  colnames(W) <- COUNTRIES

  diag(W) <- 0
  rs <- rowSums(W, na.rm=TRUE)
  if (any(!is.finite(rs) | rs <= 0)) stopf("At least one weight row has nonpositive sum.")
  W <- W / rs

  msg("[weights] loaded and row-normalized; max |rowSum-1| = %.3g",
      max(abs(rowSums(W)-1)))
  W
}

read_global_series <- function(path, value_candidates, exact=NULL, label="global") {
  if (!file.exists(path)) stopf("%s file not found: %s", label, path)
  d <- read.csv(path, check.names=FALSE)
  dc <- detect_date_col(d, label)
  qid <- quarter_id(d[[dc]])

  vc <- NA_character_
  if (!is.null(exact) && exact %in% names(d)) vc <- exact
  if (is.na(vc)) vc <- find_col_exact_norm(d, value_candidates)
  if (is.na(vc)) {
    # Fallback: choose numeric-like non-date column with most finite observations.
    candidates <- setdiff(names(d), dc)
    rates <- vapply(candidates, function(nm) {
      z <- suppressWarnings(as.numeric(d[[nm]]))
      mean(is.finite(z))
    }, numeric(1))
    if (length(rates) && max(rates, na.rm=TRUE) >= 0.8) vc <- candidates[which.max(rates)]
  }
  if (is.na(vc)) stopf("Could not identify value column in %s.", path)

  out <- data.frame(qid=qid, value=suppressWarnings(as.numeric(d[[vc]])))
  out <- out[!is.na(out$qid) & is.finite(out$value), ]
  out <- aggregate(value ~ qid, data=out, FUN=mean)
  msg("[%s] value column: %s; %d quarters", label, vc, nrow(out))
  out
}

# ------------------------- MODEL HELPERS -----------------------

lag_matrix <- function(X, lag) {
  nr <- nrow(X)
  out <- matrix(NA_real_, nr, ncol(X))
  if (lag < nr) out[(lag+1L):nr, ] <- X[1L:(nr-lag), , drop=FALSE]
  colnames(out) <- paste0(colnames(X), "_L", lag)
  out
}

companion_radius <- function(B_domestic, k, p) {
  # B_domestic: list of p coefficient matrices k x k for domestic lags.
  if (p == 1L) return(max(Mod(eigen(B_domestic[[1]], only.values=TRUE)$values)))

  top <- do.call(cbind, B_domestic)
  bottom <- cbind(diag(k*(p-1L)), matrix(0, k*(p-1L), k))
  C <- rbind(top, bottom)
  max(Mod(eigen(C, only.values=TRUE)$values))
}

fit_country <- function(cc, p, panel_wide, foreign_wide, globals) {
  Xi <- as.matrix(panel_wide[[cc]][, VARS, drop=FALSE])
  Xs <- as.matrix(foreign_wide[[cc]][, VARS, drop=FALSE])
  qid <- panel_wide[[cc]]$qid

  n <- nrow(Xi)
  if (n != nrow(Xs)) stopf("%s: domestic/foreign row mismatch.", cc)

  # Same p=1/p=2 comparison sample: always start after 2 domestic lags.
  start_idx <- max(P_CANDIDATES) + 1L
  idx <- seq.int(start_idx, n)

  Y <- Xi[idx, , drop=FALSE]
  Z <- matrix(1, length(idx), 1)
  colnames(Z) <- "const"

  dom_col_ranges <- vector("list", p)
  cursor <- 1L

  for (L in seq_len(p)) {
    XL <- lag_matrix(Xi, L)[idx, , drop=FALSE]
    colnames(XL) <- paste0("D_", VARS, "_L", L)
    st <- ncol(Z)+1L
    Z <- cbind(Z, XL)
    en <- ncol(Z)
    dom_col_ranges[[L]] <- st:en
  }

  # q_i = 1: contemporaneous and one lag of foreign variables.
  Xs0 <- Xs[idx, , drop=FALSE]
  colnames(Xs0) <- paste0("F_", VARS, "_L0")
  Xs1 <- lag_matrix(Xs, 1L)[idx, , drop=FALSE]
  colnames(Xs1) <- paste0("F_", VARS, "_L1")
  Z <- cbind(Z, Xs0, Xs1)

  # GPR and Brent: contemporaneous + one lag.
  G <- globals[match(qid, globals$qid), c("gpr","oil")]
  G0 <- as.matrix(G[idx, , drop=FALSE])
  colnames(G0) <- c("GPR_L0","OIL_L0")
  G1 <- lag_matrix(as.matrix(G), 1L)[idx, , drop=FALSE]
  colnames(G1) <- c("GPR_L1","OIL_L1")
  Z <- cbind(Z, G0, G1)

  keep <- complete.cases(Y) & complete.cases(Z)
  Y <- Y[keep, , drop=FALSE]
  Z <- Z[keep, , drop=FALSE]

  if (nrow(Y) < ncol(Z) + 10L) {
    stopf("%s p=%d: too few complete observations (%d) for %d regressors.",
          cc, p, nrow(Y), ncol(Z))
  }

  fit <- lm.fit(x=Z, y=Y)
  B <- fit$coefficients
  E <- fit$residuals

  Tn <- nrow(E)
  K  <- ncol(Y)
  m  <- ncol(Z)

  Sigma <- crossprod(E) / Tn
  detS <- determinant(Sigma, logarithm=TRUE)
  if (detS$sign <= 0) stopf("%s p=%d: residual covariance is singular.", cc, p)
  logdet <- as.numeric(detS$modulus)

  # Gaussian system IC up to constants common across candidates.
  npar <- K*m
  aic <- Tn*logdet + 2*npar
  bic <- Tn*logdet + log(Tn)*npar

  domestic_mats <- lapply(seq_len(p), function(L) {
    # lm.fit coefficient matrix: regressors x equations.
    t(B[dom_col_ranges[[L]], , drop=FALSE])
  })
  rho <- companion_radius(domestic_mats, K, p)

  lb <- do.call(rbind, lapply(seq_len(K), function(j) {
    bt <- Box.test(E[,j], lag=min(LB_LAG, max(1L, floor(Tn/5))),
                   type="Ljung-Box", fitdf=0)
    data.frame(
      country=cc, p=p, equation=VARS[j],
      lb_lag=min(LB_LAG, max(1L, floor(Tn/5))),
      statistic=unname(bt$statistic),
      p_value=bt$p.value,
      serial_corr_reject_5pct=(bt$p.value < LB_ALPHA)
    )
  }))

  diag <- data.frame(
    country=cc, p=p, nobs=Tn, regressors=m,
    AIC=aic, BIC=bic, spectral_radius=rho,
    stable=(is.finite(rho) && rho < STABILITY_CUTOFF),
    borderline=(is.finite(rho) && rho >= BORDERLINE_RHO && rho < STABILITY_CUTOFF),
    lb_rejections_5pct=sum(lb$serial_corr_reject_5pct),
    lb_all_pass_5pct=all(!lb$serial_corr_reject_5pct)
  )

  list(diag=diag, lb=lb)
}

# ---------------------------- MAIN -----------------------------

msg("Reading macro panel...")
macro_raw <- read_macro_table(MACRO_PATH, MACRO_SHEET)
macro_long <- extract_macro_long(macro_raw)

q_start <- quarter_id(SAMPLE_START)
q_end   <- quarter_id(SAMPLE_END)
macro_long <- macro_long[
  macro_long$qid >= q_start & macro_long$qid <= q_end,
]

W <- read_weights(WEIGHT_PATH)

gpr <- read_global_series(
  GPR_PATH,
  value_candidates=c("LN_GPR_QMEAN","GPR","LN_GPR","GPR_QMEAN"),
  exact=GPR_COL_EXACT,
  label="GPR"
)
names(gpr)[2] <- "gpr"

oil <- read_global_series(
  OIL_PATH,
  value_candidates=c("BRENT_LOG","LOG_BRENT","LN_BRENT","BRENT","OIL","LOG_PRICE","VALUE"),
  label="OIL"
)
names(oil)[2] <- "oil"

globals <- merge(gpr, oil, by="qid", all=FALSE)
globals <- globals[globals$qid >= q_start & globals$qid <= q_end, ]

# Common quarter set across all countries + globals.
country_q <- lapply(COUNTRIES, function(cc) {
  z <- macro_long[macro_long$country == cc, ]
  z$qid[complete.cases(z[, VARS, drop=FALSE])]
})
common_q <- Reduce(intersect, c(country_q, list(globals$qid)))
common_q <- sort(common_q)

if (length(common_q) < 30L) stopf("Common complete sample is too short: %d quarters.", length(common_q))
msg("Common complete sample: %s to %s (%d quarters)",
    quarter_label(min(common_q)), quarter_label(max(common_q)), length(common_q))

# Build aligned domestic country panels.
panel <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
for (cc in COUNTRIES) {
  z <- macro_long[macro_long$country == cc & macro_long$qid %in% common_q,
                  c("qid",VARS)]
  z <- z[match(common_q, z$qid), ]
  if (any(is.na(z$qid)) || any(!complete.cases(z[,VARS,drop=FALSE]))) {
    stopf("%s is not complete after common-sample alignment.", cc)
  }
  panel[[cc]] <- z
}
globals <- globals[match(common_q, globals$qid), ]

# Build weighted foreign variables.
foreign <- setNames(vector("list", length(COUNTRIES)), COUNTRIES)
for (i in seq_along(COUNTRIES)) {
  cc <- COUNTRIES[i]
  F <- matrix(0, length(common_q), length(VARS))
  colnames(F) <- VARS
  for (j in seq_along(COUNTRIES)) {
    F <- F + W[i,j] * as.matrix(panel[[COUNTRIES[j]]][,VARS,drop=FALSE])
  }
  foreign[[cc]] <- data.frame(qid=common_q, F, check.names=FALSE)
}

all_diag <- list()
all_lb <- list()
kk <- 1L

for (cc in COUNTRIES) {
  msg("------------------------------------------------------------")
  msg("Country: %s", cc)
  for (p in P_CANDIDATES) {
    msg("  fitting p=%d, q=%d", p, Q_FIXED)
    ans <- fit_country(cc, p, panel, foreign, globals)
    all_diag[[kk]] <- ans$diag
    all_lb[[kk]] <- ans$lb
    kk <- kk + 1L
  }
}

diag_df <- do.call(rbind, all_diag)
lb_df   <- do.call(rbind, all_lb)

# ----------------------- RECOMMENDATION ------------------------

recommend_one <- function(d) {
  d <- d[order(d$p), ]
  d1 <- d[d$p==1L, ]
  d2 <- d[d$p==2L, ]

  if (!nrow(d1) || !nrow(d2)) stopf("Missing p=1 or p=2 diagnostic for %s.", d$country[1])

  # Hard stability rule first.
  if (d1$stable && !d2$stable) {
    return(data.frame(country=d$country[1], selected_p=1L,
                      reason="p=2 unstable; p=1 stable"))
  }
  if (!d1$stable && d2$stable) {
    return(data.frame(country=d$country[1], selected_p=2L,
                      reason="p=1 unstable; p=2 stable"))
  }
  if (!d1$stable && !d2$stable) {
    # Still report the less explosive option; flag for model redesign.
    pp <- ifelse(d1$spectral_radius <= d2$spectral_radius, 1L, 2L)
    return(data.frame(country=d$country[1], selected_p=pp,
                      reason="WARNING: both p=1 and p=2 locally unstable; selected lower spectral radius"))
  }

  # If p=1 has residual autocorrelation but p=2 clears all equations, prefer p=2.
  if (!d1$lb_all_pass_5pct && d2$lb_all_pass_5pct) {
    return(data.frame(country=d$country[1], selected_p=2L,
                      reason="p=2 removes Ljung-Box rejections present under p=1"))
  }

  # Otherwise use BIC as the primary parsimony criterion.
  if (d1$BIC <= d2$BIC) {
    reason <- if (d1$borderline) "BIC prefers p=1, but p=1 is borderline stable; inspect before final TVP estimation" else
      "both acceptable; BIC prefers parsimonious p=1"
    return(data.frame(country=d$country[1], selected_p=1L, reason=reason))
  } else {
    reason <- if (d2$borderline) "BIC prefers p=2, but p=2 is borderline stable; inspect before final TVP estimation" else
      "both acceptable; BIC prefers p=2"
    return(data.frame(country=d$country[1], selected_p=2L, reason=reason))
  }
}

rec <- do.call(rbind, lapply(split(diag_df, diag_df$country), recommend_one))
rec <- rec[match(COUNTRIES, rec$country), ]

# Add comparison fields.
rec$AIC_p1 <- diag_df$AIC[match(paste(rec$country,1), paste(diag_df$country,diag_df$p))]
rec$AIC_p2 <- diag_df$AIC[match(paste(rec$country,2), paste(diag_df$country,diag_df$p))]
rec$BIC_p1 <- diag_df$BIC[match(paste(rec$country,1), paste(diag_df$country,diag_df$p))]
rec$BIC_p2 <- diag_df$BIC[match(paste(rec$country,2), paste(diag_df$country,diag_df$p))]
rec$rho_p1 <- diag_df$spectral_radius[match(paste(rec$country,1), paste(diag_df$country,diag_df$p))]
rec$rho_p2 <- diag_df$spectral_radius[match(paste(rec$country,2), paste(diag_df$country,diag_df$p))]
rec$LB_reject_p1 <- diag_df$lb_rejections_5pct[match(paste(rec$country,1), paste(diag_df$country,diag_df$p))]
rec$LB_reject_p2 <- diag_df$lb_rejections_5pct[match(paste(rec$country,2), paste(diag_df$country,diag_df$p))]

write.csv(diag_df,
          file.path(OUT_DIR, "01_country_p1_p2_model_diagnostics.csv"),
          row.names=FALSE)
write.csv(lb_df,
          file.path(OUT_DIR, "02_equation_residual_serial_correlation.csv"),
          row.names=FALSE)
write.csv(rec,
          file.path(OUT_DIR, "03_country_lag_recommendation.csv"),
          row.names=FALSE)

lag_vec <- setNames(rec$selected_p, rec$country)
lag_line <- paste(sprintf("%s=%dL", names(lag_vec), lag_vec), collapse=", ")
writeLines(
  c(
    "# Auto-generated country-specific domestic lag vector",
    sprintf("COUNTRY_P <- c(%s)", lag_line)
  ),
  file.path(OUT_DIR, "04_selected_country_lag_vector.R")
)

readme <- c(
  "Country-specific lag selection: GPR + Brent GVAR / TVP-GVAR",
  "",
  sprintf("Requested sample window: %s to %s", SAMPLE_START, SAMPLE_END),
  sprintf("Actual common complete sample: %s to %s (%d quarters)",
          quarter_label(min(common_q)), quarter_label(max(common_q)), length(common_q)),
  "Countries: AU BR CA CH CN EA UK JP KR NO SG TR US ZA",
  "Domestic variables: y dp r de deq",
  "Candidate domestic lags: p_i = 1, 2",
  "Foreign-variable lag order: q_i = 1 (current + one lag included)",
  "Global variables: GPR and Brent oil, current + one lag",
  "",
  "Selection logic:",
  "1. Reject a candidate if local domestic companion spectral radius >= 1.",
  "2. If p=1 has Ljung-Box residual rejections and p=2 clears them all, choose p=2.",
  "3. Otherwise choose the stable candidate with lower BIC.",
  "4. A spectral radius >= 0.98 but < 1 is flagged borderline.",
  "5. If both candidates are unstable, the lower-radius choice is reported with WARNING.",
  "",
  "Important:",
  "This is a conditional country-model pre-estimation diagnostic.",
  "Final TVP-GVAR posterior draws still require global/posterior stability checks.",
  "p=1 and p=2 are compared on the same common sample (trimmed to allow p=2)."
)
writeLines(readme, file.path(OUT_DIR, "README_country_specific_lag_selection.txt"))

msg("============================================================")
msg("DONE")
msg("Selected lag vector:")
print(lag_vec)
msg("Outputs written to: %s", OUT_DIR)
