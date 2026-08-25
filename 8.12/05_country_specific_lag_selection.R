#!/usr/bin/env Rscript

# ============================================================
# Country-specific lag selection for GPR + Oil TVP-GVAR
# Purpose:
#   Compare p_i = 1 vs 2 country-by-country, holding q_i = 1.
#   Diagnostics:
#     1) system AIC / BIC
#     2) equation-level Ljung-Box residual serial correlation
#     3) local VAR companion spectral radius
#     4) transparent recommended p_i
#
# IMPORTANT:
#   - This is a PRE-ESTIMATION / SPECIFICATION diagnostic.
#   - It does NOT replace the final TVP posterior stability check.
#   - p=1 and p=2 are estimated on the SAME common sample
#     (trimmed for max p=2) so AIC/BIC and residual tests are comparable.
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
STABILITY_CUTOFF <- 1.0
BORDERLINE_RHO   <- 0.98

# Current validated global variables
GPR_COL_EXACT <- "LN_GPR_QMEAN"

# If your wide macro file contains GDP_DLOG rather than y/log-level,
# the script reconstructs indexed log GDP exactly as:
#   log_index[t] = log(100) + cumsum(GDP_DLOG/100)
# with the first available in-sample log index normalized to log(100).
GDP_DLOG_DIVISOR <- 100

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------ BASIC HELPERS ------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg   <- function(...) cat(sprintf(...), "\n")

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

quarter_id <- function(x) {
  # Returns integer 4*year + quarter, useful for sorting and filtering.
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    yr <- lt$year + 1900L
    q  <- (lt$mon %/% 3L) + 1L
    return(4L * yr + q)
  }

  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))

  # 2000Q1, 2000-Q1, 2000 Q1, etc.
  m <- regexec("^([0-9]{4})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L
  if (any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    q  <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L * yr + q
  }

  # Date-like fallback.
  need <- is.na(out)
  if (any(need)) {
    suppressWarnings({
      d <- as.Date(sx[need])
    })
    ok2 <- !is.na(d)
    if (any(ok2)) {
      lt <- as.POSIXlt(d[ok2])
      yr <- lt$year + 1900L
      q  <- (lt$mon %/% 3L) + 1L
      idx <- which(need)[ok2]
      out[idx] <- 4L * yr + q
    }
  }

  out
}

quarter_label <- function(qid) {
  yr <- (qid - 1L) %/% 4L
  q  <- qid - 4L * yr
  sprintf("%dQ%d", yr, q)
}

detect_date_col <- function(df) {
  nn <- norm_name(names(df))
  preferred <- c("quarter","date","time","period","qtr")
  hit <- match(preferred, nn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(names(df)[hit[1]])

  # Try any column that parses mostly as quarters.
  rates <- vapply(df, function(z) mean(!is.na(quarter_id(z))), numeric(1))
  if (max(rates, na.rm = TRUE) >= 0.8) {
    return(names(df)[which.max(rates)])
  }
  stopf("Could not detect a quarter/date column. Available columns: %s",
        paste(names(df), collapse = ", "))
}

read_table_auto <- function(path, sheet = 1) {
  if (!file.exists(path)) stopf("File not found: %s", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx","xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stopf("Package 'readxl' is required for Excel input. Install with install.packages('readxl').")
    }
    return(as.data.frame(readxl::read_excel(path, sheet = sheet), check.names = FALSE))
  }
  if (ext == "csv") {
    return(read.csv(path, check.names = FALSE))
  }
  if (ext == "rds") return(readRDS(path))
  stopf("Unsupported file type: %s", path)
}

find_col_exact_norm <- function(df, candidates) {
  nn <- norm_name(names(df))
  cc <- norm_name(candidates)
  idx <- match(cc, nn, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) names(df)[idx[1]] else NA_character_
}

# --------------------- MACRO PANEL ADAPTER ---------------------

# Aliases allow the script to read common versions of your processed macro file.
ALIASES <- list(
  y = c("y", "gdploglevel", "gdplog", "reallogdp", "gdp"),
  y_dlog = c("gdpdlog", "ydlog", "realgdpdlog"),
  dp = c("dp", "cpidlog", "inflation", "cpiinflation"),
  r  = c("r", "rate", "interestrate", "shorttermrate", "shortrate"),
  de = c("de", "reerdlog", "exchangeratedlog", "fxdlog", "erdlog"),
  deq = c("deq", "eqdlog", "equitydlog", "stockdlog", "stockret", "equityret")
)

extract_macro_long <- function(raw, countries = COUNTRIES) {
  dcol <- detect_date_col(raw)
  qid  <- quarter_id(raw[[dcol]])
  if (mean(!is.na(qid)) < 0.8) stopf("Macro date column '%s' could not be parsed reliably.", dcol)

  # ---------- Long format ----------
  country_col <- find_col_exact_norm(raw, c("country","economy","unit","iso","code"))
  if (!is.na(country_col)) {
    value_cols <- sapply(c("y","dp","r","de","deq"), function(v) {
      find_col_exact_norm(raw, ALIASES[[v]])
    })

    if (all(!is.na(value_cols))) {
      out <- data.frame(
        qid = qid,
        country = toupper(trimws(as.character(raw[[country_col]]))),
        y   = as.numeric(raw[[value_cols["y"]]]),
        dp  = as.numeric(raw[[value_cols["dp"]]]),
        r   = as.numeric(raw[[value_cols["r"]]]),
        de  = as.numeric(raw[[value_cols["de"]]]),
        deq = as.numeric(raw[[value_cols["deq"]]])
      )
      out <- out[out$country %in% countries & !is.na(out$qid), ]
      return(out)
    }
  }

  # ---------- Wide format ----------
  nn <- norm_name(names(raw))
  out_list <- vector("list", length(countries))
  names(out_list) <- countries

  for (cc in countries) {
    find_country_var <- function(alias_vec) {
      cand1 <- paste0(cc, "_", alias_vec)
      cand2 <- paste0(alias_vec, "_", cc)
      cands <- c(cand1, cand2, paste0(cc, alias_vec), paste0(alias_vec, cc))
      idx <- match(norm_name(cands), nn, nomatch = 0L)
      idx <- idx[idx > 0L]
      if (length(idx)) names(raw)[idx[1]] else NA_character_
    }

    cy   <- find_country_var(ALIASES$y)
    cyd  <- find_country_var(ALIASES$y_dlog)
    cdp  <- find_country_var(ALIASES$dp)
    cr   <- find_country_var(ALIASES$r)
    cde  <- find_country_var(ALIASES$de)
    cdeq <- find_country_var(ALIASES$deq)

    missing_non_gdp <- c(dp=cdp, r=cr, de=cde, deq=cdeq)
    if (any(is.na(missing_non_gdp))) {
      stopf(
        paste0(
          "Could not map all 5 macro variables for %s.\n",
          "Mapped: y=%s, y_dlog=%s, dp=%s, r=%s, de=%s, deq=%s\n",
          "Edit ALIASES near the top if your column names differ.\n",
          "Available columns begin: %s"
        ),
        cc, cy, cyd, cdp, cr, cde, cdeq,
        paste(head(names(raw), 40), collapse = ", ")
      )
    }

    if (!is.na(cy)) {
      yval <- as.numeric(raw[[cy]])
    } else if (!is.na(cyd)) {
      # Reconstruct log-level index from validated GDP_DLOG.
      gd <- as.numeric(raw[[cyd]]) / GDP_DLOG_DIVISOR
      yval <- rep(NA_real_, length(gd))
      ok <- which(!is.na(qid) & !is.na(gd))
      if (!length(ok)) stopf("No usable GDP_DLOG observations for %s.", cc)
      ok <- ok[order(qid[ok])]
      yval[ok[1]] <- log(100)
      if (length(ok) >= 2L) {
        for (k in 2:length(ok)) {
          cur <- ok[k]
          prv <- ok[k-1]
          # Only cumulate over consecutive observed rows in sorted sample.
          yval[cur] <- yval[prv] + gd[cur]
        }
      }
    } else {
      stopf("Could not find y/log-GDP or GDP_DLOG for %s.", cc)
    }

    out_list[[cc]] <- data.frame(
      qid = qid,
      country = cc,
      y   = yval,
      dp  = as.numeric(raw[[cdp]]),
      r   = as.numeric(raw[[cr]]),
      de  = as.numeric(raw[[cde]]),
      deq = as.numeric(raw[[cdeq]])
    )
  }

  do.call(rbind, out_list)
}

# ------------------------- WEIGHT MATRIX ------------------------

read_weight_matrix <- function(path, countries = COUNTRIES) {
  wraw <- read_table_auto(path)
  nn <- norm_name(names(wraw))

  # Long form: origin / partner / weight
  ocol <- find_col_exact_norm(wraw, c("origin","from","reporter","countryi","i"))
  pcol <- find_col_exact_norm(wraw, c("partner","to","counterparty","countryj","j"))
  wcol <- find_col_exact_norm(wraw, c("weight","w","tradeweight"))

  if (!is.na(ocol) && !is.na(pcol) && !is.na(wcol)) {
    W <- matrix(0, length(countries), length(countries),
                dimnames = list(countries, countries))
    oo <- toupper(trimws(as.character(wraw[[ocol]])))
    pp <- toupper(trimws(as.character(wraw[[pcol]])))
    ww <- as.numeric(wraw[[wcol]])
    for (k in seq_len(nrow(wraw))) {
      if (oo[k] %in% countries && pp[k] %in% countries && is.finite(ww[k])) {
        W[oo[k], pp[k]] <- ww[k]
      }
    }
  } else {
    # Wide square matrix. Find row-country column.
    first_vals <- toupper(trimws(as.character(wraw[[1]])))
    if (sum(first_vals %in% countries) >= length(countries) - 2L) {
      row_country <- first_vals
      dat <- wraw[, -1, drop = FALSE]
    } else {
      stopf("Weight file format not recognized. Expected long origin-partner-weight or square matrix.")
    }

    col_map <- match(countries, toupper(names(dat)))
    if (anyNA(col_map)) {
      # normalized fallback
      col_map <- match(norm_name(countries), norm_name(names(dat)))
    }
    if (anyNA(col_map)) {
      stopf("Could not map all country columns in weight matrix. Columns: %s",
            paste(names(dat), collapse = ", "))
    }

    W <- matrix(NA_real_, length(countries), length(countries),
                dimnames = list(countries, countries))
    for (cc in countries) {
      ri <- which(row_country == cc)
      if (length(ri) != 1L) stopf("Weight matrix row for %s is missing or duplicated.", cc)
      W[cc, ] <- as.numeric(dat[ri, col_map])
    }
  }

  diag(W) <- 0

  rs <- rowSums(W, na.rm = TRUE)
  if (any(!is.finite(rs) | rs <= 0)) stopf("Invalid weight row sums.")
  W <- W / rs

  if (max(abs(rowSums(W) - 1)) > 1e-8) stopf("Weight normalization failed.")
  W
}

# ---------------------- GLOBAL SERIES READER --------------------

read_global_series <- function(path, value_candidates, exact = NULL, label = "global") {
  df <- read_table_auto(path)
  dcol <- detect_date_col(df)
  qid <- quarter_id(df[[dcol]])

  vcol <- NA_character_
  if (!is.null(exact)) {
    vcol <- find_col_exact_norm(df, exact)
  }
  if (is.na(vcol)) {
    vcol <- find_col_exact_norm(df, value_candidates)
  }

  if (is.na(vcol)) {
    # fuzzy fallback
    nn <- norm_name(names(df))
    keep <- setdiff(seq_along(nn), match(norm_name(dcol), nn))
    num <- vapply(df, is.numeric, logical(1))
    keep <- keep[num[keep]]
    pattern <- paste(norm_name(value_candidates), collapse = "|")
    hit <- keep[grepl(pattern, nn[keep])]
    if (length(hit)) vcol <- names(df)[hit[1]]
  }

  if (is.na(vcol)) {
    stopf("Could not detect %s value column. Available: %s",
          label, paste(names(df), collapse = ", "))
  }

  out <- data.frame(qid = qid, value = as.numeric(df[[vcol]]))
  out <- out[!is.na(out$qid), ]
  out <- out[order(out$qid), ]
  out <- out[!duplicated(out$qid), ]
  names(out)[2] <- label
  attr(out, "value_column") <- vcol
  out
}

# ---------------------- PANEL CONSTRUCTION ----------------------

build_country_frames <- function(panel_long, W, gpr, oil) {
  # Restrict sample.
  lo <- quarter_id(SAMPLE_START)
  hi <- quarter_id(SAMPLE_END)
  panel_long <- panel_long[panel_long$qid >= lo & panel_long$qid <= hi, ]

  qs <- sort(unique(panel_long$qid))
  expected <- seq(lo, hi)
  if (!all(expected %in% qs)) {
    miss <- quarter_label(setdiff(expected, qs))
    stopf("Macro panel is missing model quarters: %s", paste(miss, collapse = ", "))
  }

  # country x time x variable arrays
  frames <- vector("list", length(COUNTRIES))
  names(frames) <- COUNTRIES

  # Build matrices T x N for each domestic variable.
  Xmats <- list()
  for (v in VARS) {
    M <- matrix(NA_real_, length(expected), length(COUNTRIES),
                dimnames = list(quarter_label(expected), COUNTRIES))
    for (j in seq_along(COUNTRIES)) {
      cc <- COUNTRIES[j]
      d <- panel_long[panel_long$country == cc, c("qid", v)]
      m <- match(expected, d$qid)
      M[, j] <- as.numeric(d[[v]][m])
    }
    if (anyNA(M)) {
      bad <- which(is.na(M), arr.ind = TRUE)
      ex <- head(bad, 10)
      txt <- apply(ex, 1, function(z) sprintf("%s:%s",
                    rownames(M)[z[1]], colnames(M)[z[2]]))
      stopf("Missing macro values after sample alignment for %s, e.g. %s",
            v, paste(txt, collapse = ", "))
    }
    Xmats[[v]] <- M
  }

  mg <- match(expected, gpr$qid)
  mo <- match(expected, oil$qid)
  if (anyNA(mg)) stopf("GPR is missing quarters in model sample.")
  if (anyNA(mo)) stopf("Oil is missing quarters in model sample.")
  gv <- gpr$gpr[mg]
  ov <- oil$oil[mo]
  if (anyNA(gv) || anyNA(ov)) stopf("GPR/Oil has NA in model sample.")

  for (i in seq_along(COUNTRIES)) {
    cc <- COUNTRIES[i]
    wi <- W[cc, COUNTRIES]

    df <- data.frame(
      qid = expected,
      quarter = quarter_label(expected),
      y   = Xmats$y[, cc],
      dp  = Xmats$dp[, cc],
      r   = Xmats$r[, cc],
      de  = Xmats$de[, cc],
      deq = Xmats$deq[, cc],
      stringsAsFactors = FALSE
    )

    for (v in VARS) {
      df[[paste0("star_", v)]] <- as.numeric(Xmats[[v]] %*% wi)
    }
    df$gpr <- gv
    df$oil <- ov

    frames[[cc]] <- df
  }
  frames
}

lag_vec <- function(x, k = 1L) {
  if (k == 0L) return(x)
  c(rep(NA, k), head(x, -k))
}

make_design <- function(df, p, common_trim_p = 2L) {
  # Current Wex: 5 star variables + GPR + OIL
  wex <- c(paste0("star_", VARS), "gpr", "oil")

  Z <- data.frame(Intercept = rep(1, nrow(df)))

  # Domestic lags
  for (L in seq_len(p)) {
    for (v in VARS) Z[[paste0(v, "_lag", L)]] <- lag_vec(df[[v]], L)
  }

  # Current foreign/global and q=1 lag
  for (v in wex) Z[[v]] <- df[[v]]
  for (v in wex) Z[[paste0(v, "_lag1")]] <- lag_vec(df[[v]], 1L)

  Y <- as.matrix(df[, VARS, drop = FALSE])

  # Fair p=1 vs p=2 comparison: force common maximum-lag sample.
  keep <- seq_len(nrow(df)) > common_trim_p
  keep <- keep & complete.cases(Z) & complete.cases(Y)

  list(
    X = as.matrix(Z[keep, , drop = FALSE]),
    Y = Y[keep, , drop = FALSE],
    quarter = df$quarter[keep]
  )
}

# ------------------------- ESTIMATION ---------------------------

fit_varx_ols <- function(design, p, country) {
  X <- design$X
  Y <- design$Y
  Tn <- nrow(X)
  m  <- ncol(Y)
  kx <- ncol(X)

  if (Tn <= kx + 5L) stopf("%s p=%d has too few observations (%d) for %d regressors.",
                            country, p, Tn, kx)

  qrX <- qr(X)
  if (qrX$rank < kx) {
    return(list(
      ok = FALSE,
      country = country,
      p = p,
      n = Tn,
      k_reg = kx,
      rank = qrX$rank,
      reason = "rank_deficient"
    ))
  }

  B <- qr.coef(qrX, Y)   # regressors x equations
  E <- Y - X %*% B
  colnames(E) <- VARS

  # Gaussian system likelihood using residual covariance.
  Sigma <- crossprod(E) / Tn
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 1e-12)) {
    # numerical guard; should rarely matter
    ev <- pmax(ev, 1e-12)
  }
  logdet <- sum(log(ev))
  logLik <- -0.5 * Tn * (m * (1 + log(2*pi)) + logdet)

  k_beta  <- m * kx
  k_sigma <- m * (m + 1L) / 2L
  k_total <- k_beta + k_sigma

  AIC <- -2 * logLik + 2 * k_total
  BIC <- -2 * logLik + log(Tn) * k_total

  # Local domestic companion spectral radius.
  A <- vector("list", p)
  for (L in seq_len(p)) {
    AL <- matrix(0, m, m, dimnames = list(VARS, VARS))
    for (j in seq_along(VARS)) {
      rn <- paste0(VARS[j], "_lag", L)
      ri <- match(rn, rownames(B))
      if (is.na(ri)) stopf("Internal error: coefficient %s not found.", rn)
      # B[regressor, equation] -> A[equation, lagged-variable]
      AL[, j] <- B[ri, ]
    }
    A[[L]] <- AL
  }

  if (p == 1L) {
    companion <- A[[1]]
  } else {
    top <- do.call(cbind, A)
    lower <- cbind(diag(m), matrix(0, m, m * (p - 1L)))
    if (p > 2L) {
      stopf("This diagnostic is currently written for p<=2.")
    }
    companion <- rbind(top, lower)
  }
  rho <- max(Mod(eigen(companion, only.values = TRUE)$values))

  # Equation-level Ljung-Box residual serial correlation.
  serial <- lapply(seq_len(m), function(j) {
    bt <- tryCatch(
      Box.test(E[, j], lag = LB_LAG, type = "Ljung-Box", fitdf = 0),
      error = function(e) NULL
    )
    data.frame(
      country = country,
      p = p,
      variable = VARS[j],
      n = Tn,
      test_lag = LB_LAG,
      p_value = if (is.null(bt)) NA_real_ else unname(bt$p.value),
      reject_5pct = if (is.null(bt)) NA else unname(bt$p.value < 0.05)
    )
  })
  serial <- do.call(rbind, serial)

  list(
    ok = TRUE,
    country = country,
    p = p,
    n = Tn,
    k_reg = kx,
    rank = qrX$rank,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    rho = rho,
    stable = is.finite(rho) && rho < STABILITY_CUTOFF,
    borderline = is.finite(rho) && rho >= BORDERLINE_RHO,
    residuals = E,
    coef = B,
    serial = serial,
    serial_reject_count = sum(serial$reject_5pct, na.rm = TRUE),
    serial_reject_share = mean(serial$reject_5pct, na.rm = TRUE)
  )
}

# ---------------------- RECOMMENDATION RULE ---------------------

choose_p <- function(sub) {
  # Transparent pre-specified rule:
  # 1) instability (rho >= 1) disqualifies a candidate;
  # 2) among stable candidates, choose fewer Ljung-Box rejections;
  # 3) if tied, choose lower BIC;
  # 4) if still tied, choose lower p (parsimony).
  ok <- sub[sub$ok & sub$stable, , drop = FALSE]

  if (nrow(ok) == 0L) {
    return(list(p = NA_integer_, reason = "Neither p=1 nor p=2 passes local stability (rho<1)."))
  }

  minrej <- min(ok$serial_reject_count, na.rm = TRUE)
  cand <- ok[ok$serial_reject_count == minrej, , drop = FALSE]

  if (nrow(cand) > 1L) {
    mbic <- min(cand$BIC, na.rm = TRUE)
    cand <- cand[cand$BIC == mbic, , drop = FALSE]
  }
  if (nrow(cand) > 1L) cand <- cand[which.min(cand$p), , drop = FALSE]

  pp <- as.integer(cand$p[1])
  rr <- sprintf(
    "p=%d selected: stable rho=%.6f; Ljung-Box rejects=%d/5; BIC=%.3f.",
    pp, cand$rho[1], cand$serial_reject_count[1], cand$BIC[1]
  )

  # Add why the alternative lost.
  alt <- sub[sub$p != pp, , drop = FALSE]
  if (nrow(alt)) {
    if (!isTRUE(alt$stable[1])) {
      rr <- paste0(rr, sprintf(" Alternative p=%d rejected for instability (rho=%.6f).",
                               alt$p[1], alt$rho[1]))
    } else {
      rr <- paste0(rr, sprintf(" Alternative p=%d: rejects=%d/5, BIC=%.3f, rho=%.6f.",
                               alt$p[1], alt$serial_reject_count[1],
                               alt$BIC[1], alt$rho[1]))
    }
  }

  list(p = pp, reason = rr)
}

# ---------------------------- RUN -------------------------------

msg("Reading macro panel...")
macro_raw <- read_table_auto(MACRO_PATH, MACRO_SHEET)
panel <- extract_macro_long(macro_raw)

msg("Reading trade weights...")
W <- read_weight_matrix(WEIGHT_PATH)

msg("Reading GPR...")
gpr <- read_global_series(
  GPR_PATH,
  value_candidates = c("lngprqmean","lngpr","gpr"),
  exact = GPR_COL_EXACT,
  label = "gpr"
)

msg("Reading Brent Oil...")
oil <- read_global_series(
  OIL_PATH,
  value_candidates = c("lnbrent","logbrent","brentlog","oil","brent"),
  exact = NULL,
  label = "oil"
)

msg("Building country VARX* frames...")
frames <- build_country_frames(panel, W, gpr, oil)

results <- list()
serial_all <- list()
idx <- 1L
sidx <- 1L

for (cc in COUNTRIES) {
  msg("Estimating %s: p=1 vs p=2 ...", cc)
  for (p in P_CANDIDATES) {
    des <- make_design(frames[[cc]], p = p, common_trim_p = max(P_CANDIDATES))
    fit <- fit_varx_ols(des, p, cc)

    if (!fit$ok) {
      results[[idx]] <- data.frame(
        country = cc, p = p, ok = FALSE,
        n = fit$n, k_reg = fit$k_reg, design_rank = fit$rank,
        logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
        rho = NA_real_, stable = FALSE, borderline = NA,
        serial_reject_count = NA_integer_,
        serial_reject_share = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      results[[idx]] <- data.frame(
        country = cc, p = p, ok = TRUE,
        n = fit$n, k_reg = fit$k_reg, design_rank = fit$rank,
        logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC,
        rho = fit$rho, stable = fit$stable, borderline = fit$borderline,
        serial_reject_count = fit$serial_reject_count,
        serial_reject_share = fit$serial_reject_share,
        stringsAsFactors = FALSE
      )
      serial_all[[sidx]] <- fit$serial
      sidx <- sidx + 1L
    }
    idx <- idx + 1L
  }
}

summary_df <- do.call(rbind, results)
serial_df  <- do.call(rbind, serial_all)

# AIC/BIC choices and diagnostic recommendation.
rec_list <- lapply(COUNTRIES, function(cc) {
  sub <- summary_df[summary_df$country == cc, , drop = FALSE]
  stable_sub <- sub[sub$ok & sub$stable, , drop = FALSE]

  p_aic <- if (nrow(stable_sub)) stable_sub$p[which.min(stable_sub$AIC)] else NA_integer_
  p_bic <- if (nrow(stable_sub)) stable_sub$p[which.min(stable_sub$BIC)] else NA_integer_

  ch <- choose_p(sub)

  p1 <- sub[sub$p == 1L, , drop = FALSE]
  p2 <- sub[sub$p == 2L, , drop = FALSE]

  data.frame(
    country = cc,
    p_by_AIC_stable = p_aic,
    p_by_BIC_stable = p_bic,
    recommended_p = ch$p,
    recommendation_reason = ch$reason,
    p1_BIC = p1$BIC,
    p2_BIC = p2$BIC,
    delta_BIC_p2_minus_p1 = p2$BIC - p1$BIC,
    p1_serial_rejects = p1$serial_reject_count,
    p2_serial_rejects = p2$serial_reject_count,
    p1_rho = p1$rho,
    p2_rho = p2$rho,
    p1_stable = p1$stable,
    p2_stable = p2$stable,
    stringsAsFactors = FALSE
  )
})
rec_df <- do.call(rbind, rec_list)

# Save machine-readable lag vector.
lag_vec_out <- setNames(rec_df$recommended_p, rec_df$country)
lag_vec_text <- paste0(
  "COUNTRY_P <- c(",
  paste(sprintf("%s=%sL", names(lag_vec_out),
                ifelse(is.na(lag_vec_out), "NA", lag_vec_out)),
        collapse = ", "),
  ")\n"
)

# --------------------------- OUTPUTS ----------------------------

write.csv(
  summary_df,
  file.path(OUT_DIR, "01_country_p1_p2_model_diagnostics.csv"),
  row.names = FALSE
)

write.csv(
  serial_df,
  file.path(OUT_DIR, "02_equation_residual_serial_correlation.csv"),
  row.names = FALSE
)

write.csv(
  rec_df,
  file.path(OUT_DIR, "03_country_lag_recommendation.csv"),
  row.names = FALSE
)

writeLines(
  lag_vec_text,
  file.path(OUT_DIR, "04_selected_country_lag_vector.R")
)

# Compact audit/readme.
readme <- c(
  "Country-specific lag selection diagnostic",
  "========================================",
  sprintf("Model sample: %s - %s", SAMPLE_START, SAMPLE_END),
  sprintf("Countries: %s", paste(COUNTRIES, collapse = ", ")),
  "Candidate domestic lags: p_i in {1,2}",
  sprintf("Foreign/global lag: q_i = %d", Q_FIXED),
  sprintf("Residual serial-correlation test: Ljung-Box lag %d, 5%% level", LB_LAG),
  sprintf("Local stability rule: domestic companion spectral radius < %.3f", STABILITY_CUTOFF),
  sprintf("Borderline stability flag: rho >= %.3f", BORDERLINE_RHO),
  "",
  "IMPORTANT COMPARABILITY RULE:",
  "p=1 and p=2 are both estimated on the common p=2-trimmed sample.",
  "Therefore both candidates should have the same n (expected ~100 here).",
  "",
  "Recommendation rule (pre-specified):",
  "1. Any rho >= 1 candidate is disqualified.",
  "2. Among stable candidates, choose the p with fewer rejected equation-level Ljung-Box tests.",
  "3. If tied, choose lower BIC.",
  "4. If still tied, choose lower p for parsimony.",
  "",
  "This is a pre-estimation VARX* diagnostic. Final TVP-GVAR posterior stability,",
  "MCMC convergence, weak-exogeneity and global residual-dependence diagnostics",
  "must still be run after the country-specific lag vector is inserted.",
  "",
  "Selected lag vector:",
  lag_vec_text
)
writeLines(readme, file.path(OUT_DIR, "README_country_specific_lag_selection.txt"))

msg("")
msg("DONE.")
msg("Outputs written to: %s", OUT_DIR)
msg("Recommended lag vector:")
cat(lag_vec_text)

# Fail loudly if any country has no stable candidate.
if (anyNA(rec_df$recommended_p)) {
  bad <- rec_df$country[is.na(rec_df$recommended_p)]
  warning(sprintf("No stable p candidate for: %s. Inspect diagnostics before full TVP-GVAR.",
                  paste(bad, collapse = ", ")))
}
