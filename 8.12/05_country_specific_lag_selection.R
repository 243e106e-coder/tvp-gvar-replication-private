#!/usr/bin/env Rscript

# ============================================================
# Country-specific lag selection for GPR + Oil TVP-GVAR (v2)
# Robust input adapter + p_i = 1 vs 2 diagnostic
#
# Keeps the original research design:
#   - countries: 14
#   - domestic variables: y, dp, r, de, deq
#   - p_i in {1,2}
#   - q_i = 1
#   - AIC / BIC
#   - equation-level Ljung-Box
#   - local companion spectral radius
#
# Main v2 changes:
#   1) safe quarter/date parser: non-date text can no longer crash as.Date()
#   2) robust Excel reader for one-row OR two-row/grouped headers
#   3) explicit macro-input audit before estimation
#   4) clear country-by-country variable mapping in Actions log
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

GPR_COL_EXACT <- "LN_GPR_QMEAN"
GDP_DLOG_DIVISOR <- 100

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------ BASIC HELPERS ------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg   <- function(...) cat(sprintf(...), "\n")

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", ifelse(is.na(x), "", x)))
}

trim_chr <- function(x) trimws(as.character(x))

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

safe_date_parse <- function(x) {
  sx <- trimws(as.character(x))
  out <- rep(as.Date(NA), length(sx))

  # ISO / slash / dot formats with explicit format.
  pats <- list(
    list(re = "^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", fmt = "%Y-%m-%d"),
    list(re = "^[0-9]{4}/[0-9]{1,2}/[0-9]{1,2}$", fmt = "%Y/%m/%d"),
    list(re = "^[0-9]{4}\\.[0-9]{1,2}\\.[0-9]{1,2}$", fmt = "%Y.%m.%d")
  )

  for (pp in pats) {
    ii <- which(is.na(out) & grepl(pp$re, sx))
    if (length(ii)) out[ii] <- suppressWarnings(as.Date(sx[ii], format = pp$fmt))
  }

  out
}

quarter_id <- function(x) {
  # Returns integer 4*year + quarter.

  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    yr <- lt$year + 1900L
    q  <- (lt$mon %/% 3L) + 1L
    return(4L * yr + q)
  }

  # Excel serial dates: only accept a plausible date-serial range.
  if (is.numeric(x)) {
    xx <- as.numeric(x)
    finite <- xx[is.finite(xx)]
    if (length(finite) && mean(finite >= 20000 & finite <= 70000) >= 0.8) {
      d <- as.Date(xx, origin = "1899-12-30")
      lt <- as.POSIXlt(d)
      yr <- lt$year + 1900L
      q  <- (lt$mon %/% 3L) + 1L
      return(4L * yr + q)
    }
  }

  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))

  # 2000Q1, 2000-Q1, 2000 Q1, 2000/Q1 etc.
  m <- regexec("^([0-9]{4})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L

  if (any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    q  <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L * yr + q
  }

  # Date-like fallback. Crucially: never call as.Date() on arbitrary text.
  need <- is.na(out)
  if (any(need)) {
    d <- safe_date_parse(sx[need])
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

quarter_parse_rate <- function(z) {
  q <- tryCatch(quarter_id(z), error = function(e) rep(NA_integer_, length(z)))
  if (!length(q)) return(0)
  mean(!is.na(q))
}

detect_date_col <- function(df, label = "table") {
  if (!ncol(df)) stopf("%s has zero columns.", label)

  nn <- norm_name(names(df))
  preferred <- c("quarter","date","time","period","qtr","季度","日期")

  # Exact normalized-name preference.
  hit <- match(norm_name(preferred), nn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) {
    dcol <- names(df)[hit[1]]
    rate <- quarter_parse_rate(df[[dcol]])
    if (rate >= 0.5) {
      msg("[%s] date column by name: %s (parse rate %.1f%%)",
          label, dcol, 100 * rate)
      return(dcol)
    }
  }

  # Safe parse-rate search; arbitrary character columns can no longer crash.
  rates <- vapply(df, quarter_parse_rate, numeric(1))
  rates[!is.finite(rates)] <- 0

  best <- which.max(rates)
  if (length(best) && rates[best] >= 0.8) {
    msg("[%s] date column by parse rate: %s (%.1f%%)",
        label, names(df)[best], 100 * rates[best])
    return(names(df)[best])
  }

  ord <- order(rates, decreasing = TRUE)
  show <- head(ord, min(8L, length(ord)))
  diag_txt <- paste(sprintf("%s=%.2f", names(df)[show], rates[show]), collapse = ", ")
  stopf("Could not detect a reliable quarter/date column in %s. Best parse rates: %s",
        label, diag_txt)
}

find_col_exact_norm <- function(df, candidates) {
  nn <- norm_name(names(df))
  cc <- norm_name(candidates)
  idx <- match(cc, nn, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) names(df)[idx[1]] else NA_character_
}

# --------------------- INPUT TABLE READERS ---------------------

read_table_auto <- function(path, sheet = 1) {
  if (!file.exists(path)) stopf("File not found: %s", path)

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx","xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stopf("Package 'readxl' is required for Excel input.")
    }
    return(as.data.frame(
      readxl::read_excel(path, sheet = sheet, .name_repair = "unique"),
      check.names = FALSE
    ))
  }

  if (ext == "csv") {
    return(read.csv(path, check.names = FALSE))
  }

  if (ext == "rds") return(readRDS(path))

  stopf("Unsupported file type: %s", path)
}

# Aliases for macro input.
ALIASES <- list(
  y = c("y", "gdploglevel", "gdplog", "reallogdp", "gdp"),
  y_dlog = c("gdpdlog", "ydlog", "realgdpdlog", "gdp_dlog"),
  dp = c("dp", "cpidlog", "inflation", "cpiinflation", "cpi_dlog"),
  r  = c("r", "rate", "interestrate", "shorttermrate", "shortrate",
         "short_rate", "interest_rate"),
  de = c("de", "reerdlog", "exchangeratedlog", "fxdlog", "erdlog",
         "reer_dlog", "exchange_rate_dlog"),
  deq = c("deq", "eqdlog", "equitydlog", "stockdlog", "stockret",
          "equityret", "eq_dlog", "equity_dlog")
)

macro_mapping_score <- function(nms) {
  nn <- norm_name(nms)
  score <- 0L
  for (cc in COUNTRIES) {
    for (vv in c(ALIASES$y, ALIASES$y_dlog, ALIASES$dp, ALIASES$r,
                 ALIASES$de, ALIASES$deq)) {
      candidates <- norm_name(c(
        paste0(cc, "_", vv), paste0(vv, "_", cc),
        paste0(cc, vv), paste0(vv, cc)
      ))
      if (any(nn %in% candidates)) score <- score + 1L
    }
  }
  score
}

build_two_row_header <- function(mat, h1, h2) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (h1 > nr || h2 > nr || h1 == h2) return(NULL)

  a <- trim_chr(unlist(mat[h1, , drop = TRUE], use.names = FALSE))
  b <- trim_chr(unlist(mat[h2, , drop = TRUE], use.names = FALSE))

  # Group labels in row h1 (e.g., AU then five blanks) are filled right.
  af <- fill_right(a)

  new_names <- character(nc)
  for (j in seq_len(nc)) {
    aj <- ifelse(is.na(af[j]), "", af[j])
    bj <- ifelse(is.na(b[j]), "", b[j])

    # Date/quarter column often has only one meaningful header.
    if (j == 1L && nzchar(bj) && norm_name(bj) %in% norm_name(c("quarter","date","time","period","qtr"))) {
      new_names[j] <- bj
    } else if (j == 1L && nzchar(aj) && norm_name(aj) %in% norm_name(c("quarter","date","time","period","qtr"))) {
      new_names[j] <- aj
    } else if (nzchar(aj) && nzchar(bj) && norm_name(aj) != norm_name(bj)) {
      new_names[j] <- paste0(aj, "_", bj)
    } else if (nzchar(bj)) {
      new_names[j] <- bj
    } else if (nzchar(aj)) {
      new_names[j] <- aj
    } else {
      new_names[j] <- paste0("V", j)
    }
  }

  dat_start <- max(h1, h2) + 1L
  if (dat_start > nr) return(NULL)

  dat <- mat[dat_start:nr, , drop = FALSE]
  names(dat) <- make.unique(new_names, sep = "_dup")
  as.data.frame(dat, check.names = FALSE)
}

read_macro_table <- function(path, sheet = 1) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stopf("Package 'readxl' is required for Excel input.")
  }

  # First try the ordinary one-header-row interpretation.
  std <- as.data.frame(
    readxl::read_excel(path, sheet = sheet, .name_repair = "unique"),
    check.names = FALSE
  )

  std_score <- macro_mapping_score(names(std))
  std_date <- tryCatch(detect_date_col(std, "macro-standard"), error = function(e) NA_character_)

  msg("[macro] standard-header mapping score: %d", std_score)

  # If the ordinary header already maps many country-variable columns, use it.
  if (!is.na(std_date) && std_score >= 20L) {
    msg("[macro] using standard one-row Excel header.")
    return(std)
  }

  # Otherwise inspect raw cells and try two-row/grouped headers.
  raw <- as.data.frame(
    readxl::read_excel(
      path, sheet = sheet, col_names = FALSE,
      .name_repair = "minimal", guess_max = 200
    ),
    check.names = FALSE
  )

  if (nrow(raw) < 3L) {
    msg("[macro] raw sheet too short for two-row header recovery; using standard read.")
    return(std)
  }

  max_header <- min(5L, nrow(raw) - 1L)
  candidates <- list()
  k <- 1L

  for (h1 in seq_len(max_header)) {
    for (h2 in seq_len(max_header)) {
      if (h1 >= h2) next

      cand <- build_two_row_header(raw, h1, h2)
      if (is.null(cand)) next

      score <- macro_mapping_score(names(cand))
      dcol <- tryCatch(detect_date_col(cand, sprintf("macro-h%d-h%d", h1, h2)),
                       error = function(e) NA_character_)
      drate <- if (!is.na(dcol)) quarter_parse_rate(cand[[dcol]]) else 0

      candidates[[k]] <- list(
        h1 = h1, h2 = h2, data = cand,
        score = score, dcol = dcol, drate = drate
      )
      k <- k + 1L
    }
  }

  if (!length(candidates)) {
    msg("[macro] no usable two-row header candidate; using standard read.")
    return(std)
  }

  total_score <- vapply(
    candidates,
    function(z) z$score + ifelse(z$drate >= 0.8, 1000, 0),
    numeric(1)
  )
  best <- candidates[[which.max(total_score)]]

  msg("[macro] best grouped-header candidate: rows %d + %d; mapping score=%d; date=%s; parse=%.1f%%",
      best$h1, best$h2, best$score,
      ifelse(is.na(best$dcol), "<none>", best$dcol),
      100 * best$drate)

  # Require meaningful improvement over standard interpretation.
  if (!is.na(best$dcol) && best$score > std_score) {
    msg("[macro] using reconstructed two-row/grouped Excel header.")
    return(best$data)
  }

  msg("[macro] reconstructed header did not improve mapping; using standard read.")
  std
}

# --------------------- MACRO PANEL ADAPTER ---------------------

extract_macro_long <- function(raw, countries = COUNTRIES) {
  dcol <- detect_date_col(raw, "macro")
  qid  <- quarter_id(raw[[dcol]])
  prate <- mean(!is.na(qid))

  if (prate < 0.8) {
    stopf("Macro date column '%s' parsed only %.1f%% of rows.",
          dcol, 100 * prate)
  }

  msg("[macro] selected date column: %s", dcol)

  # ---------- Long format ----------
  country_col <- find_col_exact_norm(raw, c("country","economy","unit","iso","code"))

  if (!is.na(country_col)) {
    value_cols <- sapply(c("y","dp","r","de","deq"), function(v) {
      find_col_exact_norm(raw, ALIASES[[v]])
    })

    if (all(!is.na(value_cols))) {
      msg("[macro] detected long format; country column=%s", country_col)

      out <- data.frame(
        qid = qid,
        country = toupper(trimws(as.character(raw[[country_col]]))),
        y   = suppressWarnings(as.numeric(raw[[value_cols["y"]]])),
        dp  = suppressWarnings(as.numeric(raw[[value_cols["dp"]]])),
        r   = suppressWarnings(as.numeric(raw[[value_cols["r"]]])),
        de  = suppressWarnings(as.numeric(raw[[value_cols["de"]]])),
        deq = suppressWarnings(as.numeric(raw[[value_cols["deq"]]])),
        stringsAsFactors = FALSE
      )

      out <- out[out$country %in% countries & !is.na(out$qid), ]
      return(out)
    }
  }

  # ---------- Wide format ----------
  nn <- norm_name(names(raw))
  out_list <- vector("list", length(countries))
  names(out_list) <- countries

  mapping_rows <- list()

  for (cc in countries) {
    find_country_var <- function(alias_vec) {
      cands <- c(
        paste0(cc, "_", alias_vec),
        paste0(alias_vec, "_", cc),
        paste0(cc, alias_vec),
        paste0(alias_vec, cc)
      )
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

    mapping_rows[[cc]] <- data.frame(
      country = cc, y = cy, y_dlog = cyd, dp = cdp,
      r = cr, de = cde, deq = cdeq,
      stringsAsFactors = FALSE
    )

    msg("[macro-map] %s: y=%s | y_dlog=%s | dp=%s | r=%s | de=%s | deq=%s",
        cc,
        ifelse(is.na(cy), "<none>", cy),
        ifelse(is.na(cyd), "<none>", cyd),
        ifelse(is.na(cdp), "<none>", cdp),
        ifelse(is.na(cr), "<none>", cr),
        ifelse(is.na(cde), "<none>", cde),
        ifelse(is.na(cdeq), "<none>", cdeq))

    missing_non_gdp <- c(dp = cdp, r = cr, de = cde, deq = cdeq)

    if (any(is.na(missing_non_gdp))) {
      stopf(
        paste0(
          "Could not map all macro variables for %s.\n",
          "Mapped: y=%s, y_dlog=%s, dp=%s, r=%s, de=%s, deq=%s\n",
          "First available columns: %s"
        ),
        cc, cy, cyd, cdp, cr, cde, cdeq,
        paste(head(names(raw), 80), collapse = ", ")
      )
    }

    if (!is.na(cy)) {
      yval <- suppressWarnings(as.numeric(raw[[cy]]))
    } else if (!is.na(cyd)) {
      gd <- suppressWarnings(as.numeric(raw[[cyd]])) / GDP_DLOG_DIVISOR
      yval <- rep(NA_real_, length(gd))

      ok <- which(!is.na(qid) & !is.na(gd))
      if (!length(ok)) stopf("No usable GDP_DLOG observations for %s.", cc)

      ok <- ok[order(qid[ok])]
      yval[ok[1]] <- log(100)

      if (length(ok) >= 2L) {
        for (kk in 2:length(ok)) {
          cur <- ok[kk]
          prv <- ok[kk - 1L]
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
      dp  = suppressWarnings(as.numeric(raw[[cdp]])),
      r   = suppressWarnings(as.numeric(raw[[cr]])),
      de  = suppressWarnings(as.numeric(raw[[cde]])),
      deq = suppressWarnings(as.numeric(raw[[cdeq]])),
      stringsAsFactors = FALSE
    )
  }

  mapping_df <- do.call(rbind, mapping_rows)
  write.csv(mapping_df,
            file.path(OUT_DIR, "00_macro_column_mapping_audit.csv"),
            row.names = FALSE)

  do.call(rbind, out_list)
}

# ------------------------- WEIGHT MATRIX ------------------------

read_weight_matrix <- function(path, countries = COUNTRIES) {
  wraw <- read_table_auto(path)
  nn <- norm_name(names(wraw))

  ocol <- find_col_exact_norm(wraw, c("origin","from","reporter","countryi","i"))
  pcol <- find_col_exact_norm(wraw, c("partner","to","counterparty","countryj","j"))
  wcol <- find_col_exact_norm(wraw, c("weight","w","tradeweight"))

  if (!is.na(ocol) && !is.na(pcol) && !is.na(wcol)) {
    W <- matrix(0, length(countries), length(countries),
                dimnames = list(countries, countries))

    oo <- toupper(trimws(as.character(wraw[[ocol]])))
    pp <- toupper(trimws(as.character(wraw[[pcol]])))
    ww <- suppressWarnings(as.numeric(wraw[[wcol]]))

    for (k in seq_len(nrow(wraw))) {
      if (oo[k] %in% countries && pp[k] %in% countries && is.finite(ww[k])) {
        W[oo[k], pp[k]] <- ww[k]
      }
    }
  } else {
    first_vals <- toupper(trimws(as.character(wraw[[1]])))

    if (sum(first_vals %in% countries) >= length(countries) - 2L) {
      row_country <- first_vals
      dat <- wraw[, -1, drop = FALSE]
    } else {
      stopf("Weight file format not recognized.")
    }

    col_map <- match(countries, toupper(names(dat)))
    if (anyNA(col_map)) col_map <- match(norm_name(countries), norm_name(names(dat)))

    if (anyNA(col_map)) {
      stopf("Could not map all country columns in weight matrix. Columns: %s",
            paste(names(dat), collapse = ", "))
    }

    W <- matrix(NA_real_, length(countries), length(countries),
                dimnames = list(countries, countries))

    for (cc in countries) {
      ri <- which(row_country == cc)
      if (length(ri) != 1L) stopf("Weight matrix row for %s is missing or duplicated.", cc)
      W[cc, ] <- suppressWarnings(as.numeric(dat[ri, col_map]))
    }
  }

  diag(W) <- 0
  rs <- rowSums(W, na.rm = TRUE)

  if (any(!is.finite(rs) | rs <= 0)) stopf("Invalid weight row sums.")

  W <- W / rs

  if (max(abs(rowSums(W) - 1)) > 1e-8) {
    stopf("Weight normalization failed.")
  }

  W
}

# ---------------------- GLOBAL SERIES READER --------------------

read_global_series <- function(path, value_candidates, exact = NULL, label = "global") {
  df <- read_table_auto(path)
  dcol <- detect_date_col(df, label)
  qid <- quarter_id(df[[dcol]])

  vcol <- NA_character_

  if (!is.null(exact)) vcol <- find_col_exact_norm(df, exact)
  if (is.na(vcol)) vcol <- find_col_exact_norm(df, value_candidates)

  if (is.na(vcol)) {
    nn <- norm_name(names(df))
    didx <- match(norm_name(dcol), nn)
    keep <- setdiff(seq_along(nn), didx)
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

  msg("[%s] value column: %s", label, vcol)

  out <- data.frame(
    qid = qid,
    value = suppressWarnings(as.numeric(df[[vcol]]))
  )

  out <- out[!is.na(out$qid), ]
  out <- out[order(out$qid), ]
  out <- out[!duplicated(out$qid), ]

  names(out)[2] <- label
  attr(out, "value_column") <- vcol
  out
}

# ---------------------- PANEL CONSTRUCTION ----------------------

build_country_frames <- function(panel_long, W, gpr, oil) {
  lo <- quarter_id(SAMPLE_START)
  hi <- quarter_id(SAMPLE_END)

  panel_long <- panel_long[
    !is.na(panel_long$qid) &
      panel_long$qid >= lo &
      panel_long$qid <= hi,
    , drop = FALSE
  ]

  qs <- sort(unique(panel_long$qid))
  expected <- seq(lo, hi)

  if (!all(expected %in% qs)) {
    miss <- quarter_label(setdiff(expected, qs))
    stopf("Macro panel is missing model quarters: %s", paste(miss, collapse = ", "))
  }

  frames <- vector("list", length(COUNTRIES))
  names(frames) <- COUNTRIES

  Xmats <- list()

  for (v in VARS) {
    M <- matrix(
      NA_real_, length(expected), length(COUNTRIES),
      dimnames = list(quarter_label(expected), COUNTRIES)
    )

    for (j in seq_along(COUNTRIES)) {
      cc <- COUNTRIES[j]
      d <- panel_long[panel_long$country == cc, c("qid", v)]
      m <- match(expected, d$qid)
      M[, j] <- as.numeric(d[[v]][m])
    }

    if (anyNA(M)) {
      bad <- which(is.na(M), arr.ind = TRUE)
      ex <- head(bad, 12)
      txt <- apply(ex, 1, function(z) {
        sprintf("%s:%s", rownames(M)[z[1]], colnames(M)[z[2]])
      })
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
  wex <- c(paste0("star_", VARS), "gpr", "oil")

  Z <- data.frame(Intercept = rep(1, nrow(df)))

  for (L in seq_len(p)) {
    for (v in VARS) {
      Z[[paste0(v, "_lag", L)]] <- lag_vec(df[[v]], L)
    }
  }

  for (v in wex) Z[[v]] <- df[[v]]
  for (v in wex) Z[[paste0(v, "_lag1")]] <- lag_vec(df[[v]], 1L)

  Y <- as.matrix(df[, VARS, drop = FALSE])

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

  if (Tn <= kx + 5L) {
    stopf("%s p=%d has too few observations (%d) for %d regressors.",
          country, p, Tn, kx)
  }

  qrX <- qr(X)

  if (qrX$rank < kx) {
    return(list(
      ok = FALSE, country = country, p = p,
      n = Tn, k_reg = kx, rank = qrX$rank,
      reason = "rank_deficient"
    ))
  }

  B <- qr.coef(qrX, Y)
  E <- Y - X %*% B
  colnames(E) <- VARS

  Sigma <- crossprod(E) / Tn
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 1e-12)

  logdet <- sum(log(ev))
  logLik <- -0.5 * Tn * (m * (1 + log(2*pi)) + logdet)

  k_beta  <- m * kx
  k_sigma <- m * (m + 1L) / 2L
  k_total <- k_beta + k_sigma

  AIC <- -2 * logLik + 2 * k_total
  BIC <- -2 * logLik + log(Tn) * k_total

  A <- vector("list", p)

  for (L in seq_len(p)) {
    AL <- matrix(0, m, m, dimnames = list(VARS, VARS))

    for (j in seq_along(VARS)) {
      rn <- paste0(VARS[j], "_lag", L)
      ri <- match(rn, rownames(B))
      if (is.na(ri)) stopf("Internal error: coefficient %s not found.", rn)
      AL[, j] <- B[ri, ]
    }

    A[[L]] <- AL
  }

  if (p == 1L) {
    companion <- A[[1]]
  } else if (p == 2L) {
    top <- do.call(cbind, A)
    lower <- cbind(diag(m), matrix(0, m, m))
    companion <- rbind(top, lower)
  } else {
    stopf("This diagnostic currently supports p<=2.")
  }

  rho <- max(Mod(eigen(companion, only.values = TRUE)$values))

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
      reject_5pct = if (is.null(bt)) NA else unname(bt$p.value < 0.05),
      stringsAsFactors = FALSE
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
  ok <- sub[sub$ok & sub$stable, , drop = FALSE]

  if (nrow(ok) == 0L) {
    return(list(
      p = NA_integer_,
      reason = "Neither p=1 nor p=2 passes local stability (rho<1)."
    ))
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

  alt <- sub[sub$p != pp, , drop = FALSE]

  if (nrow(alt)) {
    if (!isTRUE(alt$stable[1])) {
      rr <- paste0(
        rr,
        sprintf(" Alternative p=%d rejected for instability (rho=%.6f).",
                alt$p[1], alt$rho[1])
      )
    } else {
      rr <- paste0(
        rr,
        sprintf(" Alternative p=%d: rejects=%d/5, BIC=%.3f, rho=%.6f.",
                alt$p[1], alt$serial_reject_count[1],
                alt$BIC[1], alt$rho[1])
      )
    }
  }

  list(p = pp, reason = rr)
}

# ---------------------------- RUN -------------------------------

msg("Reading macro panel...")
macro_raw <- read_macro_table(MACRO_PATH, MACRO_SHEET)

msg("[macro] dimensions after header normalization: %d rows x %d columns",
    nrow(macro_raw), ncol(macro_raw))
msg("[macro] first columns: %s",
    paste(head(names(macro_raw), 30), collapse = " | "))

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
    des <- make_design(
      frames[[cc]],
      p = p,
      common_trim_p = max(P_CANDIDATES)
    )

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

if (length(serial_all)) {
  serial_df <- do.call(rbind, serial_all)
} else {
  serial_df <- data.frame()
}

rec_list <- lapply(COUNTRIES, function(cc) {
  sub <- summary_df[summary_df$country == cc, , drop = FALSE]
  stable_sub <- sub[sub$ok & sub$stable, , drop = FALSE]

  p_aic <- if (nrow(stable_sub)) {
    stable_sub$p[which.min(stable_sub$AIC)]
  } else NA_integer_

  p_bic <- if (nrow(stable_sub)) {
    stable_sub$p[which.min(stable_sub$BIC)]
  } else NA_integer_

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

lag_vec_out <- setNames(rec_df$recommended_p, rec_df$country)

lag_vec_text <- paste0(
  "COUNTRY_P <- c(",
  paste(
    sprintf(
      "%s=%sL",
      names(lag_vec_out),
      ifelse(is.na(lag_vec_out), "NA", lag_vec_out)
    ),
    collapse = ", "
  ),
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

readme <- c(
  "Country-specific lag selection diagnostic (v2)",
  "==============================================",
  sprintf("Model sample: %s - %s", SAMPLE_START, SAMPLE_END),
  sprintf("Countries: %s", paste(COUNTRIES, collapse = ", ")),
  "Candidate domestic lags: p_i in {1,2}",
  sprintf("Foreign/global lag: q_i = %d", Q_FIXED),
  sprintf("Residual serial-correlation test: Ljung-Box lag %d, 5%% level", LB_LAG),
  sprintf("Local stability rule: domestic companion spectral radius < %.3f", STABILITY_CUTOFF),
  sprintf("Borderline stability flag: rho >= %.3f", BORDERLINE_RHO),
  "",
  "V2 INPUT SAFETY:",
  "- non-date character columns are never passed blindly to as.Date()",
  "- one-row and two-row/grouped Excel headers are both attempted",
  "- macro date-column and country-variable mappings are printed to the Actions log",
  "- 00_macro_column_mapping_audit.csv records the resolved wide-format mapping",
  "",
  "IMPORTANT COMPARABILITY RULE:",
  "p=1 and p=2 are both estimated on the common p=2-trimmed sample.",
  "",
  "Recommendation rule:",
  "1. Any rho >= 1 candidate is disqualified.",
  "2. Among stable candidates, choose the p with fewer rejected equation-level Ljung-Box tests.",
  "3. If tied, choose lower BIC.",
  "4. If still tied, choose lower p for parsimony.",
  "",
  "This remains a pre-estimation VARX* diagnostic.",
  "Final TVP-GVAR posterior stability, MCMC convergence, weak-exogeneity",
  "and global residual-dependence diagnostics still need to be run.",
  "",
  "Selected lag vector:",
  lag_vec_text
)

writeLines(
  readme,
  file.path(OUT_DIR, "README_country_specific_lag_selection.txt")
)

msg("")
msg("DONE.")
msg("Outputs written to: %s", OUT_DIR)
msg("Recommended lag vector:")
cat(lag_vec_text)

if (anyNA(rec_df$recommended_p)) {
  bad <- rec_df$country[is.na(rec_df$recommended_p)]
  warning(sprintf(
    "No stable p candidate for: %s. Inspect diagnostics before full TVP-GVAR.",
    paste(bad, collapse = ", ")
  ))
}
