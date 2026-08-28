#!/usr/bin/env Rscript

# ============================================================
# CN + ZA local instability diagnostic (v1.1)
# JP as borderline control
#
# Purpose:
#   Diagnose why local companion spectral radius is >= 1
#   for CN / ZA under p=1 and p=2, without changing:
#   - macro data
#   - trade weights
#   - q=1
#   - GPR / Brent specification
#   - sample
#
# Inputs:
#   8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx
#   8.12/Trade_Weights_14_Economies_2000_2014.csv
#   8.12/gpr_quarterly_processed.csv
#   8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv
#
# Outputs:
#   8.12/cn_za_instability_diagnostic/
#     01_local_stability_summary.csv
#     02_companion_eigenvalues.csv
#     03_domestic_lag_coefficient_norms.csv
#     04_equation_persistence_summary.csv
#     05_variable_univariate_persistence.csv
#     06_counterfactual_zero_block_stability.csv
#     07_suspect_ranking.csv
#     README_cn_za_instability_diagnostic.txt
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

MACRO_PATH  <- "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"
MACRO_SHEET <- "MODEL_WIDE_ALL"
WEIGHT_PATH <- "8.12/Trade_Weights_14_Economies_2000_2014.csv"
GPR_PATH    <- "8.12/gpr_quarterly_processed.csv"
OIL_PATH    <- "8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"

OUT_DIR <- "8.12/cn_za_instability_diagnostic"

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
TARGETS <- c("CN","ZA","JP")
VARS <- c("y","dp","r","de","deq")

SAMPLE_START <- "2000Q2"
SAMPLE_END   <- "2025Q3"

P_SET <- c(1L,2L)
Q_FIXED <- 1L

GPR_COL_EXACT <- "LN_GPR_QMEAN"
GDP_DLOG_DIVISOR <- 100

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
msg <- function(...) cat(sprintf(...), "\n")

quarter_id <- function(x) {
  sx <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))
  m <- regexec("^([0-9]{4})[^0-9]*Q([1-4])$", sx)
  mm <- regmatches(sx, m)
  ok <- lengths(mm) == 3L
  if(any(ok)) {
    yr <- as.integer(vapply(mm[ok], `[`, character(1), 2))
    qq <- as.integer(vapply(mm[ok], `[`, character(1), 3))
    out[ok] <- 4L * yr + qq
  }
  out
}

quarter_label <- function(qid) {
  yr <- (qid - 1L) %/% 4L
  qq <- qid - 4L * yr
  sprintf("%dQ%d", yr, qq)
}

norm_col <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]+", "_", trimws(as.character(x))))
  gsub("^_+|_+$", "", z)
}

read_macro <- function() {
  if(!requireNamespace("readxl", quietly = TRUE)) stopf("Package readxl is required")
  if(!file.exists(MACRO_PATH)) stopf("Macro workbook not found")

  d <- as.data.frame(
    readxl::read_excel(
      MACRO_PATH,
      sheet = MACRO_SHEET,
      col_names = TRUE,
      .name_repair = "minimal",
      guess_max = 2000
    ),
    check.names = FALSE
  )

  nms <- names(d)
  nn <- norm_col(nms)

  qr <- vapply(d, function(z) mean(!is.na(quarter_id(z))), numeric(1))
  qr[!is.finite(qr)] <- 0
  date_col <- which.max(qr)
  if(length(date_col) != 1L || qr[date_col] < 0.80)
    stopf("Could not detect quarter column in %s", MACRO_SHEET)

  qid <- quarter_id(d[[date_col]])

  out <- list()
  for(cc in COUNTRIES) {
    find_one <- function(patterns) {
      pat_norm <- toupper(patterns)
      hits <- integer()
      for(j in seq_along(nn)) {
        z <- nn[j]
        has_cc <- grepl(paste0("^", cc, "_"), z) ||
                  grepl(paste0("_", cc, "$"), z)
        if(!has_cc) next
        if(any(vapply(pat_norm, function(p)
          grepl(paste0("(^|_)", p, "($|_)"), z), logical(1)))) {
          hits <- c(hits, j)
        }
      }
      if(length(hits)) hits[1] else NA_integer_
    }

    jy  <- find_one(c("GDP_DLOG","GDPDLOG"))
    jdp <- find_one(c("CPI_DLOG","CPIDLOG"))
    jr  <- find_one(c("RATE_LEVEL","RATELEVEL"))
    jde <- find_one(c("REER_DLOG","REERDLOG"))
    jeq <- find_one(c("EQ_RETURN","EQRETURN"))

    if(any(is.na(c(jy,jdp,jr,jde,jeq))))
      stopf("Incomplete macro mapping for %s", cc)

    gd <- suppressWarnings(as.numeric(d[[jy]])) / GDP_DLOG_DIVISOR
    y <- rep(NA_real_, length(gd))
    valid <- is.finite(gd) & !is.na(qid)
    if(any(valid)) {
      idx <- which(valid)
      grp <- cumsum(c(1L, diff(idx) > 1L))
      runs <- split(idx, grp)
      for(ii in runs) y[ii] <- log(100) + cumsum(gd[ii])
    }

    out[[cc]] <- data.frame(
      qid = qid,
      country = cc,
      y = y,
      dp = suppressWarnings(as.numeric(d[[jdp]])),
      r  = suppressWarnings(as.numeric(d[[jr]])),
      de = suppressWarnings(as.numeric(d[[jde]])),
      deq = suppressWarnings(as.numeric(d[[jeq]]))
    )
  }

  ans <- do.call(rbind, out)
  ans <- ans[!is.na(ans$qid), ]
  ans
}

read_weights <- function() {
  if(!file.exists(WEIGHT_PATH)) stopf("Weights file not found")
  w0 <- read.csv(WEIGHT_PATH, check.names = FALSE)
  first <- toupper(trimws(as.character(w0[[1]])))

  if(all(COUNTRIES %in% first)) {
    rownames(w0) <- first
    w0 <- w0[, -1, drop = FALSE]
  }

  cn <- toupper(trimws(names(w0)))
  if(!all(COUNTRIES %in% cn)) stopf("Weight columns incomplete")
  w0 <- w0[, match(COUNTRIES, cn), drop = FALSE]

  if(is.null(rownames(w0)) || !all(COUNTRIES %in% toupper(rownames(w0)))) {
    if(nrow(w0) != length(COUNTRIES)) stopf("Weights must have 14 rows")
    rownames(w0) <- COUNTRIES
  } else {
    w0 <- w0[match(COUNTRIES, toupper(rownames(w0))), , drop = FALSE]
  }

  W <- as.matrix(data.frame(lapply(w0, as.numeric), check.names = FALSE))
  rownames(W) <- COUNTRIES
  colnames(W) <- COUNTRIES
  diag(W) <- 0

  rs <- rowSums(W, na.rm = TRUE)
  if(any(!is.finite(rs) | rs <= 0)) stopf("Invalid weight row sum")
  W / rs
}

read_global_series <- function(path, exact = NULL, label = "global") {
  if(!file.exists(path)) stopf("%s file not found: %s", label, path)
  d <- read.csv(path, check.names = FALSE)

  q_rates <- vapply(d, function(z) mean(!is.na(quarter_id(z))), numeric(1))
  q_rates[!is.finite(q_rates)] <- 0
  dc <- which.max(q_rates)
  if(q_rates[dc] < 0.50) stopf("Could not detect quarter column in %s", label)

  qid <- quarter_id(d[[dc]])

  vc <- NA_integer_
  if(!is.null(exact) && exact %in% names(d)) vc <- match(exact, names(d))

  if(is.na(vc)) {
    cand <- setdiff(seq_along(d), dc)
    rr <- vapply(cand, function(j) mean(is.finite(suppressWarnings(as.numeric(d[[j]])))), numeric(1))
    if(length(rr) && max(rr) >= 0.80) vc <- cand[which.max(rr)]
  }

  if(is.na(vc)) stopf("Could not identify value column in %s", label)

  z <- data.frame(
    qid = qid,
    value = suppressWarnings(as.numeric(d[[vc]]))
  )
  z <- z[!is.na(z$qid) & is.finite(z$value), ]
  aggregate(value ~ qid, data = z, FUN = mean)
}

lag_vec <- function(x, L) {
  out <- rep(NA_real_, length(x))
  if(L < length(x)) out[(L+1L):length(x)] <- x[1L:(length(x)-L)]
  out
}

make_panel <- function(macro, W, gpr, oil) {
  q0 <- quarter_id(SAMPLE_START)
  q1 <- quarter_id(SAMPLE_END)

  byc <- split(macro, macro$country)

  all_q <- Reduce(intersect, lapply(byc[COUNTRIES], function(z) z$qid))
  all_q <- sort(all_q[all_q >= q0 & all_q <= q1])

  if(length(all_q) < 20L) stopf("Common sample too short")

  arr <- array(
    NA_real_,
    dim = c(length(all_q), length(COUNTRIES), length(VARS)),
    dimnames = list(quarter_label(all_q), COUNTRIES, VARS)
  )

  for(cc in COUNTRIES) {
    z <- byc[[cc]]
    m <- match(all_q, z$qid)
    arr[, cc, ] <- as.matrix(z[m, VARS])
  }

  star <- array(
    NA_real_,
    dim = dim(arr),
    dimnames = dimnames(arr)
  )

  for(i in seq_along(COUNTRIES)) {
    cc <- COUNTRIES[i]
    for(v in seq_along(VARS)) {
      star[, i, v] <- as.numeric(arr[, , v] %*% W[cc, ])
    }
  }

  gg <- gpr$value[match(all_q, gpr$qid)]
  oo <- oil$value[match(all_q, oil$qid)]

  list(qid = all_q, X = arr, Xstar = star, gpr = gg, oil = oo)
}

fit_local_varx <- function(panel, cc, p) {
  i <- match(cc, COUNTRIES)
  Y <- panel$X[, i, , drop = FALSE][,1,]
  Xs <- panel$Xstar[, i, , drop = FALSE][,1,]

  Tn <- nrow(Y)
  maxlag <- max(p, Q_FIXED, 1L)

  rows <- (maxlag + 1L):Tn

  D <- data.frame(const = rep(1, length(rows)))

  # Domestic lags
  for(L in seq_len(p)) {
    for(v in seq_along(VARS)) {
      D[[paste0(VARS[v], "_L", L)]] <- lag_vec(Y[,v], L)[rows]
    }
  }

  # Foreign q=1
  for(v in seq_along(VARS)) {
    D[[paste0(VARS[v], "_star_L1")]] <- lag_vec(Xs[,v], 1L)[rows]
  }

  # Globals current + lag 1
  D$gpr_0 <- panel$gpr[rows]
  D$gpr_L1 <- lag_vec(panel$gpr, 1L)[rows]
  D$oil_0 <- panel$oil[rows]
  D$oil_L1 <- lag_vec(panel$oil, 1L)[rows]

  YY <- Y[rows, , drop = FALSE]

  ok <- complete.cases(D) & complete.cases(YY)
  D <- D[ok, , drop = FALSE]
  YY <- YY[ok, , drop = FALSE]

  Xmat <- as.matrix(D)
  B <- matrix(NA_real_, nrow = ncol(Xmat), ncol = ncol(YY),
              dimnames = list(colnames(Xmat), VARS))

  res <- matrix(NA_real_, nrow = nrow(YY), ncol = ncol(YY),
                dimnames = list(NULL, VARS))

  for(j in seq_along(VARS)) {
    fit <- lm.fit(Xmat, YY[,j])
    B[,j] <- fit$coefficients
    res[,j] <- fit$residuals
  }

  # A_l: equation rows x lagged variable columns
  A <- vector("list", p)
  for(L in seq_len(p)) {
    Am <- matrix(0, length(VARS), length(VARS),
                 dimnames = list(VARS, VARS))
    for(eq in VARS) {
      for(v in VARS) {
        rn <- paste0(v, "_L", L)
        Am[eq, v] <- B[rn, eq]
      }
    }
    A[[L]] <- Am
  }

  list(
    B = B,
    residuals = res,
    A = A,
    n = nrow(YY),
    Y = YY,
    X = Xmat,
    rows = rows[ok]
  )
}

companion_matrix <- function(A) {
  p <- length(A)
  k <- nrow(A[[1]])
  if(p == 1L) return(A[[1]])

  C <- matrix(0, k*p, k*p)
  C[1:k, 1:(k*p)] <- do.call(cbind, A)
  C[(k+1):(k*p), 1:(k*(p-1))] <- diag(k*(p-1))
  C
}

spectral_radius <- function(A) {
  C <- companion_matrix(A)
  ev <- eigen(C, only.values = TRUE)$values
  max(Mod(ev))
}

safe_ar1 <- function(x) {
  x <- as.numeric(x)
  y <- x[-1]
  z <- x[-length(x)]
  ok <- is.finite(y) & is.finite(z)
  if(sum(ok) < 10L) return(c(phi=NA, n=sum(ok)))
  fit <- lm(y[ok] ~ z[ok])
  c(phi = unname(coef(fit)[2]), n = sum(ok))
}

zero_block_counterfactuals <- function(fit, cc, p) {
  base_rho <- spectral_radius(fit$A)
  rows <- list()

  add_row <- function(label, Anew, block_type, lag = NA_integer_, variable = NA_character_) {
    rho <- spectral_radius(Anew)
    data.frame(
      country = cc,
      p = p,
      block_type = block_type,
      lag = lag,
      variable = variable,
      counterfactual = label,
      rho_base = base_rho,
      rho_cf = rho,
      delta_rho = rho - base_rho,
      improvement = base_rho - rho,
      stringsAsFactors = FALSE
    )
  }

  # Zero whole lag matrices
  for(L in seq_len(p)) {
    A2 <- lapply(fit$A, function(x) x)
    A2[[L]][,] <- 0
    rows[[length(rows)+1L]] <- add_row(
      paste0("zero_all_domestic_L", L), A2, "whole_lag", lag = L
    )
  }

  # Zero each lagged-variable column across all equations
  for(L in seq_len(p)) {
    for(v in VARS) {
      A2 <- lapply(fit$A, function(x) x)
      A2[[L]][, v] <- 0
      rows[[length(rows)+1L]] <- add_row(
        paste0("zero_", v, "_L", L, "_all_eq"),
        A2, "variable_column", lag = L, variable = v
      )
    }
  }

  # Zero each equation row at each lag
  for(L in seq_len(p)) {
    for(eq in VARS) {
      A2 <- lapply(fit$A, function(x) x)
      A2[[L]][eq, ] <- 0
      rows[[length(rows)+1L]] <- add_row(
        paste0("zero_equation_", eq, "_L", L),
        A2, "equation_row", lag = L, variable = eq
      )
    }
  }

  do.call(rbind, rows)
}

macro <- read_macro()
W <- read_weights()
gpr <- read_global_series(GPR_PATH, exact = GPR_COL_EXACT, label = "GPR")
oil <- read_global_series(OIL_PATH, label = "Brent")
panel <- make_panel(macro, W, gpr, oil)

msg("Common sample: %s to %s (%d quarters)",
    quarter_label(min(panel$qid)),
    quarter_label(max(panel$qid)),
    length(panel$qid))

stability_rows <- list()
eig_rows <- list()
coef_rows <- list()
eq_rows <- list()
univ_rows <- list()
cf_rows <- list()

for(cc in TARGETS) {
  for(p in P_SET) {
    msg("Fitting %s p=%d ...", cc, p)
    fit <- fit_local_varx(panel, cc, p)
    C <- companion_matrix(fit$A)
    ev <- eigen(C, only.values = TRUE)$values
    ord <- order(Mod(ev), decreasing = TRUE)

    stability_rows[[length(stability_rows)+1L]] <- data.frame(
      country = cc,
      p = p,
      nobs = fit$n,
      spectral_radius = max(Mod(ev)),
      stable = max(Mod(ev)) < 1,
      borderline = max(Mod(ev)) >= 0.98 & max(Mod(ev)) < 1,
      stringsAsFactors = FALSE
    )

    for(rank in seq_along(ord)) {
      e <- ev[ord[rank]]
      eig_rows[[length(eig_rows)+1L]] <- data.frame(
        country = cc,
        p = p,
        rank = rank,
        real = Re(e),
        imag = Im(e),
        modulus = Mod(e),
        angle = Arg(e),
        stringsAsFactors = FALSE
      )
    }

    for(L in seq_len(p)) {
      A <- fit$A[[L]]

      for(v in VARS) {
        coef_rows[[length(coef_rows)+1L]] <- data.frame(
          country = cc,
          p = p,
          lag = L,
          lagged_variable = v,
          column_l1 = sum(abs(A[,v]), na.rm = TRUE),
          column_l2 = sqrt(sum(A[,v]^2, na.rm = TRUE)),
          column_maxabs = max(abs(A[,v]), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }

      for(eq in VARS) {
        eq_rows[[length(eq_rows)+1L]] <- data.frame(
          country = cc,
          p = p,
          lag = L,
          equation = eq,
          own_lag_coef = A[eq,eq],
          row_l1 = sum(abs(A[eq,]), na.rm = TRUE),
          row_l2 = sqrt(sum(A[eq,]^2, na.rm = TRUE)),
          row_maxabs = max(abs(A[eq,]), na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }

    cf_rows[[length(cf_rows)+1L]] <- zero_block_counterfactuals(fit, cc, p)
  }

  # Univariate persistence from the common panel data
  i <- match(cc, COUNTRIES)
  for(v in VARS) {
    x <- panel$X[,i,v]
    ar <- safe_ar1(x)
    univ_rows[[length(univ_rows)+1L]] <- data.frame(
      country = cc,
      variable = v,
      ar1 = ar["phi"],
      n = ar["n"],
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

stability <- do.call(rbind, stability_rows)
eigs <- do.call(rbind, eig_rows)
coefnorm <- do.call(rbind, coef_rows)
eqsum <- do.call(rbind, eq_rows)
univ <- do.call(rbind, univ_rows)
cf <- do.call(rbind, cf_rows)

# Suspect ranking:
# positive improvement means zeroing the block lowers rho.
sus <- cf[cf$block_type == "variable_column", ]
sus <- sus[order(sus$country, sus$p, -sus$improvement), ]
sus$rank_within_model <- ave(
  -sus$improvement,
  interaction(sus$country, sus$p),
  FUN = function(x) rank(x, ties.method = "first")
)
sus <- sus[sus$rank_within_model <= 5, ]

write.csv(stability,
          file.path(OUT_DIR, "01_local_stability_summary.csv"),
          row.names = FALSE)

write.csv(eigs,
          file.path(OUT_DIR, "02_companion_eigenvalues.csv"),
          row.names = FALSE)

write.csv(coefnorm,
          file.path(OUT_DIR, "03_domestic_lag_coefficient_norms.csv"),
          row.names = FALSE)

write.csv(eqsum,
          file.path(OUT_DIR, "04_equation_persistence_summary.csv"),
          row.names = FALSE)

write.csv(univ,
          file.path(OUT_DIR, "05_variable_univariate_persistence.csv"),
          row.names = FALSE)

write.csv(cf,
          file.path(OUT_DIR, "06_counterfactual_zero_block_stability.csv"),
          row.names = FALSE)

write.csv(sus,
          file.path(OUT_DIR, "07_suspect_ranking.csv"),
          row.names = FALSE)

readme <- c(
  "CN + ZA local instability diagnostic",
  "",
  sprintf("Sample: %s to %s (%d quarters)",
          quarter_label(min(panel$qid)),
          quarter_label(max(panel$qid)),
          length(panel$qid)),
  "Targets: CN, ZA; JP is a borderline control.",
  "Domestic lags compared: p=1 and p=2.",
  "Foreign lag q=1 is unchanged.",
  "GPR and Brent specification is unchanged.",
  "",
  "Interpretation:",
  "01_local_stability_summary.csv:",
  "  baseline local companion spectral radius.",
  "",
  "02_companion_eigenvalues.csv:",
  "  full eigenvalue spectrum sorted by modulus.",
  "",
  "03_domestic_lag_coefficient_norms.csv:",
  "  magnitude of each lagged domestic variable's coefficient column.",
  "",
  "04_equation_persistence_summary.csv:",
  "  own-lag and total coefficient magnitude by equation.",
  "",
  "05_variable_univariate_persistence.csv:",
  "  simple AR(1) persistence and scale diagnostics for each variable.",
  "",
  "06_counterfactual_zero_block_stability.csv:",
  "  diagnostic only. Recomputes rho after zeroing a domestic lag block.",
  "  This is NOT a recommended model specification and is not causal.",
  "",
  "07_suspect_ranking.csv:",
  "  ranks lagged variables whose removal lowers rho the most.",
  "  Large positive 'improvement' identifies blocks most associated with",
  "  local instability and should be inspected in the raw/processed data.",
  "",
  "Do not change weights, oil, identification, or lag structure based only",
  "on this diagnostic. First determine whether the instability is driven",
  "by data transformation, extreme persistence, or a particular equation."
)

writeLines(readme,
           file.path(OUT_DIR, "README_cn_za_instability_diagnostic.txt"))

msg("")
msg("=== Stability summary ===")
print(stability, row.names = FALSE)

msg("")
msg("=== Top suspect variable blocks ===")
print(sus[, c("country","p","lag","variable","rho_base","rho_cf","improvement","rank_within_model")],
      row.names = FALSE)

msg("")
msg("Done. Outputs: %s", OUT_DIR)
