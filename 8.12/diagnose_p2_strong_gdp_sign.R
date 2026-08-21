#!/usr/bin/env Rscript

# Diagnose why +10% global GPR produces positive GDP-growth IRFs in the
# existing P2_strong posterior. NO MCMC re-estimation is performed.

suppressPackageStartupMessages(library(ggplot2))

src <- Sys.getenv("TVPGVAR_GDP_DIAG_SOURCE_DIR", "source-artifact")
out <- Sys.getenv("TVPGVAR_GDP_DIAG_OUT", "gdp-sign-diagnostic")
dates <- trimws(strsplit(
  Sys.getenv("TVPGVAR_GDP_DIAG_DATES",
             "2003Q1,2008Q3,2014Q3,2020Q1,2022Q1,2023Q4"),
  ",", fixed = TRUE
)[[1]])
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_GDP_DIAG_SHOCK_PCT", "10"))
near_unit <- as.numeric(Sys.getenv("TVPGVAR_NEAR_UNIT_THRESHOLD", "0.98"))

dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)

find_one <- function(name) {
  z <- list.files(src, pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"),
                  recursive = TRUE, full.names = TRUE)
  if (length(z) != 1L) stop("Expected exactly one ", name, "; found ", length(z))
  z[[1]]
}

readc <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

load(find_one("predDens_gpr_structural.rda"))
parms <- readc(find_one("stability_profile_parameters.csv"))
stab_source <- readc(find_one("stability_summary_selected_dates.csv"))

stopifnot(exists("predDens"), exists("Data.setup"))
if (nrow(parms) != 1L || parms$profile[[1]] != "P2_strong") {
  stop("Wrong artifact: expected P2_strong.")
}
p <- as.integer(parms$p[[1]])
q <- as.integer(parms$q[[1]])
if (p != 2L || q != 1L) stop("Expected GVAR(2,1), found p=", p, ", q=", q)

xglobal <- Data.setup$bigx
countries <- as.character(Data.setup$countries)
quarters <- as.character(Data.setup$quarters)
Wlist <- lapply(predDens, function(z) z$W)
Alist <- lapply(predDens, function(z) z$ALPHA)
Slist <- lapply(predDens, function(z) z$SIGMApost)

us_i <- match("US", countries)
gpr_i <- match("US_gpr", colnames(xglobal))
if (is.na(us_i) || is.na(gpr_i)) stop("US / US_gpr not found.")

ndraw <- dim(Alist[[1]])[4]
irf_dates <- quarters[-seq_len(p)]
date_idx <- match(dates, irf_dates)
if (anyNA(date_idx)) stop("Some requested dates are not in the posterior path.")

split_alpha <- function(A) {
  tags <- dimnames(A)[[2]]
  get_block <- function(tag) {
    ii <- which(tags == tag)
    if (!length(ii)) stop("Missing coefficient block: ", tag)
    A[, ii, , drop = FALSE]
  }
  lag_tags <- grep("^Ylag[0-9]+$", unique(tags), value = TRUE)
  lag_n <- as.integer(sub("^Ylag", "", lag_tags))
  pp <- max(lag_n)
  list(
    L0 = get_block("Wex"),
    Theta = lapply(seq_len(pp), function(k) get_block(paste0("Ylag", k))),
    Llag = lapply(seq_len(pp), function(k) {
      tag <- paste0("Wexlag", k)
      if (any(tags == tag)) get_block(tag) else NULL
    })
  )
}

align_lag <- function(B, Wi, k, cc) {
  B <- as.matrix(B)
  need <- nrow(Wi) - k
  if (ncol(B) == need) return(B)
  if (cc != "US" && need == 6L && ncol(B) == 5L) {
    return(cbind(B, rep(0, nrow(B))))
  }
  stop("Foreign-lag dimension mismatch for ", cc)
}

stab_cov <- function(S) {
  S <- (as.matrix(S) + t(as.matrix(S))) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  mn <- min(ev)
  if (mn <= 1e-10) S <- S + diag(abs(mn) + 1e-10, nrow(S))
  S
}

rho_companion <- function(F) {
  K <- nrow(F[[1]])
  pp <- length(F)
  if (pp == 1L) return(max(Mod(eigen(F[[1]], only.values = TRUE)$values)))
  top <- do.call(cbind, F)
  lower <- cbind(diag(K * (pp - 1L)),
                 matrix(0, K * (pp - 1L), K))
  max(Mod(eigen(rbind(top, lower), only.values = TRUE)$values))
}

build_state <- function(tt, Adraw, Sdraw) {
  Grow <- vector("list", length(countries))
  Hcountry <- vector("list", length(countries))
  S <- vector("list", length(countries))
  direct_gpr_y <- setNames(rep(NA_real_, length(countries)), countries)

  for (i in seq_along(countries)) {
    cc <- countries[[i]]
    aa <- split_alpha(Adraw[[i]])
    Wi <- as.matrix(Wlist[[i]])
    k <- dim(aa$L0)[3]

    lambda0 <- t(aa$L0[tt, , ])
    A0 <- cbind(diag(k), -lambda0)
    Grow[[i]] <- A0 %*% Wi

    # Non-US contemporaneous Wex order is:
    # foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq, global_gpr.
    # Domestic equation order begins with y.
    if (cc != "US") {
      if (ncol(lambda0) != 6L) stop("Expected 6 current Wex columns for ", cc)
      direct_gpr_y[[cc]] <- lambda0[1, 6]
    }

    Hcountry[[i]] <- lapply(seq_along(aa$Theta), function(kL) {
      theta <- t(aa$Theta[[kL]][tt, , ])
      if (is.null(aa$Llag[[kL]])) {
        lagfull <- matrix(0, nrow = k, ncol = nrow(Wi) - k)
      } else {
        lagfull <- align_lag(t(aa$Llag[[kL]][tt, , ]), Wi, k, cc)
      }
      cbind(theta, lagfull) %*% Wi
    })
    S[[i]] <- Sdraw[[i]][tt, , ]
  }

  G <- do.call(rbind, Grow)
  H <- lapply(seq_len(length(Hcountry[[1]])), function(kL) {
    do.call(rbind, lapply(Hcountry, function(z) z[[kL]]))
  })
  F <- lapply(H, function(h) solve(G, h))

  list(G = G, F = F, S = S, rho = rho_companion(F),
       direct_gpr_y = direct_gpr_y)
}

make_shocks <- function(Sus) {
  Sus <- stab_cov(Sus)

  # Current ordering: GPR, y, dp, r, de, deq
  L <- t(chol(Sus))
  current <- L[, 1]

  # Pure GPR reduced-form innovation.
  pure <- rep(0, 6)
  pure[1] <- L[1, 1]

  # Diagnostic alternative: y, dp, GPR, r, de, deq
  perm <- c(2, 3, 1, 4, 5, 6)
  Lp <- t(chol(Sus[perm, perm, drop = FALSE]))
  gpos <- which(perm == 1)
  sp <- Lp[, gpos]
  slow <- numeric(6)
  slow[perm] <- sp

  list(current = current, pure = pure, slow = slow, L = L)
}

impact_from_us_shock <- function(G, Sblocks, usshock) {
  u <- unlist(lapply(seq_along(Sblocks), function(i) {
    if (i == us_i) usshock else rep(0, nrow(Sblocks[[i]]))
  }), use.names = FALSE)
  imp <- as.numeric(solve(G, u))
  target <- log1p(shock_pct / 100)
  imp <- imp * target / imp[gpr_i]
  names(imp) <- colnames(xglobal)
  imp
}

propagate <- function(F, impact, H = 4L) {
  K <- length(impact)
  ans <- matrix(0, K, H + 1L,
                dimnames = list(names(impact), as.character(0:H)))
  ans[, 1] <- impact
  for (h in seq_len(H)) {
    for (lag in seq_len(min(length(F), h))) {
      ans[, h + 1] <- ans[, h + 1] +
        F[[lag]] %*% ans[, h - lag + 1]
    }
  }
  ans
}

rows_irf <- list()
rows_coef <- list()
rows_us <- list()
rows_stab <- list()
ni <- nc <- nu <- ns <- 0L

for (dd in seq_along(dates)) {
  d <- dates[[dd]]
  tt <- date_idx[[dd]]
  stable_n <- 0L
  near_n <- 0L
  message("Date: ", d)

  for (draw in seq_len(ndraw)) {
    Adraw <- lapply(Alist, function(z) z[, , , draw])
    Sdraw <- lapply(Slist, function(z) z[, , , draw])
    st <- build_state(tt, Adraw, Sdraw)

    if (!is.finite(st$rho) || st$rho >= 1) next
    stable_n <- stable_n + 1L
    if (st$rho >= near_unit) near_n <- near_n + 1L

    sh <- make_shocks(st$S[[us_i]])

    nu <- nu + 1L
    rows_us[[nu]] <- data.frame(
      date = d, draw = draw, rho = st$rho,
      US_y_loading_per_GPR = sh$L[2, 1] / sh$L[1, 1],
      stringsAsFactors = FALSE
    )

    for (cc in setdiff(countries, "US")) {
      nc <- nc + 1L
      rows_coef[[nc]] <- data.frame(
        date = d, draw = draw, country = cc, rho = st$rho,
        contemporaneous_GPRt_to_y = st$direct_gpr_y[[cc]],
        stringsAsFactors = FALSE
      )
    }

    shocks <- list(
      current_GPR_first_Cholesky = sh$current,
      pure_GPR_residual = sh$pure,
      slow_real_y_dp_before_GPR = sh$slow
    )

    for (mode in names(shocks)) {
      imp <- impact_from_us_shock(st$G, st$S, shocks[[mode]])
      path <- propagate(st$F, imp, H = 4L)

      for (cc in countries) {
        yn <- paste0(cc, "_y")
        for (h in 0:4) {
          ni <- ni + 1L
          rows_irf[[ni]] <- data.frame(
            date = d, draw = draw, country = cc,
            mode = mode, horizon = h, rho = st$rho,
            y_response = path[yn, h + 1],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  ns <- ns + 1L
  rows_stab[[ns]] <- data.frame(
    date = d, draws = ndraw,
    stable_draws_recomputed = stable_n,
    near_unit_draws_recomputed = near_n,
    stringsAsFactors = FALSE
  )
}

irf_raw <- do.call(rbind, rows_irf)
coef_raw <- do.call(rbind, rows_coef)
us_raw <- do.call(rbind, rows_us)
stab_re <- do.call(rbind, rows_stab)

# Reproduce source stable counts before interpreting anything.
stab_source$date <- as.character(stab_source$date)
check <- merge(
  stab_re,
  stab_source[, c("date", "stable_draws", "near_unit_draws")],
  by = "date", all.x = TRUE, sort = FALSE
)
check$stable_match <- check$stable_draws_recomputed == as.integer(check$stable_draws)
check$near_unit_match <- check$near_unit_draws_recomputed == as.integer(check$near_unit_draws)
write.csv(check, file.path(out, "01_stability_reconstruction_check.csv"),
          row.names = FALSE)
if (any(!check$stable_match)) {
  stop("Recomputed stable-draw counts do not match the source run.")
}

qv <- function(x, prob) as.numeric(quantile(x, prob, names = FALSE, type = 8))

summ <- function(df, value, groups) {
  key <- interaction(df[groups], drop = TRUE, lex.order = TRUE)
  z <- lapply(split(df, key), function(a) {
    x <- a[[value]]
    x <- x[is.finite(x)]
    first <- a[1, groups, drop = FALSE]
    cbind(first, data.frame(
      draws = length(x),
      median = median(x),
      low68 = qv(x, 0.16),
      high68 = qv(x, 0.84),
      low90 = qv(x, 0.05),
      high90 = qv(x, 0.95),
      positive_share = mean(x > 0),
      stringsAsFactors = FALSE
    ))
  })
  ans <- do.call(rbind, z)
  rownames(ans) <- NULL
  ans
}

gdp <- summ(irf_raw, "y_response",
            c("date", "country", "mode", "horizon"))
gdp_h0 <- gdp[gdp$horizon == 0, ]
coef <- summ(coef_raw, "contemporaneous_GPRt_to_y",
             c("date", "country"))
usload <- summ(us_raw, "US_y_loading_per_GPR", c("date"))

write.csv(gdp, file.path(out, "02_GDP_IRF_identification_comparison.csv"),
          row.names = FALSE)
write.csv(gdp_h0, file.path(out, "03_GDP_h0_identification_comparison.csv"),
          row.names = FALSE)
write.csv(coef, file.path(out, "04_nonUS_current_GPRt_to_GDP_coefficients.csv"),
          row.names = FALSE)
write.csv(usload, file.path(out, "05_US_GPRfirst_direct_GDP_loading.csv"),
          row.names = FALSE)

# Compact h=0 decomposition table.
pick <- function(mode, prefix) {
  z <- gdp_h0[gdp_h0$mode == mode,
              c("date", "country", "median", "positive_share")]
  names(z)[3:4] <- paste0(prefix, c("_median", "_positive_share"))
  z
}
dec <- merge(
  pick("current_GPR_first_Cholesky", "current"),
  pick("pure_GPR_residual", "pure"),
  by = c("date", "country")
)
dec <- merge(
  dec,
  pick("slow_real_y_dp_before_GPR", "slowreal"),
  by = c("date", "country")
)
dec$pure_minus_current <- dec$pure_median - dec$current_median
dec$slowreal_minus_current <- dec$slowreal_median - dec$current_median
write.csv(dec, file.path(out, "06_GDP_h0_source_decomposition.csv"),
          row.names = FALSE)

# Share of non-US countries whose median direct contemporaneous GPR->GDP
# coefficient is positive at each selected date.
share_rows <- lapply(split(coef, coef$date), function(z) {
  data.frame(
    date = z$date[[1]],
    nonUS_countries = nrow(z),
    positive_median_direct_coef = sum(z$median > 0),
    share_positive_median_direct_coef = mean(z$median > 0)
  )
})
shares <- do.call(rbind, share_rows)
rownames(shares) <- NULL
write.csv(shares, file.path(out, "07_direct_GPR_GDP_positive_share.csv"),
          row.names = FALSE)

# Plot 1: identification comparison for key countries, h=0..4.
key <- intersect(c("US", "JP", "CN", "KR", "EA", "SG"), countries)
pd <- gdp[gdp$country %in% key, ]
pd$country <- factor(pd$country, levels = key)
pd$date <- factor(pd$date, levels = dates)

p1 <- ggplot(pd, aes(horizon, median, linetype = mode, group = mode)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_grid(date ~ country, scales = "free_y") +
  labs(
    title = "GDP-growth IRF identification diagnostic",
    subtitle = paste0("+", shock_pct,
                      "% GPR shock; same P2_strong posterior; stable draws only"),
    x = "Horizon (quarters)", y = "GDP-growth response",
    linetype = "Shock identification"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(out, "plots", "GDP_identification_key_countries.png"),
       p1, width = 18, height = 14, dpi = 180)

# Plot 2: direct non-US contemporaneous GPR_t coefficient in GDP equation.
coef$date <- factor(coef$date, levels = dates)
p2 <- ggplot(coef, aes(country, median, ymin = low68, ymax = high68)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(linewidth = 0.35) +
  facet_wrap(~date, ncol = 2, scales = "free_y") +
  labs(
    title = "Non-US GDP equation: contemporaneous global GPR_t coefficient",
    subtitle = "Posterior median and 68% interval; globally stable draws only",
    x = NULL, y = "Coefficient on global GPR_t"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out, "plots", "NonUS_GPRt_to_GDP_coefficients.png"),
       p2, width = 15, height = 10, dpi = 180)

readme <- c(
  "P2_strong GDP SIGN DIAGNOSTIC",
  "=============================",
  "",
  "NO MCMC re-estimation was performed.",
  paste0("Source profile: ", parms$profile[[1]]),
  paste0("Specification: GVAR(", p, ",", q, ")"),
  paste0("Shock: +", shock_pct, "% global GPR"),
  "",
  "Three impact identifications use the SAME posterior coefficients/dynamics:",
  "1) current_GPR_first_Cholesky: current model.",
  "2) pure_GPR_residual: only the US GPR reduced-form residual is shocked.",
  "3) slow_real_y_dp_before_GPR: covariance-ordering diagnostic y,dp,GPR,r,de,deq.",
  "",
  "Interpretation:",
  "- If current GDP is positive but pure-GPR GDP becomes near zero/negative,",
  "  the current US Cholesky covariance loading is a major source.",
  "- If pure-GPR GDP stays positive outside the US, inspect file 04:",
  "  positive contemporaneous global GPR_t coefficients in non-US GDP equations",
  "  are a major source of the h=0 positive response.",
  "- If slow-real materially changes US but not foreign GDP, changing the US",
  "  recursive order alone will not solve the cross-country GDP sign pattern.",
  "",
  "Do not force GDP negative. This diagnostic is for choosing a defensible",
  "timing/identification specification before any re-estimation."
)
writeLines(readme, file.path(out, "README_GDP_SIGN_DIAGNOSTIC.txt"))

cat("\n===== Stability reconstruction =====\n")
print(check)
cat("\n===== US direct GPR-first GDP loading =====\n")
print(usload)
cat("\n===== Non-US positive direct GPR->GDP coefficient share =====\n")
print(shares)
cat("\nGDP sign diagnostic complete: ", normalizePath(out), "\n", sep = "")
