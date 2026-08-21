# =============================================================================
# VIX-aware GVAR data preparation
#
# Global ordering:
#   US:     GPR, VIX, y, dp, r, de, deq
#   non-US: y, dp, r, de, deq
#
# Non-US foreign/global block:
#   foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq,
#   global_gpr, global_vix
# =============================================================================

prepare_gvar_data <- function(data_file = "data/model_input.csv",
                              weight_file = "data/trade_weights.csv") {

  countries <- c(
    "AU","BR","CA","CH","CN","EA","UK",
    "JP","KR","NO","SG","TR","US","ZA"
  )

  macro_vars <- c("y", "dp", "r", "de", "deq")

  dat <- read.csv(
    data_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Quarter" %in% names(dat)) {
    stop("model_input.csv must contain Quarter.")
  }

  if (anyDuplicated(dat$Quarter)) {
    stop("Duplicate quarters in model_input.csv.")
  }

  expected_cols <- unlist(
    lapply(countries, function(cc) {
      own <- paste0(cc, "_", macro_vars)

      if (cc == "US") {
        c("US_gpr", "US_vix", own)
      } else {
        own
      }
    }),
    use.names = FALSE
  )

  missing_cols <- setdiff(expected_cols, names(dat))

  if (length(missing_cols)) {
    stop(
      "model_input.csv is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dat <- dat[, c("Quarter", expected_cols), drop = FALSE]

  xglobal <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(xglobal) <- "double"

  if (any(!is.finite(xglobal))) {
    stop("model_input.csv contains NA/NaN/Inf.")
  }

  quarters <- as.character(dat$Quarter)

  wdat <- read.csv(
    weight_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Country" %in% names(wdat)) {
    stop("trade_weights.csv must contain Country.")
  }

  rownames(wdat) <- as.character(wdat$Country)

  missing_rows <- setdiff(countries, rownames(wdat))
  missing_wcols <- setdiff(countries, names(wdat))

  if (length(missing_rows) || length(missing_wcols)) {
    stop("Incomplete trade matrix in trade_weights.csv.")
  }

  Wtrade <- as.matrix(
    wdat[countries, countries, drop = FALSE]
  )
  storage.mode(Wtrade) <- "double"

  if (any(!is.finite(Wtrade)) || any(Wtrade < -1e-12)) {
    stop("Invalid trade weights.")
  }

  diag(Wtrade) <- 0
  rs <- rowSums(Wtrade)

  if (any(rs <= 0)) {
    stop("Non-positive trade row sum.")
  }

  Wtrade <- Wtrade / rs

  # gW[[i]] maps the 72-variable global vector into 12 local rows:
  #
  # US:
  #   [GPR, VIX, y, dp, r, de, deq |
  #    foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq]
  #
  # non-US:
  #   [y, dp, r, de, deq |
  #    foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq,
  #    global_gpr, global_vix]
  gW <- vector("list", length(countries))
  names(gW) <- countries

  K <- ncol(xglobal)

  for (i in seq_along(countries)) {

    cc <- countries[[i]]

    own_names <- paste0(cc, "_", macro_vars)

    if (cc == "US") {
      own_names <- c(
        "US_gpr",
        "US_vix",
        own_names
      )
    }

    own_idx <- match(
      own_names,
      colnames(xglobal)
    )

    if (anyNA(own_idx)) {
      stop(
        "Could not locate all own variables for ",
        cc
      )
    }

    own_select <- matrix(
      0,
      nrow = length(own_idx),
      ncol = K
    )

    own_select[
      cbind(seq_along(own_idx), own_idx)
    ] <- 1

    rownames(own_select) <- own_names

    foreign <- matrix(
      0,
      nrow = length(macro_vars),
      ncol = K
    )

    rownames(foreign) <- paste0(
      "foreign_",
      macro_vars
    )

    for (v in seq_along(macro_vars)) {

      for (partner in countries[countries != cc]) {

        target <- paste0(
          partner,
          "_",
          macro_vars[[v]]
        )

        j <- match(
          target,
          colnames(xglobal)
        )

        if (is.na(j)) {
          stop(
            "Foreign variable not found: ",
            target
          )
        }

        foreign[v, j] <- Wtrade[cc, partner]
      }

      s <- sum(foreign[v, ])

      if (s <= 0) {
        stop(
          "Foreign weights sum to zero for ",
          cc,
          " / ",
          macro_vars[[v]]
        )
      }

      foreign[v, ] <- foreign[v, ] / s
    }

    if (cc != "US") {

      gpr_idx <- match(
        "US_gpr",
        colnames(xglobal)
      )

      vix_idx <- match(
        "US_vix",
        colnames(xglobal)
      )

      if (is.na(gpr_idx) || is.na(vix_idx)) {
        stop("US_gpr or US_vix not found in global vector.")
      }

      gpr_row <- matrix(
        0,
        nrow = 1,
        ncol = K
      )

      gpr_row[1, gpr_idx] <- 1
      rownames(gpr_row) <- "global_gpr"

      vix_row <- matrix(
        0,
        nrow = 1,
        ncol = K
      )

      vix_row[1, vix_idx] <- 1
      rownames(vix_row) <- "global_vix"

      foreign <- rbind(
        foreign,
        gpr_row,
        vix_row
      )
    }

    gW[[i]] <- rbind(
      own_select,
      foreign
    )
  }

  endogenous_total <- sum(
    vapply(
      countries,
      function(cc) {
        if (cc == "US") 7L else 5L
      },
      integer(1)
    )
  )

  if (endogenous_total != ncol(xglobal)) {
    stop(
      "Global dimension mismatch: expected ",
      endogenous_total,
      ", got ",
      ncol(xglobal)
    )
  }

  if (endogenous_total != 72L) {
    stop(
      "GPR+VIX experiment should contain exactly 72 global variables."
    )
  }

  local_rows <- vapply(
    gW,
    nrow,
    integer(1)
  )

  if (any(local_rows != 12L)) {
    stop(
      "Every GPR+VIX country mapping must have 12 local rows. Found: ",
      paste(
        names(local_rows),
        local_rows,
        sep = "=",
        collapse = ", "
      )
    )
  }

  list(
    bigx = xglobal,
    gW = gW,
    countries = countries,
    quarters = quarters,
    new.data = xglobal,
    Wtrade = Wtrade,
    macro_vars = macro_vars,
    global_risk_vars = c(
      "US_gpr",
      "US_vix"
    )
  )
}
