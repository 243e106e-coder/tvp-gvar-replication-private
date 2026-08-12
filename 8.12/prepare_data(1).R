prepare_gvar_data <- function(data_file = "data/model_input.csv",
                              weight_file = "data/trade_weights.csv") {
  countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  macro_vars <- c("y","dp","de","r","deq")

  dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE)

  if (!"Quarter" %in% names(dat)) {
    stop("model_input.csv must contain a Quarter column.")
  }
  if (anyDuplicated(dat$Quarter)) {
    stop("model_input.csv contains duplicate quarters.")
  }

  required_macro <- unlist(
    lapply(countries, function(cc) paste0(cc, "_", macro_vars)),
    use.names = FALSE
  )
  required <- c(required_macro, "US_gpr")
  missing_data <- setdiff(required, names(dat))

  if (length(missing_data)) {
    stop(
      "model_input.csv is missing required columns: ",
      paste(missing_data, collapse = ", ")
    )
  }

  # Deterministic global ordering. US_gpr belongs to the US endogenous block.
  ordered_cols <- unlist(lapply(countries, function(cc) {
    own <- paste0(cc, "_", macro_vars)
    if (cc == "US") c(own, "US_gpr") else own
  }), use.names = FALSE)

  dat <- dat[, c("Quarter", ordered_cols), drop = FALSE]

  xglobal <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(xglobal) <- "double"

  if (any(!is.finite(xglobal))) {
    stop("model_input.csv contains NA/NaN/Inf values in the estimation sample.")
  }

  quarters <- as.character(dat$Quarter)

  wdat <- read.csv(weight_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!ncol(wdat)) stop("trade_weights.csv is empty.")

  row_label <- names(wdat)[1]
  rownames(wdat) <- as.character(wdat[[row_label]])
  wdat[[row_label]] <- NULL

  missing_rows <- setdiff(countries, rownames(wdat))
  missing_cols <- setdiff(countries, names(wdat))

  if (length(missing_rows) || length(missing_cols)) {
    stop(
      "trade_weights.csv does not contain the complete 14-economy matrix. ",
      "Missing rows: ", paste(missing_rows, collapse = ", "),
      "; missing columns: ", paste(missing_cols, collapse = ", ")
    )
  }

  Wtrade <- as.matrix(wdat[countries, countries, drop = FALSE])
  storage.mode(Wtrade) <- "double"

  if (any(!is.finite(Wtrade))) {
    stop("trade_weights.csv contains NA/NaN/Inf values.")
  }
  if (any(Wtrade < -1e-12)) {
    stop("trade_weights.csv contains negative weights.")
  }

  # Home-country weights must be excluded.
  diag(Wtrade) <- 0
  rs <- rowSums(Wtrade)
  if (any(rs <= 0)) {
    stop("Trade-weight matrix contains a non-positive row sum.")
  }
  Wtrade <- Wtrade / rs

  # gW[[i]] maps the global vector into:
  # [own variables | trade-weighted foreign macro variables | global GPR if non-US].
  #
  # US:     6 own variables (5 macro + GPR) + 5 foreign macro variables
  # non-US: 5 own variables + 5 foreign macro variables + US_gpr
  #
  # This keeps US_gpr endogenous in the dominant US block while it is
  # weakly exogenous in each non-US country model.
  gW <- vector("list", length(countries))
  names(gW) <- countries
  K <- ncol(xglobal)

  for (i in seq_along(countries)) {
    cc <- countries[i]

    own_names <- paste0(cc, "_", macro_vars)
    if (cc == "US") own_names <- c(own_names, "US_gpr")

    own_idx <- match(own_names, colnames(xglobal))
    if (anyNA(own_idx)) {
      stop("Could not locate all own variables for ", cc)
    }

    own_select <- matrix(0, length(own_idx), K)
    own_select[cbind(seq_along(own_idx), own_idx)] <- 1

    foreign <- matrix(0, length(macro_vars), K)
    rownames(foreign) <- paste0("foreign_", macro_vars)

    for (v in seq_along(macro_vars)) {
      for (partner in countries[countries != cc]) {
        target <- paste0(partner, "_", macro_vars[v])
        j <- match(target, colnames(xglobal))
        if (is.na(j)) stop("Could not locate foreign variable: ", target)
        foreign[v, j] <- Wtrade[cc, partner]
      }

      s <- sum(foreign[v, ])
      if (s <= 0) {
        stop("Foreign weights sum to zero for ", cc, " / ", macro_vars[v])
      }
      foreign[v, ] <- foreign[v, ] / s
    }

    if (cc != "US") {
      gpr_idx <- match("US_gpr", colnames(xglobal))
      if (is.na(gpr_idx)) stop("US_gpr is missing from the global data vector.")

      gpr_row <- matrix(0, 1, K)
      gpr_row[1, gpr_idx] <- 1
      rownames(gpr_row) <- "global_gpr"
      foreign <- rbind(foreign, gpr_row)
    }

    gW[[i]] <- rbind(own_select, foreign)
  }

  # Stacking check: the sum of country endogenous dimensions must equal
  # the number of variables in the global vector.
  endogenous_total <- sum(vapply(countries, function(cc) {
    if (cc == "US") 6L else 5L
  }, integer(1)))

  if (endogenous_total != ncol(xglobal)) {
    stop(
      "Global dimension mismatch: expected ", endogenous_total,
      " endogenous variables but model_input.csv has ", ncol(xglobal), "."
    )
  }

  list(
    bigx = xglobal,
    gW = gW,
    countries = countries,
    quarters = quarters,
    new.data = xglobal,
    Wtrade = Wtrade
  )
}
