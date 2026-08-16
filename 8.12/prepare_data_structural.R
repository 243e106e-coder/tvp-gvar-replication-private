prepare_gvar_data <- function(data_file = "data/model_input.csv",
                              weight_file = "data/trade_weights.csv") {
  countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
  macro_vars <- c("y", "dp", "r", "de", "deq")

  dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Quarter" %in% names(dat)) stop("model_input.csv must contain Quarter.")
  if (anyDuplicated(dat$Quarter)) stop("Duplicate quarters in model_input.csv.")

  expected_cols <- unlist(lapply(countries, function(cc) {
    own <- paste0(cc, "_", macro_vars)
    if (cc == "US") c("US_gpr", own) else own
  }), use.names = FALSE)
  missing_cols <- setdiff(expected_cols, names(dat))
  if (length(missing_cols)) stop("model_input.csv is missing: ", paste(missing_cols, collapse = ", "))

  # Enforce the exact global ordering used by the structural identification.
  dat <- dat[, c("Quarter", expected_cols), drop = FALSE]
  xglobal <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(xglobal) <- "double"
  if (any(!is.finite(xglobal))) stop("model_input.csv contains NA/NaN/Inf.")

  quarters <- as.character(dat$Quarter)

  wdat <- read.csv(weight_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Country" %in% names(wdat)) stop("trade_weights.csv must contain Country.")
  rownames(wdat) <- as.character(wdat$Country)
  missing_rows <- setdiff(countries, rownames(wdat))
  missing_wcols <- setdiff(countries, names(wdat))
  if (length(missing_rows) || length(missing_wcols)) stop("Incomplete trade matrix in trade_weights.csv.")

  Wtrade <- as.matrix(wdat[countries, countries, drop = FALSE])
  storage.mode(Wtrade) <- "double"
  if (any(!is.finite(Wtrade)) || any(Wtrade < -1e-12)) stop("Invalid trade weights.")
  diag(Wtrade) <- 0
  rs <- rowSums(Wtrade)
  if (any(rs <= 0)) stop("Non-positive trade row sum.")
  Wtrade <- Wtrade / rs

  # gW[[i]] maps the 71-variable global vector into:
  # US:     [GPR,y,dp,r,de,deq | foreign y,dp,r,de,deq]
  # non-US: [y,dp,r,de,deq | foreign y,dp,r,de,deq | global GPR]
  gW <- vector("list", length(countries))
  names(gW) <- countries
  K <- ncol(xglobal)

  for (i in seq_along(countries)) {
    cc <- countries[i]
    own_names <- paste0(cc, "_", macro_vars)
    if (cc == "US") own_names <- c("US_gpr", own_names)

    own_idx <- match(own_names, colnames(xglobal))
    if (anyNA(own_idx)) stop("Could not locate all own variables for ", cc)
    own_select <- matrix(0, length(own_idx), K)
    own_select[cbind(seq_along(own_idx), own_idx)] <- 1
    rownames(own_select) <- own_names

    foreign <- matrix(0, length(macro_vars), K)
    rownames(foreign) <- paste0("foreign_", macro_vars)
    for (v in seq_along(macro_vars)) {
      for (partner in countries[countries != cc]) {
        target <- paste0(partner, "_", macro_vars[v])
        j <- match(target, colnames(xglobal))
        if (is.na(j)) stop("Foreign variable not found: ", target)
        foreign[v, j] <- Wtrade[cc, partner]
      }
      s <- sum(foreign[v, ])
      if (s <= 0) stop("Foreign weights sum to zero for ", cc, " / ", macro_vars[v])
      foreign[v, ] <- foreign[v, ] / s
    }

    if (cc != "US") {
      gpr_idx <- match("US_gpr", colnames(xglobal))
      gpr_row <- matrix(0, 1, K)
      gpr_row[1, gpr_idx] <- 1
      rownames(gpr_row) <- "global_gpr"
      foreign <- rbind(foreign, gpr_row)
    }

    gW[[i]] <- rbind(own_select, foreign)
  }

  endogenous_total <- sum(vapply(countries, function(cc) if (cc == "US") 6L else 5L, integer(1)))
  if (endogenous_total != ncol(xglobal)) {
    stop("Global dimension mismatch: expected ", endogenous_total, ", got ", ncol(xglobal))
  }
  if (any(vapply(gW, nrow, integer(1)) != 11L)) stop("Every country mapping must have 11 local rows.")

  list(
    bigx = xglobal,
    gW = gW,
    countries = countries,
    quarters = quarters,
    new.data = xglobal,
    Wtrade = Wtrade,
    macro_vars = macro_vars
  )
}
