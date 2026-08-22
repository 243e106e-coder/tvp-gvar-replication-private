# =============================================================================
# Separate dominant-unit GPR + VIX TVP-GVAR data mapping
#
# Country units (14):
#   y_i,t = [y, dp, r, de, deq]'
#
# Country foreign/global block:
#   x_i,t* = [foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq,
#             global_gpr, global_vix]'
#
# Separate dominant unit:
#   z_t = [GL_gpr, GL_vix]'
#
# Dominant-unit weak-exogeneity architecture:
#   - every country may load on GL_gpr and GL_vix;
#   - the dominant unit has NO foreign-country regressors;
#   - hence country macro variables cannot feed contemporaneously into GPR/VIX.
#
# This mirrors the dominant-unit architecture of Blagov, Dirks & Funke (2025),
# adapted to the user's [GPR,VIX] global-risk block instead of [GPR,HWWI].
# =============================================================================

prepare_gvar_data <- function(
    data_file = "data/model_input_dominant.csv",
    weight_file = "data/trade_weights.csv") {

  country_units <- c(
    "AU","BR","CA","CH","CN","EA","UK",
    "JP","KR","NO","SG","TR","US","ZA"
  )
  dominant_unit <- "GL"
  units <- c(country_units, dominant_unit)
  macro_vars <- c("y", "dp", "r", "de", "deq")
  dominant_vars <- c("GL_gpr", "GL_vix")

  dat <- read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Quarter" %in% names(dat)) stop("Dominant input must contain Quarter.")
  if (anyDuplicated(dat$Quarter)) stop("Duplicate quarters in dominant input.")

  macro_cols <- unlist(
    lapply(country_units, function(cc) paste0(cc, "_", macro_vars)),
    use.names = FALSE
  )
  expected_cols <- c(macro_cols, dominant_vars)
  missing_cols <- setdiff(expected_cols, names(dat))
  if (length(missing_cols)) {
    stop("Dominant input is missing: ", paste(missing_cols, collapse = ", "))
  }

  dat <- dat[, c("Quarter", expected_cols), drop = FALSE]
  xglobal <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(xglobal) <- "double"
  if (any(!is.finite(xglobal))) stop("Dominant input contains NA/NaN/Inf.")
  if (ncol(xglobal) != 72L) stop("Expected 72 global variables, found ", ncol(xglobal))

  quarters <- as.character(dat$Quarter)

  wdat <- read.csv(weight_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Country" %in% names(wdat)) stop("Trade-weight file must contain Country.")
  rownames(wdat) <- as.character(wdat$Country)

  if (length(setdiff(country_units, rownames(wdat))) ||
      length(setdiff(country_units, names(wdat)))) {
    stop("Trade matrix does not contain all 14 country units.")
  }

  Wtrade <- as.matrix(wdat[country_units, country_units, drop = FALSE])
  storage.mode(Wtrade) <- "double"
  if (any(!is.finite(Wtrade)) || any(Wtrade < -1e-12)) stop("Invalid trade weights.")
  diag(Wtrade) <- 0
  rs <- rowSums(Wtrade)
  if (any(rs <= 0)) stop("Non-positive trade-weight row sum.")
  Wtrade <- Wtrade / rs

  K <- ncol(xglobal)
  gpr_idx <- match("GL_gpr", colnames(xglobal))
  vix_idx <- match("GL_vix", colnames(xglobal))
  if (is.na(gpr_idx) || is.na(vix_idx)) stop("GL_gpr/GL_vix not found.")

  gW <- vector("list", length(units))
  names(gW) <- units

  # ---------------------------------------------------------------------------
  # 14 country blocks: 5 own + 5 trade-weighted foreign + 2 global risk rows
  # ---------------------------------------------------------------------------
  for (cc in country_units) {
    own_names <- paste0(cc, "_", macro_vars)
    own_idx <- match(own_names, colnames(xglobal))
    if (anyNA(own_idx)) stop("Could not locate own variables for ", cc)

    own_select <- matrix(0, nrow = 5L, ncol = K)
    own_select[cbind(seq_len(5L), own_idx)] <- 1
    rownames(own_select) <- own_names

    foreign <- matrix(0, nrow = 5L, ncol = K)
    rownames(foreign) <- paste0("foreign_", macro_vars)

    for (v in seq_along(macro_vars)) {
      for (partner in country_units[country_units != cc]) {
        j <- match(paste0(partner, "_", macro_vars[[v]]), colnames(xglobal))
        if (is.na(j)) stop("Foreign variable missing for ", cc, ": ", partner)
        foreign[v, j] <- Wtrade[cc, partner]
      }
      s <- sum(foreign[v, ])
      if (s <= 0) stop("Foreign weights sum to zero for ", cc, "/", macro_vars[[v]])
      foreign[v, ] <- foreign[v, ] / s
    }

    global_gpr <- matrix(0, nrow = 1L, ncol = K)
    global_vix <- matrix(0, nrow = 1L, ncol = K)
    global_gpr[1, gpr_idx] <- 1
    global_vix[1, vix_idx] <- 1
    rownames(global_gpr) <- "global_gpr"
    rownames(global_vix) <- "global_vix"

    gW[[cc]] <- rbind(own_select, foreign, global_gpr, global_vix)
  }

  # ---------------------------------------------------------------------------
  # Dominant block: own [GPR,VIX] only. No country macro variables enter here.
  # ---------------------------------------------------------------------------
  dominant_select <- matrix(0, nrow = 2L, ncol = K)
  dominant_select[1, gpr_idx] <- 1
  dominant_select[2, vix_idx] <- 1
  rownames(dominant_select) <- dominant_vars
  gW[[dominant_unit]] <- dominant_select

  # ---------------------------------------------------------------------------
  # Hard architecture checks
  # ---------------------------------------------------------------------------
  for (cc in country_units) {
    Wi <- gW[[cc]]
    if (nrow(Wi) != 12L) stop("Country mapping must have 12 rows for ", cc)
    if (!identical(rownames(Wi)[6:12], c(
      "foreign_y","foreign_dp","foreign_r","foreign_de","foreign_deq",
      "global_gpr","global_vix"
    ))) stop("Country Wex row order mismatch for ", cc)
  }

  Wgl <- gW[[dominant_unit]]
  if (nrow(Wgl) != 2L) stop("Dominant mapping must have exactly two own rows.")
  if (any(abs(Wgl[, -c(gpr_idx, vix_idx), drop = FALSE]) > 1e-14)) {
    stop("Dominant unit incorrectly loads on country variables.")
  }

  list(
    bigx = xglobal,
    gW = gW,
    countries = units,
    country_units = country_units,
    dominant_unit = dominant_unit,
    quarters = quarters,
    new.data = xglobal,
    Wtrade = Wtrade,
    macro_vars = macro_vars,
    dominant_vars = dominant_vars,
    specification = paste(
      "14 country blocks + separate dominant unit [GPR,VIX];",
      "countries receive current/lagged GPR,VIX; dominant unit receives no country variables"
    )
  )
}
