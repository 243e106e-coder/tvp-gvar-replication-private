# =============================================================================
# GPR + VIX US-only transmission overrides
#
# Source AFTER:
#   R/gpr_structural_irf.R
#   R/gpr_vix_overrides.R
#
# Purpose:
#   Keep the existing 7x7 US recursive structural GPR shock implementation,
#   but remove direct non-US global-GPR regressors from BOTH current and lagged
#   foreign/global blocks.
#
# non-US Wex order:
#   1 foreign_y
#   2 foreign_dp
#   3 foreign_r
#   4 foreign_de
#   5 foreign_deq
#   6 global_vix
#
# US Wex order:
#   5 foreign macro variables
# =============================================================================

tvpgvar_wex_lag <- function(
    Wex,
    lag_order = 1L,
    country_name,
    mode = c("current_only", "current_and_lag")) {

  # mode is accepted only for API compatibility with the existing estimator.
  # In this experiment GPR is absent from every non-US Wex block by design.
  mode <- match.arg(mode)
  lag_order <- as.integer(lag_order)

  if (
    length(lag_order) != 1L ||
      is.na(lag_order) ||
      lag_order != 1L
  ) {
    stop(
      "The GPR+VIX US-only experiment requires foreign lag q=1."
    )
  }

  L <- mlag(
    Wex,
    lag_order
  )

  cc <- as.character(country_name)

  expected <- if (identical(cc, "US")) 5L else 6L

  if (ncol(L) != expected) {
    stop(
      "Unexpected Wex lag width for ",
      cc,
      ": found ",
      ncol(L),
      ", expected ",
      expected,
      "."
    )
  }

  L
}


align_lag_block <- function(
    lag_block,
    W_i,
    k_i,
    country_name) {

  lag_block <- as.matrix(lag_block)
  expected_foreign <- nrow(W_i) - k_i

  if (ncol(lag_block) != expected_foreign) {
    stop(
      "Lag-block/W mismatch for ",
      country_name,
      ": lag columns=",
      ncol(lag_block),
      ", expected foreign columns=",
      expected_foreign,
      ", W rows=",
      nrow(W_i),
      ", own variables=",
      k_i
    )
  }

  lag_block
}
