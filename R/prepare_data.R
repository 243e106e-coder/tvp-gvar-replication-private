prepare_gvar_data <- function(data_file = "data/model_input.csv",
                              weight_file = "data/trade_weights.csv") {
  dat <- read.csv(data_file, check.names = FALSE)
  quarters <- dat$Quarter
  xglobal <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(xglobal) <- "double"

  countries <- c("CN","US","JP","KR","SG","UK","EA","CH","CA","AU")
  macro_vars <- c("y","dp","de","r","deq")
  Wtrade <- as.matrix(read.csv(weight_file, row.names = 1, check.names = FALSE))
  Wtrade <- Wtrade[countries, countries, drop = FALSE]

  # gW[[i]] maps global variables into country i's own variables followed by
  # foreign aggregates. GPR is endogenous in the US (dominant) block and is
  # weakly exogenous in every non-US country model.
  gW <- vector("list", length(countries))
  names(gW) <- countries
  K <- ncol(xglobal)
  for (i in seq_along(countries)) {
    cc <- countries[i]
    own_idx <- which(substr(colnames(xglobal), 1, 2) == cc)
    own_select <- matrix(0, length(own_idx), K)
    own_select[cbind(seq_along(own_idx), own_idx)] <- 1

    foreign <- matrix(0, length(macro_vars), K)
    for (v in seq_along(macro_vars)) {
      for (partner in countries[countries != cc]) {
        j <- match(paste0(partner, "_", macro_vars[v]), colnames(xglobal))
        foreign[v, j] <- Wtrade[cc, partner]
      }
      # Renormalize after excluding the home country (diagonal is already 0).
      s <- sum(foreign[v, ])
      if (s > 0) foreign[v, ] <- foreign[v, ] / s
    }
    if (cc != "US") {
      gpr_row <- matrix(0, 1, K)
      gpr_row[1, match("US_gpr", colnames(xglobal))] <- 1
      foreign <- rbind(foreign, gpr_row)
    }
    gW[[i]] <- rbind(own_select, foreign)
  }

  list(bigx = xglobal, gW = gW, countries = countries,
       quarters = quarters, new.data = xglobal)
}

