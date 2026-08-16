#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

options(stringsAsFactors = FALSE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
# Recursive macro ordering for structural identification.
macro_vars <- c("y", "dp", "r", "de", "deq")

resolve_file <- function(env_name, candidates) {
  x <- Sys.getenv(env_name, "")
  if (nzchar(x)) {
    if (!file.exists(x)) stop(env_name, " points to a missing file: ", x)
    return(x)
  }
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Could not locate ", env_name, ". Tried: ", paste(candidates, collapse = ", "))
  hit[[1]]
}

macro_file <- resolve_file(
  "TVPGVAR_MACRO_XLSX",
  c(
    "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx",
    "TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx",
    "TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成(2).xlsx"
  )
)
weight_file <- resolve_file(
  "TVPGVAR_WEIGHT_XLSX",
  c(
    "8.12/Trade_Weights_14_Economies_2000_2012(2).xlsx",
    "Trade_Weights_14_Economies_2000_2012(2).xlsx",
    "Trade_Weights_14_Economies_2000_2012(2)(2).xlsx"
  )
)
gpr_file <- resolve_file(
  "TVPGVAR_GPR_CSV",
  c("data/gpr_quarterly_processed.csv", "GPR_quarterly_processed_2000Q1_2026Q2.csv")
)

to_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
quarter_index <- function(q) {
  q <- as.character(q)
  if (any(!grepl("^[12][0-9]{3}Q[1-4]$", q))) stop("Invalid quarter label detected.")
  as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6))
}

# -----------------------------------------------------------------------------
# 1. Macro data: use the already validated balanced 14-country model sheet.
# -----------------------------------------------------------------------------
macro <- suppressWarnings(read_excel(macro_file, sheet = "MODEL_COMPLETE_14C", .name_repair = "minimal"))
macro <- as.data.frame(macro, check.names = FALSE, stringsAsFactors = FALSE)
names(macro) <- trimws(names(macro))

source_suffix <- c(
  y   = "GDP_DLOG",
  dp  = "CPI_DLOG",
  r   = "RATE_LEVEL",
  de  = "REER_DLOG",
  deq = "EQ_RETURN"
)

required_source <- c(
  "Quarter",
  unlist(lapply(countries, function(cc) paste0(cc, "_", unname(source_suffix))), use.names = FALSE)
)
missing_macro <- setdiff(required_source, names(macro))
if (length(missing_macro)) stop("MODEL_COMPLETE_14C is missing: ", paste(missing_macro, collapse = ", "))

mout <- data.frame(Quarter = as.character(macro[["Quarter"]]), stringsAsFactors = FALSE, check.names = FALSE)
for (cc in countries) {
  mout[[paste0(cc, "_y")]]   <- to_numeric(macro[[paste0(cc, "_GDP_DLOG")]])
  mout[[paste0(cc, "_dp")]]  <- to_numeric(macro[[paste0(cc, "_CPI_DLOG")]])
  mout[[paste0(cc, "_r")]]   <- to_numeric(macro[[paste0(cc, "_RATE_LEVEL")]])
  mout[[paste0(cc, "_de")]]  <- to_numeric(macro[[paste0(cc, "_REER_DLOG")]])
  mout[[paste0(cc, "_deq")]] <- to_numeric(macro[[paste0(cc, "_EQ_RETURN")]])
}

if (anyDuplicated(mout$Quarter)) stop("Duplicate quarters in MODEL_COMPLETE_14C.")
mout <- mout[order(quarter_index(mout$Quarter)), , drop = FALSE]

# -----------------------------------------------------------------------------
# 2. Global GPR: baseline is explicitly LN_GPR_QMEAN.
# Internal name US_gpr is retained only because the existing BVAR code selects
# the dominant US block using the two-character column prefix "US".
# Economically this is GLOBAL GPR, not country-specific US GPR.
# -----------------------------------------------------------------------------
gpr <- read.csv(gpr_file, check.names = FALSE, stringsAsFactors = FALSE)
gpr_column <- Sys.getenv("TVPGVAR_GPR_COLUMN", "LN_GPR_QMEAN")
allowed_gpr_columns <- c("LN_GPR_QMEAN", "LN_GPR_QMAX", "LN_GPRT_QMEAN", "LN_GPRA_QMEAN")
if (!gpr_column %in% allowed_gpr_columns) {
  stop("TVPGVAR_GPR_COLUMN must be one of: ", paste(allowed_gpr_columns, collapse = ", "))
}
if (!all(c("Quarter", gpr_column) %in% names(gpr))) {
  stop("Processed GPR CSV is missing required column: ", gpr_column)
}
gpr <- gpr[, c("Quarter", gpr_column)]
names(gpr)[2] <- "US_gpr"
gpr$Quarter <- as.character(gpr$Quarter)
gpr$US_gpr <- to_numeric(gpr$US_gpr)
if (anyDuplicated(gpr$Quarter)) stop("Duplicate quarters in processed GPR file.")

model <- merge(mout, gpr, by = "Quarter", all.x = TRUE, sort = FALSE)
model <- model[match(mout$Quarter, model$Quarter), , drop = FALSE]
if (any(!is.finite(model$US_gpr))) {
  stop("Global GPR does not fully cover the balanced macro sample.")
}

# -----------------------------------------------------------------------------
# 3. Deterministic global ordering.
# US block: GPR -> y -> dp -> r -> de -> deq
# Other blocks: y -> dp -> r -> de -> deq
# -----------------------------------------------------------------------------
ordered_cols <- unlist(lapply(countries, function(cc) {
  own <- paste0(cc, "_", macro_vars)
  if (cc == "US") c("US_gpr", own) else own
}), use.names = FALSE)
model <- model[, c("Quarter", ordered_cols), drop = FALSE]

if (ncol(model) != 72L) stop("Expected Quarter + 71 global variables, found ", ncol(model), " columns.")
if (any(!complete.cases(model))) stop("NA values remain in structural model input.")
idx <- quarter_index(model$Quarter)
if (any(diff(idx) != 1L)) stop("Structural model input is not a continuous quarterly sample.")

# -----------------------------------------------------------------------------
# 4. Trade weights: exact Trade_Weights sheet, diagonal=0, row-normalized.
# -----------------------------------------------------------------------------
wtab <- suppressWarnings(read_excel(weight_file, sheet = "Trade_Weights", .name_repair = "minimal"))
wtab <- as.data.frame(wtab, check.names = FALSE, stringsAsFactors = FALSE)
names(wtab) <- trimws(names(wtab))
if (ncol(wtab) < 15L) stop("Trade_Weights sheet is unexpectedly narrow.")
row_ids <- trimws(as.character(wtab[[1]]))
missing_rows <- setdiff(countries, row_ids)
missing_cols <- setdiff(countries, names(wtab))
if (length(missing_rows) || length(missing_cols)) {
  stop("Incomplete 14-country trade matrix. Missing rows: ", paste(missing_rows, collapse = ", "),
       "; missing cols: ", paste(missing_cols, collapse = ", "))
}
W <- matrix(NA_real_, length(countries), length(countries), dimnames = list(countries, countries))
for (r in countries) {
  rr <- match(r, row_ids)
  for (c in countries) W[r, c] <- to_numeric(wtab[[c]][rr])
}
if (any(!is.finite(W)) || any(W < -1e-12)) stop("Invalid trade weights.")
diag(W) <- 0
rs <- rowSums(W)
if (any(rs <= 0)) stop("A trade-weight row has non-positive sum.")
W <- W / rs

write.csv(model, "data/model_input.csv", row.names = FALSE, na = "")
write.csv(model, "data/model_input_structural.csv", row.names = FALSE, na = "")
write.csv(data.frame(Country = countries, W, check.names = FALSE),
          "data/trade_weights.csv", row.names = FALSE, na = "")

summary_lines <- c(
  "14-economy structural TVP-GVAR input summary",
  paste0("Macro file: ", macro_file),
  paste0("Trade weights: ", weight_file),
  paste0("Processed GPR: ", gpr_file),
  paste0("Selected GPR column: ", gpr_column),
  "Economic GPR concept: Global GPR (internal compatibility name: US_gpr)",
  "US recursive order: GPR -> y -> dp -> r -> de -> deq",
  "Non-US recursive order: y -> dp -> r -> de -> deq",
  paste0("Sample: ", model$Quarter[1], " - ", tail(model$Quarter, 1)),
  paste0("Observations: ", nrow(model)),
  paste0("Global variables: ", ncol(model) - 1L),
  paste0("Trade row-sum range: ", sprintf("%.12f", min(rowSums(W))), " - ", sprintf("%.12f", max(rowSums(W)))),
  paste0("Trade diagonal max abs: ", sprintf("%.12f", max(abs(diag(W)))))
)
writeLines(summary_lines, "results/build_structural_input_summary.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")
