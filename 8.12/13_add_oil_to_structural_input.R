#!/usr/bin/env Rscript

# =============================================================================
# Add quarterly IMF Brent oil price to the structural TVP-GVAR input.
#
# Expected local file:
#   8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv
#
# Expected columns:
#   quarter, date, brent_usd_per_barrel, ln_brent
#
# Transformation used in the model:
#   US_oil = log(brent_usd_per_barrel)
#
# IMPORTANT:
#   - NO differencing
#   - NO growth rate
#   - NO standardization
#   - NO deflation
#   - NO quarterly aggregation (the source file is already quarterly)
#
# The precomputed ln_brent column is used only as a validation check.
# The model value is recomputed from the raw Brent price to make the
# transformation fully reproducible.
# =============================================================================

options(stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

raw_file <- Sys.getenv(
  "TVPGVAR_OIL_RAW_CSV",
  "8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"
)

processed_file <- Sys.getenv(
  "TVPGVAR_OIL_PROCESSED_CSV",
  "data/oil_quarterly_processed.csv"
)

model_file <- Sys.getenv(
  "TVPGVAR_MODEL_INPUT",
  "data/model_input.csv"
)

if (!file.exists(model_file)) {
  stop(
    "Missing ", model_file,
    ". Run 8.12/02b_build_14e_input_gdp_loglevel.R first."
  )
}

if (!file.exists(raw_file)) {
  stop(
    "Missing local Brent file: ", raw_file, "\n",
    "Expected: 8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"
  )
}

quarter_index <- function(q) {
  q <- as.character(q)
  if (any(!grepl("^[12][0-9]{3}Q[1-4]$", q))) {
    stop("Invalid quarter label detected: ",
         paste(unique(q[!grepl("^[12][0-9]{3}Q[1-4]$", q)]), collapse = ", "))
  }
  as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6))
}

# -----------------------------------------------------------------------------
# 1. Read and validate the local IMF Brent file
# -----------------------------------------------------------------------------

oil_raw <- read.csv(
  raw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", ".", "NA")
)

names_lower <- tolower(names(oil_raw))

find_col <- function(candidates, required = TRUE) {
  ii <- match(tolower(candidates), names_lower)
  ii <- ii[!is.na(ii)]
  if (!length(ii)) {
    if (required) {
      stop("Could not find required column. Expected one of: ",
           paste(candidates, collapse = ", "))
    }
    return(NA_integer_)
  }
  ii[[1]]
}

q_col <- find_col(c("quarter", "Quarter"))
price_col <- find_col(c(
  "brent_usd_per_barrel",
  "brent",
  "oil_price",
  "price"
))
log_col <- find_col(c("ln_brent", "log_brent", "oil"), required = FALSE)

oil_q <- data.frame(
  Quarter = as.character(oil_raw[[q_col]]),
  BRENT_USD_PER_BARREL = suppressWarnings(as.numeric(oil_raw[[price_col]])),
  stringsAsFactors = FALSE
)

if (anyDuplicated(oil_q$Quarter)) {
  dup <- unique(oil_q$Quarter[duplicated(oil_q$Quarter)])
  stop("Duplicate Brent quarters: ", paste(dup, collapse = ", "))
}

qi <- quarter_index(oil_q$Quarter)
oil_q <- oil_q[order(qi), , drop = FALSE]
rownames(oil_q) <- NULL

if (any(!is.finite(oil_q$BRENT_USD_PER_BARREL))) {
  stop("Brent price contains NA/NaN/Inf.")
}
if (any(oil_q$BRENT_USD_PER_BARREL <= 0)) {
  stop("Brent price must be strictly positive before log transformation.")
}

# Only requested transformation: natural log of the quarterly price LEVEL.
oil_q$US_oil <- log(oil_q$BRENT_USD_PER_BARREL)

if (any(!is.finite(oil_q$US_oil))) {
  stop("Non-finite log Brent values were produced.")
}

# Validate the precomputed ln_brent supplied in the uploaded file.
log_validation_max_abs_error <- NA_real_
if (!is.na(log_col)) {
  uploaded_log <- suppressWarnings(as.numeric(oil_raw[[log_col]]))
  # Reorder uploaded log in the same way as oil_q.
  uploaded_log <- uploaded_log[order(qi)]
  good <- is.finite(uploaded_log)
  if (!all(good)) {
    stop("Uploaded ln_brent contains missing/non-finite values.")
  }
  log_validation_max_abs_error <- max(abs(uploaded_log - oil_q$US_oil))
  # The supplied file is rounded to about 10 decimal places.
  if (log_validation_max_abs_error > 1e-8) {
    stop(
      "Uploaded ln_brent is inconsistent with log(brent_usd_per_barrel). ",
      "Max absolute error = ", signif(log_validation_max_abs_error, 8)
    )
  }
}

write.csv(
  oil_q,
  processed_file,
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------
# 2. Align Brent exactly to the existing model sample
# -----------------------------------------------------------------------------

model <- read.csv(
  model_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Quarter" %in% names(model)) {
  stop("model_input.csv must contain Quarter.")
}
if (!"US_gpr" %in% names(model)) {
  stop("model_input.csv must contain US_gpr.")
}
if (anyDuplicated(model$Quarter)) {
  stop("Duplicate Quarter values in model_input.csv.")
}
if (!nrow(model)) {
  stop("model_input.csv is empty.")
}

model_qi <- quarter_index(model$Quarter)
if (is.unsorted(model_qi)) {
  stop("model_input.csv quarters are not sorted chronologically.")
}

model_start <- model$Quarter[[1]]
model_end <- tail(model$Quarter, 1)

# Defensive cleanup so repeat runs cannot leave VIX or duplicate oil behind.
if ("US_vix" %in% names(model)) model$US_vix <- NULL
if ("US_oil" %in% names(model)) model$US_oil <- NULL

ix <- match(model$Quarter, oil_q$Quarter)

if (anyNA(ix)) {
  missing_q <- model$Quarter[is.na(ix)]
  stop(
    "Brent data does not fully cover the model sample.\n",
    "Model sample: ", model_start, " - ", model_end, "\n",
    "Missing oil quarters: ", paste(missing_q, collapse = ", ")
  )
}

model$US_oil <- oil_q$US_oil[ix]

if (any(!is.finite(model$US_oil))) {
  stop("Non-finite US_oil remains after sample alignment.")
}

# -----------------------------------------------------------------------------
# 3. Insert Oil directly after GPR
#
# Compatibility ordering before dominant-unit extraction:
#   GPR -> Oil -> y -> dp -> r -> de -> deq
#
# The next build step moves US_gpr/US_oil into a separate global block:
#   [GL_gpr, GL_oil]
# -----------------------------------------------------------------------------

old_names <- names(model)
gpr_pos <- match("US_gpr", old_names)

if (is.na(gpr_pos)) {
  stop("US_gpr position could not be identified.")
}

before <- old_names[seq_len(gpr_pos)]
after <- if (gpr_pos < length(old_names)) {
  old_names[(gpr_pos + 1L):length(old_names)]
} else {
  character(0)
}

new_order <- c(before, "US_oil", after)
new_order <- new_order[!duplicated(new_order)]
model <- model[, new_order, drop = FALSE]

us_cols <- names(model)[startsWith(names(model), "US_")]
expected_us <- c(
  "US_gpr",
  "US_oil",
  "US_y",
  "US_dp",
  "US_r",
  "US_de",
  "US_deq"
)

if (!identical(us_cols, expected_us)) {
  stop(
    "US ordering after Oil insertion is wrong. Found: ",
    paste(us_cols, collapse = ", ")
  )
}

numeric_block <- model[, setdiff(names(model), "Quarter"), drop = FALSE]
numeric_matrix <- as.matrix(numeric_block)
storage.mode(numeric_matrix) <- "double"

if (any(!is.finite(numeric_matrix))) {
  stop("NA/NaN/Inf values remain after adding Oil.")
}

# -----------------------------------------------------------------------------
# 4. Save augmented structural input
# -----------------------------------------------------------------------------

write.csv(
  model,
  "data/model_input.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  model,
  "data/model_input_structural.csv",
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------
# 5. Reproducibility metadata
# -----------------------------------------------------------------------------

used_oil <- oil_q[ix, , drop = FALSE]

meta <- c(
  "GPR + Brent Oil structural input augmentation",
  paste0("Oil raw file: ", raw_file),
  "Series: IMF Global price of Brent Crude",
  "Units: U.S. dollars per barrel",
  "Source frequency used: quarterly",
  "Oil quarterly aggregation: NONE (input is already quarterly)",
  "Oil model transformation: natural log of quarterly Brent PRICE LEVEL only",
  "Oil differencing: NONE",
  "Oil growth rate: NONE",
  "Oil standardization: NONE",
  "Oil deflation: NONE",
  "Internal compatibility name: US_oil",
  "Final dominant-unit name after next build step: GL_oil",
  paste0(
    "Available Brent coverage: ",
    oil_q$Quarter[[1]], " - ", tail(oil_q$Quarter, 1)
  ),
  paste0(
    "Brent coverage used by model: ",
    used_oil$Quarter[[1]], " - ", tail(used_oil$Quarter, 1)
  ),
  paste0(
    "Brent price used range (USD/barrel): ",
    signif(min(used_oil$BRENT_USD_PER_BARREL), 8),
    " - ",
    signif(max(used_oil$BRENT_USD_PER_BARREL), 8)
  ),
  paste0(
    "Log Brent used range: ",
    signif(min(used_oil$US_oil), 8),
    " - ",
    signif(max(used_oil$US_oil), 8)
  ),
  paste0(
    "Uploaded ln_brent validation max abs error: ",
    ifelse(is.na(log_validation_max_abs_error),
           "not available",
           format(log_validation_max_abs_error, scientific = TRUE))
  ),
  "Recursive dominant-block order for the new experiment: GPR -> Oil",
  paste0("Model sample: ", model_start, " - ", model_end),
  paste0("Model observations: ", nrow(model)),
  paste0("Model variables after Oil addition: ", ncol(model) - 1L)
)

writeLines(meta, "results/oil_input_summary.txt")

if (file.exists("results/build_structural_input_summary.txt")) {
  cat(
    "\n\n--- BRENT OIL AUGMENTATION ---\n",
    paste(meta, collapse = "\n"),
    "\n",
    file = "results/build_structural_input_summary.txt",
    append = TRUE,
    sep = ""
  )
}

cat(paste(meta, collapse = "\n"), "\n")
