#!/usr/bin/env Rscript

# =============================================================================
# Add local official FRED VIXCLS data to the existing structural TVP-GVAR input
#
# Data file expected:
#   8.12/VIXCLS_daily_fred.csv
#
# Source:
#   FRED series VIXCLS
#   Original source: Chicago Board Options Exchange (CBOE)
#   Frequency: daily close
#
# Transformation:
#   Daily VIXCLS
#     -> remove missing / non-positive observations
#     -> quarterly arithmetic mean of available daily closes
#     -> log(quarterly mean)
#   Internal model name: US_vix
#
# IMPORTANT:
#   1) Run 8.12/02_build_14e_input_structural.R first.
#   2) This script NEVER downloads VIX from the internet.
#   3) The local raw CSV must cover every quarter in model_input.csv.
# =============================================================================

options(stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

raw_file <- Sys.getenv(
  "TVPGVAR_VIX_RAW_CSV",
  "8.12/VIXCLS_daily_fred.csv"
)

processed_file <- Sys.getenv(
  "TVPGVAR_VIX_PROCESSED_CSV",
  "data/vix_quarterly_processed.csv"
)

model_file <- Sys.getenv(
  "TVPGVAR_MODEL_INPUT",
  "data/model_input.csv"
)

fred_series <- "VIXCLS"
fred_url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=VIXCLS"

if (!file.exists(model_file)) {
  stop(
    "Missing ", model_file,
    ". Run 8.12/02_build_14e_input_structural.R first."
  )
}

if (!file.exists(raw_file)) {
  stop(
    "Missing local official VIX file: ", raw_file, "\n",
    "Expected file: 8.12/VIXCLS_daily_fred.csv\n",
    "This script intentionally does not download VIX from the internet."
  )
}

# -----------------------------------------------------------------------------
# 1. Read and validate local official FRED VIXCLS CSV
# -----------------------------------------------------------------------------

vix_raw <- read.csv(
  raw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", ".", "NA")
)

date_col <- intersect(
  c("DATE", "observation_date", "Date", "date"),
  names(vix_raw)
)

value_col <- intersect(
  c("VIXCLS", "vixcls"),
  names(vix_raw)
)

if (!length(date_col)) {
  stop(
    "Could not find a date column in ", raw_file,
    ". Expected one of: DATE, observation_date, Date, date."
  )
}

if (!length(value_col)) {
  stop(
    "Could not find VIXCLS in ", raw_file, "."
  )
}

date_col <- date_col[[1]]
value_col <- value_col[[1]]

vix_raw$DATE_PARSED <- as.Date(vix_raw[[date_col]])
vix_raw$VIXCLS_NUM <- suppressWarnings(
  as.numeric(vix_raw[[value_col]])
)

raw_date_min <- suppressWarnings(
  min(vix_raw$DATE_PARSED, na.rm = TRUE)
)
raw_date_max <- suppressWarnings(
  max(vix_raw$DATE_PARSED, na.rm = TRUE)
)

n_raw <- nrow(vix_raw)
n_missing_value <- sum(is.na(vix_raw$VIXCLS_NUM))
n_invalid_date <- sum(is.na(vix_raw$DATE_PARSED))

vix_clean <- vix_raw[
  !is.na(vix_raw$DATE_PARSED) &
    is.finite(vix_raw$VIXCLS_NUM) &
    vix_raw$VIXCLS_NUM > 0,
  ,
  drop = FALSE
]

if (!nrow(vix_clean)) {
  stop("No usable positive VIXCLS observations were found.")
}

# -----------------------------------------------------------------------------
# 2. Convert daily VIX to quarterly mean and log-transform
# -----------------------------------------------------------------------------

lt <- as.POSIXlt(vix_clean$DATE_PARSED)
yy <- lt$year + 1900L
qq <- (lt$mon %/% 3L) + 1L
vix_clean$Quarter <- paste0(yy, "Q", qq)

qsplit <- split(vix_clean$VIXCLS_NUM, vix_clean$Quarter)

vix_q <- data.frame(
  Quarter = names(qsplit),
  VIX_QMEAN = vapply(
    qsplit,
    mean,
    numeric(1),
    na.rm = TRUE
  ),
  VIX_N_DAILY = vapply(
    qsplit,
    function(z) sum(is.finite(z)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

quarter_index <- function(q) {
  q <- as.character(q)

  if (any(!grepl("^[12][0-9]{3}Q[1-4]$", q))) {
    stop("Invalid quarter label detected.")
  }

  as.integer(substr(q, 1, 4)) * 4L +
    as.integer(substr(q, 6, 6))
}

vix_q <- vix_q[
  order(quarter_index(vix_q$Quarter)),
  ,
  drop = FALSE
]

vix_q$US_vix <- log(vix_q$VIX_QMEAN)

if (any(!is.finite(vix_q$US_vix))) {
  stop("Non-finite log VIX values were produced.")
}

write.csv(
  vix_q,
  processed_file,
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------
# 3. Read structural model input and align VIX exactly to its sample
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

model_q_index <- quarter_index(model$Quarter)

if (is.unsorted(model_q_index)) {
  stop("model_input.csv quarters are not sorted chronologically.")
}

model_start <- model$Quarter[[1]]
model_end <- tail(model$Quarter, 1)

# Ensure repeat runs cannot duplicate US_vix.
if ("US_vix" %in% names(model)) {
  model$US_vix <- NULL
}

ix <- match(model$Quarter, vix_q$Quarter)

if (anyNA(ix)) {
  missing_q <- model$Quarter[is.na(ix)]

  stop(
    "Local VIX data does not fully cover the structural model sample.\n",
    "Model sample: ", model_start, " - ", model_end, "\n",
    "Missing VIX quarters: ",
    paste(missing_q, collapse = ", ")
  )
}

model$US_vix <- vix_q$US_vix[ix]

if (any(!is.finite(model$US_vix))) {
  stop("Non-finite US_vix remains after sample alignment.")
}

# -----------------------------------------------------------------------------
# 4. Put VIX directly after GPR
#
# US recursive order:
#   GPR -> VIX -> y -> dp -> r -> de -> deq
# -----------------------------------------------------------------------------

old_names <- names(model)
gpr_pos <- match("US_gpr", old_names)

if (is.na(gpr_pos)) {
  stop("US_gpr position could not be identified.")
}

new_order <- c(
  old_names[seq_len(gpr_pos)],
  "US_vix",
  old_names[(gpr_pos + 1L):length(old_names)]
)

new_order <- new_order[!duplicated(new_order)]
model <- model[, new_order, drop = FALSE]

us_cols <- names(model)[substr(names(model), 1, 2) == "US"]

expected_us <- c(
  "US_gpr",
  "US_vix",
  "US_y",
  "US_dp",
  "US_r",
  "US_de",
  "US_deq"
)

if (!identical(us_cols, expected_us)) {
  stop(
    "US ordering after VIX insertion is wrong. Found: ",
    paste(us_cols, collapse = ", ")
  )
}

if (any(!complete.cases(model))) {
  stop("NA values remain after adding VIX.")
}

# -----------------------------------------------------------------------------
# 5. Save augmented structural input
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
# 6. Reproducibility metadata
# -----------------------------------------------------------------------------

used_vix <- vix_q[ix, , drop = FALSE]

meta <- c(
  "GPR + VIX structural input augmentation",
  paste0("VIX raw file: ", raw_file),
  paste0("VIX FRED series: ", fred_series),
  paste0("VIX source URL: ", fred_url),
  "Original VIX source: Chicago Board Options Exchange (CBOE)",
  "VIX source frequency: daily close",
  paste0(
    "Raw VIX file date coverage: ",
    format(raw_date_min),
    " - ",
    format(raw_date_max)
  ),
  paste0("Raw VIX rows: ", n_raw),
  paste0("Raw VIX missing-value rows: ", n_missing_value),
  paste0("Raw VIX invalid-date rows: ", n_invalid_date),
  "VIX missing daily observations: dropped before quarterly aggregation",
  "VIX quarterly aggregation: arithmetic mean of available positive daily closes",
  "VIX model transformation: log(quarterly mean)",
  "Internal model name: US_vix",
  paste0(
    "Quarterly VIX coverage available: ",
    vix_q$Quarter[[1]],
    " - ",
    tail(vix_q$Quarter, 1)
  ),
  paste0(
    "Quarterly VIX coverage used by model: ",
    used_vix$Quarter[[1]],
    " - ",
    tail(used_vix$Quarter, 1)
  ),
  paste0(
    "Minimum daily observations in a used quarter: ",
    min(used_vix$VIX_N_DAILY)
  ),
  paste0(
    "Maximum daily observations in a used quarter: ",
    max(used_vix$VIX_N_DAILY)
  ),
  paste0(
    "Median daily observations in a used quarter: ",
    median(used_vix$VIX_N_DAILY)
  ),
  "US recursive order: GPR -> VIX -> y -> dp -> r -> de -> deq",
  paste0("Model sample: ", model_start, " - ", model_end),
  paste0("Model observations: ", nrow(model)),
  paste0(
    "Global variables after VIX addition: ",
    ncol(model) - 1L
  )
)

writeLines(
  meta,
  "results/vix_input_summary.txt"
)

if (file.exists("results/build_structural_input_summary.txt")) {
  cat(
    "\n\n--- VIX AUGMENTATION ---\n",
    paste(meta, collapse = "\n"),
    "\n",
    file = "results/build_structural_input_summary.txt",
    append = TRUE,
    sep = ""
  )
}

cat(
  paste(meta, collapse = "\n"),
  "\n"
)
