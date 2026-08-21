#!/usr/bin/env Rscript

# =============================================================================
# Add official VIX data to the existing structural TVP-GVAR input
#
# Source:
#   FRED series VIXCLS
#   Original source: Chicago Board Options Exchange (CBOE)
#   Frequency: daily close
#
# Transformation used here:
#   Daily VIX -> quarterly arithmetic mean -> log(quarterly mean)
#   Internal model name: US_vix
#
# IMPORTANT:
#   Run 8.12/02_build_14e_input_structural.R FIRST.
#   This script then adds VIX to data/model_input.csv without changing
#   any macro data, GPR data, or trade weights.
# =============================================================================

options(stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

fred_url <- Sys.getenv(
  "TVPGVAR_VIX_URL",
  "https://fred.stlouisfed.org/graph/fredgraph.csv?id=VIXCLS"
)

raw_file <- Sys.getenv(
  "TVPGVAR_VIX_RAW_CSV",
  "data/VIXCLS_daily_fred.csv"
)

processed_file <- Sys.getenv(
  "TVPGVAR_VIX_PROCESSED_CSV",
  "data/vix_quarterly_processed.csv"
)

model_file <- Sys.getenv(
  "TVPGVAR_MODEL_INPUT",
  "data/model_input.csv"
)

if (!file.exists(model_file)) {
  stop(
    "Missing ", model_file,
    ". Run 8.12/02_build_14e_input_structural.R first."
  )
}

# Allow a user-supplied local VIX file. Otherwise download the official FRED CSV.
if (!file.exists(raw_file) || identical(Sys.getenv("TVPGVAR_REFRESH_VIX", "1"), "1")) {
  cat("Downloading official VIXCLS from FRED:\n", fred_url, "\n", sep = "")
  ok <- tryCatch({
    utils::download.file(
      fred_url,
      destfile = raw_file,
      mode = "wb",
      quiet = FALSE,
      method = "libcurl"
    )
    TRUE
  }, error = function(e) {
    message("libcurl download failed: ", conditionMessage(e))
    FALSE
  })

  if (!ok || !file.exists(raw_file)) {
    stop(
      "Could not download VIXCLS from FRED. ",
      "You may provide a local official FRED CSV via TVPGVAR_VIX_RAW_CSV."
    )
  }
}

vix_raw <- read.csv(
  raw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", ".", "NA")
)

date_col <- intersect(c("DATE", "observation_date", "Date", "date"), names(vix_raw))
value_col <- intersect(c("VIXCLS", "vixcls"), names(vix_raw))

if (!length(date_col)) stop("Could not find the date column in the VIX CSV.")
if (!length(value_col)) stop("Could not find VIXCLS in the VIX CSV.")

date_col <- date_col[[1]]
value_col <- value_col[[1]]

vix_raw$DATE_PARSED <- as.Date(vix_raw[[date_col]])
vix_raw$VIXCLS_NUM <- suppressWarnings(as.numeric(vix_raw[[value_col]]))

vix_raw <- vix_raw[
  !is.na(vix_raw$DATE_PARSED) &
    is.finite(vix_raw$VIXCLS_NUM) &
    vix_raw$VIXCLS_NUM > 0,
  ,
  drop = FALSE
]

if (!nrow(vix_raw)) stop("No usable positive VIX observations were found.")

lt <- as.POSIXlt(vix_raw$DATE_PARSED)
yy <- lt$year + 1900L
qq <- (lt$mon %/% 3L) + 1L
vix_raw$Quarter <- paste0(yy, "Q", qq)

qsplit <- split(vix_raw$VIXCLS_NUM, vix_raw$Quarter)

vix_q <- data.frame(
  Quarter = names(qsplit),
  VIX_QMEAN = vapply(qsplit, mean, numeric(1), na.rm = TRUE),
  VIX_N_DAILY = vapply(qsplit, function(z) sum(is.finite(z)), integer(1)),
  stringsAsFactors = FALSE
)

quarter_index <- function(q) {
  q <- as.character(q)
  if (any(!grepl("^[12][0-9]{3}Q[1-4]$", q))) {
    stop("Invalid quarter label detected.")
  }
  as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6))
}

vix_q <- vix_q[order(quarter_index(vix_q$Quarter)), , drop = FALSE]
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

model <- read.csv(
  model_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Quarter" %in% names(model)) stop("model_input.csv must contain Quarter.")
if (!"US_gpr" %in% names(model)) stop("model_input.csv must contain US_gpr.")
if (anyDuplicated(model$Quarter)) stop("Duplicate Quarter values in model_input.csv.")

# Ensure repeat runs cannot accidentally duplicate US_vix.
if ("US_vix" %in% names(model)) {
  model$US_vix <- NULL
}

ix <- match(model$Quarter, vix_q$Quarter)
if (anyNA(ix)) {
  stop(
    "VIX does not fully cover the existing structural sample. Missing quarters: ",
    paste(model$Quarter[is.na(ix)], collapse = ", ")
  )
}

model$US_vix <- vix_q$US_vix[ix]

if (any(!is.finite(model$US_vix))) {
  stop("Non-finite US_vix remains after sample alignment.")
}

# Put VIX directly after GPR so the US recursive order is:
#   GPR -> VIX -> y -> dp -> r -> de -> deq
old_names <- names(model)
gpr_pos <- match("US_gpr", old_names)
new_order <- c(
  old_names[seq_len(gpr_pos)],
  "US_vix",
  old_names[(gpr_pos + 1L):length(old_names)]
)
new_order <- new_order[!duplicated(new_order)]
model <- model[, new_order, drop = FALSE]

us_cols <- names(model)[substr(names(model), 1, 2) == "US"]
expected_us <- c(
  "US_gpr", "US_vix",
  "US_y", "US_dp", "US_r", "US_de", "US_deq"
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

write.csv(model, "data/model_input.csv", row.names = FALSE, na = "")
write.csv(model, "data/model_input_structural.csv", row.names = FALSE, na = "")

meta <- c(
  "GPR + VIX structural input augmentation",
  paste0("VIX FRED series: VIXCLS"),
  paste0("VIX source URL: ", fred_url),
  "Original VIX source: Chicago Board Options Exchange (CBOE)",
  "VIX source frequency: daily close",
  "VIX quarterly aggregation: arithmetic mean of available daily closes",
  "VIX model transformation: log(quarterly mean)",
  "Internal model name: US_vix",
  "US recursive order: GPR -> VIX -> y -> dp -> r -> de -> deq",
  paste0("Model sample: ", model$Quarter[[1]], " - ", tail(model$Quarter, 1)),
  paste0("Model observations: ", nrow(model)),
  paste0("Global variables after VIX addition: ", ncol(model) - 1L)
)

writeLines(meta, "results/vix_input_summary.txt")

# Append rather than overwrite the existing builder summary.
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

cat(paste(meta, collapse = "\n"), "\n")
