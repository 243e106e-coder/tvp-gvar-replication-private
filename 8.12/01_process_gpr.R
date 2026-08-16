#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

options(stringsAsFactors = FALSE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

resolve_file <- function(env_name, candidates) {
  x <- Sys.getenv(env_name, "")
  if (nzchar(x)) {
    if (!file.exists(x)) stop(env_name, " points to a missing file: ", x)
    return(x)
  }
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not locate input for ", env_name, ". Tried: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

gpr_file <- resolve_file(
  "TVPGVAR_GPR_RAW",
  c(
    "data_gpr_export.xls",
    "data_gpr_export (1).xls",
    "8.12/data_gpr_export.xls",
    "8.12/data_gpr_export (1).xls"
  )
)

raw <- suppressWarnings(read_excel(gpr_file, sheet = "Sheet1", .name_repair = "minimal"))
raw <- as.data.frame(raw, check.names = FALSE, stringsAsFactors = FALSE)
names(raw) <- trimws(names(raw))

required <- c("month", "GPR", "GPRT", "GPRA")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) {
  stop("Raw GPR file is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))

  z <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(z)))
  fmts <- c("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%Y-%m", "%Y/%m")
  for (fmt in fmts) {
    miss <- is.na(out) & nzchar(z)
    if (!any(miss)) break
    out[miss] <- suppressWarnings(as.Date(z[miss], format = fmt))
  }
  out
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

dat <- data.frame(
  Date = to_date(raw[["month"]]),
  GPR = num(raw[["GPR"]]),
  GPRT = num(raw[["GPRT"]]),
  GPRA = num(raw[["GPRA"]]),
  stringsAsFactors = FALSE
)

dat <- dat[!is.na(dat$Date) & is.finite(dat$GPR), , drop = FALSE]
if (!nrow(dat)) stop("No usable Recent GPR observations were found.")
if (any(dat$GPR <= 0, na.rm = TRUE)) stop("GPR contains non-positive values; log transform is not defined.")
if (any(dat$GPRT <= 0, na.rm = TRUE)) stop("GPRT contains non-positive values; log transform is not defined.")
if (any(dat$GPRA <= 0, na.rm = TRUE)) stop("GPRA contains non-positive values; log transform is not defined.")

month_num <- as.integer(format(dat$Date, "%m"))
year_num <- as.integer(format(dat$Date, "%Y"))
dat$Quarter <- paste0(year_num, "Q", ((month_num - 1L) %/% 3L) + 1L)

quarter_index <- function(q) {
  as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6))
}

parts <- split(dat, dat$Quarter)
qtab <- do.call(rbind, lapply(names(parts), function(q) {
  z <- parts[[q]]
  data.frame(
    Quarter = q,
    Months = nrow(z),
    GPR_QMEAN = mean(z$GPR, na.rm = TRUE),
    GPR_QMAX = max(z$GPR, na.rm = TRUE),
    GPRT_QMEAN = mean(z$GPRT, na.rm = TRUE),
    GPRA_QMEAN = mean(z$GPRA, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(qtab) <- NULL
qtab <- qtab[order(quarter_index(qtab$Quarter)), , drop = FALSE]

qtab$LN_GPR_QMEAN <- log(qtab$GPR_QMEAN)
qtab$LN_GPR_QMAX <- log(qtab$GPR_QMAX)
qtab$LN_GPRT_QMEAN <- log(qtab$GPRT_QMEAN)
qtab$LN_GPRA_QMEAN <- log(qtab$GPRA_QMEAN)

# Stable column order used by all downstream scripts.
qtab <- qtab[, c(
  "Quarter", "Months",
  "GPR_QMEAN", "LN_GPR_QMEAN",
  "GPR_QMAX", "LN_GPR_QMAX",
  "GPRT_QMEAN", "LN_GPRT_QMEAN",
  "GPRA_QMEAN", "LN_GPRA_QMEAN"
)]

write.csv(dat[, c("Date", "GPR", "GPRT", "GPRA", "Quarter")],
          "data/gpr_monthly_selected.csv", row.names = FALSE, na = "")
write.csv(qtab, "data/gpr_quarterly_all.csv", row.names = FALSE, na = "")

lo <- 2000L * 4L + 1L
hi <- 2026L * 4L + 2L
model <- qtab[
  quarter_index(qtab$Quarter) >= lo &
    quarter_index(qtab$Quarter) <= hi &
    qtab$Months == 3L,
  , drop = FALSE
]

if (!nrow(model)) stop("No complete quarterly GPR observations remain in 2000Q1-2026Q2.")
if (any(diff(quarter_index(model$Quarter)) != 1L)) {
  stop("Processed GPR model sample is not continuous.")
}

write.csv(model, "data/gpr_quarterly_processed.csv", row.names = FALSE, na = "")

summary_lines <- c(
  "Global GPR processing summary",
  paste0("Raw source file: ", gpr_file),
  "Source series: GPR (Recent Global Geopolitical Risk Index)",
  "Baseline transformation: monthly GPR -> quarterly arithmetic mean -> ln(GPR)",
  paste0("Monthly usable sample: ", format(min(dat$Date), "%Y-%m"), " - ", format(max(dat$Date), "%Y-%m")),
  paste0("Model export: ", model$Quarter[1], " - ", tail(model$Quarter, 1)),
  paste0("Complete quarters exported: ", nrow(model)),
  paste0("Baseline LN_GPR_QMEAN range: ", sprintf("%.6f", min(model$LN_GPR_QMEAN)), " - ", sprintf("%.6f", max(model$LN_GPR_QMEAN))),
  "Robustness columns retained: quarterly maximum GPR, GPRT, GPRA",
  "No winsorization, interpolation, forward-fill, or sign reversal was applied."
)
writeLines(summary_lines, "results/gpr_processing_summary.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")
