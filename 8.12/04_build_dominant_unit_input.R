#!/usr/bin/env Rscript

# =============================================================================
# Build a semantically clean input for the separate dominant-unit experiment.
#
# Upstream data/model_input.csv currently stores the global risk series under
# compatibility names US_gpr and US_vix.  This script moves/renames them to:
#   GL_gpr, GL_vix
# and leaves the US country block with ONLY its five macro-financial variables.
#
# No data values are altered here; this is a naming/reordering operation only.
# =============================================================================

input_file <- Sys.getenv("TVPGVAR_BASE_INPUT", "data/model_input.csv")
output_file <- Sys.getenv("TVPGVAR_DOMINANT_INPUT", "data/model_input_dominant.csv")

if (!file.exists(input_file)) stop("Missing input file: ", input_file)

dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Quarter" %in% names(dat)) stop("Input must contain Quarter.")
if (anyDuplicated(dat$Quarter)) stop("Duplicate Quarter values detected.")

required_risk <- c("US_gpr", "US_vix")
missing_risk <- setdiff(required_risk, names(dat))
if (length(missing_risk)) {
  stop("Missing upstream risk columns: ", paste(missing_risk, collapse = ", "))
}

if (any(c("GL_gpr", "GL_vix") %in% names(dat))) {
  stop("GL_gpr/GL_vix already exist; refusing ambiguous duplicate risk columns.")
}

names(dat)[match("US_gpr", names(dat))] <- "GL_gpr"
names(dat)[match("US_vix", names(dat))] <- "GL_vix"

countries <- c(
  "AU","BR","CA","CH","CN","EA","UK",
  "JP","KR","NO","SG","TR","US","ZA"
)
macro_vars <- c("y", "dp", "r", "de", "deq")
macro_cols <- unlist(
  lapply(countries, function(cc) paste0(cc, "_", macro_vars)),
  use.names = FALSE
)

missing_macro <- setdiff(macro_cols, names(dat))
if (length(missing_macro)) {
  stop("Missing macro columns: ", paste(missing_macro, collapse = ", "))
}

keep <- c("Quarter", macro_cols, "GL_gpr", "GL_vix")
out <- dat[, keep, drop = FALSE]

num <- as.matrix(out[, -1, drop = FALSE])
storage.mode(num) <- "double"
if (any(!is.finite(num))) stop("Dominant-unit input contains NA/NaN/Inf.")

if (ncol(num) != 72L) {
  stop("Expected 72 model variables (14x5 + 2 dominant), found ", ncol(num))
}

us_cols <- names(out)[startsWith(names(out), "US_")]
if (!identical(us_cols, paste0("US_", macro_vars))) {
  stop("US block is not macro-only after risk-variable extraction: ",
       paste(us_cols, collapse = ", "))
}

if (!identical(tail(names(out), 2L), c("GL_gpr", "GL_vix"))) {
  stop("Dominant-unit columns must be last and ordered GL_gpr -> GL_vix.")
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.csv(out, output_file, row.names = FALSE)

cat("Created dominant-unit input: ", output_file, "\n", sep = "")
cat("Rows: ", nrow(out), "\n", sep = "")
cat("Variables: ", ncol(out) - 1L, "\n", sep = "")
cat("US block: ", paste(us_cols, collapse = ", "), "\n", sep = "")
cat("Dominant block: GL_gpr, GL_vix\n")
