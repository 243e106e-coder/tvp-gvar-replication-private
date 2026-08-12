#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

options(stringsAsFactors = FALSE)

countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
macro_vars <- c("y","dp","de","r","deq")

macro_file <- Sys.getenv(
  "TVPGVAR_MACRO_XLSX",
  "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"
)
weight_file <- Sys.getenv(
  "TVPGVAR_WEIGHT_XLSX",
  "8.12/Trade_Weights_14_Economies_Stage1_2000_2014.xlsx"
)
gpr_excel <- Sys.getenv("TVPGVAR_GPR_XLSX", "GPR处理完.xlsx")
old_model_file <- Sys.getenv("TVPGVAR_OLD_MODEL", "data/model_input.csv")

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

if (!file.exists(macro_file)) stop("Macro workbook not found: ", macro_file)
if (!file.exists(weight_file)) stop("Trade-weight workbook not found: ", weight_file)

norm_token <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  z <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  z[is.na(z)] <- x[is.na(z)]
  toupper(gsub("[^A-Z0-9]+", "", z))
}

country_aliases <- list(
  AU = c("AU","AUS","AUSTRALIA"),
  BR = c("BR","BRA","BRAZIL"),
  CA = c("CA","CAN","CANADA"),
  CH = c("CH","CHE","SWITZERLAND","SWISS"),
  CN = c("CN","CHN","CHINA"),
  EA = c("EA","EUROAREA","EUROZONE","EMU","EA20","EA19"),
  UK = c("UK","GB","GBR","UNITEDKINGDOM","BRITAIN"),
  JP = c("JP","JPN","JAPAN"),
  KR = c("KR","KOR","KOREA","SOUTHKOREA","REPUBLICOFKOREA"),
  NO = c("NO","NOR","NORWAY"),
  SG = c("SG","SGP","SINGAPORE"),
  TR = c("TR","TUR","TURKEY","TURKIYE"),
  US = c("US","USA","UNITEDSTATES","UNITEDSTATESOFAMERICA"),
  ZA = c("ZA","ZAF","SOUTHAFRICA")
)
country_aliases <- lapply(country_aliases, norm_token)

var_aliases <- list(
  y = c(
    "Y","GDP","REALGDP","RGDP","REALGDPGROWTH","GDPGROWTH",
    "DGDP","DLGDP","DLOGGDP","OUTPUTGROWTH"
  ),
  dp = c(
    "DP","CPI","INFLATION","INFL","DCPI","DLCPI","DLOGCPI",
    "CPIGROWTH","INFLATIONRATE"
  ),
  de = c(
    "DE","REER","RER","EXCHANGERATE","FX","DREER","DLREER",
    "DLOGREER","EXCHANGERATEGROWTH"
  ),
  r = c(
    "R","RATE","INTEREST","INTERESTRATE","STIR","SHORTTERM",
    "SHORTTERMINTEREST","SHORTTERMINTERESTRATE","POLICYRATE"
  ),
  deq = c(
    "DEQ","EQ","EQUITY","EQUITYPRICE","STOCK","STOCKPRICE",
    "STOCKINDEX","DSP","DLEQ","DLOGEQ","STOCKRETURN","EQUITYRETURN"
  )
)
var_aliases <- lapply(var_aliases, norm_token)

quarter_aliases <- norm_token(c("Quarter","Period","Time","Date","QuarterDate"))
year_aliases <- norm_token(c("Year","YYYY"))
month_aliases <- norm_token(c("Month","MM"))
country_col_aliases <- norm_token(c("Country","Economy","Area","Code","ISO","REF_AREA"))
partner_col_aliases <- norm_token(c("Partner","Counterpart","PartnerCountry"))
weight_col_aliases <- norm_token(c("Weight","TradeWeight","W","Share"))

canon_country <- function(x) {
  z <- norm_token(x)
  out <- rep(NA_character_, length(z))
  for (cc in names(country_aliases)) {
    out[z %in% country_aliases[[cc]]] <- cc
  }
  out
}

canon_var <- function(x) {
  z <- norm_token(x)
  out <- rep(NA_character_, length(z))
  for (v in names(var_aliases)) {
    out[z %in% var_aliases[[v]]] <- v
  }
  out
}

canon_macro_header <- function(x) {
  z <- norm_token(x)
  out <- rep(NA_character_, length(z))

  for (ii in seq_along(z)) {
    if (!nzchar(z[ii])) next

    hits <- character()
    for (cc in countries) {
      for (ca in country_aliases[[cc]]) {
        for (v in macro_vars) {
          for (va in var_aliases[[v]]) {
            if (z[ii] %in% c(paste0(ca, va), paste0(va, ca))) {
              hits <- c(hits, paste0(cc, "_", v))
            }
          }
        }
      }
    }
    hits <- unique(hits)
    if (length(hits) == 1) out[ii] <- hits
  }
  out
}

quarter_index <- function(q) {
  m <- regexec("^([12][0-9]{3})Q([1-4])$", q)
  p <- regmatches(q, m)
  vapply(p, function(z) {
    if (length(z) != 3) return(NA_integer_)
    as.integer(z[2]) * 4L + as.integer(z[3])
  }, integer(1))
}

date_to_quarter <- function(d) {
  y <- as.integer(format(d, "%Y"))
  m <- as.integer(format(d, "%m"))
  paste0(y, "Q", ((m - 1L) %/% 3L) + 1L)
}

as_quarter <- function(x, year = NULL, month = NULL) {
  if (!is.null(year) && !is.null(month)) {
    yy <- suppressWarnings(as.integer(year))
    mm <- suppressWarnings(as.integer(month))
    return(ifelse(
      is.na(yy) | is.na(mm) | mm < 1 | mm > 12,
      NA_character_,
      paste0(yy, "Q", ((mm - 1L) %/% 3L) + 1L)
    ))
  }

  if (!is.null(year)) {
    # If x itself is quarter number 1..4.
    yy <- suppressWarnings(as.integer(year))
    qq <- suppressWarnings(as.integer(x))
    if (sum(!is.na(qq) & qq %in% 1:4) >= max(1, floor(length(qq) * 0.8))) {
      return(ifelse(
        is.na(yy) | is.na(qq) | !(qq %in% 1:4),
        NA_character_,
        paste0(yy, "Q", qq)
      ))
    }
  }

  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    return(date_to_quarter(as.Date(x)))
  }

  if (is.numeric(x)) {
    # Excel serial dates are usually > 10,000.
    if (sum(x > 10000, na.rm = TRUE) >= max(1, floor(sum(!is.na(x)) * 0.8))) {
      d <- as.Date(x, origin = "1899-12-30")
      return(date_to_quarter(d))
    }
  }

  s <- toupper(trimws(as.character(x)))
  s <- gsub("\\s+", "", s)

  out <- rep(NA_character_, length(s))

  rx1 <- regexec("^([12][0-9]{3})[-_/]?Q([1-4])$", s)
  rr1 <- regmatches(s, rx1)
  for (i in seq_along(rr1)) {
    if (length(rr1[[i]]) == 3) out[i] <- paste0(rr1[[i]][2], "Q", rr1[[i]][3])
  }

  rx2 <- regexec("^Q([1-4])[-_/]?([12][0-9]{3})$", s)
  rr2 <- regmatches(s, rx2)
  for (i in seq_along(rr2)) {
    if (is.na(out[i]) && length(rr2[[i]]) == 3) {
      out[i] <- paste0(rr2[[i]][3], "Q", rr2[[i]][2])
    }
  }

  missing <- which(is.na(out) & nzchar(s))
  if (length(missing)) {
    date_formats <- c("%Y-%m-%d","%Y/%m/%d","%d/%m/%Y","%m/%d/%Y","%Y-%m","%Y/%m")
    for (fmt in date_formats) {
      d <- suppressWarnings(as.Date(s[missing], format = fmt))
      ok <- !is.na(d)
      if (any(ok)) {
        out[missing[ok]] <- date_to_quarter(d[ok])
        missing <- which(is.na(out) & nzchar(s))
        if (!length(missing)) break
      }
    }
  }

  out
}

row_score <- function(row_values) {
  z <- norm_token(row_values)
  sum(!is.na(canon_country(row_values))) +
    sum(z %in% quarter_aliases) +
    sum(z %in% year_aliases) +
    sum(z %in% country_col_aliases) +
    sum(!is.na(canon_var(row_values)))
}

read_sheet_guess_header <- function(path, sheet) {
  preview <- suppressWarnings(
    read_excel(path, sheet = sheet, col_names = FALSE, n_max = 30)
  )
  if (!nrow(preview)) return(data.frame())

  scores <- vapply(seq_len(nrow(preview)), function(i) {
    row_score(unlist(preview[i, ], use.names = FALSE))
  }, numeric(1))

  best <- which.max(scores)
  if (!length(best) || scores[best] < 1) best <- 1L

  tab <- suppressWarnings(
    read_excel(path, sheet = sheet, skip = best - 1L, .name_repair = "unique")
  )
  tab <- as.data.frame(tab, check.names = FALSE)

  # Drop completely empty rows/columns.
  if (nrow(tab)) {
    keep_rows <- apply(tab, 1, function(r) any(!is.na(r) & trimws(as.character(r)) != ""))
    tab <- tab[keep_rows, , drop = FALSE]
  }
  if (ncol(tab)) {
    keep_cols <- vapply(tab, function(v) any(!is.na(v) & trimws(as.character(v)) != ""), logical(1))
    tab <- tab[, keep_cols, drop = FALSE]
  }

  tab
}

find_named_col <- function(nms, aliases) {
  z <- norm_token(nms)
  which(z %in% aliases)[1]
}

extract_quarter <- function(tab) {
  nms <- names(tab)
  qcol <- find_named_col(nms, quarter_aliases)
  ycol <- find_named_col(nms, year_aliases)
  mcol <- find_named_col(nms, month_aliases)

  if (!is.na(qcol)) {
    if (!is.na(ycol) && qcol != ycol) {
      q_try <- as_quarter(tab[[qcol]], year = tab[[ycol]])
      if (sum(!is.na(q_try)) >= max(1, floor(nrow(tab) * 0.5))) return(q_try)
    }
    return(as_quarter(tab[[qcol]]))
  }

  # Fall back to any column that looks mostly like quarterly/date data.
  for (j in seq_along(tab)) {
    q <- if (!is.na(ycol) && !is.na(mcol) && j == mcol) {
      as_quarter(tab[[j]], year = tab[[ycol]], month = tab[[mcol]])
    } else {
      as_quarter(tab[[j]])
    }
    if (sum(!is.na(q)) >= max(3, floor(nrow(tab) * 0.7))) return(q)
  }

  rep(NA_character_, nrow(tab))
}

to_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  z <- gsub(",", "", trimws(as.character(x)))
  z[z %in% c("", "NA", "N/A", "..", "-", "—")] <- NA_character_
  suppressWarnings(as.numeric(z))
}

extract_macro_ready_wide <- function(path) {
  req <- unlist(lapply(countries, function(cc) paste0(cc, "_", macro_vars)), use.names = FALSE)
  best <- NULL
  best_count <- -1

  for (sh in excel_sheets(path)) {
    tab <- read_sheet_guess_header(path, sh)
    if (!nrow(tab) || !ncol(tab)) next

    mapped <- canon_macro_header(names(tab))
    count <- sum(!is.na(mapped) & mapped %in% req)

    if (count > best_count) {
      best_count <- count
      best <- list(sheet = sh, tab = tab, mapped = mapped)
    }
  }

  if (is.null(best) || best_count < length(req)) return(NULL)

  q <- extract_quarter(best$tab)
  if (sum(!is.na(q)) < 3) return(NULL)

  out <- data.frame(Quarter = q, stringsAsFactors = FALSE)

  for (target in req) {
    idx <- which(best$mapped == target)
    if (length(idx) != 1) return(NULL)
    out[[target]] <- to_numeric(best$tab[[idx]])
  }

  message("Macro layout detected: ready wide table on sheet '", best$sheet, "'.")
  out
}

extract_macro_long <- function(path) {
  req <- unlist(lapply(countries, function(cc) paste0(cc, "_", macro_vars)), use.names = FALSE)

  for (sh in excel_sheets(path)) {
    tab <- read_sheet_guess_header(path, sh)
    if (!nrow(tab) || ncol(tab) < 7) next

    ccol <- find_named_col(names(tab), country_col_aliases)
    if (is.na(ccol)) next

    cc <- canon_country(tab[[ccol]])
    if (sum(!is.na(cc)) < 10) next

    var_map <- canon_var(names(tab))
    if (!all(macro_vars %in% var_map)) next

    q <- extract_quarter(tab)
    if (sum(!is.na(q)) < 3) next

    out <- data.frame(Quarter = sort(unique(q[!is.na(q)])), stringsAsFactors = FALSE)

    for (country in countries) {
      for (v in macro_vars) {
        target <- paste0(country, "_", v)
        vidx <- which(var_map == v)[1]
        sub <- data.frame(
          Quarter = q[cc == country],
          value = to_numeric(tab[[vidx]][cc == country]),
          stringsAsFactors = FALSE
        )
        sub <- sub[!is.na(sub$Quarter), , drop = FALSE]
        if (anyDuplicated(sub$Quarter)) {
          stop("Duplicate country-quarter observations in macro workbook: ", country, " / ", v)
        }
        out[[target]] <- sub$value[match(out$Quarter, sub$Quarter)]
      }
    }

    if (all(req %in% names(out))) {
      message("Macro layout detected: long table on sheet '", sh, "'.")
      return(out)
    }
  }

  NULL
}

sheet_to_var <- function(sheet_name) {
  z <- norm_token(sheet_name)

  # Exact aliases first.
  v <- canon_var(sheet_name)
  if (!is.na(v)) return(v)

  patterns <- list(
    y = c("GDP","REALGDP","RGDP","OUTPUT"),
    dp = c("CPI","INFLATION","INFL"),
    de = c("REER","EXCHANGE","FX","RER"),
    r = c("INTEREST","SHORTTERM","RATE","STIR"),
    deq = c("EQUITY","STOCK","SHAREPRICE","DEQ")
  )
  for (vv in names(patterns)) {
    if (any(vapply(patterns[[vv]], function(p) grepl(p, z, fixed = TRUE), logical(1)))) {
      return(vv)
    }
  }
  NA_character_
}

extract_macro_separate_sheets <- function(path) {
  pieces <- list()

  for (sh in excel_sheets(path)) {
    v <- sheet_to_var(sh)
    if (is.na(v)) next

    tab <- read_sheet_guess_header(path, sh)
    if (!nrow(tab)) next

    q <- extract_quarter(tab)
    if (sum(!is.na(q)) < 3) next

    cc_map <- canon_country(names(tab))
    if (sum(!is.na(cc_map)) < 10) next

    tmp <- data.frame(Quarter = q, stringsAsFactors = FALSE)
    for (cc in countries) {
      idx <- which(cc_map == cc)
      if (length(idx) == 1) tmp[[paste0(cc, "_", v)]] <- to_numeric(tab[[idx]])
    }

    pieces[[v]] <- tmp
  }

  if (!all(macro_vars %in% names(pieces))) return(NULL)

  out <- pieces[[macro_vars[1]]]
  for (v in macro_vars[-1]) {
    out <- merge(out, pieces[[v]], by = "Quarter", all = TRUE, sort = FALSE)
  }

  req <- unlist(lapply(countries, function(cc) paste0(cc, "_", macro_vars)), use.names = FALSE)
  if (!all(req %in% names(out))) return(NULL)

  message("Macro layout detected: five variable sheets.")
  out
}


# -------------------------------------------------------------------------
# DIRECT MACRO IMPORT
# The workbook already provides a model-ready balanced sheet:
# MODEL_COMPLETE_14C.
# Do not guess workbook layouts; read this sheet directly and map the
# source variable names to the names expected by the TVP-GVAR code.
# -------------------------------------------------------------------------
extract_macro_direct <- function(path,
                                 sheet = "MODEL_COMPLETE_14C") {
  if (!file.exists(path)) stop("Macro workbook not found: ", path)

  sheets <- excel_sheets(path)
  if (!sheet %in% sheets) {
    stop(
      "Required macro sheet '", sheet, "' was not found. Available sheets: ",
      paste(sheets, collapse = ", ")
    )
  }

  tab <- suppressWarnings(
    read_excel(path, sheet = sheet, .name_repair = "minimal")
  )
  tab <- as.data.frame(tab, check.names = FALSE, stringsAsFactors = FALSE)
  names(tab) <- trimws(names(tab))

  source_suffix <- c(
    y   = "GDP_DLOG",
    dp  = "CPI_DLOG",
    de  = "REER_DLOG",
    r   = "RATE_LEVEL",
    deq = "EQ_RETURN"
  )

  required_source <- c(
    "Quarter",
    unlist(lapply(countries, function(cc) {
      paste0(cc, "_", unname(source_suffix))
    }), use.names = FALSE)
  )

  missing_cols <- setdiff(required_source, names(tab))
  if (length(missing_cols)) {
    stop(
      "MODEL_COMPLETE_14C is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  out <- data.frame(
    Quarter = as.character(tab[["Quarter"]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Explicit mapping:
  # GDP_DLOG   -> y
  # CPI_DLOG   -> dp
  # REER_DLOG  -> de
  # RATE_LEVEL -> r
  # EQ_RETURN  -> deq
  for (cc in countries) {
    out[[paste0(cc, "_y")]]   <- to_numeric(tab[[paste0(cc, "_GDP_DLOG")]])
    out[[paste0(cc, "_dp")]]  <- to_numeric(tab[[paste0(cc, "_CPI_DLOG")]])
    out[[paste0(cc, "_de")]]  <- to_numeric(tab[[paste0(cc, "_REER_DLOG")]])
    out[[paste0(cc, "_r")]]   <- to_numeric(tab[[paste0(cc, "_RATE_LEVEL")]])
    out[[paste0(cc, "_deq")]] <- to_numeric(tab[[paste0(cc, "_EQ_RETURN")]])
  }

  if (anyDuplicated(out$Quarter)) {
    dup <- unique(out$Quarter[duplicated(out$Quarter)])
    stop(
      "Duplicate quarter(s) in MODEL_COMPLETE_14C: ",
      paste(dup, collapse = ", ")
    )
  }

  required_model <- unlist(
    lapply(countries, function(cc) paste0(cc, "_", macro_vars)),
    use.names = FALSE
  )

  if (!all(required_model %in% names(out))) {
    stop("Direct macro mapping failed to produce all 70 model variables.")
  }

  message(
    "Macro layout: direct import from '", sheet, "'. ",
    "Rows=", nrow(out), "; macro columns=", length(required_model), "."
  )

  out
}

extract_macro <- function(path) {
  out <- extract_macro_ready_wide(path)
  if (!is.null(out)) return(out)

  out <- extract_macro_long(path)
  if (!is.null(out)) return(out)

  out <- extract_macro_separate_sheets(path)
  if (!is.null(out)) return(out)

  stop(
    "Could not recognize the macro workbook layout. Expected either: ",
    "(1) one wide sheet with columns such as AU_y/AU_dp/...; ",
    "(2) one long sheet with Quarter, Country and y/dp/de/r/deq; or ",
    "(3) separate GDP/CPI/REER/rate/equity sheets with country columns."
  )
}

extract_weight_matrix <- function(path) {
  best <- NULL
  best_score <- -Inf

  for (sh in excel_sheets(path)) {
    tab <- read_sheet_guess_header(path, sh)
    if (!nrow(tab) || ncol(tab) < 10) next

    col_cc <- canon_country(names(tab))
    if (!all(countries %in% col_cc)) next

    row_scores <- vapply(seq_along(tab), function(j) {
      sum(countries %in% canon_country(tab[[j]]))
    }, numeric(1))
    row_col <- which.max(row_scores)
    if (!length(row_col) || row_scores[row_col] < length(countries)) next

    row_cc <- canon_country(tab[[row_col]])
    if (!all(countries %in% row_cc)) next

    mat <- matrix(NA_real_, length(countries), length(countries),
                  dimnames = list(countries, countries))

    ok <- TRUE
    for (r in countries) {
      ridx <- which(row_cc == r)
      if (length(ridx) != 1) { ok <- FALSE; break }
      for (c in countries) {
        cidx <- which(col_cc == c)
        if (length(cidx) != 1) { ok <- FALSE; break }
        mat[r, c] <- to_numeric(tab[[cidx]][ridx])
      }
    }
    if (!ok || any(!is.finite(mat))) next
    if (any(mat < -1e-12)) next

    # Prefer matrices that already look like weights and have small diagonal.
    diag_pen <- mean(abs(diag(mat)))
    rs <- rowSums(mat)
    row_pen <- if (all(rs > 0)) mean(abs(rs / mean(rs) - 1)) else 1e6
    weight_like <- if (all(rs > 0)) mean(abs(rs - 1)) else 1e6
    score <- -diag_pen - min(weight_like, row_pen)

    if (score > best_score) {
      best_score <- score
      best <- list(sheet = sh, mat = mat)
    }
  }

  if (is.null(best)) {
    # Long-format fallback: Reporter/Country, Partner, Weight
    for (sh in excel_sheets(path)) {
      tab <- read_sheet_guess_header(path, sh)
      if (!nrow(tab) || ncol(tab) < 3) next

      ccol <- find_named_col(names(tab), country_col_aliases)
      pcol <- find_named_col(names(tab), partner_col_aliases)
      wcol <- find_named_col(names(tab), weight_col_aliases)
      if (any(is.na(c(ccol, pcol, wcol)))) next

      rr <- canon_country(tab[[ccol]])
      pp <- canon_country(tab[[pcol]])
      ww <- to_numeric(tab[[wcol]])

      mat <- matrix(0, length(countries), length(countries),
                    dimnames = list(countries, countries))
      for (i in seq_along(ww)) {
        if (!is.na(rr[i]) && !is.na(pp[i]) && is.finite(ww[i])) {
          mat[rr[i], pp[i]] <- mat[rr[i], pp[i]] + ww[i]
        }
      }
      if (all(rowSums(mat) > 0)) {
        best <- list(sheet = sh, mat = mat)
        break
      }
    }
  }

  if (is.null(best)) {
    stop("Could not locate a complete 14-economy trade matrix in: ", path)
  }

  mat <- best$mat
  if (any(mat < -1e-12)) stop("Trade matrix contains negative values.")

  # Enforce standard GVAR convention: no home-country weight.
  diag(mat) <- 0
  rs <- rowSums(mat)
  if (any(rs <= 0)) stop("At least one trade-weight row sums to zero after removing the diagonal.")
  mat <- mat / rs

  message("Trade matrix detected on sheet '", best$sheet, "'.")
  mat
}

extract_gpr_excel <- function(path) {
  if (!file.exists(path)) return(NULL)

  best <- NULL
  best_priority <- -Inf

  for (sh in excel_sheets(path)) {
    tab <- read_sheet_guess_header(path, sh)
    if (!nrow(tab) || ncol(tab) < 2) next

    nms_norm <- norm_token(names(tab))
    gpr_cols <- which(grepl("GPR", nms_norm))
    if (!length(gpr_cols)) next

    q <- extract_quarter(tab)
    if (sum(!is.na(q)) < 3) next

    for (j in gpr_cols) {
      vals <- to_numeric(tab[[j]])
      ok <- !is.na(q) & is.finite(vals)
      if (sum(ok) < 3) next

      nm <- nms_norm[j]
      explicit_log <- grepl("LOG|LN", nm) || nm %in% c("USGPR","GPRLOG","LOGGPR","LNGPR")
      priority <- if (explicit_log) 10 else 1

      tmp <- aggregate(vals[ok], by = list(Quarter = q[ok]), FUN = mean, na.rm = TRUE)
      names(tmp)[2] <- "value"

      # README/model convention is log quarterly GPR.
      if (!explicit_log) {
        med <- median(tmp$value, na.rm = TRUE)
        if (is.finite(med) && med > 20) {
          if (any(tmp$value <= 0)) next
          tmp$value <- log(tmp$value)
        }
      }

      if (priority > best_priority) {
        best_priority <- priority
        best <- list(sheet = sh, data = tmp)
      }
    }
  }

  if (!is.null(best)) {
    names(best$data)[2] <- "US_gpr"
    message("GPR detected in '", path, "', sheet '", best$sheet, "'.")
    return(best$data)
  }

  NULL
}

extract_gpr_old_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  old <- try(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(old, "try-error")) return(NULL)
  if (!all(c("Quarter", "US_gpr") %in% names(old))) return(NULL)

  data.frame(
    Quarter = as.character(old$Quarter),
    US_gpr = to_numeric(old$US_gpr),
    stringsAsFactors = FALSE
  )
}

combine_gpr_sources <- function(excel_gpr, old_gpr) {
  if (is.null(excel_gpr) && is.null(old_gpr)) {
    stop(
      "No usable GPR source found. Keep either GPR处理完.xlsx or an old ",
      "data/model_input.csv containing Quarter and US_gpr."
    )
  }

  if (is.null(excel_gpr)) return(old_gpr)
  if (is.null(old_gpr)) return(excel_gpr)

  # Prefer the dedicated GPR workbook; use the old model only to fill holes.
  q <- sort(unique(c(excel_gpr$Quarter, old_gpr$Quarter)))
  out <- data.frame(Quarter = q, US_gpr = NA_real_, stringsAsFactors = FALSE)

  i1 <- match(out$Quarter, excel_gpr$Quarter)
  i2 <- match(out$Quarter, old_gpr$Quarter)

  ok1 <- !is.na(i1)
  out$US_gpr[ok1] <- excel_gpr$US_gpr[i1[ok1]]

  ok2 <- is.na(out$US_gpr) & !is.na(i2)
  out$US_gpr[ok2] <- old_gpr$US_gpr[i2[ok2]]

  # Diagnostic only: large overlap discrepancies usually indicate a transform mismatch.
  common <- intersect(excel_gpr$Quarter, old_gpr$Quarter)
  if (length(common) >= 8) {
    a <- excel_gpr$US_gpr[match(common, excel_gpr$Quarter)]
    b <- old_gpr$US_gpr[match(common, old_gpr$Quarter)]
    mad <- median(abs(a - b), na.rm = TRUE)
    if (is.finite(mad) && mad > 0.25) {
      warning(
        "Dedicated GPR workbook and old US_gpr differ materially on overlapping quarters ",
        "(median absolute difference = ", round(mad, 3), "). ",
        "The dedicated GPR workbook was preferred."
      )
    }
  }

  out
}

validate_and_order_macro <- function(macro) {
  req <- unlist(lapply(countries, function(cc) paste0(cc, "_", macro_vars)), use.names = FALSE)
  if (!all(c("Quarter", req) %in% names(macro))) {
    stop("Macro extraction did not produce all 70 required country-variable columns.")
  }

  macro <- macro[, c("Quarter", req), drop = FALSE]
  macro <- macro[!is.na(macro$Quarter), , drop = FALSE]

  if (anyDuplicated(macro$Quarter)) {
    dup <- unique(macro$Quarter[duplicated(macro$Quarter)])
    stop("Duplicate quarters in macro data: ", paste(dup, collapse = ", "))
  }

  idx <- quarter_index(macro$Quarter)
  if (anyNA(idx)) stop("Unrecognized quarter labels remain in macro data.")
  macro <- macro[order(idx), , drop = FALSE]

  # Coerce values.
  for (nm in req) macro[[nm]] <- to_numeric(macro[[nm]])

  macro
}

macro <- validate_and_order_macro(extract_macro_direct(macro_file, "MODEL_COMPLETE_14C"))
W <- extract_weight_matrix(weight_file)

gpr1 <- extract_gpr_excel(gpr_excel)
gpr2 <- extract_gpr_old_csv(old_model_file)
gpr <- combine_gpr_sources(gpr1, gpr2)

gpr <- gpr[!is.na(gpr$Quarter), , drop = FALSE]
gpr <- gpr[!duplicated(gpr$Quarter), , drop = FALSE]

model <- merge(macro, gpr, by = "Quarter", all.x = TRUE, sort = FALSE)
model <- model[match(macro$Quarter, model$Quarter), , drop = FALSE]

# GPR must cover the macro sample. Do not silently forward-fill or invent values.
missing_gpr <- model$Quarter[!is.finite(model$US_gpr)]
if (length(missing_gpr)) {
  stop(
    "US_gpr does not cover the full macro sample. Missing quarter(s): ",
    paste(missing_gpr, collapse = ", "),
    ". Update only the missing GPR tail; do not forward-fill it."
  )
}

ordered_cols <- unlist(lapply(countries, function(cc) {
  own <- paste0(cc, "_", macro_vars)
  if (cc == "US") c(own, "US_gpr") else own
}), use.names = FALSE)

model <- model[, c("Quarter", ordered_cols), drop = FALSE]

# Allow missing values only at the very beginning (typically created by log differences).
complete <- complete.cases(model)
if (!any(complete)) stop("No complete quarter remains after merging all variables.")

first_complete <- which(complete)[1]
if (first_complete > 1) {
  message(
    "Dropping ", first_complete - 1L,
    " leading incomplete quarter(s): ",
    paste(model$Quarter[seq_len(first_complete - 1L)], collapse = ", ")
  )
  model <- model[first_complete:nrow(model), , drop = FALSE]
}

if (any(!complete.cases(model))) {
  bad_rows <- which(!complete.cases(model))
  bad_q <- model$Quarter[bad_rows]
  bad_cols <- unique(unlist(lapply(bad_rows, function(i) {
    names(model)[which(is.na(model[i, , drop = TRUE]))]
  })))
  stop(
    "Internal/trailing missing values remain. Quarters: ",
    paste(bad_q, collapse = ", "),
    "; columns: ", paste(bad_cols, collapse = ", ")
  )
}

# Require a continuous quarterly sample.
idx <- quarter_index(model$Quarter)
if (any(diff(idx) != 1L)) {
  gap_at <- which(diff(idx) != 1L)[1]
  stop(
    "Quarterly sample is not continuous between ",
    model$Quarter[gap_at], " and ", model$Quarter[gap_at + 1L], "."
  )
}

write.csv(model, "data/model_input.csv", row.names = FALSE, na = "")

wout <- data.frame(Country = countries, W, check.names = FALSE)
write.csv(wout, "data/trade_weights.csv", row.names = FALSE, na = "")

summary_lines <- c(
  "TVP-GVAR 14-economy input build",
  paste0("Macro workbook: ", macro_file),
  paste0("Weight workbook: ", weight_file),
  paste0("Countries: ", paste(countries, collapse = ", ")),
  paste0("Variables: ", paste(macro_vars, collapse = ", ")),
  paste0("Final sample: ", model$Quarter[1], " - ", tail(model$Quarter, 1)),
  paste0("Observations: ", nrow(model)),
  paste0("Global variables: ", ncol(model) - 1L),
  paste0("Trade row-sum range: ", sprintf("%.12f", min(rowSums(W))),
         " - ", sprintf("%.12f", max(rowSums(W)))),
  paste0("Trade diagonal max abs: ", sprintf("%.12f", max(abs(diag(W))))),
  paste0("US_gpr range: ", sprintf("%.6f", min(model$US_gpr)),
         " - ", sprintf("%.6f", max(model$US_gpr)))
)
writeLines(summary_lines, "results/build_input_summary.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")

message("Wrote data/model_input.csv")
message("Wrote data/trade_weights.csv")
message("Input validation passed.")
