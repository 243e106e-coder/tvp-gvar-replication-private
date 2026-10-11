#!/usr/bin/env Rscript
# Read-only diagnostics for saved independent GPR/VIX/EPU/OIL IRF summaries.
# Source: 243e106e-coder/tvp-gvar-replication-private, run 37872702188.
# No sourcing estimation code, posterior loading, MCMC, or input rewriting.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[1L] else Sys.getenv('DIAG_INPUT', 'results/independent_shocks')
output <- if (length(args) >= 2L) args[2L] else Sys.getenv('DIAG_OUTPUT', 'results/four_shock_irf_diagnostics')
run_id <- Sys.getenv('SOURCE_RUN_ID', '37872702188')
if (!dir.exists(input)) stop('Input directory does not exist: ', input)
inpath <- normalizePath(input, mustWork = TRUE)
outpath <- normalizePath(output, mustWork = FALSE)
if (outpath == inpath || startsWith(paste0(outpath, '/'), paste0(inpath, '/'))) stop('Output must be outside input directory')
if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) stop('Output directory must be empty: ', output)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
shocks <- c('GPR', 'VIX', 'EPU', 'OIL')
countries <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
variables <- c(y = 'GDP (log level)', dp = 'Inflation', r = 'Interest rate', de = 'Exchange rate', deq = 'Equity price')
colors <- c(GPR = '#2166AC', VIX = '#B2182B', EPU = '#1B7837', OIL = '#762A83')
find_one <- function(risk, name, required = TRUE) {
  paths <- list.files(input, recursive = TRUE, full.names = TRUE)
  paths <- paths[basename(paths) == name]
  paths <- paths[vapply(strsplit(gsub('\\\\', '/', paths), '/'), function(parts) any(toupper(parts) %in% c(risk, paste0('INDEPENDENT-', risk, '-', run_id))), logical(1))]
  if (!length(paths) && !required) return(NULL)
  if (length(paths) != 1L) stop(risk, ': expected exactly one ', name, '; found ', length(paths))
  paths
}
read_table <- function(path, columns) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(columns %in% names(x))) stop('Missing columns in ', path, ': ', paste(setdiff(columns, names(x)), collapse = ', '))
  x
}
check_numeric <- function(x, cols, label) {
  for (col in cols) if (!is.numeric(x[[col]]) || any(!is.finite(x[[col]]))) stop(label, ': invalid numeric column ', col)
}
manifest <- list(); tables <- list(); stability <- list(); notes <- list()
for (risk in shocks) {
  sf <- find_one(risk, 'stability.csv')
  s <- read_table(sf, c('shock','date','total','stable','stable_share'))
  if (!nrow(s) || anyNA(s) || any(toupper(s$shock) != risk) || anyDuplicated(s$date) || any(!grepl('^[0-9]{4}Q[1-4]$', s$date))) stop('Invalid stability keys: ', risk)
  check_numeric(s, c('total','stable','stable_share'), risk)
  if (any(s$total <= 0 | s$stable < 0 | s$stable > s$total | s$total != floor(s$total) | s$stable != floor(s$stable)) || any(abs(s$stable_share - s$stable/s$total) > 1e-8)) stop('Inconsistent stability counts: ', risk)
  s$shock <- risk
  s$excluded_draws <- s$total - s$stable
  s$flag <- ifelse(s$stable == 0, 'NO_STABLE_DRAWS', ifelse(s$stable_share < .5, 'MAJORITY_EXCLUDED', ifelse(s$stable_share < .9, 'SOME_EXCLUDED', 'HIGH_RETENTION')))
  stability[[risk]] <- s
  f <- find_one(risk, 'irf_stable_only.csv', required = any(s$stable > 0))
  pf <- find_one(risk, 'provenance.txt')
  notes[[risk]] <- readLines(pf, warn = FALSE)
  if (!any(grepl('Only posterior draws with global spectral radius < 1 retained', notes[[risk]], fixed = TRUE))) stop('Stable-only provenance not confirmed: ', risk)
  files <- c(sf, pf, f)
  manifest[[risk]] <- data.frame(shock = risk, path = files, md5 = unname(tools::md5sum(files)))
  if (is.null(f)) next
  t <- read_table(f, c('shock','date','variable','horizon','median','low90','high90','stable_draws'))
  if (!nrow(t) || anyNA(t[,c('shock','date','variable')]) || any(toupper(t$shock) != risk)) stop('Invalid IRF keys: ', risk)
  check_numeric(t, c('horizon','median','low90','high90','stable_draws'), risk)
  if (any(t$horizon < 0 | t$horizon != floor(t$horizon)) || any(t$low90 > t$median | t$median > t$high90)) stop('Invalid horizons or interval order: ', risk)
  if (anyDuplicated(t[,c('date','variable','horizon')])) stop('Duplicate IRF keys: ', risk)
  idx <- match(t$date, s$date)
  if (anyNA(idx) || any(t$stable_draws != s$stable[idx]) || any(t$stable_draws <= 0)) stop('IRF/stability mismatch: ', risk)
  macro <- unlist(lapply(countries, function(cc) paste0(cc, '_', names(variables))))
  allowed <- c(macro, 'GL_shock')
  if (any(!t$variable %in% allowed)) stop('Unknown IRF variables: ', risk)
  for (d in s$date[s$stable > 0]) {
    td <- t[t$date == d, ]
    hs <- sort(unique(td$horizon))
    if (!length(hs) || !identical(as.numeric(hs), as.numeric(0:max(hs)))) stop('Missing horizons: ', risk, '/', d)
    expected <- expand.grid(variable = allowed, horizon = hs, stringsAsFactors = FALSE)
    if (!setequal(paste(td$variable, td$horizon), paste(expected$variable, expected$horizon))) stop('Incomplete 71-variable IRF grid: ', risk, '/', d)
  }
  t$shock <- risk
  t <- t[t$variable %in% macro, ]
  t$country <- sub('_.*$', '', t$variable)
  t$measure <- sub('^[^_]+_', '', t$variable)
  t$sign90 <- ifelse(t$low90 > 0, 'positive', ifelse(t$high90 < 0, 'negative', 'includes_zero'))
  tables[[risk]] <- t
}
# Require comparable saved model sample, lag/settings, and shock normalization.
for (prefix in c('Sample:', 'p:', 'IRF normalization:')) {
  values <- vapply(notes, function(z) {
    line <- z[startsWith(z, prefix)]
    if (length(line) != 1L) stop('Missing or ambiguous provenance: ', prefix)
    line
  }, character(1))
  if (length(unique(values)) != 1L) stop('Models have differing provenance: ', prefix)
}
ss <- do.call(rbind, stability); rownames(ss) <- NULL
mf <- do.call(rbind, manifest); rownames(mf) <- NULL
write.csv(ss, file.path(output, 'stability_diagnostics.csv'), row.names = FALSE)
write.csv(mf, file.path(output, 'source_manifest.csv'), row.names = FALSE)
if (!length(tables)) stop('All four shocks have zero stable draws; stability diagnostics saved, IRF comparison unavailable')
ir <- do.call(rbind, tables); rownames(ir) <- NULL
write.csv(ir, file.path(output, 'irf_comparison_90.csv'), row.names = FALSE)
# Coverage is explicit; absent observations are never filled with zero.
coverage <- expand.grid(shock = shocks, date = sort(unique(ss$date)), stringsAsFactors = FALSE)
coverage$stable_draws <- ss$stable[match(paste(coverage$shock, coverage$date), paste(ss$shock, ss$date))]
coverage$status <- ifelse(is.na(coverage$stable_draws), 'DATE_NOT_SAVED', ifelse(coverage$stable_draws == 0, 'NO_STABLE_DRAWS', 'AVAILABLE'))
coverage$max_horizon <- vapply(seq_len(nrow(coverage)), function(i) {
  z <- ir$horizon[ir$shock == coverage$shock[i] & ir$date == coverage$date[i]]
  if (length(z)) max(z) else NA_real_
}, numeric(1))
write.csv(coverage, file.path(output, 'coverage.csv'), row.names = FALSE)
# Summary extrema are descriptive; marginal intervals are not simultaneous bands.
keys <- interaction(ir$shock, ir$date, ir$country, ir$measure, drop = TRUE)
summary <- do.call(rbind, lapply(split(ir, keys), function(z) {
  z <- z[order(z$horizon), ]; post <- z[z$horizon > 0, ]
  peak <- which.max(abs(z$median))
  data.frame(shock = z$shock[1], date = z$date[1], country = z$country[1], measure = z$measure[1], impact = z$median[1], terminal_horizon = tail(z$horizon,1), terminal_median = tail(z$median,1), max_abs_median = z$median[peak], peak_horizon = z$horizon[peak], positive_horizons90 = sum(post$sign90 == 'positive'), negative_horizons90 = sum(post$sign90 == 'negative'), zero_included_horizons90 = sum(post$sign90 == 'includes_zero'), stable_draws = z$stable_draws[1])
}))
rownames(summary) <- NULL
write.csv(summary, file.path(output, 'country_variable_summary.csv'), row.names = FALSE)
gdp <- ir[ir$measure == 'y' & ir$horizon > 0, ]
gdpcounts <- do.call(rbind, lapply(split(gdp, interaction(gdp$shock,gdp$date,gdp$horizon,drop = TRUE)), function(z) data.frame(shock=z$shock[1], date=z$date[1], horizon=z$horizon[1], countries=nrow(z), positive_median=sum(z$median>0), negative_median=sum(z$median<0), positive90=sum(z$sign90=='positive'), negative90=sum(z$sign90=='negative'), includes_zero90=sum(z$sign90=='includes_zero'))))
write.csv(gdpcounts, file.path(output, 'gdp_sign_counts.csv'), row.names = FALSE)
pdf(file.path(output, 'four_shock_irf_comparison_90.pdf'), width = 13, height = 10, onefile = TRUE)
for (d in sort(unique(ss$date))) for (v in names(variables)) {
  par(mfrow = c(4,4), mar = c(2.6,3.2,2,1), oma = c(1,1,3,1))
  for (cc in countries) {
    z <- ir[ir$date == d & ir$measure == v & ir$country == cc, ]
    if (!nrow(z)) {plot.new(); title(cc); text(.5,.5,'No saved stable IRF'); next}
    ylim <- range(c(0,z$low90,z$high90)); if (diff(ylim) == 0) ylim <- ylim + c(-1,1)*1e-6
    xmax <- max(z$horizon); if(xmax == 0) xmax <- 1
    plot(NA, xlim = c(0,xmax), ylim = ylim, xlab = 'Horizon (quarters)', ylab = 'Saved model units', main = cc)
    for (risk in shocks) {
      zz <- z[z$shock == risk, ]; zz <- zz[order(zz$horizon), ]
      if (nrow(zz)) polygon(c(zz$horizon, rev(zz$horizon)), c(zz$low90,rev(zz$high90)), col=adjustcolor(colors[risk],alpha.f=.12), border=NA)
    }
    abline(h=0, col='grey40', lty=2)
    for (risk in shocks) {
      zz <- z[z$shock == risk, ]; zz <- zz[order(zz$horizon), ]
      if(nrow(zz)) lines(zz$horizon, zz$median, col=colors[risk], lwd=1.5)
    }
  }
  plot.new(); legend('center', legend=shocks, col=colors, lwd=2, bty='n')
  plot.new(); text(.5,.6,'Median + pointwise 90% interval'); text(.5,.4,'Conditional on stable draws', cex=.9)
  mtext(paste(d, variables[v]), outer=TRUE, side=3, line=1, cex=1.2)
}
dev.off()
pdf(file.path(output,'stability_retention.pdf'), width=11, height=7)
par(mar=c(7,4,3,1))
barplot(ss$stable_share, names.arg=paste(ss$shock,ss$date), las=2, ylim=c(0,1), col=colors[ss$shock], ylab='Stable draws / saved draws', main='Global spectral radius < 1: saved retention shares', cex.names=.7)
abline(h=c(.5,.9), lty=2, col='grey50'); dev.off()
report <- c('# Four independent shock IRF diagnostics', paste('Source repository: 243e106e-coder/tvp-gvar-replication-private; run:',run_id), '',
'Inputs: original irf_stable_only.csv, stability.csv, provenance.txt. No MCMC or original-data changes.',
'Variables: GDP y, inflation dp, interest rate r, exchange rate de, equity price deq; 14 economies.',
'90% equal-tailed pointwise credible intervals are the saved 5th/95th percentiles conditional on retained stable draws. They are not full-posterior or simultaneous intervals.',
'GDP is a log-level IRF: do not cumulate. All plots preserve original model units; no unsupported scaling to percentage points.',
'Exchange-rate appreciation interpretation must follow the original series definition, not the variable name alone.',
'Four independent models: comparing medians does not test whether differences between shocks are statistically credible. A common +10% log-level shock is not a common one-standard-deviation shock.',
'Stability diagnostics reuse saved counts for spectral radius < 1. Excluded draws may be unstable or nonfinite. No eigenvalue distribution or MCMC convergence diagnostics can be recovered from these summary CSVs.',
'Flags: MAJORITY_EXCLUDED means retention below 50%; SOME_EXCLUDED below 90%. These are descriptive flags, not formal validity thresholds.',
'Missing dates and zero stable draws are reported explicitly; no imputation. GDP counts exclude impact horizon 0.', '', '## Original provenance')
for (risk in shocks) report <- c(report, paste('###',risk), notes[[risk]], '')
writeLines(report, file.path(output,'README_diagnostics.md'))
if (!identical(unname(tools::md5sum(mf$path)), mf$md5)) stop('Input checksum changed during diagnostics')
writeLines(capture.output(sessionInfo()), file.path(output,'sessionInfo.txt'))
cat('Completed read-only diagnostics:', normalizePath(output), '\n')
