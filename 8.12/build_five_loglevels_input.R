#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(readxl))
source_file <- '8.12/TVPGVAR_14C_5VAR_LOG_LEVELS.xlsx'
stopifnot(file.exists(source_file))
countries <- c('AU','BR','CA','CH','CN','EA','UK','JP','KR','NO','SG','TR','US','ZA')
vars <- c('y','p','r','q','eq'); mapped <- c('y','dp','r','de','deq')
x <- as.data.frame(readxl::read_excel(source_file,sheet='MODEL_BALANCED',.name_repair='minimal'),check.names=FALSE)
expected <- c('Quarter',unlist(lapply(countries,function(cc) paste0(cc,'_',vars))))
if(!identical(names(x),expected)) stop('Unexpected workbook schema: ',paste(names(x),collapse=','))
x$Quarter <- as.character(x$Quarter)
if(anyDuplicated(x$Quarter)) stop('Duplicate quarters')
qi <- as.integer(substr(x$Quarter,1,4))*4L+as.integer(substr(x$Quarter,6,6))
if(any(!grepl('^[12][0-9]{3}Q[1-4]$',x$Quarter))||any(diff(qi)!=1)) stop('Invalid/nonconsecutive quarters')
for(n in names(x)[-1]) {x[[n]] <- suppressWarnings(as.numeric(x[[n]])); if(any(!is.finite(x[[n]]))) stop('Missing/nonfinite: ',n)}
out <- x['Quarter']
for(cc in countries) for(j in seq_along(vars)) out[[paste0(cc,'_',mapped[j])]] <- x[[paste0(cc,'_',vars[j])]]
dir.create('data/loglevels',recursive=TRUE,showWarnings=FALSE)
write.csv(out,'data/loglevels/model_input.csv',row.names=FALSE,na='')
write.csv(x,'data/loglevels/model_semantic.csv',row.names=FALSE,na='')
writeLines(c('CANDIDATE LOG LEVELS -- not yet a validated stationary VAR or cointegrated GVAR',
  'Legacy internal variable names only: y=log GDP, dp=log CPI LEVEL, r=annual rate (%), de=log REER LEVEL, deq=log equity PRICE LEVEL.',
  'Do not interpret dp as inflation, de as REER growth, or deq as a return. Do not cumulate their IRFs.',
  paste('Sample',head(x$Quarter,1),tail(x$Quarter,1),'n=',nrow(x)),
  'No overwrites of data/model_input.csv, no NA replacement.'),'data/loglevels/MAPPING.txt')
cat('PASS build: ',nrow(out),'quarters; 14 countries; 5 macro variables\n')
