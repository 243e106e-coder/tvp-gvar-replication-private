#!/usr/bin/env Rscript

# ============================================================
# Country-specific lag selection for GPR + Brent GVAR / TVP-GVAR
# v3: robust macro Excel adapter + persistent diagnostics
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

MACRO_PATH  <- "8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx"
MACRO_SHEET <- 1
WEIGHT_PATH <- "8.12/Trade_Weights_14_Economies_2000_2014.csv"
GPR_PATH    <- "8.12/gpr_quarterly_processed.csv"
OIL_PATH    <- "8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv"
OUT_DIR     <- "8.12/country_specific_lag_selection"

COUNTRIES <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
VARS <- c("y","dp","r","de","deq")
SAMPLE_START <- "2000Q2"
SAMPLE_END   <- "2025Q3"
P_CANDIDATES <- c(1L,2L)
Q_FIXED <- 1L
LB_LAG <- 4L
LB_ALPHA <- 0.05
STABILITY_CUTOFF <- 1.0
BORDERLINE_RHO <- 0.98
GPR_COL_EXACT <- "LN_GPR_QMEAN"
GDP_DLOG_DIVISOR <- 100

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...), call.=FALSE)
msg <- function(...) cat(sprintf(...), "\n")
norm_name <- function(x) tolower(gsub("[^a-z0-9]+", "", ifelse(is.na(x),"",as.character(x))))
trim_chr <- function(x) trimws(as.character(x))

quarter_label <- function(qid) {
  yr <- (qid-1L)%/%4L; qq <- qid-4L*yr
  sprintf("%dQ%d",yr,qq)
}

safe_date_parse <- function(x) {
  sx <- trimws(as.character(x)); out <- rep(as.Date(NA),length(sx))
  specs <- list(
    c("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$","%Y-%m-%d"),
    c("^[0-9]{4}/[0-9]{1,2}/[0-9]{1,2}$","%Y/%m/%d"),
    c("^[0-9]{4}\\.[0-9]{1,2}\\.[0-9]{1,2}$","%Y.%m.%d")
  )
  for (sp in specs) {
    ii <- which(is.na(out) & grepl(sp[1],sx))
    if(length(ii)) out[ii] <- suppressWarnings(as.Date(sx[ii],format=sp[2]))
  }
  ii <- which(is.na(out) & grepl("^[0-9]{4}[-/][0-9]{1,2}$",sx))
  if(length(ii)) out[ii] <- suppressWarnings(as.Date(paste0(gsub("/","-",sx[ii]),"-01"),format="%Y-%m-%d"))
  out
}

quarter_id <- function(x) {
  if(inherits(x,"Date") || inherits(x,"POSIXt")) {
    lt <- as.POSIXlt(x); return(4L*(lt$year+1900L)+(lt$mon%/%3L)+1L)
  }
  if(is.numeric(x)) {
    xx <- as.numeric(x); ff <- xx[is.finite(xx)]
    if(length(ff) && mean(ff>=20000 & ff<=70000)>=0.80) {
      d <- as.Date(xx,origin="1899-12-30"); lt <- as.POSIXlt(d)
      return(4L*(lt$year+1900L)+(lt$mon%/%3L)+1L)
    }
  }
  sx <- toupper(trimws(as.character(x))); out <- rep(NA_integer_,length(sx))
  m <- regexec("^([0-9]{4})[^0-9]*Q([1-4])$",sx); mm <- regmatches(sx,m)
  ok <- lengths(mm)==3L
  if(any(ok)) {
    yr <- as.integer(vapply(mm[ok],`[`,character(1),2)); qq <- as.integer(vapply(mm[ok],`[`,character(1),3))
    out[ok] <- 4L*yr+qq
  }
  need <- which(is.na(out))
  if(length(need)) {
    d <- safe_date_parse(sx[need]); ok2 <- !is.na(d)
    if(any(ok2)) { lt <- as.POSIXlt(d[ok2]); out[need[ok2]] <- 4L*(lt$year+1900L)+(lt$mon%/%3L)+1L }
  }
  out
}

quarter_parse_rate <- function(z) {
  q <- tryCatch(quarter_id(z),error=function(e) rep(NA_integer_,length(z)))
  if(!length(q)) 0 else mean(!is.na(q))
}

detect_date_col <- function(df,label="table") {
  nn <- norm_name(names(df)); pref <- norm_name(c("quarter","date","time","period","qtr","季度","日期","时间"))
  hit <- which(nn %in% pref)
  if(length(hit)) for(i in hit) if(quarter_parse_rate(df[[i]])>=0.50) return(names(df)[i])
  rr <- vapply(df,quarter_parse_rate,numeric(1)); rr[!is.finite(rr)] <- 0
  if(length(rr) && max(rr)>=0.80) return(names(df)[which.max(rr)])
  stopf("Could not detect date/quarter column in %s",label)
}

find_col_exact_norm <- function(df,candidates) {
  idx <- match(norm_name(candidates),norm_name(names(df)),nomatch=0L); idx <- idx[idx>0L]
  if(length(idx)) names(df)[idx[1]] else NA_character_
}

COUNTRY_ALIASES <- list(
  AU=c("AU","AUS","AUSTRALIA","澳大利亚"), BR=c("BR","BRA","BRAZIL","巴西"),
  CA=c("CA","CAN","CANADA","加拿大"), CH=c("CH","CHE","SWITZERLAND","SWISS","瑞士"),
  CN=c("CN","CHN","CHINA","中国"), EA=c("EA","EUROAREA","EUROZONE","EMU","欧元区"),
  UK=c("UK","GB","GBR","UNITEDKINGDOM","BRITAIN","英国"), JP=c("JP","JPN","JAPAN","日本"),
  KR=c("KR","KOR","KOREA","SOUTHKOREA","韩国"), NO=c("NO","NOR","NORWAY","挪威"),
  SG=c("SG","SGP","SINGAPORE","新加坡"), TR=c("TR","TUR","TURKEY","TURKIYE","土耳其"),
  US=c("US","USA","UNITEDSTATES","UNITEDSTATESOFAMERICA","美国"), ZA=c("ZA","ZAF","SOUTHAFRICA","南非")
)

VAR_ALIASES <- list(
  y=c("y","gdploglevel","gdplog","reallogdp","gdp","lngdp","loggdp","realgdp"),
  y_dlog=c("gdpdlog","ydlog","realgdpdlog","dloggdp","gdpgrowth","gdpchange"),
  dp=c("dp","cpidlog","inflation","cpiinflation","dlogcpi","cpichange","cpi"),
  r=c("r","rate","interestrate","shorttermrate","shortrate","policy rate","money market rate","interest"),
  de=c("de","reerdlog","exchangeratedlog","fxdlog","erdlog","dlogreer","reerchange","reer"),
  deq=c("deq","eqdlog","equitydlog","stockdlog","stockret","equityret","dlogequity","equitychange","stockreturn","equity","stock")
)

fill_right <- function(x) {
  x <- trim_chr(x); x[x=="" | is.na(x)] <- NA_character_; last <- NA_character_
  for(i in seq_along(x)) { if(!is.na(x[i])) last <- x[i] else if(!is.na(last)) x[i] <- last }
  x
}

match_country_token <- function(x) {
  z <- norm_name(x)
  for(cc in COUNTRIES) if(any(z %in% norm_name(COUNTRY_ALIASES[[cc]]))) return(cc)
  NA_character_
}

match_var_token <- function(x) {
  z <- norm_name(x)
  if(!nzchar(z)) return(NA_character_)
  ord <- c("y_dlog","dp","deq","de","r","y")
  for(v in ord) {
    aa <- norm_name(VAR_ALIASES[[v]])
    if(z %in% aa) return(v)
    if(any(vapply(aa,function(a) nzchar(a) && grepl(a,z,fixed=TRUE),logical(1)))) return(v)
  }
  NA_character_
}

save_header_preview <- function(raw,n=8L) {
  p <- raw[seq_len(min(n,nrow(raw))),,drop=FALSE]
  names(p) <- paste0("col_",seq_len(ncol(p)))
  write.csv(p,file.path(OUT_DIR,"00_macro_header_preview.csv"),row.names=FALSE,na="")
}

build_macro_mapping <- function(raw,max_header=6L) {
  H <- min(max_header,nrow(raw)-1L)
  n <- ncol(raw)
  best <- NULL; best_score <- -Inf

  for(h in seq_len(H)) {
    headmat <- raw[seq_len(h),,drop=FALSE]
    dat <- raw[(h+1L):nrow(raw),,drop=FALSE]
    date_rates <- vapply(dat,quarter_parse_rate,numeric(1)); date_rates[!is.finite(date_rates)] <- 0
    date_col <- if(length(date_rates)) which.max(date_rates) else NA_integer_
    date_rate <- if(!is.na(date_col)) date_rates[date_col] else 0

    cands <- matrix(NA_character_,h,n)
    for(r in seq_len(h)) {
      fr <- fill_right(unlist(headmat[r,,drop=TRUE],use.names=FALSE))
      for(j in seq_len(n)) cands[r,j] <- match_country_token(fr[j])
    }
    country <- rep(NA_character_,n)
    for(j in seq_len(n)) {
      z <- na.omit(cands[,j]); if(length(z)) country[j] <- z[1]
    }

    variable <- rep(NA_character_,n)
    for(j in seq_len(n)) {
      cells <- trim_chr(unlist(headmat[,j,drop=FALSE],use.names=FALSE))
      hits <- na.omit(vapply(cells,match_var_token,character(1)))
      if(length(hits)) variable[j] <- hits[1]
    }

    tab <- table(country[!is.na(country)])
    five_blocks <- length(tab)>=10L && all(tab==5L)
    if(five_blocks) {
      posmap <- rep(NA_character_,5)
      for(cc in names(tab)) {
        jj <- which(country==cc)
        if(length(jj)==5L) for(k in 1:5) if(is.na(posmap[k]) && !is.na(variable[jj[k]])) posmap[k] <- variable[jj[k]]
      }
      if(sum(!is.na(posmap))==5L && length(unique(posmap))==5L) {
        for(cc in names(tab)) {
          jj <- which(country==cc)
          for(k in 1:5) if(is.na(variable[jj[k]])) variable[jj[k]] <- posmap[k]
        }
      }
    }

    audit <- data.frame(column=seq_len(n),country=country,variable=variable,stringsAsFactors=FALSE)
    required_ok <- 0L
    for(cc in COUNTRIES) {
      vv <- variable[country==cc]
      required_ok <- required_ok + as.integer(any(vv %in% c("y","y_dlog"))) +
        sum(vapply(c("dp","r","de","deq"),function(v) any(vv==v),logical(1)))
    }
    score <- required_ok + 10*date_rate
    if(score>best_score) best <- list(h=h,date_col=date_col,date_rate=date_rate,audit=audit,score=score,required_ok=required_ok)
    best_score <- max(best_score,score)
  }
  best
}

read_macro_table <- function(path,sheet=1) {
  if(!file.exists(path)) stopf("Macro file not found: %s",path)
  if(!requireNamespace("readxl",quietly=TRUE)) stopf("Package readxl is required")
  raw <- as.data.frame(readxl::read_excel(path,sheet=sheet,col_names=FALSE,.name_repair="minimal",guess_max=1000),check.names=FALSE)
  save_header_preview(raw,8L)
  best <- build_macro_mapping(raw,6L)
  if(is.null(best)) stopf("Could not inspect macro header")

  audit <- best$audit
  audit$header_rows <- best$h
  audit$date_candidate <- audit$column==best$date_col
  write.csv(audit,file.path(OUT_DIR,"00_macro_mapping_audit.csv"),row.names=FALSE,na="")

  msg("[macro] best header depth=%d; mapped required pairs=%d/70; date parse=%.1f%%",best$h,best$required_ok,100*best$date_rate)
  print(audit[!is.na(audit$country) | !is.na(audit$variable) | audit$date_candidate,],row.names=FALSE)

  missing <- character()
  for(cc in COUNTRIES) {
    vv <- audit$variable[audit$country==cc]
    if(!any(vv %in% c("y","y_dlog"))) missing <- c(missing,paste0(cc,":y/y_dlog"))
    for(v in c("dp","r","de","deq")) if(!any(vv==v)) missing <- c(missing,paste0(cc,":",v))
  }
  if(best$date_rate<0.80) missing <- c(missing,"DATE_COLUMN")
  if(length(missing)) {
    writeLines(c("Macro mapping incomplete.",paste("Missing:",paste(missing,collapse=", "))),file.path(OUT_DIR,"00_macro_mapping_failure.txt"))
    stopf("Macro header mapping incomplete (%d/70). Missing: %s. Diagnostic files were saved in %s",best$required_ok,paste(missing,collapse=", "),OUT_DIR)
  }

  dat <- raw[(best$h+1L):nrow(raw),,drop=FALSE]
  qid <- quarter_id(dat[[best$date_col]])
  list(data=dat,qid=qid,audit=audit)
}

extract_macro_long <- function(obj) {
  raw <- obj$data; qid <- obj$qid; audit <- obj$audit; out <- list()
  for(cc in COUNTRIES) {
    aa <- audit[audit$country==cc,]
    pick <- function(v) { j <- aa$column[aa$variable==v]; if(length(j)) j[1] else NA_integer_ }
    jy <- pick("y"); jyd <- pick("y_dlog"); jdp <- pick("dp"); jr <- pick("r"); jde <- pick("de"); jeq <- pick("deq")
    if(is.na(jy) && is.na(jyd)) stopf("No GDP variable for %s",cc)
    if(any(is.na(c(jdp,jr,jde,jeq)))) stopf("Incomplete macro mapping for %s",cc)

    if(!is.na(jy)) {
      y <- suppressWarnings(as.numeric(raw[[jy]]))
    } else {
      gd <- suppressWarnings(as.numeric(raw[[jyd]]))/GDP_DLOG_DIVISOR
      y <- rep(NA_real_,length(gd))
      valid <- is.finite(gd) & !is.na(qid)
      if(any(valid)) {
        idx <- which(valid)
        grp <- cumsum(c(1L, diff(idx) > 1L))
        runs <- split(idx, grp)
        for(ii in runs) y[ii] <- log(100)+cumsum(gd[ii])
      }
    }

    out[[cc]] <- data.frame(
      qid=qid,country=cc,y=y,
      dp=suppressWarnings(as.numeric(raw[[jdp]])),
      r=suppressWarnings(as.numeric(raw[[jr]])),
      de=suppressWarnings(as.numeric(raw[[jde]])),
      deq=suppressWarnings(as.numeric(raw[[jeq]]))
    )
    msg("[macro] %s mapped columns: y=%s dp=%d r=%d de=%d deq=%d",cc,ifelse(!is.na(jy),jy,jyd),jdp,jr,jde,jeq)
  }
  ans <- do.call(rbind,out)
  ans[!is.na(ans$qid),]
}

read_weights <- function(path) {
  if(!file.exists(path)) stopf("Weights file not found: %s",path)
  w0 <- read.csv(path,check.names=FALSE)
  first <- toupper(trimws(as.character(w0[[1]])))
  if(all(COUNTRIES %in% first)) { rownames(w0)<-first; w0<-w0[,-1,drop=FALSE] }
  cn <- toupper(trimws(names(w0)))
  if(!all(COUNTRIES %in% cn)) stopf("Weight columns do not contain all 14 economies")
  w0 <- w0[,match(COUNTRIES,cn),drop=FALSE]
  if(is.null(rownames(w0)) || !all(COUNTRIES %in% toupper(rownames(w0)))) {
    if(nrow(w0)!=length(COUNTRIES)) stopf("Weights must have 14 rows")
    rownames(w0)<-COUNTRIES
  } else w0 <- w0[match(COUNTRIES,toupper(rownames(w0))),,drop=FALSE]

  W <- as.matrix(data.frame(lapply(w0,as.numeric),check.names=FALSE))
  rownames(W)<-COUNTRIES; colnames(W)<-COUNTRIES
  diag(W)<-0
  rs<-rowSums(W,na.rm=TRUE)
  if(any(!is.finite(rs)|rs<=0)) stopf("Nonpositive weight row sum")
  W/rs
}

read_global_series <- function(path,value_candidates,exact=NULL,label="global") {
  if(!file.exists(path)) stopf("%s file not found: %s",label,path)
  d <- read.csv(path,check.names=FALSE)
  dc <- detect_date_col(d,label)
  qid <- quarter_id(d[[dc]])
  vc <- if(!is.null(exact) && exact %in% names(d)) exact else find_col_exact_norm(d,value_candidates)
  if(is.na(vc)) {
    cc <- setdiff(names(d),dc)
    rr <- vapply(cc,function(nm) mean(is.finite(suppressWarnings(as.numeric(d[[nm]])))),numeric(1))
    if(length(rr)&&max(rr)>=0.8) vc <- cc[which.max(rr)]
  }
  if(is.na(vc)) stopf("Could not identify value column in %s",path)
  z <- data.frame(qid=qid,value=suppressWarnings(as.numeric(d[[vc]])))
  z <- z[!is.na(z$qid)&is.finite(z$value),]
  aggregate(value~qid,data=z,FUN=mean)
}

lag_matrix <- function(X,L) {
  z <- matrix(NA_real_,nrow(X),ncol(X))
  if(L<nrow(X)) z[(L+1L):nrow(X),] <- X[1L:(nrow(X)-L),,drop=FALSE]
  z
}

companion_radius <- function(Bdom,k,p) {
  if(p==1L) return(max(Mod(eigen(Bdom[[1]],only.values=TRUE)$values)))
  C <- rbind(do.call(cbind,Bdom),cbind(diag(k*(p-1L)),matrix(0,k*(p-1L),k)))
  max(Mod(eigen(C,only.values=TRUE)$values))
}

fit_country <- function(cc,p,panel,foreign,globals) {
  Xi <- as.matrix(panel[[cc]][,VARS])
  Xs <- as.matrix(foreign[[cc]][,VARS])
  qid <- panel[[cc]]$qid
  n <- nrow(Xi)

  idx <- seq.int(max(P_CANDIDATES)+1L,n)
  Y <- Xi[idx,,drop=FALSE]
  Z <- matrix(1,length(idx),1)
  dom <- vector("list",p)

  for(L in seq_len(p)) {
    XL <- lag_matrix(Xi,L)[idx,,drop=FALSE]
    st <- ncol(Z)+1L
    Z <- cbind(Z,XL)
    dom[[L]] <- st:ncol(Z)
  }

  Z <- cbind(Z,Xs[idx,,drop=FALSE],lag_matrix(Xs,1L)[idx,,drop=FALSE])
  G <- as.matrix(globals[match(qid,globals$qid),c("gpr","oil")])
  Z <- cbind(Z,G[idx,,drop=FALSE],lag_matrix(G,1L)[idx,,drop=FALSE])

  keep <- complete.cases(Y)&complete.cases(Z)
  Y<-Y[keep,,drop=FALSE]; Z<-Z[keep,,drop=FALSE]
  if(nrow(Y)<ncol(Z)+10L) stopf("%s p=%d too few complete observations",cc,p)

  fit <- lm.fit(Z,Y)
  B<-fit$coefficients; E<-fit$residuals
  Tn<-nrow(E); K<-ncol(Y); m<-ncol(Z)
  S<-crossprod(E)/Tn
  ds<-determinant(S,logarithm=TRUE)
  if(ds$sign<=0) stopf("%s p=%d residual covariance singular",cc,p)
  logdet<-as.numeric(ds$modulus)
  npar<-K*m
  aic<-Tn*logdet+2*npar
  bic<-Tn*logdet+log(Tn)*npar

  mats <- lapply(seq_len(p),function(L) t(B[dom[[L]],,drop=FALSE]))
  rho<-companion_radius(mats,K,p)

  lb <- do.call(rbind,lapply(seq_len(K),function(j) {
    Lb<-min(LB_LAG,max(1L,floor(Tn/5)))
    bt<-Box.test(E[,j],lag=Lb,type="Ljung-Box",fitdf=0)
    data.frame(country=cc,p=p,equation=VARS[j],lb_lag=Lb,
               statistic=unname(bt$statistic),p_value=bt$p.value,
               serial_corr_reject_5pct=bt$p.value<LB_ALPHA)
  }))

  dg <- data.frame(
    country=cc,p=p,nobs=Tn,regressors=m,AIC=aic,BIC=bic,
    spectral_radius=rho,
    stable=is.finite(rho)&&rho<STABILITY_CUTOFF,
    borderline=is.finite(rho)&&rho>=BORDERLINE_RHO&&rho<STABILITY_CUTOFF,
    lb_rejections_5pct=sum(lb$serial_corr_reject_5pct),
    lb_all_pass_5pct=all(!lb$serial_corr_reject_5pct)
  )
  list(diag=dg,lb=lb)
}

msg("Reading macro panel...")
macro_obj <- read_macro_table(MACRO_PATH,MACRO_SHEET)
macro_long <- extract_macro_long(macro_obj)

q_start<-quarter_id(SAMPLE_START)
q_end<-quarter_id(SAMPLE_END)
macro_long<-macro_long[macro_long$qid>=q_start & macro_long$qid<=q_end,]

W<-read_weights(WEIGHT_PATH)

gpr<-read_global_series(
  GPR_PATH,
  c("LN_GPR_QMEAN","GPR","LN_GPR","GPR_QMEAN"),
  GPR_COL_EXACT,
  "GPR"
)
names(gpr)[2]<-"gpr"

oil<-read_global_series(
  OIL_PATH,
  c("BRENT_LOG","LOG_BRENT","LN_BRENT","BRENT","OIL","LOG_PRICE","VALUE"),
  label="OIL"
)
names(oil)[2]<-"oil"

globals<-merge(gpr,oil,by="qid")
globals<-globals[globals$qid>=q_start & globals$qid<=q_end,]

country_q<-lapply(COUNTRIES,function(cc){
  z<-macro_long[macro_long$country==cc,]
  z$qid[complete.cases(z[,VARS,drop=FALSE])]
})
common_q<-sort(Reduce(intersect,c(country_q,list(globals$qid))))
if(length(common_q)<30L) stopf("Common complete sample too short: %d quarters",length(common_q))
msg("Common complete sample: %s to %s (%d quarters)",
    quarter_label(min(common_q)),quarter_label(max(common_q)),length(common_q))

panel<-setNames(vector("list",length(COUNTRIES)),COUNTRIES)
for(cc in COUNTRIES) {
  z<-macro_long[macro_long$country==cc & macro_long$qid %in% common_q,c("qid",VARS)]
  z<-z[match(common_q,z$qid),]
  if(any(is.na(z$qid))||any(!complete.cases(z[,VARS,drop=FALSE])))
    stopf("%s incomplete after alignment",cc)
  panel[[cc]]<-z
}
globals<-globals[match(common_q,globals$qid),]

foreign<-setNames(vector("list",length(COUNTRIES)),COUNTRIES)
for(i in seq_along(COUNTRIES)) {
  F<-matrix(0,length(common_q),length(VARS))
  colnames(F)<-VARS
  for(j in seq_along(COUNTRIES))
    F<-F+W[i,j]*as.matrix(panel[[COUNTRIES[j]]][,VARS,drop=FALSE])
  foreign[[COUNTRIES[i]]]<-data.frame(qid=common_q,F,check.names=FALSE)
}

all_diag<-list()
all_lb<-list()
kk<-1L
for(cc in COUNTRIES) for(p in P_CANDIDATES) {
  msg("Fitting %s p=%d q=%d",cc,p,Q_FIXED)
  ans<-fit_country(cc,p,panel,foreign,globals)
  all_diag[[kk]]<-ans$diag
  all_lb[[kk]]<-ans$lb
  kk<-kk+1L
}

diag_df<-do.call(rbind,all_diag)
lb_df<-do.call(rbind,all_lb)

recommend_one<-function(d) {
  d1<-d[d$p==1L,]
  d2<-d[d$p==2L,]
  cc<-d$country[1]

  if(d1$stable&&!d2$stable)
    return(data.frame(country=cc,selected_p=1L,reason="p=2 unstable; p=1 stable"))
  if(!d1$stable&&d2$stable)
    return(data.frame(country=cc,selected_p=2L,reason="p=1 unstable; p=2 stable"))
  if(!d1$stable&&!d2$stable)
    return(data.frame(
      country=cc,
      selected_p=ifelse(d1$spectral_radius<=d2$spectral_radius,1L,2L),
      reason="WARNING: both locally unstable; lower spectral radius selected"
    ))

  if(!d1$lb_all_pass_5pct&&d2$lb_all_pass_5pct)
    return(data.frame(country=cc,selected_p=2L,
                      reason="p=2 removes Ljung-Box rejections present under p=1"))

  if(d1$BIC<=d2$BIC) {
    data.frame(
      country=cc,selected_p=1L,
      reason=ifelse(d1$borderline,
                    "BIC prefers p=1 but stability is borderline",
                    "both acceptable; BIC prefers parsimonious p=1")
    )
  } else {
    data.frame(
      country=cc,selected_p=2L,
      reason=ifelse(d2$borderline,
                    "BIC prefers p=2 but stability is borderline",
                    "both acceptable; BIC prefers p=2")
    )
  }
}

rec<-do.call(rbind,lapply(split(diag_df,diag_df$country),recommend_one))
rec<-rec[match(COUNTRIES,rec$country),]

key<-function(p) match(paste(rec$country,p),paste(diag_df$country,diag_df$p))
rec$AIC_p1<-diag_df$AIC[key(1)]
rec$AIC_p2<-diag_df$AIC[key(2)]
rec$BIC_p1<-diag_df$BIC[key(1)]
rec$BIC_p2<-diag_df$BIC[key(2)]
rec$rho_p1<-diag_df$spectral_radius[key(1)]
rec$rho_p2<-diag_df$spectral_radius[key(2)]
rec$LB_reject_p1<-diag_df$lb_rejections_5pct[key(1)]
rec$LB_reject_p2<-diag_df$lb_rejections_5pct[key(2)]

write.csv(diag_df,file.path(OUT_DIR,"01_country_p1_p2_model_diagnostics.csv"),row.names=FALSE)
write.csv(lb_df,file.path(OUT_DIR,"02_equation_residual_serial_correlation.csv"),row.names=FALSE)
write.csv(rec,file.path(OUT_DIR,"03_country_lag_recommendation.csv"),row.names=FALSE)

lag_vec<-setNames(rec$selected_p,rec$country)
lag_line<-paste(sprintf("%s=%dL",names(lag_vec),lag_vec),collapse=", ")
writeLines(
  c("# Auto-generated country-specific domestic lag vector",
    sprintf("COUNTRY_P <- c(%s)",lag_line)),
  file.path(OUT_DIR,"04_selected_country_lag_vector.R")
)

writeLines(
  c(
    "Country-specific lag selection: GPR + Brent GVAR / TVP-GVAR",
    "",
    "v3 macro adapter: multi-row grouped Excel header is mapped column-by-column; incomplete mappings stop safely but diagnostic files remain available.",
    sprintf("Requested sample: %s to %s",SAMPLE_START,SAMPLE_END),
    sprintf("Actual common sample: %s to %s (%d quarters)",
            quarter_label(min(common_q)),quarter_label(max(common_q)),length(common_q)),
    "Candidate domestic lags: p_i=1,2; foreign lag q_i=1",
    "Selection: stability first, residual serial correlation second, BIC third.",
    "Final TVP-GVAR still requires posterior/global stability checks."
  ),
  file.path(OUT_DIR,"README_country_specific_lag_selection.txt")
)

msg("DONE. Selected lag vector:")
print(lag_vec)
msg("Outputs: %s",OUT_DIR)
