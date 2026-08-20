#!/usr/bin/env Rscript
##
## generate_fixtures.R -- regenerate the R reference fixtures for FXRegime.jl
## ===========================================================================
##
## Every number the Julia test suite checks against R is produced by this
## script and written to test/fixtures/ as plain CSV.  Committing the CSVs
## means `] test FXRegime` runs without R installed; this script documents
## exactly how each number was produced and lets anyone regenerate them.
##
## Usage:
##     Rscript test/generate_fixtures.R [outdir]
##
## Requires: fxregime, strucchange, zoo, sandwich (car is optional -- if it is
## missing the fxpegtest fixture is skipped and the rest still runs).
##
## Runtime is dominated by the two fxregimes() calls (~5-15 min total); set
##     FXREGIME_FIXTURE_CACHE=/some/path.rds
## to cache/reuse the fitted fxregimes objects across runs.  Unset (the
## default) everything is recomputed from scratch.
##
## Numeric columns are written with sprintf("%.17g") so that every value
## round-trips exactly through the CSV.  Missing values are written as "NA".
## Dates are written as ISO yyyy-mm-dd.
##
## ---------------------------------------------------------------------------

suppressMessages({
  library("fxregime")
  library("strucchange")
  library("zoo")
  library("sandwich")
})

HAVE_CAR <- requireNamespace("car", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
FIXDIR <- if (length(args) >= 1) args[1] else {
  ## default: <dir of this script>/fixtures
  self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(self) == 1) file.path(dirname(normalizePath(self)), "fixtures") else "fixtures"
}
dir.create(FIXDIR, showWarnings = FALSE, recursive = TRUE)

CACHE <- Sys.getenv("FXREGIME_FIXTURE_CACHE", "")

pkgs <- c("fxregime", "strucchange", "zoo", "sandwich")
if (HAVE_CAR) pkgs <- c(pkgs, "car")
PKGVER <- paste(c(paste(pkgs, sapply(pkgs, function(p) as.character(packageVersion(p)))),
                  R.version.string), collapse = "; ")

cat("fixtures -> ", FIXDIR, "\n", sep = "")
cat("versions:  ", PKGVER, "\n\n", sep = "")

## --------------------------------------------------------------------------
## CSV writing helpers
## --------------------------------------------------------------------------

fmt_col <- function(x) {
  if (inherits(x, "Date")) {
    y <- format(x, "%Y-%m-%d"); y[is.na(x)] <- "NA"; y
  } else if (is.numeric(x)) {
    y <- sprintf("%.17g", as.numeric(x)); y[is.na(x)] <- "NA"; y
  } else {
    y <- as.character(x); y[is.na(y)] <- "NA"; y
  }
}

manifest <- list()

add_fixture <- function(df, file, desc, call) {
  stopifnot(nrow(df) > 0L, ncol(df) > 0L)
  out <- as.data.frame(lapply(df, fmt_col), stringsAsFactors = FALSE)
  names(out) <- names(df)
  ## no field may contain a comma, otherwise the unquoted CSV is ambiguous
  if (any(grepl(",", unlist(out), fixed = TRUE)) || any(grepl(",", names(out), fixed = TRUE)))
    stop("comma in field or header of ", file)
  write.table(out, file.path(FIXDIR, file), sep = ",",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  manifest[[length(manifest) + 1L]] <<- data.frame(
    file = file, nrow = nrow(df), ncol = ncol(df),
    description = desc, r_call = call, r_packages = PKGVER,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-46s %6d x %2d\n", file, nrow(df), ncol(df)))
  invisible(NULL)
}

## zoo/matrix -> data.frame with a leading `date` (or `obs`) column
zoo_df <- function(z, datecol = "date") {
  d <- as.data.frame(coredata(as.matrix(z)), stringsAsFactors = FALSE)
  names(d) <- colnames(as.matrix(z))
  cbind(setNames(data.frame(index(z), stringsAsFactors = FALSE), datecol), d)
}

mat_df <- function(m, rowname = "row") {
  d <- as.data.frame(unclass(m), stringsAsFactors = FALSE)
  names(d) <- colnames(m)
  cbind(setNames(data.frame(rownames(m), stringsAsFactors = FALSE), rowname), d)
}

## --------------------------------------------------------------------------
## 1. Input data: fxreturns()
## --------------------------------------------------------------------------

data("FXRatesCHF", package = "fxregime")

cat("[1] fxreturns\n")

## CNY vignette: daily returns, 2005-07-25 .. 2009-07-31
cny <- fxreturns("CNY", frequency = "daily",
  start = as.Date("2005-07-25"), end = as.Date("2009-07-31"),
  other = c("USD", "JPY", "EUR", "GBP"), data = FXRatesCHF)

## INR vignette: weekly (Friday-anchored) returns, 1993-04-01 .. 2008-01-04
inr <- fxreturns("INR", frequency = "weekly",
  start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
  other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF)

add_fixture(zoo_df(cny), "fxreturns_cny_daily.csv",
  "Daily 100*diff(log(rate)) returns for CNY and the four basket currencies (CNY vignette input); index + all columns.",
  'fxreturns("CNY", frequency = "daily", start = as.Date("2005-07-25"), end = as.Date("2009-07-31"), other = c("USD","JPY","EUR","GBP"), data = FXRatesCHF)')

add_fixture(zoo_df(inr), "fxreturns_inr_weekly.csv",
  "Weekly (next-Friday anchored) 100*diff(log(rate)) returns for INR and the four basket currencies (INR vignette input); index + all columns.",
  'fxreturns("INR", frequency = "weekly", start = as.Date("1993-04-01"), end = as.Date("2008-01-04"), other = c("USD","JPY","DUR","GBP"), data = FXRatesCHF)')

## a couple of extra fxreturns variants to pin down the na.action / trim paths
cny_trim <- fxreturns("CNY", frequency = "daily",
  start = as.Date("2005-07-25"), end = as.Date("2009-07-31"),
  other = c("USD", "JPY", "EUR", "GBP"), data = FXRatesCHF, trim = TRUE)
add_fixture(zoo_df(cny_trim), "fxreturns_cny_daily_trim.csv",
  "As fxreturns_cny_daily.csv but with trim = TRUE, i.e. rows whose CNY return is outside the 1%/99% sample quantiles (type 7) are dropped.",
  'fxreturns("CNY", frequency = "daily", start = as.Date("2005-07-25"), end = as.Date("2009-07-31"), other = c("USD","JPY","EUR","GBP"), data = FXRatesCHF, trim = TRUE)')

inr_omit <- fxreturns("INR", frequency = "weekly",
  start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
  other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF, na.action = na.omit)
add_fixture(zoo_df(inr_omit), "fxreturns_inr_weekly_naomit.csv",
  "As fxreturns_inr_weekly.csv but with na.action = na.omit instead of the default na.locf.",
  'fxreturns("INR", frequency = "weekly", start = as.Date("1993-04-01"), end = as.Date("2008-01-04"), other = c("USD","JPY","DUR","GBP"), data = FXRatesCHF, na.action = na.omit)')

## --------------------------------------------------------------------------
## 2. fxlm(): coefficients, residuals, fitted, estfun, bread, summary table
## --------------------------------------------------------------------------

cat("\n[2] fxlm\n")

emit_fxlm <- function(m, tag, desc, call) {
  cf <- coef(m)                       # coef.fxlm: regressors + "(Variance)" (MLE)
  s  <- summary(m)
  sc <- s$coefficients
  ef <- estfun(m)                     # estfun.fxlm
  br <- bread(m)                      # bread.fxlm
  n  <- length(residuals(m))
  k  <- length(cf) - 1L

  add_fixture(data.frame(term = names(cf), estimate = as.numeric(cf), stringsAsFactors = FALSE),
    sprintf("fxlm_%s_coef.csv", tag),
    paste0(desc, ": coef(), i.e. regression coefficients followed by the MLE variance `(Variance)` = mean(residuals^2)."),
    paste0("coef(", call, ")"))

  add_fixture(data.frame(term = rownames(sc), estimate = sc[, 1], std_error = sc[, 2],
                         t_value = sc[, 3], p_value = sc[, 4], stringsAsFactors = FALSE),
    sprintf("fxlm_%s_summary.csv", tag),
    paste0(desc, ": summary() coefficient table (estimate, std. error, t value, Pr(>|t|)); the OLS/df-adjusted covariance, t distribution with n-k df."),
    paste0("summary(", call, ")$coefficients"))

  add_fixture(data.frame(
      statistic = c("nobs", "nreg", "df_residual", "sigma", "sigma2_mle",
                    "r_squared", "adj_r_squared", "logLik", "rss"),
      value = c(n, k, df.residual(m), s$sigma, as.numeric(cf["(Variance)"]),
                s$r.squared, s$adj.r.squared, as.numeric(logLik(m)),
                sum(residuals(m)^2)),
      stringsAsFactors = FALSE),
    sprintf("fxlm_%s_stats.csv", tag),
    paste0(desc, ": scalar model summaries (n, #regressors incl. intercept, residual df, df-adjusted sigma, MLE variance, centred R^2, adjusted R^2, Gaussian logLik, residual sum of squares)."),
    paste0("summary(", call, ")"))

  add_fixture(data.frame(date = index(m), residual = as.numeric(residuals(m)),
                         fitted = as.numeric(fitted(m)), stringsAsFactors = FALSE),
    sprintf("fxlm_%s_resid_fitted.csv", tag),
    paste0(desc, ": per-observation residuals and fitted values."),
    paste0("residuals(", call, "); fitted(", call, ")"))

  add_fixture(zoo_df(ef), sprintf("fxlm_%s_estfun.csv", tag),
    paste0(desc, ": estfun.fxlm() -- cbind(residual * X / sigma2, (residual^2 - sigma2) / (2*sigma2^2)) with sigma2 the MLE variance."),
    paste0("estfun(", call, ")"))

  add_fixture(mat_df(br), sprintf("fxlm_%s_bread.csv", tag),
    paste0(desc, ": bread.fxlm() -- the (k+1)x(k+1) inverse expected-Hessian including the beta/variance cross terms."),
    paste0("bread(", call, ")"))

  add_fixture(mat_df(vcov(m)), sprintf("fxlm_%s_vcov.csv", tag),
    paste0(desc, ": OLS (df-adjusted) covariance matrix of the regression coefficients."),
    paste0("vcov(", call, ")"))

  ## vcovHC() must see a plain "lm" -- estfun.fxlm has the extra (Variance) column
  m_lm <- m; class(m_lm) <- "lm"
  add_fixture(mat_df(vcovHC(m_lm, type = "HC3")), sprintf("fxlm_%s_vcovHC3.csv", tag),
    paste0(desc, ": HC3 sandwich covariance matrix of the regression coefficients."),
    paste0("sandwich::vcovHC(`class<-`(", call, ', "lm"), type = "HC3")'))
}

cny_hist_data <- window(cny, end = as.Date("2005-10-31"))
cny_lm <- fxlm(CNY ~ USD + JPY + EUR + GBP, data = cny_hist_data)
emit_fxlm(cny_lm, "cny_hist",
  "CNY history model (daily returns 2005-07-26 .. 2005-10-31, i.e. the pre-revaluation window of the CNY vignette)",
  'fxlm(CNY ~ USD + JPY + EUR + GBP, data = window(cny, end = as.Date("2005-10-31")))')

inr_lm <- fxlm(INR ~ USD + JPY + DUR + GBP, data = inr)
emit_fxlm(inr_lm, "inr_full",
  "INR full-sample model (weekly returns 1993 .. 2008, INR vignette)",
  "fxlm(INR ~ USD + JPY + DUR + GBP, data = inr)")

if (HAVE_CAR) {
  emit_pegtest <- function(m, tag, desc, call) {
    pt <- fxpegtest(m)
    ## fxpegtest()'s default peg, recomputed exactly as fxregime does it
    cc <- coef(m); cc <- cc[-c(1, length(cc))]
    peg <- names(cc)[which.max(abs(cc))]
    ## row 1 of car::linearHypothesis()'s anova is the RESTRICTED model
    ## (Res.Df = n - 1 here), row 2 the unrestricted one (Res.Df = n - k)
    add_fixture(data.frame(
        statistic = c("peg", "res_df_restricted", "res_df_unrestricted",
                      "rss_restricted", "rss_unrestricted", "df", "F", "p_value"),
        value = c(peg, sprintf("%.17g", c(pt$Res.Df[1], pt$Res.Df[2], pt$RSS[1],
                                          pt$RSS[2], pt$Df[2], pt$F[2], pt$`Pr(>F)`[2]))),
        stringsAsFactors = FALSE),
      sprintf("fxpegtest_%s.csv", tag),
      paste0(desc, ": Wald/F test of `peg coefficient = 1, all other slopes = 0`, with peg = the slope with the largest |coefficient|. `df` is the number of restrictions (numerator df), `res_df_unrestricted` the denominator df."),
      paste0("fxpegtest(", call, ")"))
  }
  emit_pegtest(inr_lm, "inr_full", "INR full-sample model",
    "fxlm(INR ~ USD + JPY + DUR + GBP, data = inr)")
  emit_pegtest(cny_lm, "cny_hist", "CNY history model",
    'fxlm(CNY ~ USD + JPY + EUR + GBP, data = window(cny, end = as.Date("2005-10-31")))')
} else {
  cat("  (car not installed: fxpegtest fixtures skipped)\n")
}

## --------------------------------------------------------------------------
## 3. recresid()
## --------------------------------------------------------------------------

cat("\n[3] recresid\n")

emit_recresid <- function(formula, data, tag, desc, call) {
  df <- as.data.frame(data)
  X <- model.matrix(formula, data = df)
  attr(X, "assign") <- NULL
  y <- model.response(model.frame(formula, data = df))
  rr <- recresid(X, y)
  k <- ncol(X)
  add_fixture(data.frame(obs = seq.int(k + 1L, nrow(X)), recresid = as.numeric(rr),
                         stringsAsFactors = FALSE),
    sprintf("recresid_%s.csv", tag),
    paste0(desc, ": standardised recursive residuals; `obs` is the 1-based observation index (k+1 ... n) the residual belongs to."),
    call)
}

emit_recresid(CNY ~ USD + JPY + EUR + GBP, cny, "cny_daily",
  "CNY daily returns, full sample",
  "recresid(model.matrix(CNY ~ USD + JPY + EUR + GBP, as.data.frame(cny)), as.data.frame(cny)$CNY)")
emit_recresid(INR ~ USD + JPY + DUR + GBP, inr, "inr_weekly",
  "INR weekly returns, full sample",
  "recresid(model.matrix(INR ~ USD + JPY + DUR + GBP, as.data.frame(inr)), as.data.frame(inr)$INR)")

## --------------------------------------------------------------------------
## 4. RSS triangle (full) on a small subset
## --------------------------------------------------------------------------

cat("\n[4] RSS triangle\n")

emit_triangle <- function(formula, data, h, tag, desc, call) {
  bp <- breakpoints(formula, data = as.data.frame(data), h = h)
  tri <- bp$RSS.triang
  n <- bp$nobs
  starts <- unlist(lapply(seq_along(tri), function(i) rep.int(i, length(tri[[i]]))))
  ends   <- unlist(lapply(seq_along(tri), function(i) seq.int(i, i + length(tri[[i]]) - 1L)))
  add_fixture(data.frame(start = starts, end = ends, nseg = ends - starts + 1L,
                         rss = unlist(tri), stringsAsFactors = FALSE),
    sprintf("rss_triangle_%s.csv", tag),
    paste0(desc, ": the complete RSS triangle. Row (start, end) is the residual sum of squares of the OLS fit on observations start..end, accumulated from the recursive residuals; NA where the segment has no more observations than regressors."),
    call)
  add_fixture(data.frame(statistic = c("nobs", "nreg", "h", "ntriang"),
                         value = c(n, bp$nreg, h, length(tri)), stringsAsFactors = FALSE),
    sprintf("rss_triangle_%s_meta.csv", tag),
    paste0(desc, ": shape of the RSS triangle (number of observations, number of regressors incl. intercept, minimal segment size h, number of triangle rows = n-h+1)."),
    call)
}

emit_triangle(CNY ~ USD + JPY + EUR + GBP, cny[1:120, ], 20, "cny120_h20",
  "First 120 observations of the CNY daily return series, h = 20",
  "breakpoints(CNY ~ USD + JPY + EUR + GBP, data = as.data.frame(cny[1:120,]), h = 20)$RSS.triang")

## --------------------------------------------------------------------------
## 5. fxregimes(): breakpoints, NLL/BIC/LWZ, selection, confint, refit
## --------------------------------------------------------------------------

cat("\n[5] fxregimes (this is the slow part)\n")

fit_regimes <- function() {
  list(cny = fxregimes(CNY ~ USD + JPY + EUR + GBP, data = cny, h = 20, breaks = 10),
       inr = fxregimes(INR ~ USD + JPY + DUR + GBP, data = inr, h = 20, breaks = 10))
}

if (nzchar(CACHE) && file.exists(CACHE)) {
  cat("  (loading cached fxregimes objects from ", CACHE, ")\n", sep = "")
  regs <- readRDS(CACHE)
} else {
  t0 <- Sys.time()
  regs <- fit_regimes()
  cat("  fitted in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n", sep = "")
  if (nzchar(CACHE)) saveRDS(regs, CACHE)
}
cny_reg <- regs$cny
inr_reg <- regs$inr

MAXM <- 10L

emit_regimes <- function(reg, tag, desc, call, maxm = MAXM, h = 20) {
  n <- reg$nobs
  npar <- reg$npar
  idx <- index(reg$data)

  ## breakpoints for every m = 1 ... maxm
  rows <- do.call(rbind, lapply(1:maxm, function(m) {
    b <- breakpoints(reg, breaks = m)$breakpoints
    data.frame(m = m, k = seq_along(b), breakpoint = b, breakdate = idx[b],
               stringsAsFactors = FALSE)
  }))
  add_fixture(rows, sprintf("fxregimes_%s_breakpoints.csv", tag),
    paste0(desc, ": optimal breakpoints for every number of breaks m = 1..", maxm,
           ". `breakpoint` is the 1-based index of the LAST observation of the segment ending there; `breakdate` is index(data)[breakpoint]."),
    paste0("sapply(1:", maxm, ", function(m) breakpoints(", call, ", breaks = m)$breakpoints)"))

  ## objective / information criteria for every m = 0 ... maxm
  ic <- do.call(rbind, lapply(0:maxm, function(m) {
    b <- breakpoints(reg, breaks = m)
    data.frame(m = m, df = npar * (m + 1L) + m, nll = b$objective,
               bic = AIC(b, k = log(n)), lwz = AIC(b, k = 0.299 * log(n)^2.1),
               stringsAsFactors = FALSE)
  }))
  add_fixture(ic, sprintf("fxregimes_%s_ic.csv", tag),
    paste0(desc, ": negative Gaussian log-likelihood and the two information criteria for every m = 0..", maxm,
           ". df = npar*(m+1)+m with npar = ncol(X)+1; bic = 2*nll + df*log(n); lwz = 2*nll + df*0.299*log(n)^2.1."),
    paste0("summary(", call, ")$objective"))

  sel_bic <- ic$m[which.min(ic$bic)]
  sel_lwz <- ic$m[which.min(ic$lwz)]
  add_fixture(data.frame(
      statistic = c("nobs", "npar", "h", "maxbreaks", "m_selected_BIC", "m_selected_LWZ",
                    "m_default", "objective_default"),
      value = c(n, npar, h, maxm, sel_bic, sel_lwz,
                length(reg$breakpoints[!is.na(reg$breakpoints)]),
                breakpoints(reg)$objective),
      stringsAsFactors = FALSE),
    sprintf("fxregimes_%s_selection.csv", tag),
    paste0(desc, ": model dimensions plus the number of breaks selected by BIC and by LWZ, the number selected by the object's own criterion (LWZ, the fxregimes default) and the corresponding negative log-likelihood."),
    paste0("breakpoints(", call, ")"))

  ## confidence intervals at two levels
  for (lev in c(0.9, 0.95)) {
    ci <- withCallingHandlers(confint(reg, level = lev),
                              warning = function(w) invokeRestart("muffleWarning"))
    cm <- ci$confint
    add_fixture(data.frame(
        k = seq_len(nrow(cm)), lower = cm[, 1], breakpoint = cm[, 2], upper = cm[, 3],
        lower_date = idx[cm[, 1]], breakdate = idx[cm[, 2]], upper_date = idx[cm[, 3]],
        stringsAsFactors = FALSE),
      sprintf("confint_%s_%d.csv", tag, round(lev * 100)),
      paste0(desc, ": confidence intervals for the breakpoints of the criterion-selected partition at level ", lev,
             ". Bounds are observation indices (NA when R refuses the interval because P(argmax V <= 0) is outside [(1-level)/2, 1-(1-level)/2])."),
      paste0("confint(", call, ", level = ", lev, ")$confint"))
  }

  ## refit: one fxlm per segment
  rf <- refit(reg)
  cf <- coef(reg)
  seg <- strsplit(rownames(cf), "--", fixed = TRUE)
  add_fixture(cbind(
      data.frame(segment = seq_along(rf),
                 start_date = as.Date(sapply(seg, `[`, 1)),
                 end_date = as.Date(sapply(seg, `[`, 2)),
                 nobs = sapply(rf, function(z) length(residuals(z))),
                 stringsAsFactors = FALSE),
      setNames(as.data.frame(unclass(cf)), colnames(cf))),
    sprintf("fxregimes_%s_segment_coef.csv", tag),
    paste0(desc, ": per-segment refitted fxlm coefficients (regression coefficients plus the MLE `(Variance)`) for the criterion-selected partition, with each segment's date range and number of observations."),
    paste0("coef(", call, "); refit(", call, ")"))

  add_fixture(do.call(rbind, lapply(seq_along(rf), function(i) {
      sc <- summary(rf[[i]])$coefficients
      data.frame(segment = i, term = rownames(sc), estimate = sc[, 1], std_error = sc[, 2],
                 t_value = sc[, 3], p_value = sc[, 4], stringsAsFactors = FALSE)
    })),
    sprintf("fxregimes_%s_segment_summary.csv", tag),
    paste0(desc, ": per-segment summary() coefficient tables (estimate, std. error, t value, Pr(>|t|)) of the refitted models."),
    paste0("lapply(refit(", call, "), function(z) summary(z)$coefficients)"))

  add_fixture(data.frame(date = idx,
                         fitted = as.numeric(fitted(reg)),
                         residual = as.numeric(residuals(reg)),
                         stringsAsFactors = FALSE),
    sprintf("fxregimes_%s_fitted.csv", tag),
    paste0(desc, ": fitted values and residuals of the segmented model, concatenated across the segments of the criterion-selected partition."),
    paste0("fitted(", call, "); residuals(", call, ")"))
}

emit_regimes(cny_reg, "cny", "CNY daily returns, h = 20, breaks = 10",
  "fxregimes(CNY ~ USD + JPY + EUR + GBP, data = cny, h = 20, breaks = 10)")
emit_regimes(inr_reg, "inr", "INR weekly returns, h = 20, breaks = 10",
  "fxregimes(INR ~ USD + JPY + DUR + GBP, data = inr, h = 20, breaks = 10)")

## --------------------------------------------------------------------------
## 6. gefp() process + sctest()
## --------------------------------------------------------------------------

cat("\n[6] gefp / sctest\n")

emit_gefp <- function(m, tag, desc, call) {
  g <- gefp(m, fit = NULL)
  p <- g$process
  add_fixture(cbind(data.frame(t = 0:(NROW(p) - 1L), stringsAsFactors = FALSE), zoo_df(p)),
    sprintf("gefp_%s_process.csv", tag),
    paste0(desc, ": the decorrelated cumulative score process of gefp(x, fit = NULL). Row t = 0 is the zero start, rows t = 1..n correspond to the observations; the `date` of row 0 is index[1] - (index[2] - index[1])."),
    paste0("gefp(", call, ", fit = NULL)$process"))

  add_fixture(mat_df(g$J12), sprintf("gefp_%s_J12.csv", tag),
    paste0(desc, ": root.matrix(crossprod(estfun/sqrt(n))), the symmetric matrix square root used to decorrelate the score process."),
    paste0("gefp(", call, ", fit = NULL)$J12"))

  st <- sctest(g)
  add_fixture(data.frame(statistic = c("f_efp", "p_value", "nobs", "nreg"),
                         value = c(as.numeric(st$statistic), st$p.value, g$nobs, g$nreg),
                         stringsAsFactors = FALSE),
    sprintf("sctest_%s.csv", tag),
    paste0(desc, ": M-fluctuation test with the default double-maximum functional (maxBB): test statistic max|process| and its asymptotic p-value."),
    paste0("sctest(gefp(", call, ", fit = NULL))"))
}

emit_gefp(cny_lm, "cny_hist", "CNY history model",
  'fxlm(CNY ~ USD + JPY + EUR + GBP, data = window(cny, end = as.Date("2005-10-31")))')
emit_gefp(inr_lm, "inr_full", "INR full-sample model",
  "fxlm(INR ~ USD + JPY + DUR + GBP, data = inr)")

## --------------------------------------------------------------------------
## 7. fxmonitor() -- the CNY vignette setup
## --------------------------------------------------------------------------

cat("\n[7] fxmonitor\n")

cny_mon <- fxmonitor(CNY ~ USD + JPY + EUR + GBP,
  data = window(cny, end = as.Date("2006-05-31")),
  start = as.Date("2005-11-01"), end = 4)

MONCALL <- 'fxmonitor(CNY ~ USD + JPY + EUR + GBP, data = window(cny, end = as.Date("2006-05-31")), start = as.Date("2005-11-01"), end = 4)'

add_fixture(cbind(data.frame(t = seq_len(NROW(cny_mon$process)), stringsAsFactors = FALSE),
                  zoo_df(cny_mon$process)),
  "fxmonitor_cny_process.csv",
  "CNY monitoring fluctuation process: the decorrelated cumulative score process of the history model evaluated over history AND monitoring period (rows 1..n are the history, rows n+1.. the monitoring period).",
  paste0(MONCALL, "$process"))

add_fixture(mat_df(cny_mon$J12), "fxmonitor_cny_J12.csv",
  "CNY monitoring: root.matrix(crossprod(history scores / sqrt(n))), the decorrelating matrix square root.",
  paste0(MONCALL, "$J12"))

mon_bp <- breakpoints(cny_mon)
add_fixture(data.frame(
    statistic = c("n", "nproc", "critval", "end", "alpha", "breakpoint"),
    value = c(cny_mon$n, NROW(cny_mon$process), cny_mon$critval, 4, 0.05,
              if (is.na(mon_bp)) NA_real_ else mon_bp),
    stringsAsFactors = FALSE),
  "fxmonitor_cny_summary.csv",
  "CNY monitoring: history length n, total process length, the Bonferroni-corrected boundary constant (Table III of Zeileis et al., JAE, bilinearly interpolated at end = 4 and alpha = 1-(1-0.05)^(1/6)), and the index of the first boundary crossing.",
  paste0(MONCALL))

add_fixture(data.frame(
    statistic = c("monitor_end_of_history", "break_date", "history_start", "process_end"),
    value = c(format(cny_mon$monitor, "%Y-%m-%d"),
              if (is.na(mon_bp)) "NA" else format(breakdates(cny_mon), "%Y-%m-%d"),
              format(index(cny_mon$process)[1], "%Y-%m-%d"),
              format(index(cny_mon$process)[NROW(cny_mon$process)], "%Y-%m-%d")),
    stringsAsFactors = FALSE),
  "fxmonitor_cny_dates.csv",
  "CNY monitoring: last date of the history period, the detected break date, and the first/last date of the monitored process.",
  paste0("breakdates(", MONCALL, ")"))

## --------------------------------------------------------------------------
## 8. Small synthetic regressions with known structural breaks
## --------------------------------------------------------------------------
##
## These exist so that a failing Julia test can be localised without a
## 1000-row FX series.  The generated data itself is written out, so the
## fixtures do not depend on R's RNG being reproducible -- the RNG is only
## used when *regenerating* the data.

cat("\n[8] synthetic\n")

suppressWarnings(RNGversion("3.6.0"))

make_synth <- function(seed, n, dates, gen) {
  set.seed(seed)
  z <- gen(n)
  zoo(z, dates)
}

## synth1: single regressor, one coefficient break at observation 60
d1 <- seq(as.Date("2000-01-07"), by = "week", length.out = 120)
synth1 <- make_synth(1, 120, d1, function(n) {
  x  <- rnorm(n)
  e  <- rnorm(n, sd = 0.5)
  b0 <- c(rep(0.5, 60), rep(-0.2, 60))
  b1 <- c(rep(1.0, 60), rep(0.2, 60))
  cbind(y = b0 + b1 * x + e, x = x)
})

## synth2: two regressors, two breaks at 60 and 120, the second one a
##         pure variance break (the coefficients return to the segment-1 values)
d2 <- seq(as.Date("2000-01-07"), by = "week", length.out = 180)
synth2 <- make_synth(2, 180, d2, function(n) {
  x1 <- rnorm(n); x2 <- rnorm(n)
  b1 <- c(rep(1, 60), rep(0, 60), rep(1, 60))
  b2 <- c(rep(0, 60), rep(1, 60), rep(0, 60))
  sd <- c(rep(0.4, 60), rep(0.4, 60), rep(1.2, 60))
  cbind(y = b1 * x1 + b2 * x2 + rnorm(n, sd = sd), x1 = x1, x2 = x2)
})

emit_synth <- function(z, varname, tag, formula, h, breaks, desc, gen_call) {
  add_fixture(zoo_df(z), sprintf("synth_%s_data.csv", tag),
    paste0(desc, ": the generated data itself (the Julia tests read this file, so they do not depend on R's RNG)."),
    gen_call)

  m <- fxlm(formula, data = z)
  cf <- coef(m)
  add_fixture(data.frame(term = names(cf), estimate = as.numeric(cf), stringsAsFactors = FALSE),
    sprintf("synth_%s_fxlm_coef.csv", tag),
    paste0(desc, ": full-sample fxlm coefficients including the MLE `(Variance)`."),
    sprintf("coef(fxlm(%s, data = %s))", deparse(formula), varname))

  sc <- summary(m)$coefficients
  add_fixture(data.frame(term = rownames(sc), estimate = sc[, 1], std_error = sc[, 2],
                         t_value = sc[, 3], p_value = sc[, 4], stringsAsFactors = FALSE),
    sprintf("synth_%s_fxlm_summary.csv", tag),
    paste0(desc, ": full-sample fxlm summary() coefficient table."),
    sprintf("summary(fxlm(%s, data = %s))$coefficients", deparse(formula), varname))

  add_fixture(zoo_df(estfun(m)), sprintf("synth_%s_estfun.csv", tag),
    paste0(desc, ": estfun.fxlm of the full-sample model."),
    sprintf("estfun(fxlm(%s, data = %s))", deparse(formula), varname))

  add_fixture(mat_df(bread(m)), sprintf("synth_%s_bread.csv", tag),
    paste0(desc, ": bread.fxlm of the full-sample model."),
    sprintf("bread(fxlm(%s, data = %s))", deparse(formula), varname))

  emit_recresid(formula, z, sprintf("synth_%s", tag), desc,
    sprintf("recresid(model.matrix(%s, as.data.frame(%s)), model.response(model.frame(%s, as.data.frame(%s))))",
            deparse(formula), varname, deparse(formula), varname))

  emit_triangle(formula, z, h, sprintf("synth_%s_h%d", tag, h),
    paste0(desc, ", h = ", h),
    sprintf("breakpoints(%s, data = as.data.frame(%s), h = %d)$RSS.triang", deparse(formula), varname, h))

  reg <- fxregimes(formula, data = z, h = h, breaks = breaks)
  emit_regimes(reg, sprintf("synth_%s", tag),
    paste0(desc, ", fxregimes with h = ", h, " and breaks = ", breaks),
    sprintf("fxregimes(%s, data = %s, h = %d, breaks = %d)", deparse(formula), varname, h, breaks),
    maxm = breaks, h = h)

  emit_gefp(m, sprintf("synth_%s", tag), paste0(desc, ", full-sample model"),
    sprintf("fxlm(%s, data = %s)", deparse(formula), varname))
}

emit_synth(synth1, "synth1", "break1", y ~ x, 20, 4,
  "Synthetic series: 120 weekly observations, one coefficient break after observation 60 (intercept 0.5 -> -0.2, slope 1.0 -> 0.2, sd 0.5 throughout)",
  'set.seed(1); x <- rnorm(120); e <- rnorm(120, sd = 0.5); y <- c(rep(0.5,60),rep(-0.2,60)) + c(rep(1,60),rep(0.2,60))*x + e; zoo(cbind(y, x), seq(as.Date("2000-01-07"), by = "week", length.out = 120))')

emit_synth(synth2, "synth2", "break2", y ~ x1 + x2, 25, 4,
  "Synthetic series: 180 weekly observations, breaks after observations 60 and 120; the first is a coefficient break (x1 <-> x2), the second a pure variance break (sd 0.4 -> 1.2)",
  'set.seed(2); x1 <- rnorm(180); x2 <- rnorm(180); y <- c(rep(1,60),rep(0,60),rep(1,60))*x1 + c(rep(0,60),rep(1,60),rep(0,60))*x2 + rnorm(180, sd = c(rep(0.4,120), rep(1.2,60))); zoo(cbind(y, x1, x2), seq(as.Date("2000-01-07"), by = "week", length.out = 180))')

## --------------------------------------------------------------------------
## 9. MANIFEST
## --------------------------------------------------------------------------

cat("\n[9] MANIFEST\n")

man <- do.call(rbind, manifest)
man <- man[order(man$file), ]
write.csv(man, file.path(FIXDIR, "MANIFEST.csv"), row.names = FALSE, quote = TRUE)
cat(sprintf("  %-46s %6d x %2d\n", "MANIFEST.csv", nrow(man), ncol(man)))

sz <- sum(file.info(list.files(FIXDIR, full.names = TRUE))$size)
cat(sprintf("\n%d fixture files + MANIFEST.csv, %.2f MB total\n", nrow(man), sz / 1024^2))
