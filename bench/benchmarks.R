#!/usr/bin/env Rscript
## R side of the FXRegime.jl <-> R fxregime benchmark.
##
##   Rscript bench/benchmarks.R
##
## Writes bench/results_r.csv with one row per (task, n) and the median wall time in seconds,
## measured on exactly the same inputs the Julia benchmark uses.
##
## Note: fxregimes() is deliberately called the way the vignettes call it, i.e. without `fit`
## or `objfun`, so that gbreakpoints() takes its fast path through strucchange::breakpoints()
## and builds the RSS triangle from recursive residuals. Both implementations are therefore
## running the same O(n^2) algorithm and the comparison is a constant-factor one.

suppressMessages({
  library(fxregime)
  library(strucchange)
  library(zoo)
})

out <- file.path("bench", "results_r.csv")
data("FXRatesCHF", package = "fxregime")

## Capture the expression unevaluated once, then evaluate it repeatedly in the caller's frame.
## (Do not force() it, and do not call substitute() inside replicate() — the promise is gone by
## then and substitute() would resolve in replicate's frame, not the caller's.)
timeit <- function(expr, k = 3) {
  e  <- substitute(expr)
  pf <- parent.frame()
  eval(e, pf)                                         # warm-up
  ts <- numeric(k)
  for (i in seq_len(k)) ts[i] <- system.time(eval(e, pf))[["elapsed"]]
  median(ts)
}

rows <- list()
record <- function(task, ds, n, t) {
  rows[[length(rows) + 1L]] <<- data.frame(task = task, dataset = ds, n = n, seconds = t)
  cat(sprintf("%-22s %-10s n=%-6d %9.4f s\n", task, ds, n, t))
}

## The two published analyses -------------------------------------------------
cny <- fxreturns("CNY", frequency = "daily",
                 start = as.Date("2005-07-25"), end = as.Date("2009-07-31"),
                 other = c("USD", "JPY", "EUR", "GBP"), data = FXRatesCHF)
inr <- fxreturns("INR", frequency = "weekly",
                 start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
                 other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF)
long <- fxreturns("INR", frequency = "daily",
                  start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
                  other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF)

## 1. Data preparation --------------------------------------------------------
record("fxreturns", "CNY-daily", nrow(cny),
       timeit(fxreturns("CNY", frequency = "daily",
                        start = as.Date("2005-07-25"), end = as.Date("2009-07-31"),
                        other = c("USD", "JPY", "EUR", "GBP"), data = FXRatesCHF)))
record("fxreturns", "INR-weekly", nrow(inr),
       timeit(fxreturns("INR", frequency = "weekly",
                        start = as.Date("1993-04-01"), end = as.Date("2008-01-04"),
                        other = c("USD", "JPY", "DUR", "GBP"), data = FXRatesCHF)))

## 2. Single regression -------------------------------------------------------
record("fxlm", "CNY-daily", nrow(cny),
       timeit(fxlm(CNY ~ USD + JPY + EUR + GBP, data = cny), k = 20))

## 3. Recursive residuals -----------------------------------------------------
d <- as.data.frame(cny)
X <- model.matrix(CNY ~ USD + JPY + EUR + GBP, data = d)
y <- d$CNY
record("recresid", "CNY-daily", nrow(X), timeit(recresid(X, y), k = 5))

## 4. The O(n^2) triangle, isolated ------------------------------------------
for (nm in c("CNY-daily", "INR-weekly")) {
  s <- if (nm == "CNY-daily") cny else inr
  f <- if (nm == "CNY-daily") CNY ~ USD + JPY + EUR + GBP else INR ~ USD + JPY + DUR + GBP
  local({
    ss <- s; ff <- f
    record("rss_triangle", nm, nrow(ss),
           timeit(strucchange::breakpoints(ff, h = 20, data = as.data.frame(ss),
                                           breaks = 1)$RSS.triang, k = 2))
  })
}

## 5. Full dating -------------------------------------------------------------
record("fxregimes", "CNY-daily", nrow(cny),
       timeit(fxregimes(CNY ~ USD + JPY + EUR + GBP, data = cny, h = 20, breaks = 10), k = 1))
record("fxregimes", "INR-weekly", nrow(inr),
       timeit(fxregimes(INR ~ USD + JPY + DUR + GBP, data = inr, h = 20, breaks = 10), k = 1))

## 6. Confidence intervals ----------------------------------------------------
cny_reg <- fxregimes(CNY ~ USD + JPY + EUR + GBP, data = cny, h = 20, breaks = 10)
record("confint", "CNY-daily", nrow(cny), timeit(confint(cny_reg, level = 0.9), k = 5))

## 7. Scaling in n ------------------------------------------------------------
for (n in c(250, 500, 1000)) {
  if (n > nrow(long)) next
  local({
    s <- long[1:n, ]
    record("fxregimes-scaling", "INR-daily", n,
           timeit(fxregimes(INR ~ USD + JPY + DUR + GBP, data = s, h = 20, breaks = 10), k = 1))
  })
}

res <- do.call(rbind, rows)
write.csv(res, out, row.names = FALSE)
cat("\nR", getRversion(), "fxregime", as.character(packageVersion("fxregime")),
    "strucchange", as.character(packageVersion("strucchange")), "\n")
cat("wrote", out, "\n")
