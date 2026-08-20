"""
    FXRegime

A Julia port of the R package [**fxregime**](https://CRAN.R-project.org/package=fxregime)
(Achim Zeileis, Ajay Shah, Ila Patnaik) together with the parts of **strucchange** it depends
on: exchange rate regime analysis by dating, testing and monitoring structural change in
Frankel-Wei style regressions of currency returns.

The typical workflow is

```julia
using FXRegime, Dates

fx  = FXRatesCHF()                                  # daily rates, 25 currencies vs CHF
cny = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                frequency = :daily,
                start = Date(2005, 7, 25), stop = Date(2009, 7, 31))

m   = fxlm(cny)                                     # Frankel-Wei regression
sctest(gefp(m))                                     # is the regression stable?

reg = fxregimes(cny; h = 20, breaks = 10)           # date the regimes
breakdates(reg)
coef(reg)
confint(reg; level = 0.9)
refit(reg)                                          # one `FXLM` per regime
```

All estimates reproduce the R package numerically; see the package tests for the fixtures.
"""
module FXRegime

using Statistics
using LinearAlgebra
using Dates
using Printf
using DelimitedFiles
using Distributions

export
    # data containers
    FXSeries, FXRatesCHF, window, getcols, colindex, index, nobs,
    # returns
    fxreturns, nextfriday,
    # single-regime model
    FXLM, fxlm, coef, coefnames, residuals, fitted, npar, dof_residual, sigma, sigma2,
    r2, estfun, bread, vcov, coeftable, fstatistic, fxpegtest,
    # structural change machinery
    recresid, rss_triangle, nll_triangle, BreakpointsFull, breakpoints, breakdates,
    loglik, dof, information_criterion, bic, lwz, summarize,
    # regime dating
    FXRegimes, fxregimes, refit, confint,
    # fluctuation tests and monitoring
    GEFP, gefp, sctest, FXMonitor, fxmonitor,
    # extensions
    r2_decomposition

include("series.jl")
include("data.jl")
include("returns.jl")
include("fxlm.jl")
include("recresid.jl")
include("breakpoints.jl")
include("efp.jl")
include("monitor.jl")
include("regimes.jl")
include("r2decomp.jl")
include("confint.jl")

end # module
