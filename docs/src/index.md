```@meta
CurrentModule = FXRegime
```

# FXRegime.jl

De facto exchange rate regime classification in Julia: estimate the implicit currency basket a
managed exchange rate is anchored to, test whether that anchor is stable, date the breaks
between regimes, and monitor a live regime for change.

FXRegime.jl is a port of the R package
[**fxregime**](https://CRAN.R-project.org/package=fxregime) by Achim Zeileis, Ajay Shah and
Ila Patnaik, together with the parts of
[**strucchange**](https://CRAN.R-project.org/package=strucchange) it depends on. The test
suite asserts numerical agreement with R.

## Start with the README

The narrative documentation — what the Frankel–Wei regression is, why the residual variance is
treated as a model parameter rather than a nuisance, how the Bai–Perron dating and the
prospective monitoring procedures work, a fully worked CNY example, the R-to-Julia function
map and the list of what is deliberately not ported — lives in the repository README:

**[README.md on GitHub](https://github.com/xKDR/FXRegime.jl/blob/main/README.md)**

The rest of this page is the generated API reference. Every exported function's docstring
names the R function it ports.

```julia
using FXRegime, Dates

fx  = FXRatesCHF()                                  # daily rates, 25 currencies vs CHF
cny = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                frequency = :daily,
                start = Date(2005, 7, 25), stop = Date(2009, 7, 31))

fm  = fxlm(cny)                                     # Frankel-Wei regression
sctest(gefp(fm))                                    # is the regression stable?

reg = fxregimes(cny; h = 20, breaks = 10)           # date the regimes
breakdates(reg)
coef(reg)
confint(reg; level = 0.9)
refit(reg)                                          # one FXLM per regime

fxmonitor(cny; start = Date(2005, 11, 1))           # prospective monitoring
```

## Citation

> Zeileis A, Shah A, Patnaik I (2010). "Testing, Monitoring, and Dating Structural Changes in
> Exchange Rate Regimes." *Computational Statistics & Data Analysis*, **54**(6), 1696–1706.
> [doi:10.1016/j.csda.2009.12.005](https://doi.org/10.1016/j.csda.2009.12.005)

## API index

```@index
```

## API reference

```@autodocs
Modules = [FXRegime]
```
