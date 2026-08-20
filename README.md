# FXRegime.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xKDR.github.io/FXRegime.jl/dev)
![Build Status](https://github.com/xKDR/FXRegime.jl/actions/workflows/ci.yml/badge.svg)
![Docs](https://github.com/xKDR/FXRegime.jl/actions/workflows/documentation.yml/badge.svg)

**De facto exchange rate regime classification in Julia**: estimate the implicit currency
basket a managed exchange rate is anchored to, test whether that anchor is stable, date the
breaks between regimes, and monitor a live regime for change.

FXRegime.jl is a port of the R package
[**fxregime**](https://CRAN.R-project.org/package=fxregime) by Achim Zeileis, Ajay Shah and
Ila Patnaik, together with the parts of
[**strucchange**](https://CRAN.R-project.org/package=strucchange) it depends on. The test
suite asserts numerical agreement with R, fixture by fixture.

---

## Contents

- [What the Frankel–Wei regression is](#what-the-frankelwei-regression-is)
- [Why the variance is a parameter, not a nuisance](#why-the-variance-is-a-parameter-not-a-nuisance)
- [Dating: where did the regime change?](#dating-where-did-the-regime-change)
- [Monitoring: has the regime changed *yet*?](#monitoring-has-the-regime-changed-yet)
- [Installation](#installation)
- [Worked example: the Chinese yuan, 2005–2009](#worked-example-the-chinese-yuan-20052009)
- [Relation to the R package](#relation-to-the-r-package)
- [Tests and R fixtures](#tests-and-r-fixtures)
- [Data](#data)
- [Citation](#citation)

---

## What the Frankel–Wei regression is

Central banks announce exchange rate regimes; they do not always run them. A country may say
it targets "a basket of currencies with greater flexibility" while in practice holding a hard
peg to a single currency. The *de facto* classification problem is to recover, from observed
exchange rates alone, (i) which currencies the target is actually anchored to and with what
weights, and (ii) how tightly.

Frankel and Wei (1994) turn this into a linear regression. Suppose the monetary authority
manages its currency *T* against a basket of currencies `1, …, k` with weights
`w₁, …, w_k`, so that the value of *T* measured in units of that basket is (approximately)
constant. Pick a **numeraire currency** *N* that is not in the basket and does not co-move
with it, and express every exchange rate as a price in terms of *N*. Then taking logs and
differencing the peg relation gives the estimating equation

```
y_t = β₀ + β₁ x_{1,t} + β₂ x_{2,t} + … + β_k x_{k,t} + ε_t,     ε_t ~ N(0, σ²)
```

where

```
y_t   = 100 · Δ log (price of T in numeraire N)
x_{j,t} = 100 · Δ log (price of basket currency j in numeraire N)
```

The slope `β_j` estimates the weight of currency *j* in the implicit basket; `β₀` is the
average drift of the target against the basket (a crawling peg has `β₀ ≠ 0`); and `σ` is the
residual volatility around the basket — how much the authority lets the rate move.

Under a pure basket peg the slopes sum to one; under a hard single-currency peg one slope is
one and the others are zero. `fxpegtest` runs exactly that Wald test.

**Why the numeraire has to be independent.** The choice of *N* is not innocuous. If *N* is
itself in the basket — say you measure a USD-pegged currency against the USD — then `y_t ≡ 0`
by construction, there is nothing left to explain, and the weights are not identified. More
generally, if the numeraire co-moves with the basket, its own fluctuations are common to
both sides of the regression, the regressors are contaminated, and the estimated weights are
biased toward whatever the numeraire is tracking. The numeraire must be a currency that
floats freely and independently of the basket, so that its movements act as a source of
exogenous variation that traces out the basket weights. The bundled data set uses the Swiss
franc as numeraire for exactly this reason: CHF floats, and it is not a plausible member of
the baskets of the currencies one wants to study.

Cross-currency returns rather than levels are used because the peg relation is a statement
about *changes*: a peg constrains how the rate moves, not what number it sits at.

## Why the variance is a parameter, not a nuisance

The coefficients tell you the **anchor** of a regime. The residual standard deviation `σ`
tells you its **flexibility**. Both are policy, and both change when policy changes.

Two regimes can have identical basket weights and still be completely different regimes: a
currency that tracks the dollar to within `σ = 0.03` per cent per day is on a hard peg; the
same currency tracking the dollar with `σ = 0.26` per cent per day is on a loose band or a
managed float. Conversely, a currency can widen its band without touching its basket at all.
In the worked example below, China's mid-2006 regime change shows up almost entirely in
`σ` (0.028 → 0.106) with `β_USD` moving only from 0.999 to 0.969 — a coefficients-only
procedure would barely see it.

So FXRegime.jl treats `σ²` as a parameter of the model on the same footing as the
coefficients, not as a nuisance to be conditioned away. Three consequences run through the
whole package:

1. **The parameter vector is `[β₀, β₁, …, β_k, σ²]`,** with `σ²` last and named
   `"(Variance)"`. `coef(::FXLM)` returns all `k + 2` entries, and `σ̂² = mean(residuals²)`
   is the maximum-likelihood (not degrees-of-freedom-adjusted) estimate.

2. **Breaks are dated by likelihood, not least squares.** For a segment of `n_i`
   observations the maximised Gaussian log-likelihood is

   ```
   ℓ_i = −(n_i / 2) · ( log(2π) + log(RSS_i / n_i) + 1 )
   ```

   and the dating minimises `Σ_i (−ℓ_i)` over partitions. This is *not* the same objective as
   minimising `Σ_i RSS_i`: the log and the `n_i`-weighting mean that a segment which halves
   its residual variance is rewarded even when it contributes little to the total sum of
   squares. Least-squares dating would place breaks only where the conditional mean moves;
   likelihood dating also finds pure variance breaks — the flexibility changes that matter
   most for regime classification.

3. **The score process has a variance component.** The estimating function of an `fxlm` is

   ```
   ψ_t = [ e_t · x_t / σ² ,  (e_t² − σ²) / (2σ⁴) ]
   ```

   with `k + 2` columns, so the fluctuation tests (`gefp` / `sctest`) and the monitoring
   procedure (`fxmonitor`) track a `(k + 2)`-dimensional process in which one component is
   dedicated to the variance. In the CNY example it is that component, not any of the
   coefficient components, that crosses the boundary.

## Dating: where did the regime change?

`fxregimes` performs retrospective, Bai–Perron style dating of an arbitrary number of breaks.

**Segment likelihoods.** For every pair `(i, j)` with `j − i + 1 ≥ h` the maximised segment
log-likelihood is needed. Refitting `O(n²)` regressions directly would be `O(n³)`; instead
each row of the triangle is obtained from one pass of **recursive residuals**
(Brown–Durbin–Evans updating, `recresid`), whose squares cumulate into the residual sums of
squares of all segments starting at `i`. The triangle is built once, in parallel over its
rows, and the negative log-likelihood table follows from the formula above.

**Dynamic programming.** With `OBJ(m, i)` the best total negative log-likelihood of an
`m`-break partition of observations `1 … i`,

```
OBJ(m, i) = min over j ∈ [(m−1)·h, i−h] of  OBJ(m−1, j) + nll(j+1, i)
```

and the optimal `m`-break partition of the whole sample is recovered by minimising
`OBJ(m, i) + nll(i+1, n)`. This is exact — it is a global optimum over all admissible
partitions, not a sequential or greedy search — and costs `O(n²)` once the triangle exists.

**Minimal segment size `h`.** Every regime must contain at least `h` observations. `h` is
given either as a fraction of the sample (`h = 0.15`, the default, means `floor(0.15n)`) or
as a count (`h = 20`). It is a genuine modelling choice: it is the shortest episode you are
prepared to call a regime, and it bounds the number of breaks at `⌈n/h⌉ − 2`.

**How many breaks.** Every extra break increases the likelihood, so the number of regimes is
chosen by an information criterion. With `m` breaks the segmented model uses
`df = (k + 2)(m + 1) + m` parameters — `k + 2` per segment plus the `m` break dates — and

```
BIC(m) = 2·NLL(m) + df · log n
LWZ(m) = 2·NLL(m) + df · 0.299 · (log n)^2.1
```

LWZ (Liu, Wu and Zidek, 1997) is the default because BIC is well known to over-select breaks
in segmented regressions; both are always reported by `summarize`, and `ic = :BIC` switches
the selection. In the CNY example below BIC picks 5 breaks and LWZ picks 3.

**Confidence intervals for the break dates.** Under Bai's (1997) asymptotics the estimated
break date, rescaled by the magnitude of the parameter change, converges to the **argmax of a
two-sided Brownian motion with triangular drift** — a process built from two independent
Brownian motions `W₁`, `W₂` running left and right from the break, each pulled down by a
linear drift,

```
V(x) ∝ W₁(−x) − |x|/2       for x < 0,      V(x) ∝ W₂(x) − |x|/2      for x ≥ 0
```

with the two halves scaled relative to one another by constants `ξ, φ₁, φ₂` that depend on the
size of the parameter change and on the information on either side. The distribution function
of `argmax V` is available in closed form (this is strucchange's internal `pargmaxV`, which is
ported verbatim, log-scale normal CDFs and all). `confint` evaluates it, inverts it at
`(1 ± level)/2` by outward bracketing plus Brent root-finding, and rescales back to
observation numbers. The scaling quantities `ξ, φ₁, φ₂`
are quadratic forms in `δ = θ̂_{i+1} − θ̂_i`, the change in the **full** parameter vector
**including `(Variance)`** — consistent with treating `σ²` as a parameter — with `Q` the
inverse of the `fxlm` bread and `Ω` the score covariance (by default `Ω = Q`, i.e. the
information-matrix-equality form; pass `meat` for a sandwich version). Intervals are
asymmetric in general: a break into a high-variance regime is pinned down precisely from the
low-variance side and loosely from the other.

## Monitoring: has the regime changed *yet*?

`fxmonitor` answers a different question from `fxregimes`, and the difference is the
difference between **prospective** and **retrospective** inference.

| | `fxregimes` (dating) | `fxmonitor` (monitoring) |
|---|---|---|
| Data | the whole sample, in hand | a fixed *history* period, then observations arriving one at a time |
| Question | where were the breaks? | has a break happened by now? |
| Decision | once, after seeing everything | at every new observation |
| Output | break dates + confidence intervals | a first-crossing date, or none |
| Size control | model selection criterion | over the entire monitoring horizon, at level `α` |
| Timing | breaks placed at their most likely date | detection lags the break, by construction |

Monitoring fits an `fxlm` on the history period — a stretch you are willing to assume is one
regime — freezing `θ̂ = (β̂, σ̂²)`. Every subsequent observation contributes its score under
those frozen parameters, and the cumulative, decorrelated score path is compared with a
**linear boundary** `critval · t / n` (`n` = history length). Because the boundary must hold
for *all* future `t` up to the monitoring horizon, `critval` cannot be a pointwise normal
quantile: it comes from Table III of Zeileis, Leisch, Kleiber and Hornik (2005), bilinearly
interpolated at the horizon (`stop`, a multiple of the history length) and at the
Bonferroni-corrected level `1 − (1 − α)^(1/(k+2))`, one component of the score at a time.

The upshot is that a monitoring alarm is a *real-time* statement — it could have been made on
the day it fires, using only data available then — whereas a dated break is a *hindsight*
statement, and generally an earlier one. In the example below the monitor fires on
2006-03-27; retrospective dating on the full sample puts the break on 2006-03-14, with a 90%
interval of 2006-02-21 to 2006-03-15. The two agree, and the gap between them is the
detection delay.

## Installation

The package is not registered. From the Julia REPL:

```julia
using Pkg
Pkg.add(url = "https://github.com/xKDR/FXRegime.jl.git")
```

or in the Pkg REPL mode (`]`):

```
add https://github.com/xKDR/FXRegime.jl.git
```

Julia 1.10 or later. The only non-stdlib dependency is
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl); there is no plotting
dependency and no R is required at runtime.

## Worked example: the Chinese yuan, 2005–2009

In July 2005 China abandoned its fixed exchange rate to the US dollar and announced a move to
a basket of currencies with greater flexibility, naming USD, JPY, EUR and KRW as basket
members. This reproduces the analysis in Zeileis, Shah and Patnaik (2010) and the `CNY`
vignette of the R package: was there in fact a change, when, and of what kind?

Everything below was run as shown; the outputs are copied from the session.

### The data

```julia
using FXRegime, Dates

fx = FXRatesCHF()
```

```
FXSeries: 9819 observations × 25 series
  1971-01-04 .. 2010-02-12
  USD, JPY, DUR, EUR, DEM, GBP, AUD, BRL, CAD, CNY, DKK, HKD, INR, MYR, MXN, NOK, NZD, KRW, SEK, SGD, LKR, ZAR, TWD, THB, VEB
  1971-01-04  0.231589  82.846225  0.428514       NA  0.84377  0.096745  0.208132       NA  0.234113       NA  1.733441       NA       NA  0.714845       NA  1.652594  0.207927       NA  1.195994       NA       NA  0.166466       NA       NA       NA
  1971-01-05  0.231927  82.985829  0.429329       NA  0.845374  0.096842  0.208343       NA  0.234293       NA  1.736484       NA       NA  0.71575       NA  1.65647  0.208137       NA  1.197393       NA       NA  0.16664       NA       NA       NA
  1971-01-06  0.231949  83.005126  0.429239       NA  0.845198  0.096778  0.208212       NA  0.234407       NA  1.73588       NA       NA  0.715445       NA  1.655348  0.208007       NA  1.19718       NA       NA  0.166539       NA       NA       NA
  1971-01-07  0.232002  83.026703  0.429268       NA  0.845254  0.096817  0.208298       NA  0.235436       NA  1.736538       NA       NA  0.715727       NA  1.656172  0.208074       NA  1.198269       NA       NA  0.166601       NA       NA       NA
  ⋮
```

`fxreturns` selects the target and basket columns, restricts the window, carries the last
observation forward over missing values and returns `100 · Δlog`:

```julia
cny = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                frequency = :daily,
                start = Date(2005, 7, 25), stop = Date(2009, 7, 31))
```

```
FXSeries: 1014 observations × 5 series
  2005-07-26 .. 2009-07-31
  CNY, USD, JPY, EUR, GBP
  2005-07-26  -0.320959  -0.323426  0.310817  0.067324  -0.093981
  2005-07-27  0.197334  0.161582  0.161582  -0.179397  0.029489
  2005-07-28  0.76043  0.819613  0.855226  0.215513  0.196534
  2005-07-29  0.048065  0.07767  -0.002476  0.01171  -0.184123
  ⋮
```

Column 1 is the target currency by convention; the rest are the basket.

### The regime immediately after the announcement

Take the first three months as a history period and fit the Frankel–Wei regression:

```julia
hist = window(cny; stop = Date(2005, 10, 31))
fm = fxlm(hist)
```

```
Frankel-Wei regression (fxlm)
  68 observations, 2005-07-26 .. 2005-10-31

                     Estimate   Std. Error   t value   Pr(>|t|)
  (Intercept)       -0.004782     0.003688    -1.297      0.199
  USD                0.999653     0.008779   113.868   1.15e-74
  JPY                0.004668     0.010669     0.437      0.663
  EUR               -0.014238     0.026516    -0.537      0.593
  GBP               -0.007744     0.014568    -0.532      0.597

  Residual std. error: 0.0295329 on 63 degrees of freedom
  (Variance) [MLE]:    0.000808059
  R-squared: 0.997926,  F(4,63) = 7577.26, p = 9.824e-84
```

The USD weight is 0.9997 and nothing else is distinguishable from zero: despite the
announcement, this is a dollar peg, and a very tight one — `σ = 0.0284` per cent per day.
The formal test of "USD weight = 1, all other weights = 0" does not reject:

```julia
fxpegtest(fm)
```

```
(statistic = 0.20983842127011854, pvalue = 0.9320482013182835, df1 = 4, df2 = 63, peg = "USD", hypothesis = ["USD = 1", "JPY = 0", "EUR = 0", "GBP = 0"])
```

Over the full four-year sample, however, the model is emphatically not stable. The
generalised empirical fluctuation process and its double-maximum test:

```julia
sctest(gefp(fxlm(cny)))
```

```
(statistic = 2.222962076547519, pvalue = 0.0006121891347575792)
```

### Dating the regimes

```julia
reg = fxregimes(cny; h = 20, breaks = 10)
```

```

FX model: CNY ~ USD + JPY + EUR + GBP

Minimum LWZ partition
Number of regimes: 4
Breakdates: 2006-03-14, 2008-08-22, 2008-12-31
```

The full model-selection table — negative log-likelihood, BIC and LWZ for `m = 0 … 10`:

```julia
s = summarize(reg);
[s.m round.(s.nll, digits = 1) round.(s.BIC, digits = 1) round.(s.LWZ, digits = 1)]
```

```
11×4 Matrix{Float64}:
  0.0   -725.8  -1410.0  -1347.2
  1.0   -886.4  -1682.8  -1546.8
  2.0  -1008.6  -1878.8  -1669.6
  3.0  -1094.7  -2002.5  -1720.0
  4.0  -1141.3  -2047.3  -1691.6
  5.0  -1171.6  -2059.5  -1630.6
  6.0  -1192.5  -2052.8  -1550.6
  7.0  -1211.2  -2041.7  -1466.4
  8.0  -1223.0  -2016.8  -1368.3
  9.0  -1233.3  -1989.0  -1267.2
 10.0  -1245.1  -1964.1  -1169.1
```

The likelihood improves monotonically, BIC bottoms out at 5 breaks and LWZ — the default,
and the more conservative penalty — at 3. `summarize` also carries `s.breakpoints` and
`s.breakdates` for every `m`, so the sensitivity of the dates to the criterion can be
inspected directly.

The parameters of the four regimes (columns in the order given by `coefnames`):

```julia
coefnames(reg)
```

```
6-element Vector{String}:
 "(Intercept)"
 "USD"
 "JPY"
 "EUR"
 "GBP"
 "(Variance)"
```

```julia
round.(coef(reg), digits = 4)
```

```
4×6 Matrix{Float64}:
 -0.005   0.9994   0.0052  -0.0152   0.0068  0.0008
 -0.025   0.9694  -0.0093   0.0256  -0.0129  0.0113
 -0.0148  1.0307  -0.0265   0.0489   0.0072  0.0694
  0.0014  0.9809   0.0082  -0.0077   0.0086  0.002
```

```julia
sqrt.(coef(reg)[:, end])
```

```
4-element Vector{Float64}:
 0.02795857986249238
 0.10623383379286296
 0.2634330230239816
 0.044440068683792744
```

Read across: the anchor barely moves — USD keeps a weight near one in every regime, and no
other currency ever acquires a meaningful weight, so the announced basket never materialised.
What moves is the flexibility, `σ`: 0.028 (rigid peg) → 0.106 (a modest loosening from March
2006, with a negative intercept of −0.025 reflecting a slow crawling appreciation) → 0.263
(the autumn-2008 crisis) → 0.044 with an intercept back at zero (a return to a tight,
non-appreciating dollar peg from January 2009). This is precisely the kind of regime history
that a coefficients-only method would flatten into "pegged to USD throughout".

Confidence intervals for the three break dates:

```julia
ci = confint(reg; level = 0.9);
[ci.lower ci.breakpoints ci.upper]
```

```
3×3 Matrix{Union{Missing, Int64}}:
 143  158  159
 762  778  779
 865  866  880
```

```julia
[ci.dates.lower ci.dates.breakpoints ci.dates.upper]
```

```
3×3 Matrix{Date}:
 2006-02-21  2006-03-14  2006-03-15
 2008-07-31  2008-08-22  2008-08-25
 2008-12-30  2008-12-31  2009-01-22
```

Note the asymmetry, and that it flips: the end of a low-variance regime is dated precisely
from the right and loosely from the left (breaks 1 and 2), while the end of the
high-variance crisis regime is dated precisely from the left and loosely from the right
(break 3). A quiet regime carries a lot of information about where it stopped being quiet.

Each regime can be pulled out as a full model with `refit` — here the crisis segment:

```julia
refit(reg)[3]
```

```
Frankel-Wei regression (fxlm)
  88 observations, 2008-08-25 .. 2008-12-31

                     Estimate   Std. Error   t value   Pr(>|t|)
  (Intercept)       -0.014770     0.029756    -0.496      0.621
  USD                1.030744     0.043672    23.602   1.43e-38
  JPY               -0.026479     0.030149    -0.878      0.382
  EUR                0.048853     0.058852     0.830      0.409
  GBP                0.007187     0.035289     0.204      0.839

  Residual std. error: 0.271252 on 83 degrees of freedom
  (Variance) [MLE]:    0.069397
  R-squared: 0.95616,  F(4,83) = 452.558, p = 1.767e-55
```

### Monitoring in real time

Freeze the three-month history model and let the observations arrive:

```julia
mon = fxmonitor(cny; start = Date(2005, 11, 1))
```

```
Monitoring of FX model

History period: 2005-07-26 to 2005-10-31
Break detected: 2006-03-27
```

```julia
breakpoints(mon)
```

```
167
```

```julia
breakdates(mon)
```

```
2006-03-27
```

Which component of the process crossed its boundary is itself the answer to "what kind of
change?":

```julia
t = breakpoints(mon)
boundary = mon.critval * t / mon.n
mon.coefnames[findall(j -> abs(mon.process[t, j]) > boundary, axes(mon.process, 2))]
```

```
1-element Vector{String}:
 "(Variance)"
```

The alarm is raised on 2006-03-27, using only data available on that day, and it is the
*variance* component that crosses — none of the basket weights does. The change China made in
spring 2006 was to flexibility, not to the anchor. Dating on the full sample later places the
break on 2006-03-14, thirteen days earlier and comfortably inside the monitor's detection
delay. Prospective and retrospective inference
tell the same story about the same event; they just cannot tell it at the same time.

### An extension: whose variation is it?

Basket currencies are strongly correlated with one another, so ordinary `R²` contributions
are order-dependent and hard to read. `r2_decomposition` gives the order-independent Shapley
decomposition of `R²` across the basket, regime by regime (rows = regimes, columns = USD,
JPY, EUR, GBP). This function has **no** counterpart in the R package:

```julia
round.(r2_decomposition(reg), digits = 3)
```

```
4×4 Matrix{Float64}:
 0.689  0.084  0.115  0.109
 0.643  0.065  0.128  0.129
 0.609  0.197  0.046  0.104
 0.762  0.165  0.022  0.048
```

## Relation to the R package

FXRegime.jl is a port of **fxregime 1.0-5** (Achim Zeileis, Ajay Shah, Ila Patnaik; with Anmol
Sethy) plus the pieces of **strucchange** that fxregime builds on. It is not a
reimplementation "in the spirit of" the original: the objective is agreement to within
floating-point noise (most R comparisons in the suite are asserted at `1e-12` relative or
tighter), checked against CSV fixtures generated by R (see
[Tests and R fixtures](#tests-and-r-fixtures)). Where the two could differ, R wins — including
quirks such as R's `quantile` type 7 for trimming, `zoo::na.locf` dropping leading `NA`s,
`which(colnames(data) %in% x)` returning columns in *data* order rather than the order you
asked for, and R's `uniroot` bracketing.

Two structural differences of interface, neither of them numerical:

- **No formulas.** R writes `fxlm(CNY ~ USD + JPY + EUR + GBP, data = cny)`. Julia takes the
  model from the column order of the `FXSeries`: column 1 is the target, the rest are the
  basket, and an intercept is always included. `fxlm(cny)`.
- **No S3.** The package defines and exports its own `coef`, `residuals`, `fitted`, `vcov`,
  `confint`, `breakpoints`, … rather than committing type piracy on StatsAPI, which it does
  not depend on. If you also have StatsBase loaded, disambiguate with `FXRegime.coef`.

### Function map

| R (fxregime / strucchange / zoo) | FXRegime.jl |
|---|---|
| `data("FXRatesCHF")` | `FXRatesCHF()` |
| `fxreturns(x, other, data, frequency, start, end, na.action, trim)` | `fxreturns(target, other, data; frequency, start, stop, na_action, trim)` |
| `window()`, `index()`, `NROW()` | `window()`, `index()`, `nobs()`, `size()` |
| `fxlm(formula, data)` | `fxlm(data::FXSeries)` |
| `coef.fxlm` (incl. `(Variance)`) | `coef(::FXLM)` |
| `summary(fxlm)` | `show`, `coeftable`, `r2`, `sigma`, `sigma2`, `fstatistic`, `dof_residual` |
| `residuals`, `fitted` | `residuals`, `fitted` |
| `estfun.fxlm` | `estfun` |
| `bread.fxlm` | `bread` |
| `vcov`, `sandwich::vcovHC(., "HC3")` | `vcov(m)`, `vcov(m, :HC3)` |
| `fxpegtest` (via `car::linearHypothesis`) | `fxpegtest` |
| `strucchange::recresid` | `recresid` |
| `strucchange:::RSS.triang` | `rss_triangle` |
| `fxregime:::RSS2obj.triang` | `nll_triangle` |
| `fxregime:::gbreakpoints` (object) | `BreakpointsFull` |
| `breakpoints.gbreakpointsfull` | `breakpoints(::BreakpointsFull)` |
| `fxregimes(formula, data, h, breaks, ic)` | `fxregimes(data; h, breaks, ic)` |
| `breakpoints`, `breakdates` | `breakpoints`, `breakdates` |
| `summary.gbreakpointsfull` | `summarize` |
| `logLik`, `AIC.gbreakpoints` | `loglik`, `dof`, `information_criterion`, `bic`, `lwz` |
| `refit.fxregimes` | `refit` |
| `coef/fitted/residuals.fxregimes` | `coef(::FXRegimes)`, `fitted`, `residuals` |
| `confint.fxregimes`, `strucchange:::pargmaxV` | `confint(::FXRegimes)`, `FXRegime.pargmaxV` |
| `fxmonitor(formula, data, start, end, alpha, meat.)` | `fxmonitor(data; start, stop, alpha, meat)` |
| `breakpoints/breakdates.fxmonitor` | `breakpoints(::FXMonitor)`, `breakdates(::FXMonitor)` |
| `strucchange::gefp(x, fit = NULL)` | `gefp(::FXLM)` |
| `strucchange::sctest` (`maxBB`) | `sctest(::GEFP)` |
| `strucchange:::root.matrix` | `FXRegime.root_matrix` |
| — (no R equivalent) | `r2_decomposition` |

### What is *not* ported

Honestly and in full:

- **Every plotting method.** R's `plot.fxmonitor`, `plot.gbreakpointsfull`,
  `plot.summary.gbreakpointsfull`, `lines.fxregimes` and `lines.confint.fxregimes` have no
  Julia counterpart. FXRegime.jl takes no plotting dependency. The ingredients are all
  exported — `mon.process`, `mon.critval`, `mon.n`, `gefp(...).process`, `summarize(reg)`,
  `confint(reg)` — so the same figures are a few lines of Plots.jl or Makie.jl away, but you
  have to draw them yourself.
- **The formula interface** (`CNY ~ USD + JPY + EUR + GBP`) and the `print`/`summary` layouts
  that go with it. Julia's `show` methods convey the same information in their own format;
  they are not character-for-character reproductions of R's output.
- **`sctest` functionals other than `maxBB`.** R's strucchange offers `supLM`, `meanL2BB`,
  `maxL2BB`, `catL2BB`, `rangeBB` and more; only the double-maximum functional is
  implemented, and `sctest(p; functional = :supLM)` throws.
- **The rest of the `strucchange::efp` family** — OLS-CUSUM, RE/Rec-CUSUM, MOSUM, ME and the
  `sctest.formula` interface. Only the generalised (score-based) process `gefp` used by
  fxregime is ported.
- **General linear hypotheses on an `fxlm`.** R registers `linearHypothesis.fxlm` so that any
  `car`-style restriction can be tested; here only the pegged-regime hypothesis is available,
  via `fxpegtest`.
- **`hpc = "foreach"`.** R offers an optional `foreach` backend for the dynamic program.
  FXRegime.jl instead threads the RSS triangle natively (`Threads.@threads`); start Julia
  with `-t auto` to use it. There is no user-facing `hpc` argument.
- **`time()` methods** — use `index()`.
- **The vignettes** (`vignette("CNY")`, `vignette("INR")`). Their content is reproduced in
  this README and in the test suite, not as separate documents.
- **`fxregime`'s dependence on `zoo` semantics beyond what `FXSeries` provides.** `FXSeries`
  is a `Date`-indexed matrix with `NaN` for missing, not a general irregular time series
  type; merging, `aggregate` at other frequencies, and so on are out of scope.

## Tests and R fixtures

The test suite compares against R numerically, but does **not** need R to run: every reference
number is committed as CSV under `test/fixtures/`, written with `sprintf("%.17g")` so it
round-trips bit-exactly, and `test/fixtures/MANIFEST.csv` records the exact R call behind
every file.

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

```
Test Summary:                                                        | Pass  Total   Time
FXRegime.jl                                                          | 1592   1592  30.1s
  FXRatesCHF (fxregime::FXRatesCHF)                                  |   19     19   0.3s
  FXSeries container                                                 |   20     20   0.6s
  window (zoo::window)                                               |    8      8   0.3s
  fxreturns (fxregime::fxreturns)                                    |  161    161   0.2s
  fxlm (fxregime::fxlm)                                              |  165    165   1.9s
  recresid (strucchange::recresid)                                   |   21     21   1.0s
  RSS.triang (strucchange)                                           |   31     31   0.6s
  RSS2obj.triang (fxregime::gbreakpoints)                            |   18     18   0.0s
  breakpoints (fxregime::gbreakpoints / strucchange::extract.breaks) |  328    328   0.6s
  information criteria (fxregime:::gbreakpoints)                     |  200    200   0.4s
  fxregimes (fxregime::fxregimes)                                    |  220    220   0.8s
  confint.fxregimes (fxregime)                                       |  204    204   1.9s
  fxmonitor (fxregime::fxmonitor)                                    |   34     34   4.1s
  gefp / sctest (strucchange)                                        |  106    106   1.1s
  edge cases                                                         |   57     57   0.8s
     Testing FXRegime tests passed
```

A slower tier of brute-force cross-checks (direct least-squares refits of every segment of
every RSS triangle, and similar) is skipped by default:

```
FXREGIME_SLOW_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

### Regenerating the fixtures

To rebuild every reference number from R:

```
Rscript test/generate_fixtures.R            # writes into test/fixtures/
Rscript test/generate_fixtures.R /some/dir  # or elsewhere
```

This needs R with `fxregime`, `strucchange`, `zoo` and `sandwich` installed; `car` is
optional (without it the `fxpegtest` fixtures are skipped and everything else still runs).
Runtime is dominated by the two `fxregimes()` calls and is of the order of 5–15 minutes;

```
FXREGIME_FIXTURE_CACHE=/some/path.rds Rscript test/generate_fixtures.R
```

caches the fitted `fxregimes` objects between runs. The fixtures currently committed were
produced with fxregime 1.0-5, strucchange 1.6-0, zoo 1.9-0, sandwich 3.1-3, car 3.1-5 on R
4.3.1.

### Documentation

```
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

builds the Documenter site into `docs/build/`. Every exported function carries a docstring
naming its R equivalent.

## Data

`FXRatesCHF()` returns the data set shipped with the R package: **9819 daily observations of
25 currencies from 1971-01-04 to 2010-02-12, priced against the Swiss franc**. Missing values
(weekends, holidays, currencies not yet quoted) are `NaN`.

The source is the **US Federal Reserve's H.10 release**,
<https://www.federalreserve.gov/releases/h10/Hist/>, a US Government work in the public
domain. The copy here was exported from the R package's `FXRatesCHF` object to CSV
(`data/FXRatesCHF.csv`) so that R and Julia read identical bytes; the test suite checks the
two agree.

The `DUR` column deserves a note. It is a **synthetic continuous euro series**: from
1999-01-04 onwards it is exactly the `EUR` rate, and before that it is the Deutsche mark rate
converted to euro at a fixed conversion factor (`DEM / DUR` is constant at ≈ 1.96906 across
the whole pre-1999 sample). `DEM` and `EUR` are also present separately, each `NaN` outside
its own lifetime, whereas `DUR` alone is complete over the full 1971–2010 span. Studies of
pre-1999 baskets need a euro-like regressor that reaches back before the euro existed; the
INR analyses in the R package — and in this package's test suite — use `DUR` for exactly that
reason.

## Citation

If you use this package, please cite the paper the methods come from and the R package this
one ports:

> Zeileis A, Shah A, Patnaik I (2010). "Testing, Monitoring, and Dating Structural Changes in
> Exchange Rate Regimes." *Computational Statistics & Data Analysis*, **54**(6), 1696–1706.
> [doi:10.1016/j.csda.2009.12.005](https://doi.org/10.1016/j.csda.2009.12.005)

> Zeileis A, Shah A, Patnaik I (2010). *fxregime: Exchange Rate Regime Analysis.*
> R package. <https://CRAN.R-project.org/package=fxregime>

```bibtex
@Article{fxregime2010,
  author  = {Achim Zeileis and Ajay Shah and Ila Patnaik},
  title   = {Testing, Monitoring, and Dating Structural Changes in Exchange Rate Regimes},
  journal = {Computational Statistics \& Data Analysis},
  year    = {2010},
  volume  = {54},
  number  = {6},
  pages   = {1696--1706},
  doi     = {10.1016/j.csda.2009.12.005},
}

@Manual{fxregime-pkg,
  title  = {fxregime: Exchange Rate Regime Analysis},
  author = {Achim Zeileis and Ajay Shah and Ila Patnaik},
  note   = {R package},
  url    = {https://CRAN.R-project.org/package=fxregime},
}
```

### Further references

- Frankel JA, Wei S-J (1994). "Yen Bloc or Dollar Bloc? Exchange Rate Policies of the East
  Asian Economies." In *Macroeconomic Linkage*, University of Chicago Press, 295–333.
- Bai J (1997). "Estimation of a Change Point in Multiple Regression Models."
  *Review of Economics and Statistics*, **79**(4), 551–563.
  [doi:10.1162/003465397557132](https://doi.org/10.1162/003465397557132)
- Bai J, Perron P (2003). "Computation and Analysis of Multiple Structural Change Models."
  *Journal of Applied Econometrics*, **18**(1), 1–22.
  [doi:10.1002/jae.659](https://doi.org/10.1002/jae.659)
- Liu J, Wu S, Zidek JV (1997). "On Segmented Multivariate Regression."
  *Statistica Sinica*, **7**, 497–525.
- Zeileis A, Leisch F, Kleiber C, Hornik K (2005). "Monitoring Structural Change in Dynamic
  Econometric Models." *Journal of Applied Econometrics*, **20**(1), 99–121.
  [doi:10.1002/jae.776](https://doi.org/10.1002/jae.776)
- Zeileis A, Kleiber C, Krämer W, Hornik K (2003). "Testing and Dating of Structural Changes
  in Practice." *Computational Statistics & Data Analysis*, **44**(1–2), 109–123.
  [doi:10.1016/S0167-9473(03)00030-6](https://doi.org/10.1016/S0167-9473(03)00030-6)
- Shah A, Zeileis A, Patnaik I (2005). "What Is the New Chinese Currency Regime?"
  Research Report 23, Department of Statistics and Mathematics, WU Vienna.

## License

The Julia code in this repository is marked MIT; see [LICENSE](LICENSE).

**This is not settled yet.** FXRegime.jl is a
close port of **fxregime** and **strucchange**, both licensed GPL-2 | GPL-3. Several parts are
line-by-line translations of the R and LINPACK sources — `dqrdc2`, `pargmaxV`, R's `uniroot`
(`R_zeroin2`), the Zeileis et al. (2005) Table III critical values, and the index sets of the
dynamic program. That makes this a derivative work, and MIT redistribution of GPL-derived code
would be a licence violation.

The intended resolution is to obtain the agreement of the upstream copyright holders — Achim
Zeileis, Ajay Shah and Ila Patnaik — to release the ported portions under MIT. Until
that agreement is in writing, the licence stated above is provisional, and the alternative —
if permission is not forthcoming — is to relicense this package GPL-2 | GPL-3 to match
upstream. Publishing it or registering it in the Julia General registry under MIT before then
would be redistribution of GPL-derived code under an incompatible licence.

The bundled `FXRatesCHF` data is a US Government work in the public domain (Federal Reserve
H.10) and is not affected.
