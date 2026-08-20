## Port of fxregime/R/fxregimes.R (everything except confint.fxregimes, which lives in confint.jl).

"""
    FXRegimes

A dating of exchange rate regimes: the full triangular table of segment log-likelihoods
([`BreakpointsFull`](@ref)) together with the return series it was computed from.

The Julia analogue of R's `"fxregimes"` object (**fxregime**, `fxregimes.R`). Construct one with
[`fxregimes`](@ref); interrogate it with [`breakpoints`](@ref), [`breakdates`](@ref),
[`refit`](@ref), [`coef`](@ref), [`fitted`](@ref), [`residuals`](@ref) and [`confint`](@ref).
"""
struct FXRegimes
    bf::BreakpointsFull
    data::FXSeries
end

"""
    fxregimes(data::FXSeries; h = 0.15, breaks = nothing, ic = :LWZ) -> FXRegimes

Date the exchange rate regimes of the return series `data` by minimising the negative
log-likelihood of a segmented Frankel-Wei regression of column 1 (the target currency) on all
remaining columns. Port of `fxregime::fxregimes(formula, data, ...)`, whose default formula
regresses `colnames(data)[1]` on every other column.

* `h` — minimum segment size, a fraction of the sample when `h < 1` (giving `floor(n * h)`,
  as in R) or a number of observations. It is an error for `h` to exceed `floor(n / 2)`.
* `breaks` — largest number of breaks the dynamic program is run for, `ceil(n / h) - 2` by
  default (R's `ceiling(n/h) - 2`).
* `ic` — `:LWZ` (default) or `:BIC`, the criterion used to select the number of regimes.

```julia
fx = FXRatesCHF()
cny = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                frequency = :daily,
                start = Date(2005, 7, 25), stop = Date(2009, 7, 31))
reg = fxregimes(cny; h = 20, breaks = 10)
breakdates(reg)
```
"""
function fxregimes(data::FXSeries; h::Real = 0.15, breaks::Union{Nothing,Integer} = nothing,
                   ic::Symbol = :LWZ)
    n, k = size(data)
    k >= 2 || throw(ArgumentError("need at least a response and one regressor"))
    y = data.values[:, 1]
    X = Matrix{Float64}(undef, n, k)
    @inbounds for i in 1:n
        X[i, 1] = 1.0
    end
    @inbounds for j in 2:k, i in 1:n
        X[i, j] = data.values[i, j]
    end
    bf = BreakpointsFull(y, X, data.index; h = h, breaks = breaks, ic = ic)
    return FXRegimes(bf, data)
end

"""
    formula(r::FXRegimes) -> String

The model that was segmented, e.g. `"CNY ~ USD + JPY + EUR + GBP"`; R's `x\$formula`.
"""
formula(r::FXRegimes) = string(r.data.names[1], " ~ ", join(r.data.names[2:end], " + "))

index(r::FXRegimes) = r.data.index
nobs(r::FXRegimes) = r.bf.nobs
npar(r::FXRegimes) = r.bf.npar

"""
    breakpoints(r::FXRegimes; breaks = nothing) -> (breakpoints::Vector{Int}, objective::Float64)

Optimal partition of the sample into regimes and the value of the objective function (the
negative log-likelihood) it attains. R's `breakpoints(fxregimes_object, breaks = ...)`.

Each breakpoint is the number of the **last** observation of the segment ending there, as in R.
`breaks = m::Integer` returns the optimal `m`-break partition; `breaks = nothing` selects `m` by
minimising the information criterion `r.bf.ic`.
"""
breakpoints(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing) =
    breakpoints(r.bf; breaks = breaks)

"""
    breakdates(r::FXRegimes; breaks = nothing) -> Vector{Date}

Dates of the estimated breakpoints; R's `breakdates.fxregimes`.
"""
breakdates(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing) =
    breakdates(r.bf; breaks = breaks)

"""
    loglik(r::FXRegimes; breaks = nothing) -> Float64
    dof(r::FXRegimes; breaks = nothing) -> Int

Maximised log-likelihood of the segmented model and its degrees of freedom
(`npar * (m + 1) + m`); R's `logLik.gbreakpoints`.
"""
loglik(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing) =
    loglik(r.bf; breaks = breaks)

dof(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing) =
    dof(r.bf; breaks = breaks)

"""
    summarize(r::FXRegimes) -> NamedTuple

Negative log-likelihood, BIC, LWZ, breakpoints and breakdates of the optimal partition with
`m = 0, …, r.bf.maxbreaks` breaks. R's `summary.fxregimes`.
"""
summarize(r::FXRegimes) = summarize(r.bf)

## Segment boundaries [lo, hi] of every regime, given a choice of `breaks`.
## R does this with dates (`window(data, start = sbp[i], end = ebp[i])`); with a strictly
## increasing index that is exactly `bp[i]+1 : bp[i+1]`.
function _segments(r::FXRegimes, breaks)
    bp, _ = breakpoints(r; breaks = breaks)
    cuts = Int[0]
    for b in bp
        b > 0 && push!(cuts, b)
    end
    push!(cuts, r.bf.nobs)
    return [(cuts[i] + 1, cuts[i + 1]) for i in 1:(length(cuts) - 1)]
end

"""
    refit(r::FXRegimes; breaks = nothing) -> Vector{FXLM}

Refit the Frankel-Wei regression separately on every estimated regime. Port of
`refit.fxregimes`: segment `i` covers observations `bp[i]+1 : bp[i+1]` with
`bp = [0; breakpoints; nobs]`.
"""
function refit(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing)
    segs = _segments(r, breaks)
    return FXLM[fxlm(r.data[lo:hi, :]) for (lo, hi) in segs]
end

"""
    coef(r::FXRegimes; breaks = nothing) -> Matrix{Float64}

Regime-wise parameter estimates, one row per segment, columns
`[(Intercept), slopes…, (Variance)]` (the variance being the MLE `mean(residuals^2)`).
Port of `coef.fxregimes`, i.e. `t(sapply(refit(object), coef))`.
"""
function coef(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing)
    fits = refit(r; breaks = breaks)
    p = length(coef(fits[1]))
    out = Matrix{Float64}(undef, length(fits), p)
    @inbounds for i in eachindex(fits)
        c = coef(fits[i])
        for j in 1:p
            out[i, j] = c[j]
        end
    end
    return out
end

"""
    coefnames(r::FXRegimes) -> Vector{String}

Column names of [`coef(r)`](@ref): `["(Intercept)", regressors…, "(Variance)"]`.
"""
coefnames(r::FXRegimes) = ["(Intercept)"; r.data.names[2:end]; "(Variance)"]

"""
    fitted(r::FXRegimes; breaks = nothing) -> Vector{Float64}

Fitted values of the segmented model, concatenated over regimes; R's `fitted.fxregimes`.
"""
function fitted(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing)
    out = Vector{Float64}(undef, r.bf.nobs)
    pos = 1
    for f in refit(r; breaks = breaks)
        v = fitted(f)
        copyto!(out, pos, v, 1, length(v))
        pos += length(v)
    end
    return out
end

"""
    residuals(r::FXRegimes; breaks = nothing) -> Vector{Float64}

Residuals of the segmented model, concatenated over regimes; R's `residuals.fxregimes`.
"""
function residuals(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing)
    out = Vector{Float64}(undef, r.bf.nobs)
    pos = 1
    for f in refit(r; breaks = breaks)
        v = residuals(f)
        copyto!(out, pos, v, 1, length(v))
        pos += length(v)
    end
    return out
end

function Base.show(io::IO, ::MIME"text/plain", r::FXRegimes)
    bd = breakdates(r)
    println(io, "\nFX model: ", formula(r), "\n")
    println(io, "Minimum ", r.bf.ic, " partition")
    println(io, "Number of regimes: ", length(bd) + 1)
    print(io, "Breakdates: ", join(string.(bd), ", "))
end

Base.show(io::IO, r::FXRegimes) = show(io, MIME"text/plain"(), r)
