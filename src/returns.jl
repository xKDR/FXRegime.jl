# Port of fxregime/R/fxreturns.R

const _EPOCH = Date(1970, 1, 1)          # R's day 0
const _FRIDAY_ANCHOR = Date(1970, 1, 2)  # R's as.Date(1), a Friday

"""
    nextfriday(d::Date) -> Date

The next Friday on or after `d`. This is R's

```r
nextfri <- function(date) 7 * ceiling(as.numeric(date - 1)/7) + as.Date(1)
```

used by `fxregime::fxreturns` to anchor weekly aggregation, with `as.Date(1) == 1970-01-02`
(a Friday).
"""
function nextfriday(d::Date)
    n = Dates.value(d - _EPOCH) - 1          # as.numeric(date - 1)
    return _EPOCH + Day(7 * cld(n, 7) + 1)   # cld == ceiling division, valid for negatives too
end

"""
    aggregate_weekly(s::FXSeries) -> FXSeries

Aggregate to weekly frequency anchored on Fridays, exactly like the `aggregate()` call in
`fxregime::fxreturns`: observations are grouped by [`nextfriday`](@ref) and, **independently for
each column**, the value of the group is the last non-missing observation in it, or `NaN` if the
whole group is missing in that column.
"""
function aggregate_weekly(s::FXSeries)
    n, k = size(s)
    n == 0 && return FXSeries(Date[], s.names, Matrix{Float64}(undef, 0, k))
    anchors = Date[nextfriday(d) for d in s.index]
    # index is sorted and nextfriday is monotone, so groups are contiguous
    ngroup = 1
    @inbounds for i in 2:n
        anchors[i] != anchors[i-1] && (ngroup += 1)
    end
    outidx = Vector{Date}(undef, ngroup)
    out = fill(NaN, ngroup, k)
    g = 0
    @inbounds for i in 1:n
        if i == 1 || anchors[i] != anchors[i-1]
            g += 1
            outidx[g] = anchors[i]
        end
        for j in 1:k
            v = s.values[i, j]
            isnan(v) || (out[g, j] = v)
        end
    end
    return FXSeries(outidx, s.names, out)
end

"""
    na_locf(s::FXSeries; na_rm::Bool = true) -> FXSeries

Last-observation-carried-forward, column by column, as in `zoo::na.locf`. With `na_rm = true`
(zoo's default) the result is additionally `na.trim`med with `is.na = "all"`, i.e. leading and
trailing rows in which *every* column is missing are dropped. Leading missing values in
individual columns are therefore *not* filled and can survive into the result.
"""
function na_locf(s::FXSeries; na_rm::Bool = true)
    n, k = size(s)
    vals = copy(s.values)
    @inbounds for j in 1:k
        last = NaN
        for i in 1:n
            v = vals[i, j]
            if isnan(v)
                vals[i, j] = last
            else
                last = v
            end
        end
    end
    na_rm || return FXSeries(s.index, s.names, vals)
    allna(i) = all(isnan, view(vals, i, :))
    lo = 1
    while lo <= n && allna(lo)
        lo += 1
    end
    hi = n
    while hi >= lo && allna(hi)
        hi -= 1
    end
    rows = lo:hi
    return FXSeries(s.index[rows], s.names, vals[rows, :])
end

"""
    na_omit(s::FXSeries) -> FXSeries

Drop every row that contains a missing value, as `zoo::na.omit` does for a multivariate series.
"""
function na_omit(s::FXSeries)
    keep = [!any(isnan, view(s.values, i, :)) for i in 1:size(s, 1)]
    return FXSeries(s.index[keep], s.names, s.values[keep, :])
end

"""
    logdiff(s::FXSeries) -> FXSeries

`100 * diff(log(s))`: log returns in percent. The result has one row fewer than `s` and is
indexed by the *later* date of each pair, matching `diff.zoo`.
"""
function logdiff(s::FXSeries)
    n, k = size(s)
    n <= 1 && return FXSeries(Date[], s.names, Matrix{Float64}(undef, 0, k))
    out = Matrix{Float64}(undef, n - 1, k)
    @inbounds for j in 1:k, i in 1:(n-1)
        out[i, j] = 100.0 * (log(s.values[i+1, j]) - log(s.values[i, j]))
    end
    return FXSeries(s.index[2:end], s.names, out)
end

"""
    fxreturns(target, other, data; frequency = :weekly, start = nothing, stop = nothing,
              na_action = :locf, trim = false) -> FXSeries

Extract the exchange-rate returns needed for a Frankel-Wei regression. Port of
`fxregime::fxreturns`.

`target` is the currency to be explained, `other` the basket currencies. The columns are taken
from `data` in **`data`'s own column order** (R does `which(colnames(data) %in% x)`), with the
target first — e.g. `fxreturns("INR", ["GBP", "USD"], fx)` yields columns
`["INR", "USD", "GBP"]` because `USD` precedes `GBP` in `FXRatesCHF`.

Steps, in the order R performs them:

1. column selection, `unique([target_positions; other_positions])`;
2. `window(start, stop)`;
3. if `frequency == :weekly`, aggregation to Fridays via [`aggregate_weekly`](@ref)
   (`:daily` leaves the series alone);
4. missing-value handling: `:locf` ([`na_locf`](@ref), i.e. `zoo::na.locf`), `:omit`
   ([`na_omit`](@ref)) or `:none`;
5. `100 * diff(log(.))`;
6. optional trimming of the target column: `trim = true` means `(0.01, 0.99)`, a tuple gives the
   quantile levels directly (R's type-7 quantiles, as in `Statistics.quantile`); rows whose
   target return lies outside the two quantiles are dropped.

!!! note
    R's `rval[-wi, ]` returns an *empty* series when no row is trimmed (`x[-integer(0)]` is
    empty in R). This implementation returns the untrimmed series instead, which only differs
    for degenerate levels such as `(0.0, 1.0)`.
"""
function fxreturns(target::AbstractString, other::AbstractVector{<:AbstractString},
                   data::FXSeries;
                   frequency::Symbol = :weekly,
                   start::Union{Date,Nothing} = nothing,
                   stop::Union{Date,Nothing} = nothing,
                   na_action::Symbol = :locf,
                   trim::Union{Bool,Tuple{Real,Real},AbstractVector{<:Real}} = false)
    frequency in (:weekly, :daily) ||
        throw(ArgumentError("frequency must be :weekly or :daily, got :$frequency"))
    na_action in (:locf, :omit, :none) ||
        throw(ArgumentError("na_action must be :locf, :omit or :none, got :$na_action"))

    ## 1. column selection, in data order, target block first (R: unique(c(x, other)))
    tset = Set{String}((String(target),))
    oset = Set{String}(String.(other))
    cols = Int[]
    for j in eachindex(data.names)
        data.names[j] in tset && push!(cols, j)
    end
    for j in eachindex(data.names)
        (data.names[j] in oset && !(j in cols)) && push!(cols, j)
    end
    isempty(cols) && throw(ArgumentError("no column named \"$target\" in the data"))
    rval = FXSeries(data.index, data.names[cols], data.values[:, cols])

    ## 2. window
    if !(start === nothing && stop === nothing)
        rval = window(rval; start = start, stop = stop)
    end

    ## 3. weekly aggregation
    frequency === :weekly && (rval = aggregate_weekly(rval))

    ## 4. missing values
    if na_action === :locf
        rval = na_locf(rval)
    elseif na_action === :omit
        rval = na_omit(rval)
    end

    ## 5. returns
    rval = logdiff(rval)

    ## 6. trimming
    tr = _trim_levels(trim)
    tr === nothing && return rval
    lo, hi = tr
    x = rval.values[:, 1]
    xq1 = quantile(x, lo)
    xq2 = quantile(x, hi)
    keep = [!(x[i] < xq1 || x[i] > xq2) for i in eachindex(x)]
    all(keep) && return rval
    return FXSeries(rval.index[keep], rval.names, rval.values[keep, :])
end

fxreturns(target::AbstractString, data::FXSeries; kwargs...) =
    fxreturns(target, ["USD", "JPY", "DUR", "GBP"], data; kwargs...)

function _trim_levels(trim)
    trim === false && return nothing
    trim === true && return (0.01, 0.99)
    t = collect(Float64, trim)
    length(t) == 2 || throw(ArgumentError("trim must have two elements"))
    sort!(t)
    return (t[1], t[2])
end
