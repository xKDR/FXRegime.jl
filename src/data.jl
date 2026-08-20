const _FXRATESCHF_CACHE = Ref{Union{FXSeries,Nothing}}(nothing)

"""
    FXRatesCHF() -> FXSeries

Daily exchange rates of 25 currencies against the Swiss franc (CHF) as the unit currency,
1971-01-04 to 2010-02-12 (9819 observations). Missing values are `NaN`.

This is the `FXRatesCHF` data set of the R package **fxregime**, derived from the historical
exchange rate releases of the US Federal Reserve
(<https://www.federalreserve.gov/releases/h10/Hist/>), which are in the public domain.

`DUR` is the "German mark / euro" series: DEM returns before the introduction of the euro in
1999, adjusted to EUR rates, and EUR thereafter — it is the series used as the euro-zone
regressor in the published analyses.

```julia
fx = FXRatesCHF()
inr = fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], fx;
                frequency = :weekly, start = Date(1993, 4, 1), stop = Date(2008, 1, 4))
```
"""
function FXRatesCHF()
    cached = _FXRATESCHF_CACHE[]
    cached === nothing || return cached
    path = joinpath(dirname(@__DIR__), "data", "FXRatesCHF.csv")
    isfile(path) || error("FXRatesCHF.csv not found at $path")
    raw, header = readdlm(path, ','; header = true, quotes = true)
    colnames = String.(vec(header)[2:end])
    dates = Date.(String.(view(raw, :, 1)))
    values = Matrix{Float64}(undef, size(raw, 1), length(colnames))
    @inbounds for j in eachindex(colnames), i in axes(raw, 1)
        v = raw[i, j + 1]
        values[i, j] = v isa Real ? Float64(v) : NaN
    end
    series = FXSeries(dates, colnames, values)
    _FXRATESCHF_CACHE[] = series
    return series
end
