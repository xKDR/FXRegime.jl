## FXRegime.jl test suite.
##
## Every reference number in this suite comes from R: the CSV fixtures in `test/fixtures/` were
## produced by `test/generate_fixtures.R` running fxregime 1.0-5 / strucchange 1.6-0 /
## sandwich 3.1-3 / car 3.1-5 on R 4.3.1, printed with `sprintf("%.17g")` so that they
## round-trip bit-exactly through the CSV.  `test/fixtures/MANIFEST.csv` records the exact R
## call behind every file.  The suite itself needs neither R nor the network.

using FXRegime
using Dates
using Test
using Statistics
using LinearAlgebra
using Random

const FIXDIR = joinpath(@__DIR__, "fixtures")

## Run the extra brute-force cross-checks that are too slow for the default suite.
const SLOW = get(ENV, "FXREGIME_SLOW_TESTS", "false") == "true"

# ---------------------------------------------------------------------------------------------
# Fixture readers.  The generator guarantees that no field contains a comma, so a plain split is
# an unambiguous parse; `NA` is R's missing value.
# ---------------------------------------------------------------------------------------------

"Read a fixture CSV; returns `(header::Vector{String}, rows::Vector{Vector{String}})`."
function readfix(name::AbstractString)
    lines = readlines(joinpath(FIXDIR, name))
    isempty(lines) && error("empty fixture $name")
    hdr = String.(split(lines[1], ','))
    rows = [String.(split(l, ',')) for l in lines[2:end]]
    isempty(rows) && error("fixture $name has no data rows")
    for (i, r) in enumerate(rows)
        length(r) == length(hdr) ||
            error("fixture $name row $i has $(length(r)) fields, header has $(length(hdr))")
    end
    return (hdr, rows)
end

"Column `name` of a fixture, as raw strings."
function scol(hdr::Vector{String}, rows::Vector{Vector{String}}, name::AbstractString)
    j = findfirst(==(String(name)), hdr)
    j === nothing && error("no column \"$name\" in fixture (have: $(join(hdr, ", ")))")
    return String[r[j] for r in rows]
end

"Column `name` as `Float64`; R's `NA` becomes `NaN`."
fcol(hdr, rows, name) = Float64[s == "NA" ? NaN : parse(Float64, s) for s in scol(hdr, rows, name)]

"Column `name` as `Int`."
icol(hdr, rows, name) = Int[parse(Int, s) for s in scol(hdr, rows, name)]

"Column `name` as `Int`, with R's `NA` becoming `missing`."
mcol(hdr, rows, name) =
    Union{Missing,Int}[s == "NA" ? missing : parse(Int, s) for s in scol(hdr, rows, name)]

"Column `name` as `Date`."
dcol(hdr, rows, name) = Date[Date(s) for s in scol(hdr, rows, name)]

"Columns `names` of a fixture, side by side, as a `Matrix{Float64}`."
fmat(hdr, rows, names) = hcat((fcol(hdr, rows, c) for c in names)...)

"A two-column `statistic,value` fixture, as a `Dict{String,String}`."
function kvfix(name::AbstractString)
    hdr, rows = readfix(name)
    return Dict(zip(scol(hdr, rows, "statistic"), scol(hdr, rows, "value")))
end

kvnum(d::Dict{String,String}, key::AbstractString) = parse(Float64, d[key])

"""
    designof(s::FXSeries) -> (X, y)

Design matrix (intercept in column 1) and response of an `FXSeries`, i.e. what R's
`model.matrix(y ~ ., data)` / `model.response(model.frame(...))` produce for the Frankel-Wei
regression of column 1 on the remaining columns.
"""
designof(s::FXSeries) = (hcat(ones(size(s, 1)), s.values[:, 2:end]), s.values[:, 1])

"""
    brute_force_rss_dev(tri, X, y) -> Float64

Worst deviation between a triangular RSS table and a *direct* least-squares fit on every one of
its segments — an independent check on the recursive-residual recursion the table is built from.

The deviation of a segment is measured against `max(RSS, 1e-12 * TSS)` rather than against `RSS`
alone: a few segments of the pegged-currency samples are fitted essentially exactly (INR tracked
USD to the last digit for weeks on end, giving an RSS of ~1e-31 against a total sum of squares of
order 1), and there a pure relative comparison measures nothing but rounding noise.
"""
function brute_force_rss_dev(tri::Vector{Vector{Float64}}, X::Matrix{Float64}, y::Vector{Float64})
    n, k = size(X)
    worst = 0.0
    for i in 1:length(tri), m in (k + 1):(n - i + 1)
        rr = i:(i + m - 1)
        Xs = X[rr, :]
        ys = y[rr]
        direct = sum(abs2, ys - Xs * (Xs \ ys))
        tss = sum(abs2, ys .- mean(ys))
        denom = max(direct, 1e-12 * tss, floatmin(Float64))
        worst = max(worst, abs(tri[i][m] - direct) / denom)
    end
    return worst
end

"Load one of the synthetic-data fixtures as an `FXSeries` (column 1 is the response)."
function loadseries(name::AbstractString)
    hdr, rows = readfix(name)
    return FXSeries(dcol(hdr, rows, hdr[1]), hdr[2:end], fmat(hdr, rows, hdr[2:end]))
end

# ---------------------------------------------------------------------------------------------
# Agreement metrics.  Both treat `NaN` on both sides as agreement, so tables carrying R's `NA`
# can be compared directly.  Tests are written as `@test <metric>(julia, r) < tol` so that a
# failure reports the deviation actually achieved.
# ---------------------------------------------------------------------------------------------

"""
    reldev(x, r) -> Float64

Worst **element-wise relative** deviation `|x - r| / |r|` between the Julia result `x` and the R
reference `r` (the absolute difference is used where `r == 0`).  Use this where every reference
element is comfortably away from zero: coefficients, standard errors, p-values (whose magnitudes
span 200 orders of magnitude and must agree relatively), residual sums of squares, information
criteria.
"""
function reldev(x, r)
    xv = vec(collect(Float64, x))
    rv = vec(collect(Float64, r))
    length(xv) == length(rv) || return Inf
    worst = 0.0
    @inbounds for i in eachindex(xv)
        a, b = xv[i], rv[i]
        (isnan(a) && isnan(b)) && continue
        d = b == 0.0 ? abs(a - b) : abs(a - b) / abs(b)
        (isnan(d) || d > worst) && (worst = isnan(d) ? Inf : d)
    end
    return worst
end

"""
    scaledev(x, r) -> Float64

Worst deviation `max|x - r| / max|r|`, i.e. measured against the magnitude of the reference
array as a whole.  Use this for arrays that legitimately pass through zero, where an
element-wise relative error carries no information: residuals and fitted values, estimating
functions, the empirical fluctuation processes (whose first row is exactly `0`) and the
beta/variance cross block of `bread` (whose entries are the sample score means, ~1e-17).
"""
function scaledev(x, r)
    xv = vec(collect(Float64, x))
    rv = vec(collect(Float64, r))
    length(xv) == length(rv) || return Inf
    scale = 0.0
    @inbounds for b in rv
        isnan(b) || (abs(b) > scale && (scale = abs(b)))
    end
    scale = max(scale, floatmin(Float64))
    worst = 0.0
    @inbounds for i in eachindex(xv)
        a, b = xv[i], rv[i]
        (isnan(a) && isnan(b)) && continue
        d = abs(a - b) / scale
        (isnan(d) || d > worst) && (worst = isnan(d) ? Inf : d)
    end
    return worst
end

# ---------------------------------------------------------------------------------------------
# Shared inputs.  These are the exact data the fixtures were generated from (see MANIFEST.csv);
# the expensive datings are computed once here and reused by several test files.
# ---------------------------------------------------------------------------------------------

const FX = FXRatesCHF()

const CNY = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], FX;
                      frequency = :daily, start = Date(2005, 7, 25), stop = Date(2009, 7, 31))
const INR = fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], FX;
                      frequency = :weekly, start = Date(1993, 4, 1), stop = Date(2008, 1, 4))
const CNY_HIST = window(CNY; stop = Date(2005, 10, 31))
const CNY120 = CNY[1:120, :]

const SYNTH1 = loadseries("synth_break1_data.csv")
const SYNTH2 = loadseries("synth_break2_data.csv")

const CNY_LM = fxlm(CNY_HIST)
const INR_LM = fxlm(INR)
const SYNTH1_LM = fxlm(SYNTH1)
const SYNTH2_LM = fxlm(SYNTH2)

const CNY_REG = fxregimes(CNY; h = 20, breaks = 10)
const INR_REG = fxregimes(INR; h = 20, breaks = 10)
const SYNTH1_REG = fxregimes(SYNTH1; h = 20, breaks = 4)
const SYNTH2_REG = fxregimes(SYNTH2; h = 25, breaks = 4)

## (tag used in fixture file names, model, series, dating, maxbreaks, h)
const MODELS = (
    (tag = "cny_hist", lm = CNY_LM, data = CNY_HIST),
    (tag = "inr_full", lm = INR_LM, data = INR),
)
const DATINGS = (
    (tag = "cny", reg = CNY_REG, data = CNY, maxbreaks = 10, h = 20),
    (tag = "inr", reg = INR_REG, data = INR, maxbreaks = 10, h = 20),
    (tag = "synth_break1", reg = SYNTH1_REG, data = SYNTH1, maxbreaks = 4, h = 20),
    (tag = "synth_break2", reg = SYNTH2_REG, data = SYNTH2, maxbreaks = 4, h = 25),
)

@testset verbose = true "FXRegime.jl" begin
    include(joinpath(@__DIR__, "test_series.jl"))
    include(joinpath(@__DIR__, "test_returns.jl"))
    include(joinpath(@__DIR__, "test_fxlm.jl"))
    include(joinpath(@__DIR__, "test_recresid.jl"))
    include(joinpath(@__DIR__, "test_breakpoints.jl"))
    include(joinpath(@__DIR__, "test_regimes.jl"))
    include(joinpath(@__DIR__, "test_confint.jl"))
    include(joinpath(@__DIR__, "test_monitor.jl"))
    include(joinpath(@__DIR__, "test_efp.jl"))
    include(joinpath(@__DIR__, "test_edge.jl"))
end
