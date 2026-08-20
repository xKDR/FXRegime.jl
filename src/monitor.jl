# Monitoring an exchange rate regime — port of fxregime/R/fxmonitor.R
# (`fxmonitor`, `breakpoints.fxmonitor`, `breakdates.fxmonitor`, `print.fxmonitor`).

# Table III of Zeileis, Leisch, Kleiber & Hornik (2005), Journal of Applied Econometrics,
# "Monitoring structural change in dynamic econometric models". Reproduced verbatim from
# fxregime/R/fxmonitor.R: R fills the matrix column-wise with nrow = 7, so row `i` is
# monitoring horizon `_FXMONITOR_CV_ROW[i]` and column `j` is level `_FXMONITOR_CV_COL[j]`.
const _FXMONITOR_CV = reshape([
    1.159, 1.329, 1.430, 1.472, 1.502, 1.541, 1.567,
    1.253, 1.445, 1.544, 1.589, 1.619, 1.668, 1.688,
    1.383, 1.590, 1.695, 1.753, 1.789, 1.838, 1.860,
    1.467, 1.688, 1.793, 1.861, 1.899, 1.961, 1.964,
    1.568, 1.814, 1.939, 2.006, 2.046, 2.090, 2.128,
    1.616, 1.896, 2.022, 2.076, 2.131, 2.159, 2.219,
    1.680, 1.997, 2.103, 2.177, 2.226, 2.257, 2.311,
    1.801, 2.114, 2.217, 2.301, 2.397, 2.380, 2.454,
    1.976, 2.300, 2.423, 2.525, 2.573, 2.597, 2.650,
    2.118, 2.478, 2.599, 2.712, 2.812, 2.766, 2.888,
    2.435, 2.789, 2.973, 3.288, 3.226, 3.230, 3.401], 7, 11)

const _FXMONITOR_CV_ROW = Float64[2, 3, 4, 5, 6, 8, 10]
const _FXMONITOR_CV_COL = Float64[0.2, 0.15, 0.1, 0.075, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]

"""
    approx_rule2(x, y, xout) -> Float64

Linear interpolation of `(x, y)` at `xout` with constant extrapolation outside the range of
`x`; R's `approx(x, y, xout, rule = 2)\$y`. Like R's `approx`, the knots are sorted by `x`
first, so `x` may be given in decreasing order.
"""
function approx_rule2(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xout::Real)
    ord = sortperm(collect(Float64, x))
    xs = Float64.(x)[ord]
    ys = Float64.(y)[ord]
    v = Float64(xout)
    v <= xs[1] && return ys[1]
    v >= xs[end] && return ys[end]
    j = searchsortedlast(xs, v)
    j >= length(xs) && return ys[end]
    xs[j] == v && return ys[j]
    w = (v - xs[j]) / (xs[j + 1] - xs[j])
    return ys[j] + w * (ys[j + 1] - ys[j])
end

"""
    FXMonitor

Result of [`fxmonitor`](@ref); the Julia analogue of an R `"fxmonitor"` object.

* `process`   — `N × (k+1)` decorrelated cumulative score process over the *full* sample
  (history plus monitoring period), one column per model parameter.
* `index`     — the `N` dates of `process`.
* `n`         — number of historical observations (rows `1:n` of `process`).
* `monitor`   — last date of the history period.
* `critval`   — Bonferroni-corrected critical value; the boundary at observation `t > n` is
  `critval * t / n`.
"""
struct FXMonitor
    process::Matrix{Float64}
    index::Vector{Date}
    n::Int
    monitor::Date
    critval::Float64
    coefnames::Vector{String}
    J12::Matrix{Float64}
end

"""
    fxmonitor(data::FXSeries; start::Date, stop = 3, alpha = 0.05, meat = nothing) -> FXMonitor

Monitor an exchange rate regime for structural change. Port of `fxmonitor`
(fxregime/R/fxmonitor.R). Column 1 of `data` is the target currency, the remaining columns the
basket currencies.

`start` is the first date of the *monitoring* period; the history period is everything before
it, on which an [`fxlm`](@ref) is fitted. `stop` is R's `end` argument — the monitoring
horizon as a multiple of the history length, used to pick the critical value from Table III of
Zeileis et al. (2005, JAE). `alpha` is Bonferroni corrected across the parameters,
`alpha = 1 - (1 - alpha)^(1/ncol)`, and the two-way lookup into the table is R's
`approx(..., rule = 2)` first over the horizon, then over the level.

`meat` (R's `meat.`) may be a function of the historical `FXLM` returning the score
covariance; by default the outer product of the historical scores is used.
"""
function fxmonitor(data::FXSeries; start::Date, stop::Real = 3, alpha::Real = 0.05,
                   meat = nothing)
    hist = window(data; stop = start - Day(1))
    size(hist, 1) > 0 || throw(ArgumentError("no observations before `start`"))
    monitor = last(hist.index)
    n = size(hist, 1)

    histfm = fxlm(hist)
    sigma2 = histfm.sigma2
    beta = histfm.beta

    # full sample design matrix and response
    N, ncols = size(data)
    k = ncols - 1                          # number of regressors excluding the intercept
    y = data.values[:, 1]
    X = Matrix{Float64}(undef, N, k + 1)
    @inbounds for i in 1:N
        X[i, 1] = 1.0
        for j in 1:k
            X[i, j + 1] = data.values[i, j + 1]
        end
    end

    # estimating functions evaluated at the historical parameter estimates
    p = k + 2                              # regressors + intercept + (Variance)
    sq = sqrt(n)
    mscore = Matrix{Float64}(undef, N, p)
    @inbounds for i in 1:N
        r = y[i]
        for j in 1:(k + 1)
            r -= X[i, j] * beta[j]
        end
        for j in 1:(k + 1)
            mscore[i, j] = r * X[i, j] / sq
        end
        mscore[i, p] = (r * r - sigma2) / sq
    end

    hist_score = view(mscore, 1:n, :)
    J12 = root_matrix(meat === nothing ? transpose(hist_score) * hist_score :
                                         Matrix{Float64}(meat(histfm)))

    proc = Matrix{Float64}(undef, N, p)
    @inbounds for j in 1:p
        s = 0.0
        for i in 1:N
            s += mscore[i, j]
            proc[i, j] = s
        end
    end
    proc = proc * inv(cholesky(Symmetric(J12)))   # R: t(chol2inv(chol(J12)) %*% t(mscore))

    # Bonferroni-corrected critical value from Table III
    a = 1.0 - (1.0 - Float64(alpha))^(1 / p)
    critrow = [approx_rule2(_FXMONITOR_CV_ROW, view(_FXMONITOR_CV, :, j), stop)
               for j in axes(_FXMONITOR_CV, 2)]
    critval = approx_rule2(_FXMONITOR_CV_COL, critrow, a)

    return FXMonitor(proc, copy(data.index), n, monitor, critval,
                     copy(histfm.coefnames), J12)
end

"""
    breakpoints(m::FXMonitor) -> Union{Int,Nothing}

Observation number of the first crossing of the monitoring boundary, or `nothing` if the
process never crosses. Port of `breakpoints.fxmonitor`: the boundary at observation `t > n` is
`critval * t / n` and the crossing is detected on `maximum(abs, process[t, :])`.
"""
function breakpoints(m::FXMonitor)
    N = size(m.process, 1)
    @inbounds for t in (m.n + 1):N
        ma = 0.0
        for j in axes(m.process, 2)
            a = abs(m.process[t, j])
            a > ma && (ma = a)
        end
        ma > (t / m.n) * m.critval && return t
    end
    return nothing
end

"""
    breakdates(m::FXMonitor) -> Union{Date,Nothing}

Date of the first boundary crossing, or `nothing`. Port of `breakdates.fxmonitor`.
"""
function breakdates(m::FXMonitor)
    bp = breakpoints(m)
    return bp === nothing ? nothing : m.index[bp]
end

function Base.show(io::IO, ::MIME"text/plain", m::FXMonitor)
    bd = breakdates(m)
    println(io, "Monitoring of FX model")
    println(io)
    println(io, "History period: ", first(m.index), " to ", m.monitor)
    print(io, "Break detected: ", bd === nothing ? "none" : string(bd))
end

Base.show(io::IO, m::FXMonitor) = show(io, MIME"text/plain"(), m)
