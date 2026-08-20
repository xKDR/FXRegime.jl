## Dating structural changes — port of strucchange/R/breakpoints.R (`RSSi`, `RSS.triang`,
## `extend.RSS.table`, `extract.breaks`) and of fxregime/R/gbreakpoints.R (`RSS2obj.triang`, the
## dynamic program, and the information criteria).
##
## The expensive part is the triangular table of segment objective functions.  Computing it by
## brute force costs O(n^3 k^2); strucchange instead runs one recursive-residual pass per starting
## observation and accumulates the residual sum of squares, which is O(n^2 k^2).

"""
    rss_triangle(y, X, h) -> Vector{Vector{Float64}}

Triangular table of segment residual sums of squares: `rss[i][m]` is the RSS of the OLS regression
of `y` on `X` over the rows `i:(i + m - 1)`, and is `NaN` for `m <= size(X, 2)` (R's `NA`).  The
outer vector has `n - h + 1` entries, one for each admissible segment start.

Port of the `RSSi()` / `RSS.triang` part of `breakpoints.formula` in **strucchange**, including its
closed-form fast path for the intercept-only regression.  Row `i` is obtained from
`cumsum(recresid(X[i:n, :], y[i:n]) .^ 2)`.  The rows are computed in parallel with
`Threads.@threads`; because their cost decreases with `i`, dynamic scheduling is used.
"""
function rss_triangle(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, h::Int)
    n, k = size(X)
    length(y) == n || throw(DimensionMismatch("y has $(length(y)) entries, X has $n rows"))
    h > k || throw(ArgumentError("minimum segment size must be greater than the number of regressors"))
    h <= div(n, 2) || throw(ArgumentError("minimum segment size must be smaller than half the number of observations"))

    yv = Vector{Float64}(y)
    Xv = Matrix{Float64}(X)
    ## R: intercept_only <- isTRUE(all.equal(as.vector(X), rep(1L, n)))
    intercept_only = k == 1 && _all_equal(vec(Xv), ones(n), sqrt(eps(Float64)))

    nrows = n - h + 1
    triang = Vector{Vector{Float64}}(undef, nrows)
    Threads.@threads :dynamic for i in 1:nrows
        triang[i] = _rssi(yv, Xv, i, n, k, intercept_only)
    end
    return triang
end

function _rssi(y::Vector{Float64}, X::Matrix{Float64}, i::Int, n::Int, k::Int,
               intercept_only::Bool)
    len = n - i + 1
    rval = Vector{Float64}(undef, len)
    @inbounds for j in 1:k
        rval[j] = NaN
    end
    s = 0.0
    if intercept_only
        ## (y - cumsum(y)/(1:m))[-1] * sqrt(1 + 1/(1:(m-1)))  --- R's closed form
        cs = 0.0
        @inbounds for u in 1:(len - 1)
            cs += y[i + u - 1]
            csn = cs + y[i + u]
            e = (y[i + u] - csn / (u + 1)) * sqrt(1.0 + 1.0 / u)
            s += e * e
            rval[k + u] = s
        end
    else
        e = recresid(view(X, i:n, :), view(y, i:n))
        @inbounds for u in eachindex(e)
            s += e[u] * e[u]
            rval[k + u] = s
        end
    end
    return rval
end

"""
    nll_triangle(rss_triang, k) -> Vector{Vector{Float64}}

Convert a table of segment residual sums of squares into a table of segment *negative*
log-likelihoods of the Gaussian model,
`0.5 * m * (log(rss) + 1 - log(m) + log(2π))` for a segment of length `m`, and drop the first
`h - 1` entries of every row so that entry `t` of row `i` refers to the segment `i:(i + t + h - 2)`.

Port of `RSS2obj.triang` inside `gbreakpoints` in **fxregime**.  `k` is the number of regressors
(R's `nreg`) and is used only for validation; the minimum segment size `h` is recovered from the
shape of the table exactly as R does (`h = nobs - length(rval)`, i.e. R's local `h` is `h - 1`).
"""
function nll_triangle(rss_triang::Vector{Vector{Float64}}, k::Int)
    isempty(rss_triang) && throw(ArgumentError("empty RSS triangle"))
    n = length(rss_triang[1])
    h = n - length(rss_triang) + 1
    h > k || throw(ArgumentError("minimum segment size must be greater than the number of regressors"))
    log2pi = log(2 * pi)
    out = Vector{Vector{Float64}}(undef, length(rss_triang))
    @inbounds for i in eachindex(rss_triang)
        r = rss_triang[i]
        len = length(r) - (h - 1)
        v = Vector{Float64}(undef, len)
        for t in 1:len
            m = t + h - 1
            v[t] = 0.5 * m * (log(r[m]) + 1.0 - log(m) + log2pi)
        end
        out[i] = v
    end
    return out
end

"""
    BreakpointsFull

All partitions of a sample into segments with at least `h` observations, together with the optimal
`m`-break partition for every `m = 1, …, maxbreaks`.  The Julia analogue of R's
`"gbreakpointsfull"` object (**fxregime**, `gbreakpoints.R`).

Fields:
* `index`      — the date of every observation,
* `nobs`       — number of observations `n`,
* `npar`       — number of parameters of one segment model, `size(X, 2) + 1` (variance included),
* `h`          — minimum segment size,
* `maxbreaks`  — largest number of breaks the table was built for,
* `obj_triang` — triangular table of segment negative log-likelihoods ([`nll_triangle`](@ref)),
* `obj_table`  — R's `obj.table`: column `1` holds the candidate observations `h:(n-h)`, column
  `2m` the optimal objective with `m` breaks the last of which is at that observation, and column
  `2m-1` (for `m >= 2`) the position of the preceding break,
* `breakpoints`— the partition selected by `ic`,
* `ic`         — `:LWZ` or `:BIC`.
"""
struct BreakpointsFull
    index::Vector{Date}
    nobs::Int
    npar::Int
    h::Int
    maxbreaks::Int
    obj_triang::Vector{Vector{Float64}}
    obj_table::Matrix{Float64}
    breakpoints::Vector{Int}
    ic::Symbol
end

"""
    BreakpointsFull(y, X, index; h = 0.15, breaks = nothing, ic = :LWZ)

Compute the full triangular table of segment negative log-likelihoods for the regression of `y` on
`X` and the optimal partitions derived from it.  Port of `gbreakpoints(formula, data, ...)` in
**fxregime** for its default arguments `fit = lm`, `objfun = nlogLik` (the fast, RSS-based path).

`h` is the minimum segment size, either a fraction of the sample (`h < 1`, giving `floor(n * h)`)
or a number of observations.  `breaks` is the largest number of breaks considered, defaulting to
`ceil(n / h) - 2`.
"""
function BreakpointsFull(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real},
                         index::AbstractVector{Date};
                         h::Real = 0.15, breaks::Union{Nothing,Integer} = nothing,
                         ic::Symbol = :LWZ)
    n = size(X, 1)
    length(index) == n || throw(DimensionMismatch("index and X disagree on the number of rows"))
    ic in (:LWZ, :BIC) || throw(ArgumentError("`ic` must be :LWZ or :BIC"))
    hh = h < 1 ? floor(Float64(n) * Float64(h)) : Float64(h)
    hi = round(Int, hh)
    hi <= div(n, 2) ||
        throw(ArgumentError("minimum segment size must be smaller than half the number of observations"))
    nbreaks = breaks === nothing ? cld(n, hi) - 2 : Int(breaks)
    nbreaks >= 1 || throw(ArgumentError("number of breaks must be at least 1"))

    npar = size(X, 2) + 1
    triang = nll_triangle(rss_triangle(y, X, hi), size(X, 2))
    table = _obj_table(triang, hi, n, nbreaks)

    bf = BreakpointsFull(collect(Date, index), n, npar, hi, nbreaks, triang, table, Int[], ic)
    bp, _ = breakpoints(bf)
    return BreakpointsFull(bf.index, n, npar, hi, nbreaks, triang, table, bp, ic)
end

"""
    objective(bf::BreakpointsFull, i, j) -> Float64

Negative log-likelihood of the segment running from observation `i` to observation `j`, or `NaN`
(R's `NA`) if that segment is shorter than `bf.h`.  R's `objective(i, j)` closure.
"""
objective(bf::BreakpointsFull, i::Int, j::Int) = _objective(bf.obj_triang, bf.h, bf.nobs, i, j)

@inline function _objective(triang::Vector{Vector{Float64}}, h::Int, n::Int, i::Int, j::Int)
    (j < i + h - 1 || i < 1 || j > n) && return NaN
    return @inbounds triang[i][j - i - h + 2]
end

## R's obj.table: rows are the candidate observations h:(n-h), columns come in (break, objective)
## pairs, one pair per number of breaks.  Column 1 holds the observation numbers themselves.
function _obj_table(triang::Vector{Vector{Float64}}, h::Int, n::Int, breaks::Int)
    idx = h:(n - h)
    table = fill(NaN, length(idx), 2)
    @inbounds for (r, i) in enumerate(idx)
        table[r, 1] = i
        table[r, 2] = _objective(triang, h, n, 1, i)
    end
    return _extend_obj_table(table, breaks, triang, h, n)
end

## R's extend.obj.table(): fill in the (break, objective) pairs for m = have+1, ..., breaks.
function _extend_obj_table(table::Matrix{Float64}, breaks::Int,
                           triang::Vector{Vector{Float64}}, h::Int, n::Int)
    have = div(size(table, 2), 2)
    breaks <= have && return table
    idx = h:(n - h)
    out = fill(NaN, length(idx), 2 * breaks)
    out[:, 1:(2 * have)] .= table
    @inbounds for m in (have + 1):breaks
        prevobj = 2 * (m - 1)
        for i in (m * h):(n - h)
            best = Inf
            bestj = 0
            for j in ((m - 1) * h):(i - h)
                v = out[j - h + 1, prevobj] + _objective(triang, h, n, j + 1, i)
                (isnan(v) || !(v < best)) && continue
                best = v
                bestj = j
            end
            if bestj != 0
                out[i - h + 1, 2 * m - 1] = bestj
                out[i - h + 1, 2 * m] = best
            end
        end
    end
    return out
end

## R's extract.breaks()
function _extract_breaks(table::Matrix{Float64}, breaks::Int,
                         triang::Vector{Vector{Float64}}, h::Int, n::Int)
    2 * breaks <= size(table, 2) ||
        throw(ArgumentError("compute the objective table with enough breaks before"))
    best = Inf
    opt = 0
    @inbounds for r in axes(table, 1)
        i = Int(table[r, 1])
        v = table[r, 2 * breaks] + _objective(triang, h, n, i + 1, n)
        (isnan(v) || !(v < best)) && continue
        best = v
        opt = i
    end
    opt == 0 && return Int[]
    bp = Vector{Int}(undef, breaks)
    bp[breaks] = opt
    @inbounds for m in breaks:-1:2
        opt = Int(table[opt - h + 1, 2 * m - 1])
        bp[m - 1] = opt
    end
    return bp
end

"""
    breakpoints(bf::BreakpointsFull; breaks = nothing) -> (breakpoints, objective)

Optimal partition of the sample and the corresponding value of the objective function (the
negative log-likelihood).  Each breakpoint is the number of the **last** observation of the segment
that ends there, as in R.

With `breaks = m::Integer` the optimal `m`-break partition is returned; with `breaks = nothing`
(or a `Symbol`, `:BIC` / `:LWZ`) the number of breaks is chosen by minimising the information
criterion, R's `breakpoints.gbreakpointsfull`.
"""
function breakpoints(bf::BreakpointsFull; breaks::Union{Nothing,Integer,Symbol} = nothing)
    m = if breaks === nothing || breaks isa Symbol
        crit = breaks === nothing ? bf.ic : breaks
        crit in (:LWZ, :BIC) || throw(ArgumentError("`breaks` must be :BIC or :LWZ"))
        ics = [information_criterion(bf, mm; penalty = _ic_penalty(crit, bf.nobs))
               for mm in 0:bf.maxbreaks]
        any(isnan, ics) ? bf.maxbreaks : argmin(ics) - 1
    else
        Int(breaks)
    end
    if m < 1
        return (Int[], objective(bf, 1, bf.nobs))
    end
    table = _extend_obj_table(bf.obj_table, m, bf.obj_triang, bf.h, bf.nobs)
    bp = _extract_breaks(table, m, bf.obj_triang, bf.h, bf.nobs)
    obj = 0.0
    prev = 0
    for b in vcat(bp, bf.nobs)
        obj += objective(bf, prev + 1, b)
        prev = b
    end
    return (bp, obj)
end

"""
    breakdates(bf::BreakpointsFull; breaks = nothing) -> Vector{Date}

Dates of the optimal breakpoints, R's `breakdates.gbreakpoints`.
"""
function breakdates(bf::BreakpointsFull; breaks::Union{Nothing,Integer,Symbol} = nothing)
    bp, _ = breakpoints(bf; breaks = breaks)
    return bf.index[bp]
end

"""
    loglik(bf::BreakpointsFull; breaks = nothing) -> Float64

Maximised log-likelihood of the segmented model, i.e. minus the objective function.  R's
`logLik.gbreakpoints`.
"""
function loglik(bf::BreakpointsFull; breaks::Union{Nothing,Integer,Symbol} = nothing)
    _, obj = breakpoints(bf; breaks = breaks)
    return -obj
end

"""
    dof(bf::BreakpointsFull; breaks = nothing) -> Int

Degrees of freedom of the segmented model with `m` breaks, `npar * (m + 1) + m`; the `df`
attribute of R's `logLik.gbreakpoints`.
"""
function dof(bf::BreakpointsFull; breaks::Union{Nothing,Integer,Symbol} = nothing)
    bp, _ = breakpoints(bf; breaks = breaks)
    return _seg_dof(bf.npar, length(bp))
end

_seg_dof(npar::Int, m::Int) = npar * (m + 1) + m

_ic_penalty(ic::Symbol, n::Int) = ic === :BIC ? log(n) : 0.299 * log(n)^2.1

"""
    information_criterion(bf::BreakpointsFull, m; penalty = 2.0) -> Float64

Information criterion of the optimal `m`-break partition, `2 * nll + penalty * df` with
`df = npar * (m + 1) + m`.  R's `AIC(breakpoints(object, breaks = m), k = penalty)`.
"""
function information_criterion(bf::BreakpointsFull, m::Integer; penalty::Float64 = 2.0)
    bp, obj = breakpoints(bf; breaks = Int(m))
    return 2.0 * obj + penalty * _seg_dof(bf.npar, length(bp))
end

"""
    bic(bf::BreakpointsFull, m) -> Float64

Schwarz criterion `2 * nll + df * log(n)` of the optimal `m`-break partition.
"""
bic(bf::BreakpointsFull, m::Integer) =
    information_criterion(bf, m; penalty = _ic_penalty(:BIC, bf.nobs))

"""
    lwz(bf::BreakpointsFull, m) -> Float64

Liu-Wu-Zidek criterion `2 * nll + df * 0.299 * log(n)^2.1` of the optimal `m`-break partition.
"""
lwz(bf::BreakpointsFull, m::Integer) =
    information_criterion(bf, m; penalty = _ic_penalty(:LWZ, bf.nobs))

"""
    summarize(bf::BreakpointsFull) -> NamedTuple

Table of the optimal partitions with `m = 0, …, bf.maxbreaks` breaks: the negative log-likelihood,
BIC and LWZ of each, and the breakpoints and breakdates themselves.  R's
`summary.gbreakpointsfull`.
"""
function summarize(bf::BreakpointsFull)
    ms = collect(0:bf.maxbreaks)
    nll = Vector{Float64}(undef, length(ms))
    bicv = Vector{Float64}(undef, length(ms))
    lwzv = Vector{Float64}(undef, length(ms))
    bps = Vector{Vector{Int}}(undef, length(ms))
    bds = Vector{Vector{Date}}(undef, length(ms))
    pb = _ic_penalty(:BIC, bf.nobs)
    pl = _ic_penalty(:LWZ, bf.nobs)
    for (t, m) in enumerate(ms)
        bp, obj = breakpoints(bf; breaks = m)
        df = _seg_dof(bf.npar, length(bp))
        nll[t] = obj
        bicv[t] = 2.0 * obj + pb * df
        lwzv[t] = 2.0 * obj + pl * df
        bps[t] = bp
        bds[t] = bf.index[bp]
    end
    return (m = ms, nll = nll, BIC = bicv, LWZ = lwzv, breakpoints = bps, breakdates = bds)
end

function Base.show(io::IO, ::MIME"text/plain", bf::BreakpointsFull)
    println(io, "BreakpointsFull: ", bf.nobs, " observations, h = ", bf.h,
            ", up to ", bf.maxbreaks, " breaks")
    println(io, "Minimum ", bf.ic, " partition: ", length(bf.breakpoints) + 1, " segments")
    println(io, "Breakpoints at observation number: ", join(bf.breakpoints, " "))
    print(io, "Corresponding to breakdates: ", join(string.(bf.index[bf.breakpoints]), " "))
end

Base.show(io::IO, bf::BreakpointsFull) = show(io, MIME"text/plain"(), bf)
