# Port of fxregime/R/fxlm.R (fxlm, coef.fxlm, estfun.fxlm, bread.fxlm)
# and of fxregime/R/fxtools.R (fxpegtest).

"""
    FXLM

A fitted Frankel-Wei regression: the target currency's returns regressed on the basket
currencies' returns, with an intercept. Julia analogue of the `"fxlm"` object of the R package
**fxregime** (which is an `lm` carrying the series index).

Fields: `index`, `coefnames` (`["(Intercept)", basket..., "(Variance)"]`), `beta` (the `k`
regression coefficients), `sigma2` (the *MLE* variance `mean(res^2)`), the model matrix `X`
(intercept in column 1), the response `y`, `residuals` and `fitted`.
"""
struct FXLM
    index::Vector{Date}
    coefnames::Vector{String}
    beta::Vector{Float64}
    sigma2::Float64
    X::Matrix{Float64}
    y::Vector{Float64}
    residuals::Vector{Float64}
    fitted::Vector{Float64}
end

"""
    fxlm(y::Vector{Float64}, X::Matrix{Float64}, index; coefnames) -> FXLM

Low-level constructor: OLS of `y` on `X` (which must already contain the intercept column).
`coefnames` names the `k` columns of `X`; `"(Variance)"` is appended internally.
"""
function fxlm(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real},
              index::AbstractVector{Date};
              coefnames::AbstractVector{<:AbstractString})
    yy = Vector{Float64}(y)
    XX = Matrix{Float64}(X)
    n, k = size(XX)
    n == length(yy) || throw(DimensionMismatch("X has $n rows but y has $(length(yy))"))
    length(index) == n || throw(DimensionMismatch("index has $(length(index)) entries, need $n"))
    length(coefnames) == k || throw(DimensionMismatch("$(length(coefnames)) names for $k columns"))
    beta = XX \ yy
    fit = XX * beta
    res = yy .- fit
    s2 = zero(Float64)
    @inbounds for i in 1:n
        s2 += res[i]^2
    end
    s2 /= n
    return FXLM(collect(Date, index), [String.(coefnames); "(Variance)"], beta, s2,
                XX, yy, res, fit)
end

"""
    fxlm(data::FXSeries) -> FXLM

Fit the Frankel-Wei regression of column 1 of `data` (the target currency) on all remaining
columns, with an intercept. Port of `fxregime::fxlm(data = data)`, whose default formula is
`colnames(data)[1] ~ .`.
"""
function fxlm(data::FXSeries)
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
    return fxlm(y, X, data.index; coefnames = ["(Intercept)"; data.names[2:end]])
end

"""
    coef(m::FXLM) -> Vector{Float64}

The full parameter vector `[beta; sigma2]` of length `k + 1`, with the MLE variance
`mean(residuals^2)` last. Port of `coef.fxlm`, which appends `"(Variance)"` to `coef.lm`.
Names are in `m.coefnames`.
"""
coef(m::FXLM) = [m.beta; m.sigma2]

"""
    residuals(m::FXLM) -> Vector{Float64}

OLS residuals; R's `residuals.lm`.
"""
residuals(m::FXLM) = m.residuals

"""
    fitted(m::FXLM) -> Vector{Float64}

Fitted values; R's `fitted.lm`.
"""
fitted(m::FXLM) = m.fitted

"""
    nobs(m::FXLM) -> Int

Number of observations; R's `nobs`.
"""
nobs(m::FXLM) = length(m.y)

"""
    npar(m::FXLM) -> Int

Number of estimated parameters including the variance, i.e. `ncol(X) + 1` — R's `k` in
`fxregimes`.
"""
npar(m::FXLM) = size(m.X, 2) + 1

"""
    dof_residual(m::FXLM) -> Int

Residual degrees of freedom `n - ncol(X)`; R's `df.residual`.
"""
dof_residual(m::FXLM) = nobs(m) - size(m.X, 2)

"""
    sigma2(m::FXLM) -> Float64

The MLE error variance `mean(residuals^2)`, i.e. the `"(Variance)"` entry of [`coef`](@ref).
"""
sigma2(m::FXLM) = m.sigma2

"""
    index(m::FXLM) -> Vector{Date}

Time index of the observations used in the fit; R's `index.fxlm` / `time.fxlm`.
"""
index(m::FXLM) = m.index

"""
    r2(m::FXLM) -> Float64

Centred coefficient of determination `1 - rss/tss`, as reported by R's `summary.lm` for a model
with an intercept.
"""
function r2(m::FXLM)
    ybar = mean(m.y)
    rss = zero(Float64)
    tss = zero(Float64)
    @inbounds for i in eachindex(m.y)
        rss += m.residuals[i]^2
        tss += (m.y[i] - ybar)^2
    end
    return 1.0 - rss / tss
end

"""
    estfun(m::FXLM) -> Matrix{Float64}

The `n × (k+1)` matrix of empirical estimating functions (scores) of the Gaussian likelihood.
Port of `estfun.fxlm`:

```r
cbind(estfun.lm(x)/sigma2, (res^2 - sigma2)/(2 * sigma2^2))
```

with `estfun.lm(x) == residuals * X`, so column `j <= k` is `res * X[,j] / sigma2` and the last
column is `(res^2 - sigma2) / (2 sigma2^2)`.
"""
function estfun(m::FXLM)
    n, k = size(m.X)
    s2 = m.sigma2
    out = Matrix{Float64}(undef, n, k + 1)
    inv_s2 = 1.0 / s2
    inv_2s4 = 0.5 / (s2 * s2)
    @inbounds for j in 1:k, i in 1:n
        out[i, j] = m.residuals[i] * m.X[i, j] * inv_s2
    end
    @inbounds for i in 1:n
        out[i, k+1] = (m.residuals[i]^2 - s2) * inv_2s4
    end
    return out
end

"""
    bread(m::FXLM) -> Matrix{Float64}

The `(k+1) × (k+1)` "bread" matrix of the sandwich covariance, i.e. the inverse of the expected
negative Hessian scaled by `n`. Port of `bread.fxlm`:

```r
br  <- bread.lm(x)                # = n * solve(crossprod(X))
b1  <- solve(br)/sigma2           # = crossprod(X)/(n * sigma2)
b2  <- colMeans(estfun(x)[,-(k+1)]/sigma2)
b3  <- 0.5/sigma2^2
solve(rbind(cbind(b1, b2), c(b2, b3)))
```

`b2` is numerically zero for a model with an intercept, but it is computed and included here
exactly as R does.
"""
function bread(m::FXLM)
    n, k = size(m.X)
    s2 = m.sigma2
    A = Matrix{Float64}(undef, k + 1, k + 1)
    ## b1 = X'X / (n * sigma2)
    c1 = 1.0 / (n * s2)
    @inbounds for j in 1:k, l in 1:k
        acc = zero(Float64)
        for i in 1:n
            acc += m.X[i, l] * m.X[i, j]
        end
        A[l, j] = acc * c1
    end
    ## b2 = colMeans(res * X / sigma2^2)
    c2 = 1.0 / (n * s2 * s2)
    @inbounds for j in 1:k
        acc = zero(Float64)
        for i in 1:n
            acc += m.residuals[i] * m.X[i, j]
        end
        v = acc * c2
        A[j, k+1] = v
        A[k+1, j] = v
    end
    A[k+1, k+1] = 0.5 / (s2 * s2)
    return inv(A)
end

"""
    _xtxinv(X) -> Matrix{Float64}

`(X'X)^{-1}` computed from the pivoted QR factorisation of `X`, never forming the Gram matrix.

This mirrors what R does: `summary.lm` sets `cov.unscaled <- chol2inv(Qr\$qr[p1, p1])`, i.e. it
inverts the triangular factor `R` of the QR decomposition that `lm.fit` already computed, so that
`(X'X)^{-1} = R^{-1} R^{-T}`. Forming `X'X` and inverting it instead would square the condition
number of `X`. For exchange-rate return regressors that is harmless (`cond(X)` is typically < 10),
but a basket holding near-collinear currencies — DEM alongside EUR, or several currencies pegged
to the same anchor — is exactly the case this package is pointed at, and there it is not.
"""
function _xtxinv(X::AbstractMatrix{Float64})
    k = size(X, 2)
    F = qr(X, ColumnNorm())
    Rf = UpperTriangular(F.R)
    rank_def = any(i -> abs(Rf[i, i]) <= eps(Float64) * abs(Rf[1, 1]) * k, 1:k)
    rank_def && throw(SingularException(k))
    Rinv = Rf \ Matrix{Float64}(I, k, k)
    A = Rinv * Rinv'                     # = ((X[:, p])' X[:, p])^{-1}
    out = Matrix{Float64}(undef, k, k)
    p = F.p
    @inbounds for j in 1:k, i in 1:k     # undo the column pivoting
        out[p[i], p[j]] = A[i, j]
    end
    return out
end

"""
    vcov(m::FXLM, type::Symbol = :ols) -> Matrix{Float64}
    vcov(m::FXLM; type::Symbol = :ols)  -> Matrix{Float64}

Covariance matrix of the `k` regression coefficients (the variance parameter is not included).

* `:ols` — the usual `s^2 (X'X)^{-1}` with the **df-adjusted** `s^2 = rss/(n-k)`, as returned by
  R's `vcov.lm`;
* `:HC3` — the heteroskedasticity-consistent estimator
  `(X'X)^{-1} X' diag(u_i^2/(1-h_i)^2) X (X'X)^{-1}`, as returned by
  `sandwich::vcovHC(model, type = "HC3")`.
"""
vcov(m::FXLM, type::Symbol) = _vcov(m, type)
vcov(m::FXLM; type::Symbol = :ols) = _vcov(m, type)

function _vcov(m::FXLM, type::Symbol)
    n, k = size(m.X)
    XtXinv = _xtxinv(m.X)
    if type === :ols
        rss = zero(Float64)
        @inbounds for i in 1:n
            rss += m.residuals[i]^2
        end
        return Symmetric(XtXinv) * (rss / (n - k)) |> Matrix
    elseif type === :HC3
        ## leverages h_i = x_i' (X'X)^{-1} x_i
        omega = Vector{Float64}(undef, n)
        @inbounds for i in 1:n
            h = zero(Float64)
            for a in 1:k, b in 1:k
                h += m.X[i, a] * XtXinv[a, b] * m.X[i, b]
            end
            omega[i] = m.residuals[i]^2 / (1.0 - h)^2
        end
        meat = Matrix{Float64}(undef, k, k)
        @inbounds for a in 1:k, b in 1:k
            acc = zero(Float64)
            for i in 1:n
                acc += m.X[i, a] * omega[i] * m.X[i, b]
            end
            meat[a, b] = acc
        end
        return XtXinv * meat * XtXinv
    else
        throw(ArgumentError("unknown vcov type :$type (use :ols or :HC3)"))
    end
end

"""
    coeftable(m::FXLM; type = :ols)
      -> (names, estimate, stderror, tvalue, pvalue)

Coefficient table for the `k` regression coefficients: estimates, standard errors from
[`vcov`](@ref)`(m, type)`, `t` statistics and two-sided p-values from a `TDist(n - k)` — the
`summary.lm` table (`type = :HC3` gives the `coeftest`/HC3 version).
"""
function coeftable(m::FXLM; type::Symbol = :ols)
    V = vcov(m, type)
    k = size(m.X, 2)
    se = Float64[sqrt(V[j, j]) for j in 1:k]
    tv = m.beta ./ se
    d = TDist(dof_residual(m))
    pv = Float64[2 * ccdf(d, abs(t)) for t in tv]
    return (names = m.coefnames[1:k], estimate = copy(m.beta), stderror = se,
            tvalue = tv, pvalue = pv)
end

"""
    fstatistic(m::FXLM) -> (statistic, df1, df2, pvalue)

Overall F test that all slopes are zero, as printed by R's `summary.lm`.
"""
function fstatistic(m::FXLM)
    n, k = size(m.X)
    ybar = mean(m.y)
    rss = zero(Float64)
    tss = zero(Float64)
    @inbounds for i in 1:n
        rss += m.residuals[i]^2
        tss += (m.y[i] - ybar)^2
    end
    df1 = k - 1
    df2 = n - k
    stat = ((tss - rss) / df1) / (rss / df2)
    return (statistic = stat, df1 = df1, df2 = df2,
            pvalue = ccdf(FDist(df1, df2), stat))
end

"""
    sigma(m::FXLM) -> Float64

Residual standard error `sqrt(rss/(n-k))`, i.e. `summary.lm(.)$sigma` (note this is the
df-adjusted quantity, unlike [`sigma2`](@ref)).
"""
function sigma(m::FXLM)
    n, k = size(m.X)
    rss = zero(Float64)
    @inbounds for i in 1:n
        rss += m.residuals[i]^2
    end
    return sqrt(rss / (n - k))
end

function Base.show(io::IO, ::MIME"text/plain", m::FXLM)
    n, k = size(m.X)
    ct = coeftable(m)
    println(io, "Frankel-Wei regression (fxlm)")
    println(io, "  ", n, " observations, ", first(m.index), " .. ", last(m.index))
    println(io, "")
    @printf(io, "  %-14s %12s %12s %9s %10s\n", "", "Estimate", "Std. Error", "t value", "Pr(>|t|)")
    for j in 1:k
        @printf(io, "  %-14s %12.6f %12.6f %9.3f %10.3g\n",
                ct.names[j], ct.estimate[j], ct.stderror[j], ct.tvalue[j], ct.pvalue[j])
    end
    fs = fstatistic(m)
    println(io, "")
    @printf(io, "  Residual std. error: %.6g on %d degrees of freedom\n", sigma(m), dof_residual(m))
    @printf(io, "  (Variance) [MLE]:    %.6g\n", m.sigma2)
    @printf(io, "  R-squared: %.6g,  F(%d,%d) = %.6g, p = %.4g",
            r2(m), fs.df1, fs.df2, fs.statistic, fs.pvalue)
end

Base.show(io::IO, m::FXLM) = show(io, MIME"text/plain"(), m)

"""
    fxpegtest(m::FXLM; peg = nothing)
      -> (statistic, pvalue, df1, df2, peg, hypothesis)

Wald test for a pegged exchange-rate regime: the null is that the slope on `peg` equals 1 and all
other slopes are 0, with the intercept and the variance left unrestricted. Port of
`fxregime::fxpegtest`, which builds that hypothesis and hands it to `car::linearHypothesis`; the
statistic is the classical F test

```
F = (R b - q)' (R V R')^{-1} (R b - q) / q,     V = vcov.lm(m) = s^2 (X'X)^{-1}
```

on `q = k - 1` and `n - k` degrees of freedom, with the df-adjusted `s^2`. `peg` defaults to the
slope with the largest absolute coefficient, exactly as in R.
"""
function fxpegtest(m::FXLM; peg::Union{AbstractString,Nothing} = nothing)
    n, k = size(m.X)
    k >= 2 || throw(ArgumentError("model has no slopes to test"))
    slopes = view(m.beta, 2:k)
    slopenames = m.coefnames[2:k]
    pegname = peg === nothing ? slopenames[argmax(abs.(slopes))] : String(peg)
    jpeg = findfirst(==(pegname), slopenames)
    jpeg === nothing && throw(ArgumentError("no slope named \"$pegname\""))

    q = k - 1
    V = vcov(m, :ols)
    Vs = V[2:k, 2:k]                       # R V R' with R selecting the slopes
    d = Vector{Float64}(undef, q)
    @inbounds for j in 1:q
        d[j] = slopes[j] - (j == jpeg ? 1.0 : 0.0)
    end
    stat = dot(d, Symmetric(Vs) \ d) / q
    df2 = n - k
    return (statistic = stat, pvalue = ccdf(FDist(q, df2), stat), df1 = q, df2 = df2,
            peg = pegname,
            hypothesis = [string(slopenames[j], " = ", j == jpeg ? 1 : 0) for j in 1:q])
end
