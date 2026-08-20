## Shapley-value decomposition of R^2. This is an addition of ours: it has NO equivalent in the
## R package fxregime. Ported from `calculate_r2_decomposition` in err_functions.jl.

"""
    r2_decomposition(y::AbstractVector, X::AbstractMatrix) -> Vector{Float64}

Shapley-value decomposition of the coefficient of determination of the regression of `y` on
`X` (given **without** an intercept column; an intercept is always included in every subset
regression).

Element `j` is the average marginal increase in `R²` obtained by adding regressor `j` to a
subset `S` of the other regressors, averaged over all `2^(k-1)` subsets with the Shapley
weights `γ(|S|) = |S|! (k - |S| - 1)! / k!`. The contributions sum to the `R²` of the full
regression (up to the clamping described below) and, unlike the sequential "type I" or the
partial "type III" decompositions, do not depend on the order of the regressors — which is what
makes them usable for correlated currency baskets.

Each subset `R²` is clamped to `[0, 1]`, and a subset whose design matrix is singular
contributes `0`. Returns `NaN`s when `k > 10` (the `2^k` enumeration is refused) or when
`n < k + 1`, and zeros when `y` has (numerically) no variance.

!!! note
    This function has no counterpart in the R package **fxregime**; it is an extension provided
    by FXRegime.jl.
"""
function r2_decomposition(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real})
    yy = Vector{Float64}(y)
    XX = Matrix{Float64}(X)
    n, k = size(XX)
    n == length(yy) || throw(DimensionMismatch("X has $n rows but y has $(length(yy))"))
    k > 10 && return fill(NaN, k)
    n < k + 1 && return fill(NaN, k)

    ybar = mean(yy)
    tss = 0.0
    @inbounds for i in 1:n
        d = yy[i] - ybar
        tss += d * d
    end
    tss < 1e-9 && return zeros(Float64, k)

    Xfull = Matrix{Float64}(undef, n, k + 1)
    @inbounds for i in 1:n
        Xfull[i, 1] = 1.0
    end
    @inbounds for j in 1:k, i in 1:n
        Xfull[i, j + 1] = XX[i, j]
    end

    nsub = 1 << k
    r2cache = zeros(Float64, nsub)

    Threads.@threads :dynamic for i in 1:(nsub - 1)
        cols = Vector{Int}(undef, count_ones(i) + 1)
        cols[1] = 1
        p = 2
        @inbounds for bit in 0:(k - 1)
            if (i >> bit) & 1 == 1
                cols[p] = bit + 2
                p += 1
            end
        end
        Xsub = view(Xfull, :, cols)
        r2cache[i + 1] = try
            beta = Xsub \ yy
            res = yy - Xsub * beta
            clamp(1.0 - dot(res, res) / tss, 0.0, 1.0)
        catch
            0.0
        end
    end

    gamma = Vector{Float64}(undef, k)
    fk = Float64(factorial(k))          # k <= 10, so this is exact in Float64
    @inbounds for s in 0:(k - 1)
        gamma[s + 1] = Float64(factorial(s)) * Float64(factorial(k - s - 1)) / fk
    end

    contrib = zeros(Float64, k)
    Threads.@threads :dynamic for j in 1:k
        bitj = 1 << (j - 1)
        acc = 0.0
        @inbounds for i in 0:(nsub - 1)
            if (i & bitj) == 0
                acc += gamma[count_ones(i) + 1] * (r2cache[(i | bitj) + 1] - r2cache[i + 1])
            end
        end
        contrib[j] = acc
    end
    return contrib
end

"""
    r2_decomposition(m::FXLM) -> Vector{Float64}

Shapley `R²` contribution of every basket currency of a fitted [`fxlm`](@ref) model (the
intercept column of `m.X` is dropped and re-added internally). The elements correspond to
`m.coefnames[2:end-1]`.
"""
r2_decomposition(m::FXLM) = r2_decomposition(m.y, m.X[:, 2:end])

"""
    r2_decomposition(s::FXSeries) -> Vector{Float64}

Shapley `R²` contribution of columns `2:end` of `s` to column 1, i.e. of the basket currencies
to the target currency.
"""
r2_decomposition(s::FXSeries) = r2_decomposition(s.values[:, 1], s.values[:, 2:end])

"""
    r2_decomposition(r::FXRegimes; breaks = nothing) -> Matrix{Float64}

Shapley `R²` contributions regime by regime: one row per segment of the dating `r`, one column
per basket currency (`r.data.names[2:end]`).
"""
function r2_decomposition(r::FXRegimes; breaks::Union{Nothing,Integer,Symbol} = nothing)
    fits = refit(r; breaks = breaks)
    k = size(fits[1].X, 2) - 1
    out = Matrix{Float64}(undef, length(fits), k)
    for i in eachindex(fits)
        out[i, :] = r2_decomposition(fits[i])
    end
    return out
end
