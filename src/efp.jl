# Generalised empirical fluctuation process (M-fluctuation test) for an `FXLM`.
#
# Port of `strucchange::gefp(x, fit = NULL)` (strucchange/R/gefp.R), `strucchange::sctest.gefp`
# with the `maxBB` functional, and the "Brownian bridge"/"max" branch of
# `strucchange:::pvalue.efp` (strucchange/R/efp.R).

"""
    root_matrix(X::AbstractMatrix) -> Matrix{Float64}

Symmetric square root of a symmetric positive semi-definite matrix via its eigen
decomposition: `V * sqrt(D) * V'`. Port of `strucchange:::root.matrix`.

This is *not* a Cholesky factor — the result is symmetric, and `strucchange` relies on that
when it decorrelates the fluctuation process.
"""
function root_matrix(X::AbstractMatrix)
    A = Matrix{Float64}(X)
    if size(A, 1) == 1 && size(A, 2) == 1
        return reshape([sqrt(A[1, 1])], 1, 1)
    end
    F = eigen(Symmetric(A))
    any(<(0.0), F.values) && error("matrix is not positive semidefinite")
    V = F.vectors
    return V * Diagonal(sqrt.(F.values)) * transpose(V)
end

"""
    GEFP

Generalised empirical fluctuation process, the Julia analogue of an R `"gefp"` object.

* `process` — `(nobs + 1) × nreg` matrix; row 1 is the zero starting value, row `t + 1` is the
  decorrelated cumulative score after `t` observations.
* `index`   — `nobs + 1` dates; the leading entry is extrapolated one spacing before the first
  observation, as R's `gefp` does.
* `J12`     — symmetric square root of the outer-product score covariance.
"""
struct GEFP
    process::Matrix{Float64}
    index::Vector{Date}
    coefnames::Vector{String}
    nobs::Int
    nreg::Int
    J12::Matrix{Float64}
end

Base.size(p::GEFP) = size(p.process)

function Base.show(io::IO, ::MIME"text/plain", p::GEFP)
    println(io, "Generalised empirical fluctuation process (M-fluctuation test)")
    println(io, "  ", p.nobs, " observations, ", p.nreg, " parameters")
    println(io, "  ", first(p.index), " .. ", last(p.index))
    print(io, "  parameters: ", join(p.coefnames, ", "))
end

Base.show(io::IO, p::GEFP) = show(io, MIME"text/plain"(), p)

"""
    gefp(m::FXLM; fit_ = nothing) -> GEFP

Generalised empirical fluctuation process of a fitted [`fxlm`](@ref) model, i.e. the
decorrelated cumulative sum of its estimating functions. Matches
`strucchange::gefp(fxlm_object, fit = NULL)`.

With `psi = estfun(m)` and `process = psi / sqrt(n)`, the score covariance is
`J = process' * process`, its symmetric square root `J12 = root_matrix(J)`, and the returned
process is `cumsum([0; process]; dims = 1) * inv(J12)`. `fit_` exists only to mirror R's
`fit` argument and must be `nothing`.
"""
function gefp(m::FXLM; fit_ = nothing)
    fit_ === nothing || throw(ArgumentError("gefp(::FXLM) only supports `fit_ = nothing`"))

    psi = Matrix{Float64}(estfun(m))
    n = size(psi, 1)
    k = size(psi, 2)

    process = psi ./ sqrt(n)
    J = transpose(process) * process
    J12 = root_matrix(J)

    cp = Matrix{Float64}(undef, n + 1, k)
    @inbounds for j in 1:k
        s = 0.0
        cp[1, j] = 0.0
        for t in 1:n
            s += process[t, j]
            cp[t + 1, j] = s
        end
    end
    # R: t(chol2inv(chol(J12)) %*% t(process)); J12 is symmetric so this is cp * inv(J12)
    iJ12 = inv(cholesky(Symmetric(J12)))
    cp = cp * iJ12

    idx = m.index
    z0 = length(idx) >= 2 ? idx[1] - (idx[2] - idx[1]) : idx[1] - Day(1)

    return GEFP(cp, vcat(z0, idx), copy(m.coefnames), n, k, J12)
end

"""
    pvalue_maxBB(x::Real, k::Integer) -> Float64

Asymptotic p-value of the double-maximum functional of a `k`-dimensional Brownian bridge,
`1 - (1 + 2 * sum_{i=1}^{100} (-1)^i exp(-2 i^2 x^2))^k`, and `1` for `x < 0.1`.
The "Brownian bridge"/"max" branch of `strucchange:::pvalue.efp`.
"""
function pvalue_maxBB(x::Real, k::Integer)
    xf = Float64(x)
    xf < 0.1 && return 1.0
    p = 0.0
    @inbounds for i in 1:100
        p += exp(-2.0 * (i^2) * (xf^2)) * (isodd(i) ? -1.0 : 1.0)
    end
    return 1.0 - (1.0 + 2.0 * p)^k
end

"""
    sctest(p::GEFP; functional = :maxBB) -> (statistic = …, pvalue = …)

Structural change test based on a generalised empirical fluctuation process. Port of
`strucchange::sctest.gefp` with the default `maxBB` (double-maximum) functional: the statistic
is `maximum(abs, process)` and the p-value comes from [`pvalue_maxBB`](@ref) with `k` equal to
the number of columns of the process.
"""
function sctest(p::GEFP; functional::Symbol = :maxBB)
    functional === :maxBB ||
        throw(ArgumentError("only the :maxBB functional is implemented (got :$functional)"))
    proc = p.process
    # R's maxBB$computeStatistic drops the (all-zero) first row; irrelevant for max|.|
    stat = 0.0
    @inbounds for j in axes(proc, 2), t in 2:size(proc, 1)
        a = abs(proc[t, j])
        a > stat && (stat = a)
    end
    return (statistic = stat, pvalue = pvalue_maxBB(stat, size(proc, 2)))
end
