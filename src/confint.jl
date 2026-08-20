# Confidence intervals for the breakpoints of an `FXRegimes` segmentation.
#
# Port of `fxregime:::confint.fxregimes` (fxregime/R/fxregimes.R) together with the
# unexported helper `strucchange:::pargmaxV`.

"""
    pargmaxV(x; xi = 1.0, phi1 = 1.0, phi2 = 1.0) -> Float64

Distribution function of the argmax of the two-sided Brownian-motion-with-triangular-drift
limiting process of Bai (1997), evaluated at `x`.

Direct port of the (unexported) R function `strucchange:::pargmaxV`. `xi` is the ratio of the
two quadratic forms in the *Q* matrices, `phi1`/`phi2` the ratios of the *Omega* to the *Q*
quadratic forms on either side of the break.

Both branches (`G1` for `x < 0`, `G2` for `x >= 0`) are evaluated on `abs(x)` and use the
log-scale normal CDF (R's `pnorm(..., log.p = TRUE)`), which here is
`logcdf(Normal(), ...)`. Where R produces `NaN` from `log()` of a negative number this
function does the same rather than throwing.
"""
function pargmaxV(x::Real; xi::Real = 1.0, phi1::Real = 1.0, phi2::Real = 1.0)
    xf = Float64(x)
    xif = Float64(xi)
    phi = xif * (Float64(phi2) / Float64(phi1))^2
    return xf < 0.0 ? _pargmaxV_G1(abs(xf), xif, phi) : _pargmaxV_G2(abs(xf), xif, phi)
end

# R's log() returns NaN (with a warning) for negative arguments; Julia's throws.
@inline _rlog(v::Float64) = v < 0.0 ? NaN : log(v)

const _LOG2PI = log(2.0 * pi)

function _pargmaxV_G1(x::Float64, xi::Float64, phi::Float64)
    frac = xi / phi
    t1 = -exp(_rlog(x) / 2 - x / 8 - _LOG2PI / 2)
    t2 = -(phi / xi * (phi + 2 * xi) / (phi + xi)) *
         exp((frac * (1 + frac) * x / 2) + logcdf(Normal(), -(0.5 + frac) * sqrt(x)))
    t3 = exp(_rlog(x / 2 - 2 + ((phi + 2 * xi)^2) / ((phi + xi) * xi)) +
             logcdf(Normal(), -sqrt(x) / 2))
    return t1 + t2 + t3
end

function _pargmaxV_G2(x::Float64, xi::Float64, phi::Float64)
    frac = xi^2 / phi
    t1 = 1.0
    t2 = sqrt(frac) * exp(_rlog(x) / 2 - (frac * x) / 8 - _LOG2PI / 2)
    t3 = (xi / phi * (2 * phi + xi) / (phi + xi)) *
         exp(((phi + xi) * x / 2) + logcdf(Normal(), -(phi + xi / 2) / sqrt(phi) * sqrt(x)))
    t4 = -exp(_rlog(((2 * phi + xi)^2) / ((phi + xi) * phi) - 2 + frac * x / 2) +
              logcdf(Normal(), -sqrt(frac) * sqrt(x) / 2))
    return t1 + t2 + t3 + t4
end

"""
    _uniroot(f, lower, upper; tol = eps()^0.25, maxiter = 1000) -> Float64

Brent's `zeroin` root finder, ported from R's C routine `R_zeroin2` (which backs
`stats::uniroot`) so that roots agree with R's to the last bit for the same `f`.
Requires `f(lower)` and `f(upper)` to bracket a sign change.
"""
function _uniroot(f, lower::Float64, upper::Float64;
                  tol::Float64 = eps(Float64)^0.25, maxiter::Int = 1000)
    a = lower
    b = upper
    fa = f(a)
    fb = f(b)
    fa == 0.0 && return a
    fb == 0.0 && return b
    (fa * fb > 0.0) && throw(ArgumentError("f(lower) and f(upper) must have opposite signs"))

    c = a
    fc = fa
    for _ in 1:(maxiter + 1)
        prev_step = b - a
        if abs(fc) < abs(fb)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        end
        tol_act = 2 * eps(Float64) * abs(b) + tol / 2
        new_step = (c - b) / 2
        if abs(new_step) <= tol_act || fb == 0.0
            return b
        end
        if abs(prev_step) >= tol_act && abs(fa) > abs(fb)
            cb = c - b
            local p::Float64, q::Float64
            if a == c
                t1 = fb / fa
                p = cb * t1
                q = 1.0 - t1
            else
                q = fa / fc
                t1 = fb / fc
                t2 = fb / fa
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0))
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0)
            end
            if p > 0.0
                q = -q
            else
                p = -p
            end
            if p < (0.75 * cb * q - abs(tol_act * q) / 2) && p < abs(prev_step * q / 2)
                new_step = p / q
            end
        end
        if abs(new_step) < tol_act
            new_step = new_step > 0.0 ? tol_act : -tol_act
        end
        a = b
        fa = fb
        b += new_step
        fb = f(b)
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)
            c = a
            fc = fa
        end
    end
    return b
end

# Expanding bracket search, R's `while(pargmaxV(ub) < 1 - a2) ub <- ub + 1000`.
# R loops forever on a non-convergent case; we cap the search and raise instead.
const _CONFINT_MAX_BRACKET_STEPS = 10_000

function _expand_upper(target::Float64, xi::Float64, phi1::Float64, phi2::Float64)
    ub = 0.0
    steps = 0
    while pargmaxV(ub; xi = xi, phi1 = phi1, phi2 = phi2) < target
        ub += 1000.0
        steps += 1
        steps > _CONFINT_MAX_BRACKET_STEPS && error(
            "confint: upper bracket search did not converge (ub = $ub, xi = $xi, " *
            "phi1 = $phi1, phi2 = $phi2); pargmaxV never reached $target")
    end
    return ub
end

function _expand_lower(target::Float64, xi::Float64, phi1::Float64, phi2::Float64)
    lb = 0.0
    steps = 0
    while pargmaxV(lb; xi = xi, phi1 = phi1, phi2 = phi2) > target
        lb -= 1000.0
        steps += 1
        steps > _CONFINT_MAX_BRACKET_STEPS && error(
            "confint: lower bracket search did not converge (lb = $lb, xi = $xi, " *
            "phi1 = $phi1, phi2 = $phi2); pargmaxV never fell below $target")
    end
    return lb
end

@inline _quadform(delta::Vector{Float64}, mat::Matrix{Float64}) = dot(delta, mat * delta)

"""
    confint(r::FXRegimes; level = 0.95, breaks = nothing, meat = nothing)

Confidence intervals for the breakpoints of a fitted [`fxregimes`](@ref) segmentation.
Port of `confint.fxregimes` in fxregime/R/fxregimes.R.

Returns a named tuple

```
(lower, breakpoints, upper, dates = (lower = …, breakpoints = …, upper = …),
 level, nobs, npar)
```

where `lower`/`breakpoints`/`upper` are observation numbers (1-based, the last observation of
the segment ending at the break, as everywhere in this package) and `dates` holds the
corresponding `Date`s. Bounds are `missing` for a break whose interval R declines to compute.

The interval for break `i` is built from `delta = coef(seg_{i+1}) - coef(seg_i)` over the
**full** coefficient vector *including* the `(Variance)` parameter, with `Q = inv(bread(fit))`
and `Omega = Q` (or `meat(fit)` when `meat` is supplied) — both `(k+1) × (k+1)`. This is what
R does; using only the regressors and `X'X/n` gives different numbers.

`meat` may be a function mapping an `FXLM` to a `(k+1) × (k+1)` matrix.
"""
function confint(r::FXRegimes; level::Real = 0.95, breaks = nothing, meat = nothing)
    a2 = (1.0 - Float64(level)) / 2

    bp = collect(Int, first(breakpoints(r; breaks = breaks)))
    isempty(bp) && throw(ArgumentError("cannot compute confidence interval when `breaks = 0'"))
    any(b -> b <= 0, bp) &&
        throw(ArgumentError("cannot compute confidence interval when `breaks = 0'"))

    nbp = length(bp)
    upper = Vector{Union{Missing,Float64}}(undef, nbp)
    lower = Vector{Union{Missing,Float64}}(undef, nbp)

    rf = refit(r; breaks = breaks)
    cf = [Vector{Float64}(coef(f)) for f in rf]
    Q = [Matrix{Float64}(inv(bread(f))) for f in rf]
    Omega = meat === nothing ? Q : [Matrix{Float64}(meat(f)) for f in rf]

    for i in 1:nbp
        delta = cf[i + 1] .- cf[i]
        Oprod1 = _quadform(delta, Omega[i])
        Oprod2 = _quadform(delta, Omega[i + 1])
        Qprod1 = _quadform(delta, Q[i])
        Qprod2 = _quadform(delta, Q[i + 1])

        xi = Qprod2 / Qprod1
        phi1 = sqrt(Oprod1 / Qprod1)
        phi2 = sqrt(Oprod2 / Qprod2)

        p0 = pargmaxV(0.0; xi = xi, phi1 = phi1, phi2 = phi2)
        if isnan(p0) || p0 < a2 || p0 > (1 - a2)
            @warn "Confidence interval $i cannot be computed: P(argmax V <= 0) = $(round(p0, digits = 4))"
            upper[i] = missing
            lower[i] = missing
        else
            ub = _expand_upper(1 - a2, xi, phi1, phi2)
            lb = _expand_lower(a2, xi, phi1, phi2)

            fup = x -> pargmaxV(x; xi = xi, phi1 = phi1, phi2 = phi2) - (1 - a2)
            flo = x -> pargmaxV(x; xi = xi, phi1 = phi1, phi2 = phi2) - a2

            u = _uniroot(fup, 0.0, ub)
            l = _uniroot(flo, lb, 0.0)

            upper[i] = u * phi1^2 / Qprod1
            lower[i] = l * phi1^2 / Qprod1
        end
    end

    lo = Vector{Union{Missing,Int}}(undef, nbp)
    hi = Vector{Union{Missing,Int}}(undef, nbp)
    for i in 1:nbp
        hi[i] = ismissing(lower[i]) ? missing : bp[i] - Int(floor(lower[i]))
        lo[i] = ismissing(upper[i]) ? missing : bp[i] - Int(ceil(upper[i]))
    end

    idx = r.data.index
    todate = v -> Vector{Union{Missing,Date}}([ismissing(x) ? missing : idx[x] for x in v])

    return (lower = lo,
            breakpoints = bp,
            upper = hi,
            dates = (lower = todate(lo), breakpoints = todate(bp), upper = todate(hi)),
            level = Float64(level),
            nobs = r.bf.nobs,
            npar = r.bf.npar)
end
