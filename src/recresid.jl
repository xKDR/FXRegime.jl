## Recursive residuals — port of strucchange/R/recresid.R (`recresid_r`).
##
## The recursion is the standard Brown/Durbin/Evans updating formula, but it is checked against a
## full QR refit at every step until the two agree to `tol` (after that the cheap recursion is
## trusted).  To reproduce R's behaviour bit-for-bit we also reproduce R's `lm.fit`, i.e. the
## LINPACK routine `dqrdc2` with its "limited column pivoting" rank detection, `chol2inv` for
## `(X'X)^{-1}`, and `all.equal.numeric` for the agreement test.

# ---------------------------------------------------------------------------------------------
# R's all.equal.numeric(target, current, tolerance = tol) for plain Float64 vectors.
# countEQ = FALSE, scale = NULL: positions with target == current (or NaN target) are dropped,
# the scale is the mean |target| over the remaining positions and the reported error is the mean
# absolute difference divided by that scale (absolute difference if the scale is <= tol).
# ---------------------------------------------------------------------------------------------
function _all_equal(target::AbstractVector{Float64}, current::AbstractVector{Float64},
                    tol::Float64)
    length(target) == length(current) || return false
    # 'is.NA' value mismatch => not equal
    @inbounds for i in eachindex(target)
        isnan(target[i]) == isnan(current[i]) || return false
    end
    n = 0
    sabs = 0.0
    @inbounds for i in eachindex(target)
        (isnan(target[i]) || target[i] == current[i]) && continue
        n += 1
        sabs += abs(target[i])
    end
    n == 0 && return true                       # all(out) => TRUE
    scale = sabs / n
    if !(isfinite(scale) && scale > tol)
        scale = 1.0
    end
    xy = 0.0
    @inbounds for i in eachindex(target)
        (isnan(target[i]) || target[i] == current[i]) && continue
        xy += abs(target[i] - current[i]) / (n * scale)
    end
    return !isnan(xy) && !(xy > tol)
end

@inline function _nrm2(A::AbstractMatrix{Float64}, i0::Int, i1::Int, j::Int)
    s = 0.0
    @inbounds @simd for i in i0:i1
        s += A[i, j] * A[i, j]
    end
    return sqrt(s)
end

"""
    _dqrdc2!(A, tol, pivot, qraux, work1, work2) -> rank

In-place QR decomposition of `A` with the limited column pivoting of LINPACK's `dqrdc2`, the
routine behind R's `qr(..., LAPACK = FALSE)` and hence behind `lm.fit`.  A column is moved to the
right-hand edge (and excluded from the rank) as soon as its reduced norm has fallen below `tol`
times its original norm.  The Householder vectors are left in the lower triangle of `A` and in
`qraux` in LINPACK's packed form; the upper triangle of `A` holds R.
"""
function _dqrdc2!(A::AbstractMatrix{Float64}, tol::Float64, pivot::Vector{Int},
                  qraux::Vector{Float64}, work1::Vector{Float64}, work2::Vector{Float64})
    n, p = size(A)
    @inbounds for j in 1:p
        pivot[j] = j
        nj = _nrm2(A, 1, n, j)
        qraux[j] = nj
        work1[j] = nj
        work2[j] = nj == 0.0 ? 1.0 : nj
    end
    lup = min(n, p)
    kk = p + 1
    @inbounds for l in 1:lup
        # cycle negligible columns to the right-hand edge
        while l < kk && qraux[l] < work2[l] * tol
            for i in 1:n
                t = A[i, l]
                for j in (l + 1):p
                    A[i, j - 1] = A[i, j]
                end
                A[i, p] = t
            end
            ipv = pivot[l]; t1 = qraux[l]; t2 = work1[l]; t3 = work2[l]
            for j in (l + 1):p
                pivot[j - 1] = pivot[j]
                qraux[j - 1] = qraux[j]
                work1[j - 1] = work1[j]
                work2[j - 1] = work2[j]
            end
            pivot[p] = ipv; qraux[p] = t1; work1[p] = t2; work2[p] = t3
            kk -= 1
        end
        l == n && continue
        nrmxl = _nrm2(A, l, n, l)
        nrmxl == 0.0 && continue
        if A[l, l] != 0.0
            nrmxl = copysign(nrmxl, A[l, l])
        end
        inv_nrmxl = 1.0 / nrmxl
        @simd for i in l:n
            A[i, l] *= inv_nrmxl
        end
        A[l, l] += 1.0
        alpha = A[l, l]
        for j in (l + 1):p
            t = 0.0
            @simd for i in l:n
                t += A[i, l] * A[i, j]
            end
            t = -t / alpha
            @simd for i in l:n
                A[i, j] += t * A[i, l]
            end
            if qraux[j] != 0.0
                tt = 1.0 - (abs(A[l, j]) / qraux[j])^2
                tt = max(tt, 0.0)
                if tt < 1.0e-6
                    ## re-compute the norm exactly after a large reduction; `work2` keeps the
                    ## *original* norm, which is what the negligibility test above compares to
                    qraux[j] = _nrm2(A, l + 1, n, j)
                    work1[j] = qraux[j]
                else
                    qraux[j] = qraux[j] * sqrt(tt)
                end
            end
        end
        qraux[l] = A[l, l]
        A[l, l] = -nrmxl
    end
    return min(kk - 1, n)
end

"""
    _lmfit!(coefv, xinv, A, yv, tol, pivot, qraux, work1, work2) -> (rank, hasna)

R's `lm.fit(A, yv, tol = tol)`, restricted to what `recresid` needs.  `A` and `yv` are overwritten.
`coefv` receives the coefficients in the original column order with `NaN` in the aliased positions
(R's `NA`); `xinv` receives `Xinv0()`, i.e. the zero-padded `chol2inv` of the leading `rank` block
of R, which equals `(X'X)^{-1}` when `A` has full rank.
"""
function _lmfit!(coefv::Vector{Float64}, xinv::Matrix{Float64}, A::AbstractMatrix{Float64},
                 yv::AbstractVector{Float64}, tol::Float64, pivot::Vector{Int},
                 qraux::Vector{Float64}, work1::Vector{Float64}, work2::Vector{Float64})
    n, p = size(A)
    rank = _dqrdc2!(A, tol, pivot, qraux, work1, work2)

    # apply Q' to y (dqrsl), using the first `rank` Householder transformations
    ju = min(rank, n - 1)
    @inbounds for l in 1:ju
        vl = qraux[l]
        vl == 0.0 && continue
        saved = A[l, l]
        A[l, l] = vl
        t = 0.0
        @simd for i in l:n
            t += A[i, l] * yv[i]
        end
        t = -t / vl
        @simd for i in l:n
            yv[i] += t * A[i, l]
        end
        A[l, l] = saved
    end

    # back-solve R[1:rank, 1:rank] b = (Q'y)[1:rank]
    fill!(coefv, NaN)
    @inbounds for i in rank:-1:1
        s = yv[i]
        for j in (i + 1):rank
            s -= A[i, j] * coefv[pivot[j]]
        end
        coefv[pivot[i]] = s / A[i, i]
    end

    # chol2inv of the leading rank block, scattered back to the original column order
    fill!(xinv, 0.0)
    if rank > 0
        # Rinv = inv(UpperTriangular(R[1:rank, 1:rank])), stored in work of size rank^2
        Rinv = Matrix{Float64}(undef, rank, rank)
        fill!(Rinv, 0.0)
        @inbounds for j in 1:rank
            Rinv[j, j] = 1.0 / A[j, j]
            for i in (j - 1):-1:1
                s = 0.0
                for l in (i + 1):j
                    s += A[i, l] * Rinv[l, j]
                end
                Rinv[i, j] = -s / A[i, i]
            end
        end
        @inbounds for b in 1:rank, a in 1:rank
            s = 0.0
            for l in max(a, b):rank
                s += Rinv[a, l] * Rinv[b, l]
            end
            xinv[pivot[a], pivot[b]] = s
        end
    end
    hasna = rank < p
    return rank, hasna
end

"""
    recresid(X, y; start = size(X, 2) + 1, stop = size(X, 1), tol = sqrt(eps()) / size(X, 2))

Standardised recursive residuals of the linear regression of `y` on `X`, computed for observations
`start:stop`; the result has length `stop - start + 1`.

Port of `recresid_r` in **strucchange** (`recresid.R`), i.e. of `strucchange::recresid(x, y)` with
`engine = "R"`.  The recursion is initialised with an OLS fit on the first `start - 1` rows and
then updated by the rank-one formula, while a full QR refit is recomputed at every step until the
two agree to `tol` in the sense of R's `all.equal` (from then on only the recursion is used).
Rank-deficient fits are handled as in R: aliased coefficients are treated as zero and the
corresponding rows/columns of `(X'X)^{-1}` are zero.
"""
function recresid(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real};
                  start::Int = size(X, 2) + 1, stop::Int = size(X, 1),
                  tol::Float64 = sqrt(eps(Float64)) / size(X, 2),
                  qr_tol::Float64 = 1.0e-7)
    n = stop
    q = start - 1
    k = size(X, 2)
    start > k || throw(ArgumentError("`start` must be greater than ncol(X)"))
    start <= size(X, 1) || throw(ArgumentError("`start` must be at most nrow(X)"))
    (stop >= start && stop <= size(X, 1)) || throw(ArgumentError("invalid `stop`"))

    rval = zeros(Float64, n - q)

    # work space
    Abuf = Matrix{Float64}(undef, n, k)
    ybuf = Vector{Float64}(undef, n)
    pivot = Vector{Int}(undef, k)
    qraux = Vector{Float64}(undef, k)
    work1 = Vector{Float64}(undef, k)
    work2 = Vector{Float64}(undef, k)
    betar = Vector{Float64}(undef, k)
    coefv = Vector{Float64}(undef, k)
    X1 = Matrix{Float64}(undef, k, k)
    xr = Vector{Float64}(undef, k)
    v = Vector{Float64}(undef, k)
    w = Vector{Float64}(undef, k)

    @inline function refit!(r::Int)
        @inbounds for j in 1:k, i in 1:r
            Abuf[i, j] = Float64(X[i, j])
        end
        @inbounds for i in 1:r
            ybuf[i] = Float64(y[i])
        end
        _lmfit!(coefv, X1, view(Abuf, 1:r, 1:k), view(ybuf, 1:r), qr_tol,
                pivot, qraux, work1, work2)
    end

    ## initialise recursion: lm.fit on rows 1:q
    _, fm_hasna = refit!(q)
    @inbounds for j in 1:k
        betar[j] = isnan(coefv[j]) ? 0.0 : coefv[j]   # coef0()
    end
    @inbounds for j in 1:k
        xr[j] = Float64(X[q + 1, j])
    end
    _rr_symv!(v, X1, xr)
    fr = 1.0 + _rr_dot(xr, v)
    rval[1] = (Float64(y[q + 1]) - _rr_dot(xr, betar)) / sqrt(fr)

    check = true

    if (q + 1) < n
        @inbounds for r in (q + 2):n
            nona = !fm_hasna

            ## recursion: X1 <- X1 - X1 xr xr' X1 / fr   (v == X1 * xr)
            invfr = 1.0 / fr
            for b in 1:k
                vb = v[b] * invfr
                @simd for a in 1:k
                    X1[a, b] -= v[a] * vb
                end
            end
            ## betar <- betar + X1 xr * rval[r-q-1] * sqrt(fr)
            _rr_symv!(w, X1, xr)
            c = rval[r - q - 1] * sqrt(fr)
            @simd for a in 1:k
                betar[a] += c * w[a]
            end

            ## full QR decomposition
            if check
                _, hasna = refit!(r - 1)
                nona = nona && !any(isnan, betar) && !hasna
                if nona && _all_equal(coefv, betar, tol)
                    check = false
                end
                fm_hasna = hasna
                ## X1 has already been overwritten with Xinv0(fm) by refit!
                for j in 1:k
                    betar[j] = isnan(coefv[j]) ? 0.0 : coefv[j]
                end
            end

            ## residual
            for j in 1:k
                xr[j] = Float64(X[r, j])
            end
            _rr_symv!(v, X1, xr)
            fr = 1.0 + _rr_dot(xr, v)
            s = 0.0
            for j in 1:k
                t = xr[j] * betar[j]
                isnan(t) || (s += t)          # sum(..., na.rm = TRUE)
            end
            rval[r - q] = (Float64(y[r]) - s) / sqrt(fr)
        end
    end
    return rval
end

@inline function _rr_dot(a::Vector{Float64}, b::Vector{Float64})
    s = 0.0
    @inbounds @simd for i in eachindex(a)
        s += a[i] * b[i]
    end
    return s
end

@inline function _rr_symv!(out::Vector{Float64}, A::Matrix{Float64}, x::Vector{Float64})
    k = length(x)
    @inbounds for i in 1:k
        out[i] = 0.0
    end
    @inbounds for j in 1:k
        xj = x[j]
        xj == 0.0 && continue
        @simd for i in 1:k
            out[i] += A[i, j] * xj
        end
    end
    return out
end
