## The dating machinery — ports of strucchange's RSS.triang / extend.RSS.table / extract.breaks
## and of fxregime's gbreakpoints (RSS2obj.triang, the dynamic program, the information criteria).
## Fixtures: rss_triangle_*.csv (R's `breakpoints(f, data, h)$RSS.triang`),
##           fxregimes_*_breakpoints.csv, fxregimes_*_ic.csv, fxregimes_*_selection.csv.

"Check a full RSS triangle against R's `RSS.triang` and its `_meta` companion."
function check_triangle(tag::AbstractString, s::FXSeries, h::Int; rtol::Float64)
    hdr, rows = readfix("rss_triangle_$(tag).csv")
    meta = kvfix("rss_triangle_$(tag)_meta.csv")
    X, y = designof(s)
    n, k = size(X)
    @test n == Int(kvnum(meta, "nobs"))
    @test k == Int(kvnum(meta, "nreg"))
    @test h == Int(kvnum(meta, "h"))

    tri = rss_triangle(y, X, h)
    ## R computes one row per admissible segment start, 1 … n-h+1, and row i runs over the
    ## segment lengths 1 … n-i+1
    @test length(tri) == n - h + 1 == Int(kvnum(meta, "ntriang"))
    @test all(i -> length(tri[i]) == n - i + 1, 1:length(tri))

    st = icol(hdr, rows, "start")
    en = icol(hdr, rows, "end")
    ns = icol(hdr, rows, "nseg")
    @test ns == en .- st .+ 1
    got = Float64[tri[st[i]][ns[i]] for i in eachindex(st)]
    ref = fcol(hdr, rows, "rss")
    ## R's NA (a segment with no residual degrees of freedom) must be NaN in exactly the same
    ## places, which is what `reldev` requires of the two arrays
    @test isnan.(got) == isnan.(ref)
    @test all(i -> isnan(got[i]) == (ns[i] <= k), eachindex(got))
    ## the RSS is accumulated from the recursive residuals in both languages, but the two
    ## recursions are not evaluated in the same order, so agreement is float-error limited
    @test scaledev(got, ref) < rtol
end

@testset "RSS.triang (strucchange)" begin
    @testset "CNY, first 120 observations, h = 20" begin
        check_triangle("cny120_h20", CNY120, 20; rtol = 1e-11)
    end
    @testset "synthetic break1, h = 20" begin
        check_triangle("synth_break1_h20", SYNTH1, 20; rtol = 1e-11)
    end
    @testset "synthetic break2, h = 25" begin
        check_triangle("synth_break2_h25", SYNTH2, 25; rtol = 1e-11)
    end

    @testset "against direct OLS fits" begin
        ## Independent check of the recursion: every entry of the triangle must be the residual
        ## sum of squares of a plain least-squares fit on that segment.  Tolerance 1e-8 rather
        ## than the 1e-12 used against the R fixtures: this compares two *different* algorithms
        ## (an accumulation of squared recursive residuals against a direct QR fit), and the
        ## recursion loses several digits on long, badly conditioned segments.
        X, y = designof(CNY120)
        tri = rss_triangle(y, X, 20)
        @test brute_force_rss_dev(tri, X, y) < 1e-8
    end

    @testset "argument checks" begin
        X, y = designof(CNY120)
        @test_throws ArgumentError rss_triangle(y, X, 5)     # h must exceed ncol(X)
        @test_throws ArgumentError rss_triangle(y, X, 61)    # h must not exceed n/2
        @test_throws DimensionMismatch rss_triangle(y[1:10], X, 20)
    end
end

@testset "RSS2obj.triang (fxregime::gbreakpoints)" begin
    X, y = designof(SYNTH1)
    n, k = size(X)
    h = 20
    tri = rss_triangle(y, X, h)
    nll = nll_triangle(tri, k)

    ## RSS2obj.triang drops the first h-1 entries of every row, so entry t of row i refers to
    ## the segment i … i+t+h-2
    @test length(nll) == length(tri)
    @test all(i -> length(nll[i]) == length(tri[i]) - (h - 1), eachindex(nll))

    ## the Gaussian negative log-likelihood of the segment
    for i in (1, 5, 40, length(nll)), t in (1, min(3, length(nll[i])), length(nll[i]))
        m = t + h - 1
        @test reldev(nll[i][t], 0.5 * m * (log(tri[i][m]) + 1 - log(m) + log(2pi))) < 1e-14
    end

    ## and it agrees with the log-likelihood of an fxlm fitted on that segment
    for (i, j) in ((1, 40), (11, 60), (25, 120))
        seg = fxlm(SYNTH1[i:j, :])
        m = j - i + 1
        @test reldev(nll[i][m - h + 1],
                     0.5 * m * (log(2pi) + log(sigma2(seg)) + 1)) < 1e-12
    end

    @test_throws ArgumentError nll_triangle(tri, h)          # h must exceed nreg
end

@testset "breakpoints (fxregime::gbreakpoints / strucchange::extract.breaks)" begin
    for D in DATINGS
        @testset "$(D.tag)" begin
            bf = D.reg.bf
            n = size(D.data, 1)
            k = size(D.data, 2)
            @test bf.nobs == n
            @test bf.npar == k + 1                # regressors (incl. intercept) + the variance
            @test bf.h == D.h
            @test bf.maxbreaks == D.maxbreaks
            @test bf.index == D.data.index

            hdr, rows = readfix("fxregimes_$(D.tag)_breakpoints.csv")
            ms = icol(hdr, rows, "m")
            bps = icol(hdr, rows, "breakpoint")
            bds = dcol(hdr, rows, "breakdate")

            for m in 1:D.maxbreaks
                sel = ms .== m
                bp, obj = breakpoints(bf; breaks = m)
                ## breakpoint indices must match R exactly, they are integers
                @test bp == bps[sel]
                @test length(bp) == m
                ## breakdates must line up with R's to the day
                @test breakdates(bf; breaks = m) == bds[sel]
                @test breakdates(bf; breaks = m) == bf.index[bp]

                ## index convention: breakpoint i is the LAST observation of segment i, so
                ## every segment [bp[i-1]+1, bp[i]] holds at least h observations
                @test issorted(bp) && allunique(bp)
                @test bp[1] >= bf.h
                @test n - bp[end] >= bf.h
                @test all(diff(bp) .>= bf.h)

                ## the objective is the sum of the segment negative log-likelihoods
                cuts = vcat(0, bp, n)
                total = sum(FXRegime.objective(bf, cuts[i] + 1, cuts[i+1])
                            for i in 1:(length(cuts) - 1))
                @test reldev(obj, total) < 1e-12
                @test reldev(loglik(bf; breaks = m), -obj) < 1e-15
                @test dof(bf; breaks = m) == bf.npar * (m + 1) + m
            end
        end
    end
end

@testset "information criteria (fxregime:::gbreakpoints)" begin
    for D in DATINGS
        @testset "$(D.tag)" begin
            bf = D.reg.bf
            hdr, rows = readfix("fxregimes_$(D.tag)_ic.csv")
            ms = icol(hdr, rows, "m")
            @test ms == collect(0:D.maxbreaks)

            nll = Float64[breakpoints(bf; breaks = m)[2] for m in ms]
            @test reldev(nll, fcol(hdr, rows, "nll")) < 1e-12
            @test Int[dof(bf; breaks = m) for m in ms] == icol(hdr, rows, "df")
            @test reldev([bic(bf, m) for m in ms], fcol(hdr, rows, "bic")) < 1e-12
            @test reldev([lwz(bf, m) for m in ms], fcol(hdr, rows, "lwz")) < 1e-12

            ## BIC = 2*nll + df*log(n), LWZ = 2*nll + df*0.299*log(n)^2.1
            n = bf.nobs
            for m in ms
                df = bf.npar * (m + 1) + m
                @test reldev(bic(bf, m), 2 * nll[m+1] + df * log(n)) < 1e-13
                @test reldev(lwz(bf, m), 2 * nll[m+1] + df * 0.299 * log(n)^2.1) < 1e-13
                @test reldev(information_criterion(bf, m; penalty = 2.0),
                             2 * nll[m+1] + 2 * df) < 1e-13
            end

            ## model selection, R's breakpoints(<gbreakpointsfull>)
            sel = kvfix("fxregimes_$(D.tag)_selection.csv")
            @test bf.nobs == Int(kvnum(sel, "nobs"))
            @test bf.npar == Int(kvnum(sel, "npar"))
            @test bf.h == Int(kvnum(sel, "h"))
            @test length(breakpoints(bf; breaks = :BIC)[1]) == Int(kvnum(sel, "m_selected_BIC"))
            @test length(breakpoints(bf; breaks = :LWZ)[1]) == Int(kvnum(sel, "m_selected_LWZ"))
            bp, obj = breakpoints(bf)
            @test length(bp) == Int(kvnum(sel, "m_default"))
            @test reldev(obj, kvnum(sel, "objective_default")) < 1e-12
            ## the default criterion is LWZ, and it is the argmin of the LWZ column
            @test bf.ic === :LWZ
            @test bp == breakpoints(bf; breaks = :LWZ)[1]
            @test length(bp) == argmin([lwz(bf, m) for m in ms]) - 1
            @test bf.breakpoints == bp

            @test_throws ArgumentError breakpoints(bf; breaks = :AIC)

            ## summarize() reproduces the whole table
            sm = summarize(bf)
            @test sm.m == ms
            @test reldev(sm.nll, nll) < 1e-15
            @test reldev(sm.BIC, fcol(hdr, rows, "bic")) < 1e-12
            @test reldev(sm.LWZ, fcol(hdr, rows, "lwz")) < 1e-12
            @test sm.breakpoints[1] == Int[]
            @test sm.breakpoints[end] == breakpoints(bf; breaks = D.maxbreaks)[1]
            @test sm.breakdates[end] == bf.index[sm.breakpoints[end]]

            @test occursin("Breakpoints at observation number", sprint(show, bf))
        end
    end

    @testset "zero breaks" begin
        bf = SYNTH1_REG.bf
        bp, obj = breakpoints(bf; breaks = 0)
        @test bp == Int[]
        @test dof(bf; breaks = 0) == bf.npar
        ## with no break the objective is the negative log-likelihood of the whole sample
        m = fxlm(SYNTH1)
        @test reldev(obj, 0.5 * nobs(m) * (log(2pi) + log(sigma2(m)) + 1)) < 1e-12
        hdr, rows = readfix("fxregimes_synth_break1_ic.csv")
        @test reldev(obj, fcol(hdr, rows, "nll")[1]) < 1e-12
    end
end
