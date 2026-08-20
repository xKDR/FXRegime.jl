## Edge cases and boundary behaviour.  Where R's behaviour was checked with Rscript while
## writing these tests, the R call is quoted in a comment; where FXRegime.jl deliberately
## departs from R, the departure is named.

@testset "edge cases" begin

    @testset "h at its limits" begin
        s = CNY[1:60, :]                       # n = 60, 4 regressors + intercept
        n, k = size(s, 1), size(s, 2)

        ## h may go all the way up to floor(n/2): with n = 60 and h = 30 the only admissible
        ## one-break partition is the one that splits the sample in half
        reg = fxregimes(s; h = 30, breaks = 1)
        @test reg.bf.h == 30
        @test breakpoints(reg; breaks = 1)[1] == [30]
        @test [nobs(f) for f in refit(reg; breaks = 1)] == [30, 30]

        ## one observation more is refused, as in R's breakpoints()
        X, y = designof(s)
        @test_throws ArgumentError fxregimes(s; h = 31, breaks = 1)
        @test_throws ArgumentError rss_triangle(y, X, 31)

        ## and h must leave at least one residual degree of freedom per segment: h > ncol(X)
        @test_throws ArgumentError rss_triangle(y, X, k)          # h == ncol(X)
        tri = rss_triangle(y, X, k + 1)                            # the smallest legal h
        @test length(tri) == n - k
        @test all(isnan, tri[1][1:k])
        @test isfinite(tri[1][k+1])

        ## n < 2h is exactly the condition that is refused
        @test_throws ArgumentError fxregimes(CNY[1:39, :]; h = 20, breaks = 1)
        @test fxregimes(CNY[1:40, :]; h = 20, breaks = 1).bf.h == 20
    end

    @testset "zero breaks" begin
        bf = SYNTH2_REG.bf
        bp, obj = breakpoints(bf; breaks = 0)
        @test bp == Int[]
        @test breakdates(bf; breaks = 0) == Date[]
        @test dof(bf; breaks = 0) == bf.npar

        ## the m = 0 objective is the negative log-likelihood of the unsegmented fxlm
        m0 = fxlm(SYNTH2)
        @test reldev(obj, 0.5 * nobs(m0) * (log(2pi) + log(sigma2(m0)) + 1)) < 1e-12
        @test reldev(loglik(bf; breaks = 0), -obj) < 1e-15

        ## the segmented accessors degenerate to the single full-sample fit
        rf = refit(SYNTH2_REG; breaks = 0)
        @test length(rf) == 1
        @test nobs(rf[1]) == size(SYNTH2, 1)
        @test coef(rf[1]) == coef(m0)
        @test size(coef(SYNTH2_REG; breaks = 0)) == (1, size(SYNTH2, 2) + 1)
        @test scaledev(fitted(SYNTH2_REG; breaks = 0), fitted(m0)) < 1e-15

        ## but a confidence interval for "no break" is meaningless, as in R
        @test_throws ArgumentError confint(SYNTH2_REG; breaks = 0)
        ## and the dynamic program refuses to be built for zero breaks
        @test_throws ArgumentError fxregimes(SYNTH1; h = 20, breaks = 0)
    end

    @testset "rank-deficient design" begin
        rng = MersenneTwister(11)
        n = 40
        x = randn(rng, n)
        y = 1 .+ 2 .* x .+ 0.3 .* randn(rng, n)
        X2 = hcat(ones(n), x)
        X3 = hcat(ones(n), x, x)               # the third column is perfectly collinear

        ## R treats the aliased coefficient as NA and drops it from the prediction, so the
        ## recursion simply starts one observation later.  Checked against R:
        ##   recresid(cbind(1, x, x), y) == recresid(cbind(1, x), y)[-1]   (max diff 1.3e-15)
        r2v = recresid(X2, y)
        r3v = recresid(X3, y)
        @test length(r3v) == n - 3
        @test length(r2v) == n - 2
        @test scaledev(r3v, r2v[2:end]) < 1e-12
        @test all(isfinite, r3v)

        ## a rank-deficient RSS triangle stays finite and monotone (the RSS of a nested
        ## sequence of segments can only grow)
        tri = rss_triangle(y, X3, 10)
        @test all(isnan, tri[1][1:3])
        @test all(isfinite, tri[1][4:end])
        @test issorted(tri[1][4:end])

        ## a segment on which a regressor is constant is rank-deficient too; the dating must
        ## still run and produce a legal partition
        z = copy(x); z[1:15] .= 0.5
        w = 1 .+ 0.5 .* z .+ 0.1 .* randn(rng, n)
        s = FXSeries([Date(2000, 1, 1) + Day(7 * (i - 1)) for i in 1:n], ["y", "x"], hcat(w, z))
        reg = fxregimes(s; h = 10, breaks = 2)
        bp, obj = breakpoints(reg; breaks = 1)
        @test length(bp) == 1
        @test 10 <= bp[1] <= n - 10
        @test isfinite(obj)
    end

    @testset "all-missing column" begin
        ## EUR does not exist before 1999, so it is missing throughout this window.  R:
        ##   fxreturns("INR", other = c("USD","EUR"), frequency = "weekly",
        ##             start = "1993-04-01", end = "1995-01-01")
        ## gives a 91 x 3 series whose EUR column is NA everywhere.
        s = fxreturns("INR", ["USD", "EUR"], FX; frequency = :weekly,
                      start = Date(1993, 4, 1), stop = Date(1995, 1, 1))
        @test size(s) == (91, 3)
        @test s.names == ["INR", "USD", "EUR"]
        j = colindex(s, "EUR")
        @test count(isnan, s.values[:, j]) == 91
        @test count(isnan, s.values[:, 1:2]) == 0
        @test s.values[1, 1] ≈ 0.97726135807612025 rtol = 1e-12   # unaffected by the NA column

        ## a model that includes the missing column is all-NaN rather than silently dropping
        ## rows: FXRegime.jl has no equivalent of R's model-frame `na.action`
        m = fxlm(s)
        @test all(isnan, coef(m))
    end

    @testset "single-column input" begin
        ## R: fxreturns("INR", other = character(0), ...) returns a one-column series
        s = fxreturns("INR", String[], FX; frequency = :weekly,
                      start = Date(1993, 4, 1), stop = Date(1995, 1, 1))
        @test size(s) == (91, 1)
        @test s.names == ["INR"]
        ## but there is nothing to regress on
        @test_throws ArgumentError fxlm(s)
        @test_throws ArgumentError fxregimes(s; h = 20, breaks = 1)

        ## an intercept-only model is still a legal FXLM through the low-level constructor,
        ## and hits the closed-form path of the RSS triangle
        y = s.values[1:60, 1]
        X = reshape(ones(60), 60, 1)
        m = fxlm(y, X, s.index[1:60]; coefnames = ["(Intercept)"])
        @test coef(m)[1] ≈ mean(y) rtol = 1e-12
        @test length(coef(m)) == 2
        tri = rss_triangle(y, X, 20)
        @test isnan(tri[1][1])
        @test reldev(tri[1][60], sum(abs2, y .- mean(y))) < 1e-10
    end

    @testset "NaN handling" begin
        ## na_action = :none keeps the missing values.  R:
        ##   fxreturns("INR", other = c("USD","EUR"), start = "1998-01-01", end = "1999-06-01",
        ##             na.action = function(z) z)   ->  74 x 3, 53 NAs
        s = fxreturns("INR", ["USD", "EUR"], FX; frequency = :weekly,
                      start = Date(1998, 1, 1), stop = Date(1999, 6, 1), na_action = :none)
        @test size(s) == (74, 3)
        @test count(isnan, s.values) == 53

        ## na.omit drops every row with a missing value
        so = fxreturns("INR", ["USD", "EUR"], FX; frequency = :weekly,
                       start = Date(1998, 1, 1), stop = Date(1999, 6, 1), na_action = :omit)
        @test count(isnan, so.values) == 0
        @test size(so, 1) < size(s, 1)

        ## NaN propagates through the estimators rather than being dropped
        v = copy(CNY120.values); v[3, 2] = NaN
        mn = fxlm(FXSeries(CNY120.index, CNY120.names, v))
        @test all(isnan, coef(mn))
        @test all(isnan, recresid(designof(FXSeries(CNY120.index, CNY120.names, v))...))
    end

    @testset "trim edge cases" begin
        ## FXRegime.jl deviation from R, documented in the `fxreturns` docstring: R's
        ## `rval[-wi, ]` with `wi = integer(0)` yields an *empty* series, so R returns 0 rows
        ## for trim = c(0, 1); we return the untrimmed series instead.
        full = fxreturns("CNY", ["USD"], FX; frequency = :daily,
                         start = Date(2005, 7, 25), stop = Date(2006, 7, 31))
        t01 = fxreturns("CNY", ["USD"], FX; frequency = :daily,
                        start = Date(2005, 7, 25), stop = Date(2006, 7, 31), trim = (0.0, 1.0))
        @test size(t01) == size(full)

        ## a symmetric trim removes observations from both tails
        t10 = fxreturns("CNY", ["USD"], FX; frequency = :daily,
                        start = Date(2005, 7, 25), stop = Date(2006, 7, 31), trim = (0.1, 0.9))
        @test size(t10, 1) < size(full, 1)
        x = full.values[:, 1]
        @test t10.index == full.index[.!((x .< quantile(x, 0.1)) .| (x .> quantile(x, 0.9)))]
        @test_throws ArgumentError fxreturns("CNY", ["USD"], FX; trim = [0.1, 0.5, 0.9])
    end

    ## An independent brute-force check of the whole RSS triangle on a full-size sample: every
    ## entry must be the residual sum of squares of a plain OLS fit on that segment.  This is
    ## O(n^2) least-squares fits, far too slow for the default suite, so it only runs when
    ## FXREGIME_SLOW_TESTS=true.  The fast suite covers the same code path on the 120-observation
    ## CNY sample in test_breakpoints.jl.
    if SLOW
        @testset "brute-force RSS triangle, full INR sample (slow)" begin
            X, y = designof(INR)
            tri = rss_triangle(y, X, 20)
            ## 1e-8, not 1e-12: recursive accumulation against a direct QR fit, see the same
            ## check on the CNY sample in test_breakpoints.jl
            @test brute_force_rss_dev(tri, X, y) < 1e-8
        end
    end

end
