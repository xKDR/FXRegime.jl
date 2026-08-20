## fxlm and its methods — port of fxregime::fxlm, coef.fxlm, estfun.fxlm, bread.fxlm,
## sandwich::vcovHC and fxregime::fxpegtest.
## Fixtures: fxlm_*.csv, synth_break*_{fxlm_coef,fxlm_summary,estfun,bread}.csv, fxpegtest_*.csv.

"Check `coef`, the summary table, `estfun` and `bread` against the `<prefix>` fixtures."
function check_fxlm_core(prefix::AbstractString, m::FXLM;
                         coeffile = "$(prefix)_coef.csv", sumfile = "$(prefix)_summary.csv",
                         estfile = "$(prefix)_estfun.csv", breadfile = "$(prefix)_bread.csv")
    k = size(m.X, 2)

    @testset "coef.fxlm" begin
        hdr, rows = readfix(coeffile)
        ## the parameter vector is [intercept; slopes; (Variance)], variance last and MLE
        @test scol(hdr, rows, "term") == m.coefnames
        @test m.coefnames[end] == "(Variance)"
        @test length(coef(m)) == k + 1
        @test coef(m)[end] == sigma2(m)
        @test sigma2(m) ≈ mean(residuals(m) .^ 2) rtol = 1e-15
        @test reldev(coef(m), fcol(hdr, rows, "estimate")) < 1e-12
    end

    @testset "summary.lm coefficients" begin
        hdr, rows = readfix(sumfile)
        ct = coeftable(m)
        @test scol(hdr, rows, "term") == ct.names == m.coefnames[1:k]
        @test reldev(ct.estimate, fcol(hdr, rows, "estimate")) < 1e-12
        @test reldev(ct.stderror, fcol(hdr, rows, "std_error")) < 1e-12
        @test reldev(ct.tvalue, fcol(hdr, rows, "t_value")) < 1e-12
        ## p-values span ~200 orders of magnitude, so they are compared relatively; the
        ## incomplete-beta implementations of R and Distributions.jl differ in the last few ulp
        @test reldev(ct.pvalue, fcol(hdr, rows, "p_value")) < 1e-10
    end

    @testset "estfun.fxlm" begin
        hdr, rows = readfix(estfile)
        @test hdr[2:end] == m.coefnames
        @test dcol(hdr, rows, "date") == m.index
        E = estfun(m)
        @test size(E) == (nobs(m), k + 1)
        ## scores cross zero, so they are judged on the scale of the whole matrix
        @test scaledev(E, fmat(hdr, rows, hdr[2:end])) < 1e-12
        ## the defining identities (R: cbind(estfun.lm(x)/sigma2, (res^2 - sigma2)/(2*sigma2^2)))
        @test scaledev(E[:, 1:k], (residuals(m) .* m.X) ./ sigma2(m)) < 1e-14
        @test scaledev(E[:, k+1], (residuals(m) .^ 2 .- sigma2(m)) ./ (2 * sigma2(m)^2)) < 1e-14
        ## the scores of a fitted model sum to zero (exactly, in exact arithmetic)
        @test maximum(abs, sum(E; dims = 1)) / maximum(abs, E) < 1e-10
    end

    @testset "bread.fxlm" begin
        hdr, rows = readfix(breadfile)
        B = bread(m)
        @test size(B) == (k + 1, k + 1)
        ## the beta/variance cross block is the (numerically zero) mean score, ~1e-17, so the
        ## comparison is against the scale of the matrix rather than element by element
        @test scaledev(B, fmat(hdr, rows, hdr[2:end])) < 1e-12
        @test scaledev(B, transpose(B)) < 1e-14        # symmetric
    end
end

@testset "fxlm (fxregime::fxlm)" begin

    for M in MODELS
        @testset "$(M.tag)" begin
            m = M.lm
            k = size(m.X, 2)
            check_fxlm_core("fxlm_$(M.tag)", m)

            @testset "summary.lm scalars" begin
                st = kvfix("fxlm_$(M.tag)_stats.csv")
                @test nobs(m) == Int(kvnum(st, "nobs"))
                @test size(m.X, 2) == Int(kvnum(st, "nreg"))
                @test dof_residual(m) == Int(kvnum(st, "df_residual"))
                @test npar(m) == k + 1
                @test reldev(sigma(m), kvnum(st, "sigma")) < 1e-12
                @test reldev(sigma2(m), kvnum(st, "sigma2_mle")) < 1e-12
                @test reldev(r2(m), kvnum(st, "r_squared")) < 1e-12
                @test reldev(sum(abs2, residuals(m)), kvnum(st, "rss")) < 1e-12
                ## adjusted R^2 and the Gaussian log-likelihood, as printed by summary.lm/logLik
                adj = 1 - (1 - r2(m)) * (nobs(m) - 1) / dof_residual(m)
                @test reldev(adj, kvnum(st, "adj_r_squared")) < 1e-12
                ll = -0.5 * nobs(m) * (log(2pi) + log(sigma2(m)) + 1)
                @test reldev(ll, kvnum(st, "logLik")) < 1e-12
            end

            @testset "residuals.lm / fitted.lm" begin
                hdr, rows = readfix("fxlm_$(M.tag)_resid_fitted.csv")
                @test dcol(hdr, rows, "date") == m.index == M.data.index
                @test scaledev(residuals(m), fcol(hdr, rows, "residual")) < 1e-12
                @test scaledev(fitted(m), fcol(hdr, rows, "fitted")) < 1e-12
                @test scaledev(fitted(m) .+ residuals(m), m.y) < 1e-14
                @test m.y == M.data.values[:, 1]
                @test m.X[:, 1] == ones(nobs(m))
            end

            @testset "vcov.lm / sandwich::vcovHC" begin
                hdr, rows = readfix("fxlm_$(M.tag)_vcov.csv")
                V = vcov(m)
                @test size(V) == (k, k)
                @test reldev(V, fmat(hdr, rows, hdr[2:end])) < 1e-12
                @test V == vcov(m, :ols) == vcov(m; type = :ols)
                ## vcov.lm uses the df-adjusted variance
                @test scaledev(V, inv(m.X' * m.X) .* sigma(m)^2) < 1e-12

                hdr, rows = readfix("fxlm_$(M.tag)_vcovHC3.csv")
                @test scaledev(vcov(m, :HC3), fmat(hdr, rows, hdr[2:end])) < 1e-12
                @test_throws ArgumentError vcov(m, :HC2)
            end

            @testset "fxpegtest (fxregime::fxpegtest)" begin
                st = kvfix("fxpegtest_$(M.tag).csv")
                pt = fxpegtest(m)
                @test pt.peg == st["peg"]
                ## default peg = the slope with the largest |coefficient|
                @test pt.peg == m.coefnames[1 + argmax(abs.(m.beta[2:k]))]
                @test pt.df1 == Int(kvnum(st, "df")) == k - 1
                @test pt.df2 == Int(kvnum(st, "res_df_unrestricted")) == nobs(m) - k
                @test Int(kvnum(st, "res_df_restricted")) == pt.df2 + pt.df1
                @test reldev(pt.statistic, kvnum(st, "F")) < 1e-12
                @test reldev(pt.pvalue, kvnum(st, "p_value")) < 1e-10
                @test length(pt.hypothesis) == k - 1

                ## the restricted/unrestricted residual sums of squares R reports: the restricted
                ## model regresses y - x_peg on an intercept alone
                rss_u = sum(abs2, residuals(m))
                z = m.y .- m.X[:, colindex(M.data, pt.peg)]
                rss_r = sum(abs2, z .- mean(z))
                @test reldev(rss_u, kvnum(st, "rss_unrestricted")) < 1e-12
                @test reldev(rss_r, kvnum(st, "rss_restricted")) < 1e-12
                ## and the F statistic is the usual one built from them
                @test reldev(pt.statistic,
                             ((rss_r - rss_u) / pt.df1) / (rss_u / pt.df2)) < 1e-10

                ## naming the peg explicitly reproduces the default choice
                @test fxpegtest(m; peg = pt.peg).statistic == pt.statistic
                @test_throws ArgumentError fxpegtest(m; peg = "ZZZ")
            end

            @test occursin("Frankel-Wei", sprint(show, m))
        end
    end

    @testset "synthetic series" begin
        for (tag, m) in (("break1", SYNTH1_LM), ("break2", SYNTH2_LM))
            @testset "synth_$tag" begin
                check_fxlm_core("synth_$(tag)", m;
                                coeffile = "synth_$(tag)_fxlm_coef.csv",
                                sumfile = "synth_$(tag)_fxlm_summary.csv",
                                estfile = "synth_$(tag)_estfun.csv",
                                breadfile = "synth_$(tag)_bread.csv")
            end
        end
    end

    @testset "fxlm low-level constructor" begin
        m = CNY_LM
        m2 = fxlm(m.y, m.X, m.index; coefnames = m.coefnames[1:end-1])
        @test coef(m2) == coef(m)
        @test m2.coefnames == m.coefnames
        @test_throws DimensionMismatch fxlm(m.y[1:10], m.X, m.index; coefnames = m.coefnames[1:end-1])
        @test_throws DimensionMismatch fxlm(m.y, m.X, m.index[1:10]; coefnames = m.coefnames[1:end-1])
        @test_throws DimensionMismatch fxlm(m.y, m.X, m.index; coefnames = ["a"])
    end

    @testset "fstatistic (summary.lm F)" begin
        ## no fixture; check the identity that defines it
        m = CNY_LM
        fs = fstatistic(m)
        n, k = size(m.X)
        tss = sum(abs2, m.y .- mean(m.y))
        rss = sum(abs2, residuals(m))
        @test fs.df1 == k - 1
        @test fs.df2 == n - k
        @test reldev(fs.statistic, ((tss - rss) / (k - 1)) / (rss / (n - k))) < 1e-12
        @test 0 <= fs.pvalue <= 1
    end
end
