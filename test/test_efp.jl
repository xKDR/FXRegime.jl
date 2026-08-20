## Generalised empirical fluctuation process and the M-fluctuation test — ports of
## strucchange::gefp(x, fit = NULL) and strucchange::sctest.gefp with the maxBB functional.
## Fixtures: gefp_*_{process,J12}.csv, sctest_*.csv.

function check_gefp(tag::AbstractString, m::FXLM; rtol::Float64)
    g = gefp(m)
    k = length(m.coefnames)

    hdr, rows = readfix("gefp_$(tag)_process.csv")
    ## R's gefp process has n+1 rows, indexed t = 0 … n, and its first row is exactly zero
    @test icol(hdr, rows, "t") == collect(0:nobs(m))
    @test size(g.process) == (nobs(m) + 1, k)
    @test hdr[3:end] == g.coefnames == m.coefnames
    @test all(g.process[1, :] .== 0.0)
    ## the leading index entry is extrapolated one spacing before the first observation
    @test dcol(hdr, rows, "date") == g.index
    @test g.index[2:end] == m.index
    @test g.index[1] == m.index[1] - (m.index[2] - m.index[1])
    @test scaledev(g.process, fmat(hdr, rows, hdr[3:end])) < rtol

    hdr, rows = readfix("gefp_$(tag)_J12.csv")
    @test hdr[2:end] == m.coefnames
    @test scaledev(g.J12, fmat(hdr, rows, hdr[2:end])) < rtol
    ## J12 is the *symmetric* square root of J = crossprod(estfun)/n (not a Cholesky factor)
    @test scaledev(g.J12, transpose(g.J12)) < 1e-13
    psi = estfun(m) ./ sqrt(nobs(m))
    @test scaledev(g.J12 * g.J12, transpose(psi) * psi) < 1e-12
    @test g.nobs == nobs(m)
    @test g.nreg == k
    @test size(g) == size(g.process)

    st = kvfix("sctest_$(tag).csv")
    s = sctest(g)
    @test Int(kvnum(st, "nobs")) == nobs(m)
    @test Int(kvnum(st, "nreg")) == k
    @test reldev(s.statistic, kvnum(st, "f_efp")) < 1e-12
    ## p-values of the double-maximum functional: relative agreement across their whole range
    @test reldev(s.pvalue, kvnum(st, "p_value")) < 1e-10
    ## the statistic is max|process| (the all-zero first row cannot be the maximum)
    @test s.statistic == maximum(abs, g.process)
    @test 0 <= s.pvalue <= 1
    return g
end

@testset "gefp / sctest (strucchange)" begin
    @testset "cny_hist" begin
        g = check_gefp("cny_hist", CNY_LM; rtol = 1e-12)
        @test occursin("fluctuation process", sprint(show, g))
    end
    @testset "inr_full" begin
        check_gefp("inr_full", INR_LM; rtol = 1e-12)
    end
    @testset "synth_break1" begin
        check_gefp("synth_break1", SYNTH1_LM; rtol = 1e-12)
    end
    @testset "synth_break2" begin
        check_gefp("synth_break2", SYNTH2_LM; rtol = 1e-12)
    end

    @testset "root.matrix (strucchange:::root.matrix)" begin
        ## the symmetric square root via the eigen decomposition, not a Cholesky factor
        A = [4.0 1.0; 1.0 3.0]
        R = FXRegime.root_matrix(A)
        @test scaledev(R * R, A) < 1e-14
        @test scaledev(R, transpose(R)) < 1e-15
        @test FXRegime.root_matrix(reshape([9.0], 1, 1)) == reshape([3.0], 1, 1)
        @test_throws ErrorException FXRegime.root_matrix([1.0 0.0; 0.0 -1.0])
    end

    @testset "pvalue.efp, maxBB branch" begin
        ## 1 - (1 + 2*sum_{i=1}^{100} (-1)^i exp(-2 i^2 x^2))^k, and 1 below 0.1.
        ## Reference values from R: strucchange::maxBB$computePval(x, k)
        @test FXRegime.pvalue_maxBB(0.05, 5) == 1.0
        @test FXRegime.pvalue_maxBB(0.0999, 2) == 1.0
        @test FXRegime.pvalue_maxBB(0.0, 1) == 1.0
        for (k, x, p) in ((1, 0.5, 0.96394524366487511), (1, 1.0, 0.2699996716773545),
                          (1, 1.358, 0.050026797334447037), (1, 2.0, 0.00067092525577971962),
                          (1, 3.0, 3.0459959443618345e-08), (3, 1.0, 0.61098247511035053),
                          (3, 2.0, 0.0020114256472534287), (6, 1.0, 0.84866536532873105),
                          (6, 1.358, 0.2650325122013365), (6, 3.0, 1.8275974278392226e-07))
            @test reldev(FXRegime.pvalue_maxBB(x, k), p) < 1e-10
        end
        ## monotone decreasing in x, increasing in k
        xs = 0.2:0.1:3.0
        @test issorted([FXRegime.pvalue_maxBB(x, 6) for x in xs]; rev = true)
        @test FXRegime.pvalue_maxBB(1.5, 6) > FXRegime.pvalue_maxBB(1.5, 1)
    end

    @testset "argument checks" begin
        @test_throws ArgumentError gefp(CNY_LM; fit_ = identity)
        @test_throws ArgumentError sctest(gefp(CNY_LM); functional = :maxL2BB)
    end
end
