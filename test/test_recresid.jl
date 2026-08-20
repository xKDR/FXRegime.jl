## recresid — port of strucchange:::recresid_r (recresid.R).
## Fixtures: recresid_*.csv, produced by `recresid(model.matrix(f, d), model.response(...))`.

function check_recresid_fixture(name::AbstractString, s::FXSeries; rtol::Float64)
    hdr, rows = readfix(name)
    X, y = designof(s)
    n, k = size(X)
    rr = recresid(X, y)
    ## R returns one recursive residual per observation k+1 … n, and the fixture records the
    ## observation number each belongs to
    @test length(rr) == n - k
    @test icol(hdr, rows, "obs") == collect((k + 1):n)
    ## recursive residuals cross zero, so they are judged against the scale of the whole vector
    @test scaledev(rr, fcol(hdr, rows, "recresid")) < rtol
end

@testset "recresid (strucchange::recresid)" begin

    @testset "CNY daily" begin
        check_recresid_fixture("recresid_cny_daily.csv", CNY; rtol = 1e-12)
    end
    @testset "INR weekly" begin
        check_recresid_fixture("recresid_inr_weekly.csv", INR; rtol = 1e-12)
    end
    @testset "synthetic, one coefficient break" begin
        check_recresid_fixture("recresid_synth_break1.csv", SYNTH1; rtol = 1e-12)
    end
    @testset "synthetic, coefficient and variance break" begin
        check_recresid_fixture("recresid_synth_break2.csv", SYNTH2; rtol = 1e-12)
    end

    @testset "start / stop arguments" begin
        ## The recursive residual for observation t depends only on observations 1 … t, so
        ## restricting the range must just slice the full vector.  Verified against R:
        ##   recresid(X, y, start = s) == recresid(X, y)[(s - k):length(.)]
        ##   recresid(X, y, end   = e) == recresid(X, y)[1:(e - k)]
        X, y = designof(SYNTH1)
        n, k = size(X)
        full = recresid(X, y)
        ## the restart re-initialises the recursion from a fresh OLS fit, so the agreement is
        ## up to accumulated floating point error rather than exact
        @test scaledev(recresid(X, y; start = 30), full[(30 - k):end]) < 1e-12
        @test recresid(X, y; stop = 80) == full[1:(80 - k)]
        @test length(recresid(X, y; start = 30, stop = 80)) == 51
        @test_throws ArgumentError recresid(X, y; start = k)         # start must exceed ncol(X)
        @test_throws ArgumentError recresid(X, y; start = n + 1)
        @test_throws ArgumentError recresid(X, y; stop = n + 1)
    end

    @testset "exact fit and known values" begin
        ## an exactly fitting regression has zero recursive residuals
        t = range(0, 1; length = 60)
        X = hcat(ones(60), t, t .^ 2)
        y = X * [1.0, 2.0, -1.0]
        rr = recresid(X, y)
        @test length(rr) == 57
        @test maximum(abs, rr) < 1e-8    # absolute: the exact answer is 0

        ## intercept-only model: the standardised recursive residual has the closed form
        ## (y_t - mean(y_1..t)) * sqrt(t/(t-1)) — this is the fast path strucchange special-cases
        z = Float64[3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
        r1 = recresid(reshape(ones(10), 10, 1), z)
        closed = [(z[t] - mean(z[1:t])) * sqrt(t / (t - 1)) for t in 2:10]
        @test scaledev(r1, closed) < 1e-14
    end
end
