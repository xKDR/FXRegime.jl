## FXSeries container and the FXRatesCHF data set.
## R equivalents: the `zoo` series machinery used throughout fxregime, and
## `data("FXRatesCHF", package = "fxregime")`.

@testset "FXRatesCHF (fxregime::FXRatesCHF)" begin
    ## Reference values from R:
    ##   data("FXRatesCHF", package = "fxregime"); dim(z); colnames(z); coredata(z)[i, j]
    @test size(FX) == (9819, 25)
    @test size(FX, 1) == 9819
    @test size(FX, 2) == 25
    @test length(FX) == 9819
    @test nobs(FX) == 9819

    @test FX.index[1] == Date(1971, 1, 4)
    @test FX.index[end] == Date(2010, 2, 12)
    @test issorted(FX.index)
    @test allunique(FX.index)

    @test names(FX) == ["USD", "JPY", "DUR", "EUR", "DEM", "GBP", "AUD", "BRL", "CAD", "CNY",
                        "DKK", "HKD", "INR", "MYR", "MXN", "NOK", "NZD", "KRW", "SEK", "SGD",
                        "LKR", "ZAR", "TWD", "THB", "VEB"]

    ## exact bit patterns, R printed them with %.17g
    @test FX.values[1, 1] === 0.23158869847151461
    @test FX[:, "USD"][1] === 0.23158869847151461
    @test FX[:, "CNY"][9819] === 6.3359606824925816
    @test FX[:, "JPY"][5000] === 104.88073683381877

    ## missing values are NaN, and there are exactly as many as R has NAs
    @test isnan(FX.values[1, colindex(FX, "EUR")])
    @test count(isnan, FX.values) == 45303
    ## EUR starts on 1999-01-04 (observation 7022), as in R
    @test findfirst(!isnan, FX[:, "EUR"]) == 7022
    @test FX.index[7022] == Date(1999, 1, 4)

    ## the loader caches, so repeated calls are free and identical
    @test FXRatesCHF() === FX
end

@testset "FXSeries container" begin
    idx = [Date(2020, 1, i) for i in 1:5]
    v = Float64[1 2; 3 4; 5 6; 7 8; 9 10]
    s = FXSeries(idx, ["A", "B"], v)

    @test size(s) == (5, 2)
    @test index(s) == idx
    @test names(s) == ["A", "B"]
    @test colindex(s, "B") == 2
    @test_throws ArgumentError colindex(s, "Z")
    @test s[:, "A"] == Float64[1, 3, 5, 7, 9]
    @test s[:, 2] == Float64[2, 4, 6, 8, 10]
    @test first(s) == Date(2020, 1, 1)
    @test last(s) == Date(2020, 1, 5)

    sub = s[2:4, :]
    @test sub isa FXSeries
    @test index(sub) == idx[2:4]
    @test sub.values == v[2:4, :]

    g = getcols(s, ["B", "A"])
    @test names(g) == ["B", "A"]
    @test g.values == v[:, [2, 1]]

    ## constructor validation
    @test_throws DimensionMismatch FXSeries(idx[1:4], ["A", "B"], v)
    @test_throws DimensionMismatch FXSeries(idx, ["A"], v)
    @test_throws ArgumentError FXSeries(reverse(idx), ["A", "B"], v)

    ## missing values survive construction as NaN
    vn = copy(v); vn[3, 1] = NaN
    @test isnan(FXSeries(idx, ["A", "B"], vn).values[3, 1])

    ## show must not error
    @test occursin("FXSeries", sprint(show, s))
    @test occursin("FXSeries", sprint(show, MIME"text/plain"(), s))
end

@testset "window (zoo::window)" begin
    idx = [Date(2020, 1, i) for i in 1:5]
    s = FXSeries(idx, ["A"], reshape(Float64.(1:5), 5, 1))

    ## both bounds are inclusive, exactly like zoo::window
    @test index(window(s; start = Date(2020, 1, 2), stop = Date(2020, 1, 4))) == idx[2:4]
    @test index(window(s; start = Date(2020, 1, 2))) == idx[2:5]
    @test index(window(s; stop = Date(2020, 1, 4))) == idx[1:4]
    @test index(window(s)) == idx
    ## a bound that falls between observations
    @test index(window(s; start = Date(2019, 12, 31), stop = Date(2020, 1, 3))) == idx[1:3]
    ## an empty window is legal
    @test size(window(s; start = Date(2021, 1, 1)), 1) == 0

    ## the window used for the CNY history model in the vignette
    @test size(CNY_HIST, 1) == 68
    @test CNY_HIST.index[end] == Date(2005, 10, 31)
end
