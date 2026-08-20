## fxreturns — port of fxregime::fxreturns.
## Fixtures: fxreturns_*.csv (see MANIFEST.csv for the exact R calls).

"Compare an FXSeries against a `date, col1, col2, …` fixture."
function check_returns_fixture(name::AbstractString, s::FXSeries; rtol::Float64)
    hdr, rows = readfix(name)
    @test dcol(hdr, rows, "date") == s.index          # dates must agree to the day
    @test hdr[2:end] == s.names
    @test size(s) == (length(rows), length(hdr) - 1)
    ## returns pass through zero, so they are judged on the scale of the whole panel
    @test scaledev(s.values, fmat(hdr, rows, hdr[2:end])) < rtol
end

@testset "fxreturns (fxregime::fxreturns)" begin

    @testset "daily, CNY vs USD/JPY/EUR/GBP" begin
        ## fxreturns("CNY", frequency = "daily", start = 2005-07-25, end = 2009-07-31, ...)
        check_returns_fixture("fxreturns_cny_daily.csv", CNY; rtol = 1e-12)
        @test size(CNY) == (1014, 5)
        @test CNY.index[1] == Date(2005, 7, 26)      # one row is lost to differencing
        @test CNY.index[end] == Date(2009, 7, 31)
    end

    @testset "weekly, INR vs USD/JPY/DUR/GBP" begin
        ## fxreturns("INR", frequency = "weekly", start = 1993-04-01, end = 2008-01-04, ...)
        check_returns_fixture("fxreturns_inr_weekly.csv", INR; rtol = 1e-12)
        @test size(INR) == (770, 5)
        ## weekly aggregation is anchored on Fridays, so every index date is a Friday
        @test all(Dates.dayofweek.(INR.index) .== Dates.Friday)
    end

    @testset "trim = TRUE" begin
        ## fxreturns(..., trim = TRUE) == trim = c(0.01, 0.99) on the target column
        trimmed = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], FX; frequency = :daily,
                            start = Date(2005, 7, 25), stop = Date(2009, 7, 31), trim = true)
        check_returns_fixture("fxreturns_cny_daily_trim.csv", trimmed; rtol = 1e-12)
        @test size(trimmed, 1) == 992

        ## the surviving rows are exactly those inside R's type-7 quantiles of column 1
        x = CNY.values[:, 1]
        lo, hi = quantile(x, 0.01), quantile(x, 0.99)
        keep = .!((x .< lo) .| (x .> hi))
        @test trimmed.index == CNY.index[keep]
        ## explicit levels give the same thing, in either order (R sorts them)
        @test fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], FX; frequency = :daily,
                        start = Date(2005, 7, 25), stop = Date(2009, 7, 31),
                        trim = (0.99, 0.01)).index == trimmed.index
    end

    @testset "na.action = na.omit" begin
        omitted = fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], FX; frequency = :weekly,
                            start = Date(1993, 4, 1), stop = Date(2008, 1, 4), na_action = :omit)
        check_returns_fixture("fxreturns_inr_weekly_naomit.csv", omitted; rtol = 1e-12)
        ## on this window na.locf and na.omit happen to coincide
        @test omitted.index == INR.index
    end

    @testset "column selection order" begin
        ## R does `which(colnames(data) %in% x)`, so the columns come out in *data* order,
        ## not in the order the caller listed them: USD precedes GBP in FXRatesCHF.
        s = fxreturns("INR", ["GBP", "USD"], FX; frequency = :weekly,
                      start = Date(1993, 4, 1), stop = Date(1994, 1, 1))
        @test s.names == ["INR", "USD", "GBP"]
        ## the target is always first even though INR sits after USD/GBP in the data
        @test s.names[1] == "INR"
        ## a currency named twice is taken once (R: unique(c(x, other)))
        @test fxreturns("INR", ["INR", "USD"], FX; frequency = :weekly,
                        start = Date(1993, 4, 1), stop = Date(1994, 1, 1)).names == ["INR", "USD"]
        ## R's `which(colnames(data) %in% x)` returns integer(0) for an unknown currency, so
        ## R silently drops the target instead of failing; checked with
        ##   fxreturns("XXX", other = "USD", ...)  ->  a one-column series.
        @test size(fxreturns("XXX", ["USD"], FX; frequency = :weekly,
                             start = Date(1993, 4, 1), stop = Date(1994, 1, 1))) == (39, 1)
        @test_throws ArgumentError fxreturns("XXX", String[], FX)   # nothing left at all
        @test_throws ArgumentError fxreturns("INR", ["USD"], FX; frequency = :monthly)
        @test_throws ArgumentError fxreturns("INR", ["USD"], FX; na_action = :interpolate)
    end

    @testset "nextfriday (fxregime:::nextfri)" begin
        ## R: nextfri <- function(date) 7 * ceiling(as.numeric(date - 1)/7) + as.Date(1)
        ## reference values obtained by evaluating that expression in R
        @test nextfriday(Date(1970, 1, 1)) == Date(1970, 1, 2)   # Thu -> Fri
        @test nextfriday(Date(1970, 1, 2)) == Date(1970, 1, 2)   # Fri -> itself
        @test nextfriday(Date(1970, 1, 3)) == Date(1970, 1, 9)   # Sat -> next Fri
        @test nextfriday(Date(1970, 1, 8)) == Date(1970, 1, 9)
        @test nextfriday(Date(1969, 12, 25)) == Date(1969, 12, 26)  # before the epoch
        @test nextfriday(Date(2005, 7, 25)) == Date(2005, 7, 29)
        @test nextfriday(Date(2005, 7, 29)) == Date(2005, 7, 29)
        @test nextfriday(Date(2005, 7, 30)) == Date(2005, 8, 5)
        @test nextfriday(Date(2008, 2, 29)) == Date(2008, 2, 29)  # leap day, a Friday
        ## the anchor is always a Friday and never more than 6 days ahead
        for d in Date(1999, 1, 1):Day(1):Date(1999, 3, 1)
            f = nextfriday(d)
            @test Dates.dayofweek(f) == Dates.Friday
            @test 0 <= Dates.value(f - d) <= 6
        end
    end
end
