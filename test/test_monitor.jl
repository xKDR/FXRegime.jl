## fxmonitor — port of fxregime::fxmonitor, breakpoints.fxmonitor, breakdates.fxmonitor.
## Fixtures: fxmonitor_cny_{process,J12,summary,dates}.csv.

@testset "fxmonitor (fxregime::fxmonitor)" begin
    ## the monitoring example of the CNY vignette:
    ##   fxmonitor(CNY ~ USD + JPY + EUR + GBP, data = window(cny, end = "2006-05-31"),
    ##             start = as.Date("2005-11-01"), end = 4)
    data = window(CNY; stop = Date(2006, 5, 31))
    mon = fxmonitor(data; start = Date(2005, 11, 1), stop = 4)

    st = kvfix("fxmonitor_cny_summary.csv")
    dt = kvfix("fxmonitor_cny_dates.csv")

    @testset "history period" begin
        @test mon.n == Int(kvnum(st, "n")) == 68
        @test mon.monitor == Date(dt["monitor_end_of_history"])
        ## the history is everything strictly before `start`
        @test mon.monitor == last(window(data; stop = Date(2005, 10, 31)).index)
        @test mon.index[1] == Date(dt["history_start"])
        @test mon.index[end] == Date(dt["process_end"])
        @test mon.index == data.index
        @test mon.coefnames == CNY_LM.coefnames
    end

    @testset "process" begin
        hdr, rows = readfix("fxmonitor_cny_process.csv")
        @test icol(hdr, rows, "t") == collect(1:size(mon.process, 1))
        @test dcol(hdr, rows, "date") == mon.index
        @test hdr[3:end] == mon.coefnames
        @test size(mon.process) == (Int(kvnum(st, "nproc")), length(mon.coefnames))
        ## the process crosses zero by construction, so it is judged against its own scale
        @test scaledev(mon.process, fmat(hdr, rows, hdr[3:end])) < 1e-12

        hdr, rows = readfix("fxmonitor_cny_J12.csv")
        @test scaledev(mon.J12, fmat(hdr, rows, hdr[2:end])) < 1e-12
        @test scaledev(mon.J12, transpose(mon.J12)) < 1e-14        # symmetric square root
    end

    @testset "critical value (Zeileis et al. 2005, Table III)" begin
        @test reldev(mon.critval, kvnum(st, "critval")) < 1e-12
        ## the same lookup at other horizons and levels, from R
        for (stop, alpha, cv) in ((3, 0.05, 2.352956971853843),
                                  (10, 0.01, 3.3146038934411859),
                                  (2, 0.20, 1.6383518975935751),
                                  (12, 0.001, 3.4009999999999998))   # beyond the table: clamped
            m = fxmonitor(data; start = Date(2005, 11, 1), stop = stop, alpha = alpha)
            @test reldev(m.critval, cv) < 1e-12
        end
    end

    @testset "breakpoints.fxmonitor / breakdates.fxmonitor" begin
        bp = breakpoints(mon)
        @test bp == Int(kvnum(st, "breakpoint")) == 167
        @test breakdates(mon) == Date(dt["break_date"]) == Date(2006, 3, 27)
        @test mon.index[bp] == breakdates(mon)

        ## the boundary is critval * t/n, crossed for the first time exactly at bp
        maxabs = t -> maximum(abs, view(mon.process, t, :))
        @test maxabs(bp) > (bp / mon.n) * mon.critval
        @test all(t -> maxabs(t) <= (t / mon.n) * mon.critval, (mon.n + 1):(bp - 1))
        ## monitoring only starts after the history period
        @test bp > mon.n

        ## breakpoints found at the other critical values (from R)
        for (stop, alpha, b) in ((3, 0.05, 167), (10, 0.01, 171), (2, 0.20, 113),
                                 (12, 0.001, 171))
            m = fxmonitor(data; start = Date(2005, 11, 1), stop = stop, alpha = alpha)
            @test breakpoints(m) == b
        end
    end

    @testset "no break detected" begin
        ## stop the series inside the stable history extension: nothing to detect
        quiet = fxmonitor(window(CNY; stop = Date(2005, 12, 31));
                          start = Date(2005, 11, 1), stop = 4)
        @test breakpoints(quiet) === nothing
        @test breakdates(quiet) === nothing
        @test occursin("none", sprint(show, quiet))
    end

    @testset "argument checks" begin
        @test_throws ArgumentError fxmonitor(CNY; start = Date(2000, 1, 1))
    end

    @test occursin("Monitoring of FX model", sprint(show, mon))
end
