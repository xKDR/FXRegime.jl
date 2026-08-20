## fxregimes and its methods — port of fxregime::fxregimes, refit.fxregimes, coef.fxregimes,
## fitted.fxregimes, residuals.fxregimes.
## Fixtures: fxregimes_*_{breakpoints,segment_coef,segment_summary,fitted}.csv.

@testset "fxregimes (fxregime::fxregimes)" begin
    for D in DATINGS
        @testset "$(D.tag)" begin
            reg = D.reg
            n = size(D.data, 1)

            @testset "breakpoints / breakdates" begin
                ## the FXRegimes wrapper must delegate to its BreakpointsFull unchanged
                hdr, rows = readfix("fxregimes_$(D.tag)_breakpoints.csv")
                ms = icol(hdr, rows, "m")
                for m in 1:D.maxbreaks
                    sel = ms .== m
                    @test breakpoints(reg; breaks = m)[1] == icol(hdr, rows, "breakpoint")[sel]
                    @test breakdates(reg; breaks = m) == dcol(hdr, rows, "breakdate")[sel]
                end
                sel = kvfix("fxregimes_$(D.tag)_selection.csv")
                bp, obj = breakpoints(reg)
                @test length(bp) == Int(kvnum(sel, "m_default"))
                @test reldev(obj, kvnum(sel, "objective_default")) < 1e-12
                @test breakdates(reg) == D.data.index[bp]
                @test nobs(reg) == n
                @test npar(reg) == size(D.data, 2) + 1
                @test index(reg) == D.data.index
                @test reldev(loglik(reg), -obj) < 1e-15
                @test dof(reg) == npar(reg) * (length(bp) + 1) + length(bp)
            end

            @testset "refit.fxregimes" begin
                hdr, rows = readfix("fxregimes_$(D.tag)_segment_coef.csv")
                rf = refit(reg)
                bp, _ = breakpoints(reg)
                @test length(rf) == length(bp) + 1 == length(rows)

                ## segment i covers observations bp[i-1]+1 … bp[i]: the breakpoint is the LAST
                ## observation of its segment, exactly as in R
                cuts = vcat(0, bp, n)
                @test [nobs(f) for f in rf] == diff(cuts) == icol(hdr, rows, "nobs")
                @test [index(f)[1] for f in rf] == dcol(hdr, rows, "start_date")
                @test [index(f)[end] for f in rf] == dcol(hdr, rows, "end_date")
                for i in eachindex(rf)
                    @test index(rf[i]) == D.data.index[(cuts[i] + 1):cuts[i+1]]
                end
                ## the segments tile the sample exactly
                @test vcat([index(f) for f in rf]...) == D.data.index
                ## and each breakdate is the last date of the segment that ends there
                @test [index(rf[i])[end] for i in 1:length(bp)] == breakdates(reg)

                @testset "coef.fxregimes" begin
                    cf = coef(reg)
                    cn = hdr[5:end]
                    @test cn == coefnames(reg)
                    @test cn[end] == "(Variance)"
                    @test size(cf) == (length(rf), size(D.data, 2) + 1)
                    @test reldev(cf, fmat(hdr, rows, cn)) < 1e-11
                    ## row i is just coef() of the i-th refitted segment
                    for i in eachindex(rf)
                        @test cf[i, :] == coef(rf[i])
                    end
                end
            end

            @testset "summary.lm per regime" begin
                hdr, rows = readfix("fxregimes_$(D.tag)_segment_summary.csv")
                rf = refit(reg)
                segs = icol(hdr, rows, "segment")
                est = Float64[]; se = Float64[]; tv = Float64[]; pv = Float64[]
                terms = String[]
                for f in rf
                    ct = coeftable(f)
                    append!(terms, ct.names)
                    append!(est, ct.estimate); append!(se, ct.stderror)
                    append!(tv, ct.tvalue); append!(pv, ct.pvalue)
                end
                @test segs == vcat([fill(i, size(rf[i].X, 2)) for i in eachindex(rf)]...)
                @test scol(hdr, rows, "term") == terms
                @test reldev(est, fcol(hdr, rows, "estimate")) < 1e-11
                @test reldev(se, fcol(hdr, rows, "std_error")) < 1e-12
                @test reldev(tv, fcol(hdr, rows, "t_value")) < 1e-11
                ## p-values are compared relatively across ~200 orders of magnitude
                @test reldev(pv, fcol(hdr, rows, "p_value")) < 1e-10
            end

            @testset "fitted.fxregimes / residuals.fxregimes" begin
                hdr, rows = readfix("fxregimes_$(D.tag)_fitted.csv")
                @test dcol(hdr, rows, "date") == D.data.index
                fv = fitted(reg)
                rv = residuals(reg)
                @test length(fv) == length(rv) == n
                ## both cross zero, so they are judged on the scale of the series
                @test scaledev(fv, fcol(hdr, rows, "fitted")) < 1e-12
                @test scaledev(rv, fcol(hdr, rows, "residual")) < 1e-11
                ## they are the segment-wise values concatenated, and they add up to y
                rf = refit(reg)
                @test fv == vcat([fitted(f) for f in rf]...)
                @test rv == vcat([residuals(f) for f in rf]...)
                @test scaledev(fv .+ rv, D.data.values[:, 1]) < 1e-14
            end

            @test occursin("Minimum LWZ partition", sprint(show, reg))
            @test FXRegime.formula(reg) ==
                  string(D.data.names[1], " ~ ", join(D.data.names[2:end], " + "))
        end
    end

    @testset "h given as a fraction" begin
        ## R: h < 1 means floor(n * h) observations
        reg = fxregimes(SYNTH1; h = 1 / 6, breaks = 2)
        @test reg.bf.h == 20                              # floor(120/6)
        @test breakpoints(reg; breaks = 1)[1] == breakpoints(SYNTH1_REG; breaks = 1)[1]
    end

    @testset "default maxbreaks" begin
        ## R: breaks = ceiling(n/h) - 2
        reg = fxregimes(SYNTH1; h = 20)
        @test reg.bf.maxbreaks == cld(120, 20) - 2 == 4
    end

    @testset "ic = :BIC" begin
        reg = fxregimes(SYNTH2; h = 25, breaks = 4, ic = :BIC)
        sel = kvfix("fxregimes_synth_break2_selection.csv")
        @test reg.bf.ic === :BIC
        @test length(breakpoints(reg)[1]) == Int(kvnum(sel, "m_selected_BIC"))
        @test_throws ArgumentError fxregimes(SYNTH2; h = 25, breaks = 4, ic = :AIC)
    end
end
