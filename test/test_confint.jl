## Confidence intervals for breakpoints — port of fxregime:::confint.fxregimes and of the
## unexported helper strucchange:::pargmaxV.
## Fixtures: confint_<tag>_{90,95}.csv.

@testset "confint.fxregimes (fxregime)" begin
    for D in DATINGS
        @testset "$(D.tag)" begin
            for (level, suffix) in ((0.9, "90"), (0.95, "95"))
                @testset "level = $level" begin
                    hdr, rows = readfix("confint_$(D.tag)_$(suffix).csv")
                    ci = confint(D.reg; level = level)

                    ## the bounds are observation numbers: they must match R exactly
                    @test ci.breakpoints == icol(hdr, rows, "breakpoint")
                    @test isequal(ci.lower, mcol(hdr, rows, "lower"))
                    @test isequal(ci.upper, mcol(hdr, rows, "upper"))
                    @test ci.level == level
                    @test ci.nobs == size(D.data, 1)
                    @test ci.npar == size(D.data, 2) + 1

                    ## and the dates that go with them line up with R's to the day
                    @test isequal(ci.dates.lower, dcol(hdr, rows, "lower_date"))
                    @test isequal(ci.dates.breakpoints, dcol(hdr, rows, "breakdate"))
                    @test isequal(ci.dates.upper, dcol(hdr, rows, "upper_date"))
                    @test ci.dates.breakpoints == breakdates(D.reg)

                    ## structural sanity: the interval brackets the point estimate
                    for i in eachindex(ci.breakpoints)
                        @test ci.lower[i] <= ci.breakpoints[i] <= ci.upper[i]
                        @test 1 <= ci.lower[i]
                        @test ci.upper[i] <= size(D.data, 1)
                        @test ci.dates.lower[i] == D.data.index[ci.lower[i]]
                        @test ci.dates.upper[i] == D.data.index[ci.upper[i]]
                    end
                    @test icol(hdr, rows, "k") == collect(1:length(ci.breakpoints))
                end
            end

            ## a wider level can only widen the interval
            c90 = confint(D.reg; level = 0.9)
            c95 = confint(D.reg; level = 0.95)
            @test all(c95.lower .<= c90.lower)
            @test all(c95.upper .>= c90.upper)
        end
    end

    @testset "breaks argument" begin
        ## an explicitly requested partition gets its own intervals
        ci = confint(SYNTH2_REG; level = 0.9, breaks = 1)
        @test ci.breakpoints == breakpoints(SYNTH2_REG; breaks = 1)[1]
        @test length(ci.lower) == 1
        ## R refuses an interval when there is no break at all
        @test_throws ArgumentError confint(SYNTH1_REG; breaks = 0)
    end

    @testset "pargmaxV (strucchange:::pargmaxV)" begin
        ## reference values from Rscript -e 'strucchange:::pargmaxV(x)', printed with %.17g
        @test reldev(FXRegime.pargmaxV(-5.0), 0.092766506878279031) < 1e-12
        @test reldev(FXRegime.pargmaxV(-1.0), 0.30114608758464689) < 1e-12
        @test reldev(FXRegime.pargmaxV(-0.001), 0.49951645290936342) < 1e-12
        @test reldev(FXRegime.pargmaxV(0.0), 0.5) < 1e-15
        @test reldev(FXRegime.pargmaxV(0.001), 0.50048354709063614) < 1e-12
        @test reldev(FXRegime.pargmaxV(1.0), 0.69885391241535311) < 1e-12
        @test reldev(FXRegime.pargmaxV(5.0), 0.90723349312172097) < 1e-12
        @test reldev(FXRegime.pargmaxV(50.0), 0.99995809115668377) < 1e-12
        ## with the general xi / phi arguments
        @test reldev(FXRegime.pargmaxV(2.0; xi = 0.5, phi1 = 2.0, phi2 = 0.7),
                     0.98602650262456604) < 1e-12
        ## it is a distribution function: increasing, symmetric about 0 in the standard case
        xs = -20.0:0.5:20.0
        ps = [FXRegime.pargmaxV(x) for x in xs]
        @test issorted(ps)
        @test all(0 .<= ps .<= 1)
        @test reldev(FXRegime.pargmaxV(3.0) + FXRegime.pargmaxV(-3.0), 1.0) < 1e-12
    end

    @testset "uniroot (R_zeroin2)" begin
        ## the bisection used to invert pargmaxV must find roots to R's default tolerance
        ## R: uniroot(function(x) strucchange:::pargmaxV(x) - 0.95, c(0, 1000))$root
        ## and the mirror image on the negative side.  Brent's method is stopped on the
        ## *argument*, at R's default tol = .Machine$double.eps^0.25 ~ 1.2e-4, so the root
        ## itself is only pinned to that many digits — hence the loose tolerance here.
        f = x -> FXRegime.pargmaxV(x) - 0.95
        r = FXRegime._uniroot(f, 0.0, 1000.0)
        @test reldev(r, 7.6872777677720174) < 1e-6
        g = x -> FXRegime.pargmaxV(x) - 0.05
        @test reldev(FXRegime._uniroot(g, -1000.0, 0.0), -7.6872777677719935) < 1e-6
        @test_throws ArgumentError FXRegime._uniroot(x -> x^2 + 1, -1.0, 1.0)
    end
end
