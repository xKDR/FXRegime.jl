#!/usr/bin/env julia
# Julia side of the FXRegime.jl <-> R fxregime benchmark.
#
#   julia --project=bench bench/benchmarks.jl [--threads-note]
#
# Writes bench/results_julia.csv with one row per (task, n) and the median wall time in seconds.
# Run bench/benchmarks.R for the R side, then bench/compare.jl to join them.

using FXRegime, Dates, Statistics, Printf, DelimitedFiles

const OUT = joinpath(@__DIR__, "results_julia.csv")

"Median of `k` timings of `f()`, after one warm-up call to take JIT out of the measurement."
function timeit(f; k::Int = 5)
    f()
    ts = Float64[]
    for _ in 1:k
        push!(ts, @elapsed f())
    end
    return median(ts)
end

fx = FXRatesCHF()

## The two published analyses -------------------------------------------------
cny = fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                frequency = :daily, start = Date(2005, 7, 25), stop = Date(2009, 7, 31))
inr = fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], fx;
                frequency = :weekly, start = Date(1993, 4, 1), stop = Date(2008, 1, 4))

## A long daily series for the scaling study ---------------------------------
long = fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], fx;
                 frequency = :daily, start = Date(1993, 4, 1), stop = Date(2008, 1, 4))

rows = Vector{Any}[]
push!(rows, Any["task", "dataset", "n", "seconds"])

record!(task, ds, n, t) = (push!(rows, Any[task, ds, n, t]);
                           @printf("%-22s %-10s n=%-6d %9.4f s\n", task, ds, n, t))

## 1. Data preparation --------------------------------------------------------
record!("fxreturns", "CNY-daily", length(cny),
        timeit(() -> fxreturns("CNY", ["USD", "JPY", "EUR", "GBP"], fx;
                               frequency = :daily, start = Date(2005, 7, 25),
                               stop = Date(2009, 7, 31))))
record!("fxreturns", "INR-weekly", length(inr),
        timeit(() -> fxreturns("INR", ["USD", "JPY", "DUR", "GBP"], fx;
                               frequency = :weekly, start = Date(1993, 4, 1),
                               stop = Date(2008, 1, 4))))

## 2. Single regression -------------------------------------------------------
record!("fxlm", "CNY-daily", length(cny), timeit(() -> fxlm(cny); k = 50))

## 3. Recursive residuals -----------------------------------------------------
let m = fxlm(cny), X = m.X, y = m.y
    record!("recresid", "CNY-daily", size(X, 1), timeit(() -> recresid(X, y); k = 10))
end

## 4. The O(n^2) triangle, isolated ------------------------------------------
for (name, s) in (("CNY-daily", cny), ("INR-weekly", inr))
    m = fxlm(s)
    record!("rss_triangle", name, size(m.X, 1),
            timeit(() -> rss_triangle(m.y, m.X, 20); k = 3))
end

## 5. Full dating -------------------------------------------------------------
record!("fxregimes", "CNY-daily", length(cny),
        timeit(() -> fxregimes(cny; h = 20, breaks = 10); k = 3))
record!("fxregimes", "INR-weekly", length(inr),
        timeit(() -> fxregimes(inr; h = 20, breaks = 10); k = 3))

## 6. Confidence intervals ----------------------------------------------------
let reg = fxregimes(cny; h = 20, breaks = 10)
    record!("confint", "CNY-daily", length(cny),
            timeit(() -> confint(reg; level = 0.9); k = 10))
end

## 7. Scaling in n ------------------------------------------------------------
for n in (250, 500, 1000, 2000, 3000)
    n > length(long) && continue
    s = long[1:n, :]
    record!("fxregimes-scaling", "INR-daily", n,
            timeit(() -> fxregimes(s; h = 20, breaks = 10); k = 2))
end

open(OUT, "w") do io
    writedlm(io, rows, ',')
end
@printf("\nthreads = %d\nwrote %s\n", Threads.nthreads(), OUT)
