#!/usr/bin/env julia
# Join bench/results_julia.csv and bench/results_r.csv into a comparison table.
#
#   julia --project=bench bench/compare.jl
#
# Run bench/benchmarks.jl and bench/benchmarks.R first.

using DelimitedFiles, Printf

read_results(path) = begin
    raw, _ = readdlm(path, ','; header = true)
    Dict((string(raw[i, 1]), string(raw[i, 2]), Int(raw[i, 3])) => Float64(raw[i, 4])
         for i in axes(raw, 1))
end

jl = read_results(joinpath(@__DIR__, "results_julia.csv"))
r  = read_results(joinpath(@__DIR__, "results_r.csv"))

keys_common = sort(collect(intersect(keys(jl), keys(r))), by = k -> (k[1], k[3]))

fmt(t) = t >= 10 ? @sprintf("%.1f s", t) :
         t >= 0.1 ? @sprintf("%.3f s", t) : @sprintf("%.4f s", t)

println("| task | dataset | n | Julia | R | speed-up |")
println("|---|---|---:|---:|---:|---:|")
for k in keys_common
    task, ds, n = k
    @printf("| `%s` | %s | %d | %s | %s | %.0f× |\n", task, ds, n, fmt(jl[k]), fmt(r[k]), r[k] / jl[k])
end

only_jl = sort(collect(setdiff(keys(jl), keys(r))), by = k -> (k[1], k[3]))
if !isempty(only_jl)
    println("\nJulia only (R not run at this size):")
    for k in only_jl
        @printf("  %-20s %-11s n=%-5d %s\n", k[1], k[2], k[3], fmt(jl[k]))
    end
end
