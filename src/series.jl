"""
    FXSeries(index, names, values)

A minimal irregular time-series container: a strictly increasing vector of `Date`s, a vector of
column names, and a `Matrix{Float64}` of observations with `NaN` encoding a missing value.

This plays the role that a `zoo` series plays in the R package **fxregime**. By convention the
first column of a series passed to [`fxlm`](@ref) or [`fxregimes`](@ref) is the target currency
and the remaining columns are the basket currencies.
"""
struct FXSeries
    index::Vector{Date}
    names::Vector{String}
    values::Matrix{Float64}

    function FXSeries(index::AbstractVector{Date}, names::AbstractVector{<:AbstractString},
                      values::AbstractMatrix{<:Real})
        length(index) == size(values, 1) ||
            throw(DimensionMismatch("index has $(length(index)) entries but values has $(size(values,1)) rows"))
        length(names) == size(values, 2) ||
            throw(DimensionMismatch("$(length(names)) names given for $(size(values,2)) columns"))
        issorted(index) || throw(ArgumentError("index must be sorted in increasing order"))
        new(collect(Date, index), String.(names), Matrix{Float64}(values))
    end
end

Base.size(s::FXSeries) = size(s.values)
Base.size(s::FXSeries, d::Int) = size(s.values, d)
Base.length(s::FXSeries) = size(s.values, 1)
nobs(s::FXSeries) = size(s.values, 1)
Base.names(s::FXSeries) = s.names
index(s::FXSeries) = s.index

"""
    colindex(s::FXSeries, name) -> Int

Position of the column called `name`. Errors if it is not present.
"""
function colindex(s::FXSeries, name::AbstractString)
    j = findfirst(==(String(name)), s.names)
    j === nothing && throw(ArgumentError("no column named \"$name\" (have: $(join(s.names, ", ")))"))
    return j
end

Base.getindex(s::FXSeries, ::Colon, name::AbstractString) = s.values[:, colindex(s, name)]
Base.getindex(s::FXSeries, ::Colon, j::Integer) = s.values[:, j]
Base.getindex(s::FXSeries, i, ::Colon) = FXSeries(s.index[i], s.names, s.values[i, :])

"""
    getcols(s::FXSeries, names) -> FXSeries

Sub-series with the named columns, in the order given.
"""
getcols(s::FXSeries, cols::AbstractVector{<:AbstractString}) =
    FXSeries(s.index, String.(cols), s.values[:, [colindex(s, c) for c in cols]])

"""
    window(s::FXSeries; start = nothing, stop = nothing) -> FXSeries

Rows with `start <= date <= stop`; either bound may be `nothing`. Mirrors `zoo::window`.
"""
function window(s::FXSeries; start::Union{Date,Nothing} = nothing,
                             stop::Union{Date,Nothing} = nothing)
    keep = trues(length(s.index))
    start === nothing || (keep .&= s.index .>= start)
    stop  === nothing || (keep .&= s.index .<= stop)
    return FXSeries(s.index[keep], s.names, s.values[keep, :])
end

Base.first(s::FXSeries) = first(s.index)
Base.last(s::FXSeries) = last(s.index)

function Base.show(io::IO, ::MIME"text/plain", s::FXSeries)
    n, k = size(s)
    println(io, "FXSeries: $n observations × $k series")
    n == 0 && return
    println(io, "  ", first(s.index), " .. ", last(s.index))
    println(io, "  ", join(s.names, ", "))
    nshow = min(n, 4)
    for i in 1:nshow
        println(io, "  ", s.index[i], "  ",
                join((v -> isnan(v) ? "     NA" : rpad(round(v, digits = 6), 7)).(s.values[i, :]), "  "))
    end
    n > nshow && print(io, "  ⋮")
end

Base.show(io::IO, s::FXSeries) = show(io, MIME"text/plain"(), s)
