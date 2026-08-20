using FXRegime
using Documenter

DocMeta.setdocmeta!(FXRegime, :DocTestSetup, :(using FXRegime); recursive=true)

makedocs(;
    modules=[FXRegime],
    authors="xKDR Forum",
    repo="https://github.com/xKDR/FXRegime.jl/blob/{commit}{path}#{line}",
    sitename="$FXRegime.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xKDR.github.io/FXRegime.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xKDR/FXRegime.jl",
    target = "build",
    devbranch="main"
)
