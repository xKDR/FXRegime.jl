using FXRegime
using Documenter

DocMeta.setdocmeta!(FXRegime, :DocTestSetup, :(using FXRegime); recursive=true)

makedocs(;
    modules=[FXRegime],
    authors="xKDR Forum",
    repo=Remotes.GitHub("xKDR", "FXRegime.jl"),
    sitename="FXRegime.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xKDR.github.io/FXRegime.jl",
        edit_link="main",
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
