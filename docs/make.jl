using PAMG
using Documenter

DocMeta.setdocmeta!(PAMG, :DocTestSetup, :(using PAMG); recursive=true)

makedocs(;
    modules=[PAMG],
    authors="Quesys-tech",
    repo="https://github.com/Quesys-tech/PAMG.jl/blob/{commit}{path}#{line}",
    sitename="PAMG.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Quesys-tech.github.io/PAMG.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Quesys-tech/PAMG.jl",
    devbranch="main",
)
