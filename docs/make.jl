using Xtallography
using Documenter

DocMeta.setdocmeta!(Xtallography, :DocTestSetup, :(using Xtallography); recursive=true)

makedocs(;
    modules=[Xtallography],
    authors="Kevin Chu <kevin@velexi.com> and contributors",
    sitename="Xtallography.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://velexi-research.github.io/Xtallography.jl",
        repolink="https://github.com/velexi-research/Xtallography.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Overview" => "index.md",
        "Examples" => "examples.md",
        "Reference" => [
            "Symmetry" => "symmetry.md",
            "Unit Cell" => "unit-cell.md",
            "Utility Functions" => "functions.md",
            "Python Interface" => "python-interface.md",
        ],
        "Index" => "docs-index.md",
    ],
    checkdocs=:exports,
    warnonly=[:missing_docs, :cross_references],
)

deploydocs(; repo="github.com/velexi-research/Xtallography.jl", devbranch="main")
