using XtallographyUtils
using Documenter

DocMeta.setdocmeta!(XtallographyUtils, :DocTestSetup, :(using XtallographyUtils); recursive=true)

makedocs(;
    modules=[XtallographyUtils],
    authors="Kevin Chu <kevin@velexi.com> and contributors",
    repo="https://github.com/velexi-research/XtallographyUtils.jl/blob/{commit}{path}#{line}",
    sitename="XtallographyUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://velexi-research.github.io/XtallographyUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/velexi-research/XtallographyUtils.jl",
    devbranch="main",
)
