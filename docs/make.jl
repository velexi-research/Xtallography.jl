using XtallographyUtils
using Documenter

DocMeta.setdocmeta!(
    XtallographyUtils, :DocTestSetup, :(using XtallographyUtils); recursive=true
)

makedocs(;
    modules=[XtallographyUtils],
    authors="Kevin Chu <kevin@velexi.com> and contributors",
    sitename="XtallographyUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://velexi-research.github.io/XtallographyUtils.jl",
        repolink="https://github.com/velexi-research/XtallographyUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Types" => "types.md",
        "Functions" => "functions.md",
        "Index" => "docs-index.md",
    ],
    warnonly=[:missing_docs],
)

deploydocs(; repo="github.com/velexi-research/XtallographyUtils.jl", devbranch="main")
