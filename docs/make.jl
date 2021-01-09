using Bolt
using Documenter

makedocs(;
    modules=[Bolt],
    authors="Zack Li",
    repo="https://github.com/xzackli/Bolt.jl/blob/{commit}{path}#L{line}",
    sitename="Bolt.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/Bolt.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xzackli/Bolt.jl",
)
