
# ensure in right directory and environment
cd(dirname(@__FILE__))
using Pkg
Pkg.activate(".")
using Bolt
using Documenter
ENV["GKSwstype"] = 100 # https://github.com/jheinen/GR.jl/issues/318#issuecomment-651890036


# change to true to skip running doc code
Documenter.Utilities.Selectors.disable(::Type{Documenter.Expanders.ExampleBlocks}) = false

md_filenames = []
for filename in readdir("src")
    if endswith(filename, "ipynb")
        contents = read(`jupytext src/$filename --to md --opt notebook_metadata_filter=-jupytext,-kernelspec -o -`, String)
        # replace input cells with @example blocks so Documenter runs them
        contents = replace(contents, "```julia"=>"```@example 1")
        # change $$ ... $$ equations to ```math ... ````
        contents = replace(contents, r"\$\$(.*?)\$\$"s => s"""```math
        \g<1>
        ```""")
        # Documenter doesn't suppress outputs with ';', so add nothing return values
        contents = replace(contents, r"^(.*?;)\s*$"m => s"""\g<1>
        nothing # hide""")
        # starting a line with inline math screws up tex2jax for some reason
        contents = replace(contents, r"\* \$(.*?)\$" => s"* ``\g<1>``")
        md_filename = joinpath("src", splitext(filename)[1] * ".md")
        push!(md_filenames, md_filename)
        open(md_filename, "w") do f
            write(f, contents)
        end
    end
end

rm("src/index.md", force=true)
symlink("../../README.md", "src/index.md")


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
        "Bolt.jl" => "index.md",
        "Basic Usage" => "basic_usage.md",
        "dev.md",
        "api.md"
    ],
)

rm.(md_filenames)

deploydocs(;
    repo="github.com/xzackli/Bolt.jl",
    push_preview = true,
    forcepush = true
)
