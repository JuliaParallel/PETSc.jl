using Documenter, PETSc

#makedocs(modules = [PETSc])

#deploydocs(
#  deps = Deps.pip("mkdocs", "python-markdown-math"),
#  repo = "github.com/JuliaParallel/PETSc.jl.git",
#  julia = "0.4",
#)

#using Documenter

DocMeta.setdocmeta!(PETSc, :DocTestSetup, :(using PETSc); recursive=true)

makedocs(;
    modules=[PETSc],
    repo="https://github.com/JuliaParallel/PETSc/{commit}{path}#{line}",
    sitename="PETSc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/JuliaParallel/PETSc.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "PETSc" => Any[
            "DMStag" =>  "dmstag.md",
        ],
        "List of functions"  => "listfunctions.md"
    ],
)

#deploydocs(;
#    repo="github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl.git",
#    branch = "gh-pages",
#    target = "build",
#    devbranch = "main",
#    devurl = "dev",
#)
