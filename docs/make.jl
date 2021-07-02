using Documenter, PETSc

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
            "Mat" =>  "mat.md",
        ],
        "List of functions"  => "listfunctions.md"
    ],
)
