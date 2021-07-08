using Documenter, PETSc

makedocs(;
    modules=[PETSc],
    sitename="PETSc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "man/getting_started.md",
        "PETSc" => Any[
            "Mat" =>  "man/mat.md",
        ],
        "List of functions"  => "man/listfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaParallel/PETSc.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
)