using Documenter, PETSc

<<<<<<< HEAD
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
=======
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
>>>>>>> 30e66665e4f65f42c22d953b96ab1fb234d08fb8
    ],
)

deploydocs(;
    repo="github.com/JuliaParallel/PETSc.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
)