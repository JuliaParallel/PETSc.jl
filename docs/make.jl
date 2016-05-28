using Documenter, PETSc

makedocs(modules = [PETSc])

deploydocs(
  deps = Deps.pip("mkdocs", "python-markdown-math"),
  repo = "github.com/JuliaParallel/PETSc.jl.git",
  julia = "0.4",
)
