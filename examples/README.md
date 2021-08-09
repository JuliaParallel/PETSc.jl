# `PETSc.jl` Examples

The examples directory includes a `Project.toml` file. Since there is no
`Manifest.toml` when instantiated, this will use the latest released version of
`PETSc.jl`.

To use a specific branch, for example `main`, you should run
```sh
julia --project=examples -e "import Pkg; Pkg.add(name=\"PETSc\", rev=\"main\");"
```

If you are developing the repository you should run
```sh
julia --project=examples -e "import Pkg; Pkg.develop(path=@__DIR__);"
```
in order for your current developments to be used.
