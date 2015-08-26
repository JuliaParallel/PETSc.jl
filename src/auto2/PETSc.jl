module PETSc

generated_path = joinpath(Pkg.dir("PETSc"), "src/auto2/generated")
push!(LOAD_PATH, generated_path)
import C

include("petsc_com.jl")


end
