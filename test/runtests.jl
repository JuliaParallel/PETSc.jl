using PETSc
PETSc.PetscInitialize()
using FactCheck

import MPI
MPI.Init()

include("error.jl")
include("vec.jl")
include("mat.jl")
PETSc.C.PetscFinalize(Float64)
