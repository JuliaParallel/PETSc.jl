using PETSc
PETSc.PetscInitialize()
using FactCheck

import MPI
MPI.Init()

include("test_error.jl")
include("test_vec.jl")
include("test_mat.jl")
PETSc.C.PetscFinalize(Float64)
