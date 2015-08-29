using PETSc
using FactCheck

import MPI
MPI.Init()

include("test_vec.jl")

PETSc.C.PetscFinalize(Float64)
