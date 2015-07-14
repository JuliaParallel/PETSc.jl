
using PETSc
using FactCheck


import MPI
MPI.Init()


# write your own tests here


# define some arrays to make looping easier
vec_norms = [PETSC_NORM_1, PETSC_NORM_2, PETSC_NORM_INFINITY]

# for some reason only the first four vector format works
#vec_formats = [VECSEQ, VECMPI, VECSTANDARD, VECSHARED, VECSEQCUSP, VECMPICUSP, VECCUSP, VECSEQVIENNACL, VECMPIVIENNACL, VECVIENNACL, VECNEST]  # no pthread
vec_formats = [VECSEQ, VECMPI, VECSTANDARD, VECSHARED]

PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]);

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)

sys_size = 3

rhs = Array(1.0:3)
A_julia = [1.0 2.0 3; 4 5 7; 7 8 9]
#println("rhs = ", rhs)
#println("A_julia = ", A_julia)
x_julia = A_julia\rhs
#println("x_julia = ", x_julia)

include("test_vec.jl")
include("test_mat.jl")
include("test_ksp.jl")

FactCheck.exitstatus()
