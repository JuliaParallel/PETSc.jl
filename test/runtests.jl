
using PETSc
using FactCheck


import MPI
println("pwd = ", pwd())
println("PETSC_DIR = ", ENV["PETSC_DIR"])
MPI.Init()

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)


# write your own tests here


# define some arrays to make looping easier
#vec_norms = [PETSC_NORM_1, PETSC_NORM_2, PETSC_NORM_INFINITY]

# skip 1 norm untill Petsc 1 complex 1 norm fix is released
vec_norms = [ PETSC_NORM_2, PETSC_NORM_INFINITY]
julia_vec_norms = [ 2, Inf]
# for some reason only the first four vector format works
#vec_formats = [VECSEQ, VECMPI, VECSTANDARD, VECSHARED, VECSEQCUSP, VECMPICUSP, VECCUSP, VECSEQVIENNACL, VECMPIVIENNACL, VECVIENNACL, VECNEST]  # no pthread

if comm_size > 1  # this is a parallel run
  vec_formats = [ VECMPI, VECSTANDARD]
else  # serial run 
  vec_formats = [VECSEQ]
end

PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]);

# size of the system owned by this process (ie. local sys size)
sys_size = PetscInt(3)

# indicies of the vector owned by this process
# create these with smallest precision, so they can be promoted
tmp2 = convert(Array{Float32,1}, Array(1.0:3))
tmp = convert(Array{Float32, 2}, [1.0 2.0 3; 4 5 7; 7 8 9])

tmp3 = convert(Array{Complex64, 1}, [1.0 + 0im; 2.0 + 1.0im; 3.0 + 2.0im])
tmp4 = convert( Array{Complex64, 2}, [1.0 + 1im   2 + 2im  3 + 3im; 4 + 4im  5 + 5im 7 + 7im; 7 + 7im 8 + 8im 9 + 9im])

# Have both local and global versions of system for testing
A_julia = zeros(PetscScalar, sys_size, sys_size)
rhs = zeros(PetscScalar, sys_size)
A_julia_global = zeros(PetscScalar, comm_size*sys_size, comm_size*sys_size)
rhs_global = zeros(PetscScalar, comm_size*sys_size)
# convert to arrays of the proper Petsc type
# this facilitates testing with the different Petsc build options
if PetscScalar <: Real
  for i=1:sys_size
    rhs[i] = convert(PetscScalar, tmp2[i])
    
    for j=0:(comm_size - 1)  # global 
      idx = i + j*sys_size
      rhs_global[idx] = convert(PetscScalar, tmp2[i])
    end    

  end

  for i=1:sys_size
    for j=1:sys_size
      A_julia[i,j] = convert(PetscScalar, tmp[i,j])

      for k=0:(comm_size - 1)
        idxm = i + k*sys_size
        idxn = j + k*sys_size
        A_julia_global[idxm, idxn] = convert(PetscScalar, tmp[i,j] )
      end
    end
  end

end  # end if statement

if PetscScalar <: Complex
  for i=1:sys_size
    rhs[i] = convert(PetscScalar, tmp3[i])

    for j=0:(comm_size - 1)  # global 
      idx = i + j*sys_size
      rhs_global[idx] = convert(PetscScalar, tmp3[i])
    end    


  end

  for i=1:sys_size
    for j=1:sys_size
      A_julia[i,j] = convert(PetscScalar, tmp4[i,j])

      for k=0:(comm_size - 1)
        idxm = i + k*sys_size
        idxn = j + k*sys_size
        A_julia_global[idxm, idxn] = convert(PetscScalar, tmp4[i,j] )
      end
 
    end
  end
end
#

if comm_rank == 0
  println("rhs_global = ", rhs_global)
  println("A_julia_global = ", A_julia_global)
end

#println("rhs = ", rhs)
#println("A_julia = ", A_julia)
x_julia = A_julia\rhs
#println("x_julia = ", x_julia)

include("test_vec.jl")
#include("test_mat.jl")
#include("test_ksp.jl")

PetscFinalize()
FactCheck.exitstatus()
