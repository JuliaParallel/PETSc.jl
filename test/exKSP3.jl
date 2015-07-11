using PETSc

import MPI
MPI.Init()

PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]);
#
#   Create a vector and put values in it

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)

sys_size = 3
println("my comm rank = ", comm_rank)

b = PetscVec(comm);
PetscVecSetType(b,"mpi");
PetscVecSetSizes(b,sys_size, comm_size*sys_size);

x = PetscVec(comm);
PetscVecSetType(x,"mpi");
PetscVecSetSizes(x,sys_size, comm_size*sys_size);




rhs = zeros(sys_size)
for i=1:sys_size
  idxm = [(comm_rank)*10 + i]  # index
  PetscVecSetValues(b,idxm,[ convert(Float64,i)],PETSC_INSERT_VALUES);
  rhs[i] = i
end




PetscVecAssemblyBegin(b);
println("began vector assembly")
PetscVecAssemblyEnd(b);
println("finished vector assembly")
PetscView(b);

println("finished assembling vector")
#
#    Create a matrix


A_julia = [1.0 2.0 3; 4 5 7; 7 8 9]
println("rhs = ", rhs)
println("A_julia = ", A_julia)
x_julia = A_julia\rhs
println("x_julia = ", x_julia)

A = PetscMat(comm);
PetscMatSetType(A,"mpiaij");
PetscMatSetSizes(A,sys_size,sys_size,comm_size*sys_size,comm_size*sys_size);
PetscSetUp(A);
for i=1:sys_size
  for j = 1:sys_size
    idxm = [(comm_rank)*sys_size + i]  # row index
    idxn = [(comm_rank)*sys_size + j]  # column index
    PetscMatSetValues(A,idxm, idxn, [A_julia[i,j]],PETSC_INSERT_VALUES);
  end
end

println("finished setting matrix values")
PetscMatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
PetscMatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);
println("finished assembling matrix")
PetscView(A);


println("finished viewing matrix")

println("setting up linear solve")
# setup linear solve
ksp = KSP(comm)
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)
println("performing linear solve")
KSPSolve(ksp, b, x)
println("finished performing linear solve")
PetscView(x)

# copy solution back to Julia
x_copy = zeros(sys_size)
idx = Array(0:2)

PetscVecGetValues(x, sys_size, idx, x_copy)

println("x_copy = ", x_copy)

err_norm = norm(x_copy - x_julia)
println("err norm = ", err_norm)

PetscDestroy(b)
PetscDestroy(A)
PetscDestroy(ksp)
# remove references to all Petsc variables
# perform solution
#A_julia\rhs
#println("x = ", x)

#PetscFinalize()

