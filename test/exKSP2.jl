include("PETSc.jl")

import MPI
MPI.Init()

PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]);
#
#   Create a vector and put values in it

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)

println("my comm rank = ", comm_rank)

b = PetscVec(comm);
PetscVecSetType(b,"mpi");
PetscVecSetSizes(b,10, comm_size*10);


for i=1:10
  idxm = [(comm_rank)*10 + i]  # index
  PetscVecSetValues(b,idxm,[10.0 + i],PETSC_ADD_VALUES);
end




#PetscVecSetValues(b,Array(1.0:10.));
println("assigned first set of values")
#PetscVecSetValues(b,[1,2],[11.5,12.5],PETSC_ADD_VALUES);
println("assigned second set of values")
PetscVecAssemblyBegin(b);
println("began vector assembly")
PetscVecAssemblyEnd(b);
println("finished vector assembly")
PetscView(b);

println("finished assembling vector")
#
#    Create a matrix
A = PetscMat(comm);
PetscMatSetType(A,"mpiaij");
PetscMatSetSizes(A,10,10,comm_size*10,comm_size*10);
PetscSetUp(A);
for i=1:10
  idxm = [(comm_rank)*10 + i]  # row index
  idxn = [(comm_rank)*10 + i]  # column index
  PetscMatSetValues(A,idxm, idxn,[10.0 + i],PETSC_ADD_VALUES);
end

println("finished setting matrix values")
PetscMatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
PetscMatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);
println("finished assembling matrix")
PetscView(A);

println("finished viewing matrix")
