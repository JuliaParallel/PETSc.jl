facts("\n   ---testing KSP solvers---") do
# create vectors and matrices


  b = PetscVec(comm);
  PetscVecSetType(b, "mpi");
  PetscVecSetSizes(b,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(b)
  b_global_indices = Array(low:PetscInt(high - 1))
 

  x = PetscVec(comm);
  PetscVecSetType(x,"mpi");
  PetscVecSetSizes(x,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(x)
  x_global_indices = Array(low:PetscInt(high - 1))
 

  for i=1:sys_size
    idxm = [ b_global_indices[i] ]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)



  A = PetscMat(comm)
  PetscMatSetType(A, "mpiaij")

  PetscMatSetSizes(A,sys_size,sys_size,PetscInt(comm_size*sys_size),PetscInt(comm_size*sys_size));
  PetscSetUp(A);

  low, high = PetscMatGetOwnershipRange(A)
  mat_global_indices = Array(low:PetscInt(high - 1))
 

  for i=1:sys_size
    for j = 1:sys_size
      idxm = [ mat_global_indices[i] ]  # row index
      idxn = [ mat_global_indices[j] ]  # column index
      PetscMatSetValues(A,idxm, idxn, [A_julia[i,j]],PETSC_INSERT_VALUES);
    end
  end

  PetscMatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
  PetscMatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);



# perform solve
ksp = KSP(comm)
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)
KSPSolve(ksp, b, x)


# copy solution back to Julia
x_copy = zeros(PetscScalar, sys_size)
#idx = Array(0:2)
#idx = zeros(PetscInt, sys_size)
#for i=1:sys_size
#  idx[i] = i-1
#end

PetscVecGetValues(x, sys_size, x_global_indices, x_copy)
println("x_copy = ", x_copy)
println("x_julia = ", x_julia)
for i=1:sys_size
    @fact x_copy[i] => roughly(x_julia[i], atol=1e-14)
end

#PetscView(x, 0)

PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(A)
PetscDestroy(ksp)
end
