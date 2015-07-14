facts("testing KSP solvers") do
# create vectors and matrices


  b = PetscVec(comm);
  PetscVecSetType(b, "mpi");
  PetscVecSetSizes(b,sys_size, comm_size*sys_size);

  x = PetscVec(comm);
  PetscVecSetType(x,"mpi");
  PetscVecSetSizes(x,sys_size, comm_size*sys_size);

  for i=1:sys_size
    idxm = [(comm_rank)*sys_size + i]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)



  A = PetscMat(comm)
  PetscMatSetType(A, "mpiaij")

  PetscMatSetSizes(A,sys_size,sys_size,comm_size*sys_size,comm_size*sys_size);
  PetscSetUp(A);

  for i=1:sys_size
    for j = 1:sys_size
      idxm = [(comm_rank)*sys_size + i]  # row index
      idxn = [(comm_rank)*sys_size + j]  # column index
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
x_copy = zeros(sys_size)
idx = Array(0:2)
PetscVecGetValues(x, sys_size, idx, x_copy)

for i=1:sys_size
    @fact x_copy[i] => roughly(x_julia[i])
end

end
