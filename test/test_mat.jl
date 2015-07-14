# test Petsc matrix functions
facts("testing matrix functions") do 
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

  for i=1:sys_size
    for j=1:sys_size
      idxm = [convert(PetscInt, (comm_rank)*sys_size + i)]  # row index
      idxn = [convert(PetscInt, (comm_rank)*sys_size + j)]  # column index
      v = zeros(1,1)
      MatGetValues(A, idxm, idxn, v)
      println("i = ", i, " j = ", j, " v = ", v[1,1], " A[i,j] = ", A_julia[i,j])
      @fact v[1,1] => roughly(A_julia[i,j]) "mismatch at i=$i, j=$j"
    end
  end


  @fact PetscMatGetSize(A) => (sys_size, sys_size)

  # testing non zero exist code is all we can really do here
  @fact PetscView(A, 0) => 0
     
  @fact PetscDestroy(A) => 0 
end
