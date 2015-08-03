# test Petsc matrix functions
facts("\n   ---testing matrix functions---") do 
  A = PetscMat(comm)
  PetscMatSetType(A, "mpiaij")

  PetscMatSetSizes(A,sys_size,sys_size, PetscInt(comm_size*sys_size),PetscInt(comm_size*sys_size));
  PetscSetUp(A);


  B = PetscMat(comm)
  PetscMatSetType(B, "mpiaij")

  PetscMatSetSizes(B,sys_size,sys_size, PetscInt(comm_size*sys_size),PetscInt(comm_size*sys_size));
  PetscSetUp(B);

  x = PetscVec(comm);
  PetscVecSetType(x, VECMPI);
  PetscVecSetSizes(x,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(x)
  global_indices = Array(low:PetscInt(high - 1))
 

  for i=1:sys_size
    idxm = [global_indices[i]]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(x, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(x)
  PetscVecAssemblyEnd(x)

  xvec_julia = deepcopy(rhs)


  y = PetscVec(comm);
  PetscVecSetType(y, VECMPI);
  PetscVecSetSizes(y,sys_size, PetscInt(comm_size*sys_size));

  for i=1:sys_size
    idxm = [global_indices[i]]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(y, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(y)
  PetscVecAssemblyEnd(y)

  y_julia = deepcopy(rhs)


  low, high = PetscMatGetOwnershipRange(A)
  global_indices = Array(low:PetscInt(high - 1))
  println("comm_rank = ", comm_rank, " , global_indices = ", global_indices)

  global_row_ind = zeros(PetscInt, sys_size*sys_size)
  global_col_ind = zeros(PetscInt, sys_size*sys_size)


  # get global row and column indicies, in column major order
  pos = 1
  for i=1:sys_size  # loop over columns
    for j=1:sys_size  # loop over rows
      global_row_ind[pos] = j - 1 + low
      global_col_ind[pos] = i - 1 + low
      pos += 1
    end
  end

  println("global_row_ind = ", global_row_ind)
  println("global_col_ind = ", global_col_ind)
  B_julia = A_julia + PetscScalar(1)



  for i=1:sys_size
    for j = 1:sys_size
      idxm = [ global_indices[i] ] # row index
      idxn = [ global_indices[j] ] # column index
      PetscMatSetValues(A,idxm, idxn, [A_julia[i,j]],PETSC_INSERT_VALUES);
      PetscMatSetValues(B,idxm, idxn, [A_julia[i,j] + PetscScalar(1)],PETSC_INSERT_VALUES);
    end
  end

  B_copy = zeros(PetscScalar, sys_size, sys_size)

  PetscMatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
  PetscMatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);

  PetscMatAssemblyBegin(B,PETSC_MAT_FINAL_ASSEMBLY);
  PetscMatAssemblyEnd(B,PETSC_MAT_FINAL_ASSEMBLY);



  for i=1:sys_size
    for j=1:sys_size
      idxm = [global_indices[i] ]  # row index
      idxn = [ global_indices[j] ] # column index
      v = zeros(PetscScalar, 1,1)
      PetscMatGetValues(A, idxm, idxn, v)
#      println("i = ", i, " j = ", j, " v = ", v[1,1], " A[i,j] = ", A_julia[i,j])
      @fact v[1,1] => roughly(A_julia[i,j]) "mismatch at i=$i, j=$j"
    end
  end


  @fact PetscMatGetSize(A) => (comm_size*sys_size, comm_size*sys_size)
  println("Mat local size = ", PetscMatGetLocalSize)
  @fact PetscMatGetLocalSize(A) => (sys_size, sys_size)
  # testing non zero exist code is all we can really do here
  println("Printing A")
  @fact PetscView(A, 0) => 0
  println("Printing B")
  PetscView(B, 0)
 

  alpha = PetscScalar(2.3)


  PetscMatAXPY(B, alpha, A, DIFFERENT_NONZERO_PATTERN)
  B_julia = alpha*A_julia + B_julia
  PetscMatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size
    for j=1:sys_size
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  println("finished testing MatAXPY")

  PetscMatAYPX(B, alpha, A, DIFFERENT_NONZERO_PATTERN)
  B_julia = alpha*B_julia + A_julia
  PetscMatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size
    for j=1:sys_size
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  println("finished testing MatAYPX")

  PetscMatScale(B, alpha)
  B_julia = alpha*B_julia
  PetscMatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size
    for j=1:sys_size
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  println("finished testing MatScale")

  PetscMatShift(B, alpha)
  B_julia += alpha*eye(PetscScalar, sys_size)
  PetscMatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size
    for j=1:sys_size
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  println("finished testing PetscMatShift")

  y_copy = zeros(PetscScalar, sys_size)
  PetscMatMult(B, x, y)
  y_julia = B_julia*xvec_julia

  PetscVecGetValues(y, sys_size, global_indices, y_copy)
  for i=1:sys_size
      @fact y_julia[i] => roughly(y_copy[i])
  end


  PetscMatMultAdd(B, x, y, y)
  y_julia = y_julia + B_julia*xvec_julia
  PetscVecGetValues(y, sys_size, global_indices, y_copy)
  for i=1:sys_size
      @fact y_julia[i] => roughly(y_copy[i])
  end

  PetscMatMultTranspose(A, x, y)
  y_julia = A_julia.'*xvec_julia
  PetscVecGetValues(y, sys_size, global_indices, y_copy)
  for i=1:sys_size
      @fact y_julia[i] => roughly(y_copy[i])
  end

  PetscMatMultHermitianTranspose(A, x, y)
  y_julia = A_julia'*xvec_julia
  PetscVecGetValues(y, sys_size, global_indices, y_copy)
  for i=1:sys_size
      @fact y_julia[i] => roughly(y_copy[i])
  end








  @fact PetscDestroy(A) => 0
  PetscDestroy(B)
  PetscDestroy(y)
  PetscDestroy(x)
end
