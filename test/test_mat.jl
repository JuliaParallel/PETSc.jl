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


  D = PetscMat(comm)
  PetscMatMatMult(A, B, MAT_INITIAL_MATRIX, PetscReal(1.0), D)
  D_julia = A_julia*B_julia
  println("D_julia = ", D_julia)
  D_copy = zeros(PetscScalar, sys_size, sys_size)
  PetscMatGetValues(D, global_indices, global_indices, D_copy)

  for i=1:sys_size
    for j=1:sys_size
      @fact D_copy[j,i] => roughly(D_julia[i,j])
    end
  end

  # test matrix norms
  fnorm = PetscMatNorm(D, NORM_FROBENIUS)
  infnorm = PetscMatNorm(D, NORM_INFINITY)
  
  @fact fnorm => roughly(sqrt(comm_size)*vecnorm(D_julia))
  @fact infnorm => roughly(norm(D_julia, Inf))


  println("testing preallocation")
  nb = PetscInt(100)  # number of times/blocks to insert
  C = PetscMat(comm)
  PetscMatSetType(C, "mpiaij")
#  PetscMatSetFromOptions(C)
  println("nb = ", nb)
  println("comm_size = ", comm_size)
  println("sys_size = ", sys_size)

  PetscMatSetSizes(C, nb*sys_size, nb*sys_size, PetscInt(nb*comm_size*sys_size),PetscInt(nb*comm_size*sys_size));

  # preallocation parameters
  bs = PetscInt(1)
  dnnz = 3*ones(PetscInt, nb*sys_size)  # on diagonal (row + column owned by this process)
  onnz = zeros(PetscInt, nb*sys_size)  # no off diagonal (column not owned by this process)
  dnnzu = Array(PetscInt, 0)  # this is not a symmetric matrix, so unused
  onnzu = Array(PetscInt, 0)  # this is not a symmetric matrix, so unused

  println("dnnz = ", dnnz)
  PetscMatXAIJSetPreallocation(C, bs, dnnz, onnz, dnnzu, onnzu)

#  PetscSetUp(C);


  println("A_julia = ", A_julia)
  println("typeof(A_julia) = ", typeof(A_julia))

  A_julia_t = A_julia.'  # transpose because C is row major

  low, high = PetscMatGetOwnershipRange(C)
  idi = zeros(PetscInt, sys_size)  # row indices
  idj = zeros(PetscInt, sys_size)  # column indices

  for i=1:nb

    # get indices
    pos = 1
    for j=1:sys_size # insert block on diagonal
	idi[pos] = low + sys_size*(i - 1) + j - 1
	idj[pos] = low + sys_size*(i - 1) + j - 1
	pos += 1
    end

    # add a random component
    for i=1:sys_size
      for j=1:sys_size
	A_julia_t[i,j] += rand()
      end
    end

#    println("idi = ", idi)
#    println("idj = ", idj)
    PetscMatSetValues(C, idi, idj, A_julia_t, PETSC_INSERT_VALUES)
  end


  PetscMatAssemblyBegin(C, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(C, PETSC_MAT_FINAL_ASSEMBLY)

  matinfo = PetscMatGetInfo(C, MAT_LOCAL)

  @fact matinfo.mallocs => roughly(0.0)

#  PetscView(C, 0)

  println("finished testing preallocation")



  @fact PetscDestroy(A) => 0
  PetscDestroy(B)
  PetscDestroy(C)
  PetscDestroy(D)
  PetscDestroy(y)
  PetscDestroy(x)
end
