
facts("\n   ---Testing vector functions---") do

for j=1:length(vec_formats)

  format_j = vec_formats[j]
  println("testing vector format ", format_j)
  b = PetscVec(comm);
  PetscVecSetType(b, format_j);
  PetscVecSetSizes(b,sys_size, PetscInt(comm_size*sys_size));

  low, high = PetscVecGetOwnershipRange(b)
  global_indices = Array(low:PetscInt(high - 1))
  println("comm_rank = ", comm_rank, " , global_indices = ", global_indices)


  #=
  x = PetscVec(comm);
  PetscVecSetType(x,"mpi");
  PetscVecSetSizes(x,sys_size, comm_size*sys_size);
  =#

  for i=1:sys_size
    idxm = [global_indices[i]]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)

  # check that the vector was set/assembled correctly
  # check all the methods of copying/accesing a vector work
  b_copy = zeros(PetscScalar, sys_size)
  b2_copy = zeros(PetscScalar, sys_size)
#  idx = Array(0:2)  
#  idx = Array(PetscInt, 3)
#  for i=1:sys_size
#    idx[i] = global_indices[i]
#  end

  b_arr, ptr_arr = PetscVecGetArray(b)
  b_arr_ro, ptr_arr2 = PetscVecGetArrayRead(b)
  b2 = PetscVecDuplicate(b)
  PetscVecCopy(b, b2)
  
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  PetscVecGetValues(b, sys_size, global_indices, b2_copy)
  for i=1:sys_size
     @fact b_copy[i] => roughly(rhs[i])
     @fact b_arr[i] => roughly(rhs[i])
     @fact b_arr_ro[i] => roughly(rhs[i])
     @fact b2_copy[i] => roughly(rhs[i])
  end

  PetscVecRestoreArray(b, ptr_arr)
  PetscVecRestoreArrayRead(b, ptr_arr2)

  

  @fact PetscVecGetSize(b) => comm_size*sys_size

  for i=1:length(vec_norms)
    norm_petsc = PetscVecNorm(b, vec_norms[i])
    norm_julia = norm(rhs, julia_vec_norms[i])
    println("petsc_norm = ", norm_petsc, ", julia_norm = ", norm_julia)
    println("petsc norm type: ", vec_norms[i], " , julia norm type: ", julia_vec_norms[i])
    @fact PetscVecNorm(b, vec_norms[i]) => roughly( norm(rhs_global, julia_vec_norms[i]))
    print("\n")
  end
  
  # check for non zero exit status is all we can really do here 
  @fact PetscView(b, 0) => 0


  # test math functions
  PetscVecSqrtAbs(b)
  rhs_tmp = deepcopy(rhs)  # don't modifiy original rhs
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = sqrt(abs(rhs_tmp[i]))
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  PetscVecLog(b)
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = log(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  PetscVecExp(b)
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = exp(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  PetscVecAbs(b)
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = abs(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end
  
  
  PetscDestroy(b)
  PetscDestroy(b2)
  end


end # end fact check block



