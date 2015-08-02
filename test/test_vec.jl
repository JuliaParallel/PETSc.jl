
facts("\n   ---Testing vector functions---") do

for j=1:length(vec_formats)
#  format_j = "standard"
  format_j = vec_formats[j]
  println("testing vector format ", format_j)
  b = PetscVec(comm);
  PetscVecSetType(b, format_j);
  PetscVecSetSizes(b,sys_size, PetscInt(comm_size*sys_size));


  # create 3rd vector to store results in
  b3 = PetscVec(comm);
  PetscVecSetType(b3, format_j);
  PetscVecSetSizes(b3,sys_size, PetscInt(comm_size*sys_size));

  rhs_tmp2 = zeros(PetscScalar, sys_size)

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
  PetscVecGetValues(b2, sys_size, global_indices, b2_copy)
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

  petsc_r, petsc_idx = PetscVecMax(b)
  julia_r = maximum(real(b_copy))
  julia_idx = indmax(real(b_copy))

  println("petsc_r = ", petsc_r, ", julia_r = ", julia_r)
  println("petsc_idx = ", petsc_idx, ", julia_idx = ", julia_idx)

  @fact petsc_r => roughly(julia_r)
  @fact petsc_idx => (julia_idx - 1)

  println("finished testing VecMax")


  petsc_r, petsc_idx = PetscVecMin(b)
  julia_r = minimum(real(b_copy))
  julia_idx = indmin(real(b_copy))

  @fact petsc_r => roughly(julia_r)
  @fact petsc_idx => (julia_idx - 1)

  println("finished testing VecMin")

  PetscVecReciprocal(b)
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = 1/rhs_tmp[i]
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end
  println("finished testing VecReciprocal")

  PetscVecShift(b, PetscScalar(2.0))
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp[i] = rhs_tmp[i] + 2.0
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end
  println("finished testing VecShift")

  PetscVecPointwiseMult(b3, b, b2)
  PetscVecGetValues(b3, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp2[i] = rhs_tmp[i]*rhs[i]
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end

  println("finished testing VecPointwiseMult")

  PetscVecPointwiseDivide(b3, b, b2)
  PetscVecGetValues(b3, sys_size, global_indices, b_copy)
  for i=1:sys_size
    rhs_tmp2[i] = rhs_tmp[i]/rhs[i]
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end

  println("finished testing VecPointwiseDivide")






  # test vector multiplication, addition functions
  # use b, b2 as the Petsc vectors
  # use rhs_tmp, rhs as the julia vectors

  alpha = PetscScalar(2.3)
  beta = PetscScalar(3.5)
  gamma = PetscScalar(4.8)

  println("alpha = ", alpha, ", beta = ", beta, ", gamma = ", gamma)

  PetscView(b)
  PetscView(b2)

  PetscVecAXPY(b, alpha, b2)
  rhs_tmp += alpha*rhs
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  println("b_copy = ", b_copy)
  println("rhs_tmp = ", rhs_tmp)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing AXPY")

  #=
  PetscVecAXPBY(b, alpha, beta, b2)
  rhs_tmp = alpha*rhs + beta*rhs_tmp
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end
=#

  PetscVecAYPX(b, alpha, b2)
  rhs_tmp = alpha*rhs_tmp + rhs
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing AYPX")

  PetscVecWAXPY(b3, alpha, b, b2)
  rhs_tmp2 = alpha*rhs_tmp + rhs
  PetscVecGetValues(b3, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end

  println("finished testing WAXPY")

  # MAXPY
  scalar_arr = [alpha, beta]
  vec_arr = [b2.pobj, b3.pobj]
  PetscVecMAXPY(b, PetscInt(2), scalar_arr, vec_arr)
  rhs_tmp += alpha*rhs + beta*rhs_tmp2
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing MAXPY")

  PetscVecAXPBYPCZ(b, alpha, beta, gamma, b2, b3)
  rhs_tmp = alpha*rhs + beta*rhs_tmp2 + gamma*rhs_tmp
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing AXPBYPCZ")



  PetscVecScale(b, alpha)
  rhs_tmp *= alpha
   PetscVecGetValues(b, sys_size, global_indices, b_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing VecScale")

  PetscView(b)
  PetscView(b2)

  petsc_r = PetscVecDot(b, b2)
  julia_r = comm_size*rhs_tmp'*rhs

  println("petsc_r = ", petsc_r)
  println("julia_r = ", julia_r)

  @fact petsc_r => roughly(julia_r[1])

  println("finished testing VecDot")
#=
  petsc_r = PetscVecTDot(b, b2)
  julia_r = rhs_tmp.'*rhs

  @fact petsc_r => roughly(julia_r[1])

  println("finished testing PetscVecTDot")

  petsc_r =  PetscVecSum(b2)
  julia_r = sum(rhs_global)

  @fact petsc_r => roughly(julia_r)
  
  println("finished testing VecSum")

  PetscVecSwap(b, b2)
  b2_copy = zeros(PetscScalar, sys_size)
  PetscVecGetValues(b, sys_size, global_indices, b_copy)
  PetscVecGetValues(b2, sys_size, global_indices, b2_copy)
  for i=1:sys_size
    @fact b_copy[i] => roughly(rhs[i])
    @fact b2_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing VecSwap")

  PetscVecSwap(b, b2)  # swap back

=#

 
  
  PetscDestroy(b)
  PetscDestroy(b2)
  PetscDestroy(b3)
  end


end # end fact check block



