
facts("Testing vector function") do

for j=1:length(vec_formats)

  format_j = vec_formats[j]
  println("testing vector format ", format_j)
  b = PetscVec(comm);
  PetscVecSetType(b, format_j);
  PetscVecSetSizes(b,sys_size, comm_size*sys_size);

  #=
  x = PetscVec(comm);
  PetscVecSetType(x,"mpi");
  PetscVecSetSizes(x,sys_size, comm_size*sys_size);
  =#

  for i=1:sys_size
    idxm = [(comm_rank)*sys_size + i]   # index
    val = [ rhs[i] ]  # value
    PetscVecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)

  # check that the vector was set/assembled correctly
  b_copy = zeros(sys_size)
  idx = Array(0:2)  

  PetscVecGetValues(b, sys_size, idx, b_copy)
  for i=1:sys_size
     @fact b_copy[i] => roughly(rhs[i])
  end

  @fact PetscVecGetSize(b) => sys_size

  julia_vec_norms = [1, 2, Inf]
  for i=1:length(vec_norms)
    @fact PetscVecNorm(b, vec_norms[i]) => roughly( norm(rhs, julia_vec_norms[i]))
  end
  
  # check for non zero exit status is all we can really do here 
  @fact PetscView(b, 0) => 0
  
  PetscDestroy(b)

  end


end # end fact check block
