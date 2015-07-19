using FactCheck
using PETSc
facts("  ---Checking Petsc data types---") do

  @fact PetscScalar => Float32
  @fact PetscReal => Float32
  @fact PetscInt => Int64

end

FactCheck.exitstatus()
