using FactCheck
using PETSc

facts("testing double precision real w/ 32 bit indices") do

  @fact PetscScalar => Float64
  @fact PetscReal => Float64
  @fact PetscInt => Int32

end

FactCheck.exitstatus()
