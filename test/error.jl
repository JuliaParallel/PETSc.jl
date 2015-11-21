facts("--- Testing Error Handling ---") do
  msg = PETSc.PetscErrorMessage(73)
  @fact msg --> "Object is in wrong state"
  @fact_throws PETSc.PetscError PETSc.chk(76)
end
