facts("--- Testing Error Handling ---") do
  msg = PETSc.PetscErrorMessage(73)
  @fact msg --> "Object is in wrong state"
end
