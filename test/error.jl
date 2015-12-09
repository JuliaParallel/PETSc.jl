@testset "Error Handling{$ST}" begin
  msg = PETSc.PetscErrorMessage(73)
  @test msg == "Object is in wrong state"
  @test_throws PETSc.PetscError PETSc.chk(76)
end
