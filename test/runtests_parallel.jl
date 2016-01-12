# run tests in parallel
include("runtests_setup.jl")

for ST in PETSc.C.petsc_type
  # @testset "Scalar type $ST" begin # uncomment when nested test results can be printed
  include("vecp.jl")
  # end
end


