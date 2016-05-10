include("runtests_setup.jl")
println("testing types: ", PETSc.C.petsc_type)
for (i, ST) in enumerate(PETSc.C.petsc_type)
  if PETSc.have_petsc[i]
    println("testing datatype ", ST)
  # @testset "Scalar type $ST" begin # uncomment when nested test results can be printed
    include("error.jl")
    include("ksp.jl")
    include("vec.jl")
    include("is.jl")
    include("mat.jl")
    include("ts.jl")
  end
  # end
end

@test PETSc.petsc_sizeof(PETSc.C.PETSC_BOOL) == 4


