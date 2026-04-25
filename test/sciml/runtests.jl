using Test
using PETSc
using SciMLBase
using DiffEqBase

@testset "PETSc.jl SciML extension" begin
    ext = Base.get_extension(PETSc, :PETScSciMLExt)
    @test ext !== nothing
    # Algorithm-level test files will be added in Steps 2–8:
    #   include("test_rk.jl")
    #   include("test_rosenbrock.jl")
    #   include("test_implicit.jl")
    #   include("test_imex.jl")
    #   include("test_output.jl")
    #   include("test_callbacks.jl")
    #   include("test_integrator.jl")
end
