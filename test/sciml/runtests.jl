using Test
using Pkg

# When this file is invoked directly (`julia --project=. test/sciml/runtests.jl`)
# the active project is the package and SciMLBase is not on the load path. The
# top-level `test/runtests.jl` performs the same activation; doing it here too
# keeps the standalone script entrypoint working.
if !haskey(Pkg.project().dependencies, "SciMLBase")
    Pkg.activate(joinpath(@__DIR__, ".."))
end

using PETSc
using SciMLBase
using DiffEqBase

@testset "PETSc SciML extension" begin
    # Smoke test: extension activates with just SciMLBase + DiffEqBase loaded.
    @testset "Extension activation" begin
        ext = Base.get_extension(PETSc, :PETScSciMLExt)
        @test ext !== nothing
    end

    include("test_rk.jl")
    include("test_rosenbrock.jl")
    include("test_implicit.jl")
    include("test_imex.jl")
    include("test_output.jl")
    include("test_callbacks.jl")
    include("test_integrator.jl")
    include("test_polish.jl")
    include("test_review_fixes.jl")
end
