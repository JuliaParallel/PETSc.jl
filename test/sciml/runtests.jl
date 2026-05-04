using Test
using Pkg

# Mirror `test/runtests.jl`'s project activation so test-only deps such as
# SciMLBase are on the load path even when this file is `include`d directly
# from a session whose active project is the package itself. Note: invoking
# the file as a standalone script (`julia --project=. test/sciml/runtests.jl`)
# is *not* a fully supported entry point — the process may exit non-zero
# during PETSc/MPI teardown. Use `Pkg.test("PETSc")` or `include` the file
# from `test/runtests.jl` instead.
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
    include("test_ex51_sciml.jl")
    include("test_rosenbrock.jl")
    include("test_implicit.jl")
    include("test_imex.jl")
    include("test_output.jl")
    include("test_callbacks.jl")
    include("test_integrator.jl")
    include("test_polish.jl")
    include("test_review_fixes.jl")
end
