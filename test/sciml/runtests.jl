using Test
using PETSc
using SciMLBase

@testset "PETSc SciML extension" begin
    # Smoke test: extension activates with just SciMLBase loaded.
    @testset "Extension activation" begin
        ext = Base.get_extension(PETSc, :PETScSciMLExt)
        @test ext !== nothing
    end

    include("test_rk.jl")
    include("test_float32.jl")
    include("test_ex51_sciml.jl")
    include("test_ex16_sciml.jl")
    include("test_rosenbrock.jl")
    include("test_implicit.jl")
    include("test_imex.jl")
    include("test_output.jl")
    include("test_callbacks.jl")
    include("test_integrator.jl")
    include("test_polish.jl")
    include("test_review_fixes.jl")

    # Distributed (MPI) explicit integration runs in a separate `mpiexec`
    # subprocess. Launch it with the *currently active* project so the child
    # inherits SciMLBase (a test-only dependency not present in PETSc's main
    # `[deps]`). Skipped on Windows, where the PETSc + MPI stack is unavailable.
    if !Sys.iswindows()
        @testset "MPI explicit integration (2 ranks)" begin
            using MPI: mpiexec
            script = joinpath(@__DIR__, "mpi_sciml.jl")
            cmd = `$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$(Base.active_project()) $script`
            @test success(pipeline(cmd; stderr = stderr))
        end
    end
end
