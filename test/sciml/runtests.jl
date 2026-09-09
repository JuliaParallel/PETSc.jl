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
    include("test_api.jl")
    include("test_regressions.jl")

    # Regression guard: abandoned integrators must not crash at process exit
    # when their GC finalizers fire after `PetscFinalize` (PETSC_ERR_MPI / MPICH
    # abort). Run in a subprocess and assert a clean exit with no finalizer
    # banners. Uses the active project so the child inherits SciMLBase.
    @testset "finalizer safety at process exit" begin
        script = joinpath(@__DIR__, "finalizer_safety.jl")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $script`
        out = IOBuffer()
        ok = success(pipeline(cmd; stdout = out, stderr = out))
        log = String(take!(out))
        @test ok
        @test occursin("finalizer_safety: reached clean exit", log)
        @test !occursin("error in running finalizer", log)
        @test !occursin("after finalizing", log)
    end

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
