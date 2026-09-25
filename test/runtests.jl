using Test
using MPI: MPI, mpiexec
using PETSc, Pkg

# Make sure that all dependencies are installed also on a clean system
Pkg.instantiate()

# When set_library! has been used, petsc_library is a path string and PETSc_jll
# is not loaded.  Only import JLL-specific symbols when using the default binaries.
const USING_CUSTOM_LIB = PETSc.petsclibs[1].petsc_library isa AbstractString

import MPIPreferences  # always import to ensure local MPI API is configured

if USING_CUSTOM_LIB
    @info "Testing PETSc.jl with custom library" path=PETSc.petsclibs[1].petsc_library MPIPreferences.binary MPIPreferences.abi
else
    using PETSc_jll
    @info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi PETSc_jll.host_platform
end

# Do the MPI tests first so we do not have mpi running inside MPI
mpi_tests = ("mpivec.jl", "mpimat.jl", "ksp.jl", "dmstag.jl", "mpi_blas_threads.jl")

# PETSc_jll >= 3.25.4 has MPI-enabled Windows binaries (MicrosoftMPI), so the MPI tests run everywhere
do_mpi = true

include("init.jl")
include("blas_threads.jl")   # -blas_num_threads sizes the BLAS pool PETSc runs in
include("lib.jl")
include("vec.jl")           # autowrapped
include("mat.jl")           # autowrapped
include("options.jl")       # autowrapped
include("ksp.jl")           # autowrapped
include("pc.jl")            # high-level PC
include("mat_vec_methods.jl")  # copyto!, fill!, zero_rows!, set_option!, diagonal!, isassembled
include("snes.jl")          # autowrapped
include("snes_helpers.jl")  # small helper tests for SNES return-style wrappers
include("bang_returns.jl")   # ! functions return the object they mutate
include("solver_api.jl")     # KSP/SNES/TS/PC readers, operators, prefixes, set_from_options!
include("ts.jl")            # high-level TS interface
include("ts_long_runs.jl")  # equation type, step number, pre/post-step hooks, failure limits
include("dmda.jl")          # autowrapped
include("dmstag.jl")        # autowrapped
include("dmplex.jl")        # DMPlex tests
include("dmnetwork.jl")     # new test for DMNetwork example
include("dmshell.jl")       # new test for DMShell example
include("dmproduct.jl")     # test for DMProduct example
include("matshell.jl")      # autowrapped!
include("test_dmstag.jl")
include("dm_dimension.jl")  # dimension-correct DM returns: shape, inference, allocations
include("test_snes.jl")
include("test_audit.jl")    # leak auditor
include("test_errors.jl")   # argument validation
include("wrapper_signatures.jl")  # wrapper arguments take the abstract types
include("wrapper_quality.jl")     # every generated method infers a concrete return type; no allocations
include("wrapper_leaks.jl")       # repeated create/destroy and Get/Restore pairs do not grow memory
@testset "method ambiguities do not grow" begin
    # 130 with the regenerated wrappers (all in the high-level layer); regenerate or rename
    # without adding new ones
    @test length(detect_ambiguities(PETSc; recursive = true)) <= 130
end
include("old_test.jl")
include("low_level_viewer.jl")  # Low-level viewer convenience functions
include("low_level_ts.jl")      # Low-level TS functions
include("ts_ex51.jl")           # Regression test for repeated ex51 solves
include("ts_ex51_implicit.jl")  # Regression test for repeated implicit Gauss solves
include("ts_ex16.jl")           # Regression test for the van der Pol IMEX example
include("ts_scenarios.jl")      # High-level TS: ex16 equivalence + a damping sweep
include("low_level_is.jl")      # Low-level IS functions
include("low_level_petscsection.jl")  # Low-level PetscSection functions
include("low_level_petscsf.jl")      # Low-level PetscSF graph and communication functions
include("petscbool.jl")              # PetscBool is one byte (PETSc >= 3.24)
include("low_level_tao.jl")     # Low-level Tao functions
include("test_destroy.jl")      # destroy! guards: stale cycle, double destroy
include("handles.jl")           # destroy! on IS, AO, PF, Tao; ISColoringGetIS ownership
include("lifetimes.jl")         # arrays a Vec or Mat wraps live as long as it does
include("test_deprecations.jl") # the v0.4 names still work and warn (the only file calling them)
include("test_api_surface.jl")  # scripts/api_surface.jl --check: the register covers the surface

include("testutils.jl")

# Run helper tests for SNES and TAO
include("snes_helpers.jl")
include("tao_helpers.jl")

# Run the examples to make sure they all work
include("examples.jl")

# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
# XXX: Currently not working on windows (since we have no PETSc + MPI)
if do_mpi
    include("mpi_examples.jl")
end

# Do the MPI tests
# XXX: Currently not working on windows (since we have no PETSc + MPI)
if do_mpi
    @testset "MPI Tests" begin
        for testfile in mpi_tests
            testpath = joinpath(@__DIR__, testfile)
            cmd = `$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=. $testpath`
            @test success(pipeline(cmd, stderr = stderr))
        end
    end
end
