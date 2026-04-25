using Test
using MPI: MPI, mpiexec
using PETSc, Pkg

# Make sure that all dependencies are installed also on a clean system
Pkg.instantiate()

# When set_library! has been used, petsc_library is a path string and PETSc_jll
# is not loaded.  Only import JLL-specific symbols when using the default binaries.
const _using_custom_lib = PETSc.petsclibs[1].petsc_library isa AbstractString

if _using_custom_lib
    @info "Testing PETSc.jl with custom library" path=PETSc.petsclibs[1].petsc_library
else
    using PETSc_jll
    import MPIPreferences
    @info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi PETSc_jll.host_platform
end

# Do the MPI tests first so we do not have mpi running inside MPI
mpi_tests = ("mpivec.jl", "mpimat.jl", "ksp.jl", "dmstag.jl")

do_mpi = true
if Sys.iswindows()
    do_mpi = false
end

include("init.jl")
include("lib.jl")
include("vec.jl")           # autowrapped
include("mat.jl")           # autowrapped
include("options.jl")       # autowrapped
include("ksp.jl")           # autowrapped
include("snes.jl")          # autowrapped
include("snes_helpers.jl")  # small helper tests for SNES return-style wrappers
include("dmda.jl")          # autowrapped
include("dmstag.jl")        # autowrapped
include("dmnetwork.jl")     # new test for DMNetwork example
include("dmshell.jl")       # new test for DMShell example
include("dmproduct.jl")     # test for DMProduct example
include("matshell.jl")      # autowrapped!
include("test_dmstag.jl")
include("test_snes.jl")
include("old_test.jl")
include("low_level_viewer.jl")  # Low-level viewer convenience functions
include("low_level_ts.jl")      # Low-level TS functions
include("ts_ex51.jl")           # Regression test for repeated ex51 solves
include("ts_ex51_implicit.jl")  # Regression test for repeated implicit Gauss solves
include("ts_ex16.jl")           # Regression test for the van der Pol IMEX example
include("low_level_is.jl")      # Low-level IS functions
include("low_level_petscsection.jl")  # Low-level PetscSection functions
include("low_level_tao.jl")     # Low-level Tao functions

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
