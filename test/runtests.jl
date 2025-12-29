using Test
using MPI: MPI, mpiexec
using PETSc, PETSc_jll


import MPIPreferences
@info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi PETSc_jll.host_platform

# Do the MPI tests first so we do not have mpi running inside MPI
mpi_tests = ("mpivec.jl", "mpimat.jl", "ksp.jl", "dmstag.jl")

do_mpi = true
if Sys.iswindows()
    do_mpi = false
end




include("vec.jl")           # autowrapped
include("mat.jl")           # autowrapped
include("options.jl")       # autowrapped
include("ksp.jl")           # autowrapped
include("snes.jl")          # autowrapped
include("dmda.jl")          # autowrapped
include("dmstag.jl")        # autowrapped
include("matshell.jl")      # autowrapped!
#include("test_dmstag.jl")   # "old" dmstag tests - need to be finalized ; also needs KSP to run
include("old_test.jl")

#=

# Run the examples to make sure they all work
include("examples.jl")
=#


# Do the MPI tests
# XXX: Currently not working on windows (since we have no PETSc + MPI)
if do_mpi
    @testset "MPI Tests" begin
        for testfile in mpi_tests
            cmd = `$(mpiexec())  -n 4 $(Base.julia_cmd()) --project $testfile`
            run(cmd)
            success(pipeline(cmd, stderr = stderr))
        end
    end
end

# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
# XXX: Currently not working on windows (since we have no PETSc + MPI)
if do_mpi
#    include("mpi_examples.jl")
end

