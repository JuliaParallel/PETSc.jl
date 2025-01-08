using Test
using MPI: MPI, mpiexec
using PETSc, PETSc_jll


import MPIPreferences
@info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi PETSc_jll.host_platform

do_mpi = true
if Sys.iswindows()
    do_mpi = false
end

# Do the MPI tests first so we do not have mpi running inside MPI
# XXX: Currently not working on windows (since we have no PETSc + MPI)
#=
if do_mpi
    cmd = `$(mpiexec())  -n 4 $(Base.julia_cmd()) --project dmda.jl`
    run(cmd)
    success(pipeline(cmd, stderr = stderr))
end
=#

# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
# XXX: Currently not working on windows (since we have no PETSc + MPI)
if do_mpi
#    include("mpi_examples.jl")
end

include("options.jl")
include("vec.jl")           # not yet autowrapped!
include("mat.jl")           # not yet autowrapped!
include("matshell.jl")      # not yet autowrapped!
include("dmda.jl")          # not yet autowrapped!
include("dmstag.jl")        # mostly autowrapped
##include("test_dmstag.jl")    # "old" dmstag tests - need to be finalized ; also needs KSP to run
include("ksp.jl")            # not yet autowrapped!
include("snes.jl")          # not yet autowrapped!
include("old_test.jl")

#=

# Run the examples to make sure they all work
include("examples.jl")
=#
