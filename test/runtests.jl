using Test
using MPI


do_mpi = true
#if Sys.iswindows()
#    do_mpi = false
#end


# Do the MPI tests first so we do not have mpi running inside MPI
# XXX: Currently not working on windows, not sure why
#if do_mpi
#    cmd = `$(mpiexec())  -n 4 $(Base.julia_cmd()) --project dmda.jl`
#    run(cmd)
#    success(pipeline(cmd, stderr = stderr))
#end

# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
# XXX: Currently not working on windows reliably, not sure why
#if do_mpi
#    include("mpi_examples.jl")
#end

include("options.jl")
#include("dmda.jl")
#include("old_test.jl")
#include("test_dmstag.jl")
include("test_snes.jl")

# Run the examples to make sure they all work
include("examples.jl")
