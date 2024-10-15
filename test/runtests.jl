using Test
using MPI: MPI, mpiexec

import MPIPreferences
@info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi

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


@testset   "MPI" begin
    using MPI
    MPI.install_mpiexecjl(force=true)

    # Do a dummy `@test true`:
    # If the process errors out the testset would error out as well,
    # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
    @test true

    @info "Starting parallel tests"

    run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --threads=1 --check-bounds=yes --project=$(dirname(@__DIR__)) $(abspath("01-hello.jl"))`)

    run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --threads=1 --check-bounds=yes --project=$(dirname(@__DIR__)) $(abspath("dmda.jl"))`)
   # cmd = `$(mpiexec())  -n 4 $(Base.julia_cmd()) --project dmda.jl`

    @info "Finished parallel tests"
end


# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
# XXX: Currently not working on windows reliably, not sure why
#if do_mpi
#    include("mpi_examples.jl")
#end

#=

include("options.jl")
#include("dmda.jl")
#include("old_test.jl")
#include("test_dmstag.jl")
include("test_snes.jl")

# Run the examples to make sure they all work
include("examples.jl")

=#