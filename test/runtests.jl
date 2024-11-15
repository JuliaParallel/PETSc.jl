using Test
using MPI: mpiexec

# Do the MPI tests first so we do not have mpi running inside MPI
mpi_tests = ("mpivec.jl", "mpimat.jl", "ksp.jl", "dmstag.jl")

# XXX: We have problems with MPI and Windows with this test
if !(Sys.iswindows())
    mpi_tests = (mpi_tests..., "dmda.jl")
end

#=
for file in mpi_tests
    @testset "MPI test: $file" begin
        @test mpiexec() do mpi_cmd
            cmd =
                `$mpi_cmd -n 4 $(Base.julia_cmd()) --startup-file=no --project $file`
            success(pipeline(cmd, stderr = stderr))
        end
    end
end
=#

# run example with MPI that include comment line
# INCLUDE IN MPI TEST
##include("mpi_examples.jl")

# Run the serial tests
include("init.jl")
include("lib.jl")
include("options.jl")
include("vec.jl")
include("mat.jl")
include("matshell.jl")
include("ksp.jl")
include("snes.jl")
#include("dmda.jl")
#include("dmstag.jl")
include("test_dmstag.jl")

# Run the examples
# Exclude examples with first line comment
# EXCLUDE FROM TESTING
#include("examples.jl")
