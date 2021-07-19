using Test
using MPI

# Do the MPI tests first so we do not have mpi running inside MPI
# XXX: Currently not working on windows, not sure why
if !Sys.iswindows()
    @testset "mpi tests" begin
        @test mpiexec() do mpi_cmd
            cmd = `$mpi_cmd -n 4 $(Base.julia_cmd()) --project dmda.jl`
            success(pipeline(cmd, stderr = stderr))
        end
    end
end

# Examples with the comment
#   # INCLUDE IN MPI TEST
# will be run here
include("mpi_examples.jl")

include("options.jl")
include("dmda.jl")
include("old_test.jl")

# Run the examples to make sure they are all work
include("examples.jl")
