using Test
using MPI: mpiexec

# Do the MPI tests first so we do not have mpi running inside MPI
@testset "mpi tests" begin
    for file in ("mpivec.jl", "mpimat.jl")
        @test mpiexec() do mpi_cmd
            cmd =
                `$mpi_cmd -n 4 $(Base.julia_cmd()) --startup-file=no --project $file`
            success(pipeline(cmd, stderr = stderr))
        end
    end
end

# Run the serial tests
include("init.jl")
include("options.jl")
include("vec.jl")
include("mat.jl")
include("matshell.jl")
