using Test
using MPI: mpiexec

# Do the MPI tests first so we do not have mpi running inside MPI
for file in ("mpivec.jl", "mpimat.jl", "ksp.jl", "dmda.jl")
    @testset "MPI test: $file" begin
        @test mpiexec() do mpi_cmd
            cmd =
            `$mpi_cmd -n 4 $(Base.julia_cmd()) --startup-file=no --project $file`
            success(pipeline(cmd, stderr = stderr))
        end
    end
end

# run example with MPI that include comment line
# INCLUDE IN MPI TEST
include("mpi_examples.jl")

# Run the serial tests
include("init.jl")
include("options.jl")
include("vec.jl")
include("mat.jl")
include("matshell.jl")
include("ksp.jl")
include("dmda.jl")

# Run the examples
# Exclude examples with first line comment
# EXCLUDE FROM TESTING
include("examples.jl")
