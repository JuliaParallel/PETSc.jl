using Test
using MPI

# Do the MPI tests first so we do not have mpi running inside MPI
@testset "mpi dmda tests" begin
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    @test mpiexec() do cmd
        run(`$cmd -n 4 $(julia) --project dmda.jl`)
        true
    end
end

include("options.jl")
include("dmda.jl")
include("old_test.jl")
include("examples.jl")
