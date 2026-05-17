# INCLUDE IN MPI TEST
# MPI tests for examples/ex12.jl (run with mpiexec -n 4 by the test framework).
# Each test launches ex12.jl as an MPI subprocess and checks L² error.
#
# Run via:
#   julia --project examples/ex12_mpi_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test
using MPI

# The test framework runs this file with `mpiexec -n 4`.  Without a rank
# guard every rank would independently spawn subprocesses, executing each
# test 4×.  We detect rank via PMI_RANK (MPICH/MPItrampoline) or
# OMPI_COMM_WORLD_RANK (Open MPI) without calling MPI.Init(), which would
# block nested mpiexec calls on some implementations.
let rank = parse(Int, get(ENV, "PMI_RANK", get(ENV, "OMPI_COMM_WORLD_RANK", "0")))
    rank == 0 || exit(0)
end

const EX12 = joinpath(@__DIR__, "ex12.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

"""
    run_ex12_mpi(nranks, args...; tol=nothing) -> Bool

Run ex12.jl under `nranks` MPI ranks with the given arguments.
"""
function run_ex12_mpi(nranks::Int, args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$(MPI.mpiexec()) -n $nranks $JL $EX12 $args`;
                   stdout = out, stderr = devnull)
    ok = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true
    parse(Float64, m[1]) < tol
end

@testset "ex12 MPI" begin

  @testset "2 ranks — 3D field dirichlet (Q1 hex)" begin
    @test run_ex12_mpi(2, "-dim","3","-coeff","field","-bc","dirichlet",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","4,4,4";
                          tol = 0.1)
  end

  @testset "2 ranks — 2D nonlinear neumann (Q1)" begin
    @test run_ex12_mpi(2, "-dim","2","-coeff","nonlinear","-bc","neumann",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","4,4";
                          tol = 0.1)
  end

  @testset "2 ranks — 2D field dirichlet GAMG" begin
    @test run_ex12_mpi(2, "-dim","2","-coeff","field","-bc","dirichlet",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
                          "-pc_type","gamg","-ksp_type","cg"; tol = 0.1)
  end

end
