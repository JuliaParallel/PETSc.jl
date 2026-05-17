# INCLUDE IN MPI TEST
# MPI tests for examples/ex17.jl (run with mpiexec -n 4 by the test framework).
#
# Run via:
#   julia --project examples/ex17_mpi_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test
using MPI

let rank = parse(Int, get(ENV, "PMI_RANK", get(ENV, "OMPI_COMM_WORLD_RANK", "0")))
    rank == 0 || exit(0)
end

const EX17 = joinpath(@__DIR__, "ex17.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

function run_ex17_mpi(nranks::Int, args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$(MPI.mpiexec()) -n $nranks $JL $EX17 $args`;
                   stdout = out, stderr = devnull)
    ok = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true
    parse(Float64, m[1]) < tol
end

@testset "ex17 MPI" begin

  @testset "2 ranks — 3D elas_quad GAMG" begin
    @test run_ex17_mpi(2,
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","3,3,3",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1",
      "-pc_type","gamg","-ksp_type","cg"; tol = 0.1)
  end

  @testset "2 ranks — 2D elas_quad LU" begin
    @test run_ex17_mpi(2,
      "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1",
      "-pc_type","lu"; tol = 0.1)
  end

end
