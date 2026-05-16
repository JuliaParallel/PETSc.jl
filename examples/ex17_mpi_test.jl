# INCLUDE IN MPI TEST
# MPI tests for examples/ex17.jl (run with mpiexec -n 4 by the test framework).

using Test
using MPI

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

  @testset "4 ranks — 3D elas_quad GAMG" begin
    @test run_ex17_mpi(4,
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","4,4,4",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1",
      "-pc_type","gamg","-ksp_type","cg"; tol = 0.1)
  end

  @testset "4 ranks — 2D elas_quad LU" begin
    @test run_ex17_mpi(4,
      "-dm_plex_simplex","0","-dm_plex_box_faces","8,8",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1",
      "-pc_type","lu"; tol = 0.05)
  end

end
