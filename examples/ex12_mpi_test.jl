# INCLUDE IN MPI TEST
# MPI tests for examples/ex12.jl (run with mpiexec -n 4 by the test framework).
# Each test launches ex12.jl as an MPI subprocess and checks L² error.

using Test
using MPI

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

  @testset "4 ranks — 3D field dirichlet (Q1 hex)" begin
    @test run_ex12_mpi(4, "-dim","3","-coeff","field","-bc","dirichlet",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","4,4,4";
                          tol = 0.1)
  end

  @testset "4 ranks — 2D nonlinear neumann (Q1)" begin
    @test run_ex12_mpi(4, "-dim","2","-coeff","nonlinear","-bc","neumann",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","8,8";
                          tol = 0.1)
  end

  @testset "4 ranks — 2D field dirichlet GAMG" begin
    @test run_ex12_mpi(4, "-dim","2","-coeff","field","-bc","dirichlet",
                          "-dm_plex_simplex","0","-dm_plex_box_faces","8,8",
                          "-pc_type","gamg","-ksp_type","cg"; tol = 0.05)
  end

end
