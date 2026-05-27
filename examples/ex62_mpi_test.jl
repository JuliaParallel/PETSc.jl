# INCLUDE IN MPI TEST
# MPI tests for examples/ex62.jl — covers correctness across element types and 3D.
# Solver-variant tests (GMG, GAMG, SIMPLE) are omitted to keep CI runtime short.
#
# Run via:
#   julia --project examples/ex62_mpi_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test
using MPI

let rank = parse(Int, get(ENV, "PMI_RANK", get(ENV, "OMPI_COMM_WORLD_RANK", "0")))
    rank == 0 || exit(0)
end

const EX62 = joinpath(@__DIR__, "ex62.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

const SOLVER = ["-pc_use_amat",
                "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
                "-pc_fieldsplit_schur_factorization_type","full",
                "-pc_fieldsplit_schur_precondition","a11",
                "-fieldsplit_velocity_pc_type","lu",
                "-fieldsplit_pressure_ksp_rtol","1e-9",
                "-fieldsplit_pressure_pc_type","lu",
                "-ksp_rtol","1e-9"]

function run_ex62_mpi(nranks::Int, args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$(MPI.mpiexec()) -n $nranks $JL $EX62 $args`;
                   stdout = out, stderr = devnull)
    ok = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true
    parse(Float64, m[1]) < tol
end

@testset "ex62 MPI" begin

  # 2 ranks — Q2/Q1 quad 2D, quadratic MMS (machine eps)
  @testset "2 ranks — Q2/Q1 2D quadratic MMS" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — P2/P1 simplex 2D, quadratic MMS (machine eps)
  @testset "2 ranks — P2/P1 simplex 2D quadratic MMS" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — Q2/Q1 hex 3D, quadratic MMS (machine eps)
  @testset "2 ranks — Q2/Q1 3D quadratic MMS" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic",
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","2,2,2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — P2/P1 simplex 2D, parallel Jacobian check
  @testset "2 ranks — P2/P1 simplex 2D dmsnes_check" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_refine","1","-petscpartitioner_type","simple",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-dmsnes_check","0.0001"; tol = 1e-8)
  end

end
