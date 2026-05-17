# INCLUDE IN MPI TEST
# MPI tests for examples/ex62.jl (run with mpiexec -n 4 by the test framework).
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

  # 2 ranks — Q2/Q1 quad 2D, trig MMS
  @testset "2 ranks — Q2/Q1 2D trig MMS (L²≈0.65)" begin
    @test run_ex62_mpi(2,
      "-sol","trig","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1.0)
  end

  # 2 ranks — P2/P1 simplex 2D, quadratic MMS (machine eps)
  @testset "2 ranks — P2/P1 simplex 2D quadratic MMS" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — P2/P1-disc Crouzeix-Raviart
  @testset "2 ranks — P2/P1-disc simplex 2D Crouzeix-Raviart (machine eps)" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — Q2/Q1-disc discontinuous pressure on quads
  @testset "2 ranks — Q2/Q1-disc quad 2D discontinuous pressure (machine eps)" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — P2/P1 simplex 2D, parallel dmsnes_check
  @testset "2 ranks — P2/P1 simplex 2D dmsnes_check" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_refine","2","-petscpartitioner_type","simple",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-dmsnes_check","0.0001"; tol = 1e-8)
  end

  # 2 ranks — Q2/Q1 hex 3D, quadratic MMS (machine eps)
  @testset "2 ranks — Q2/Q1 3D quadratic MMS" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic",
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","3,3,3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # 2 ranks — P2/P1 SIMPLE preconditioner
  @testset "2 ranks — P2/P1 simplex 2D SIMPLE precond (L²≈6e-7)" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1","-petscds_jac_pre","0",
      "-snes_error_if_not_converged",
      "-ksp_type","fgmres","-ksp_gmres_restart","100","-ksp_rtol","1.0e-9","-ksp_error_if_not_converged",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","schur","-pc_fieldsplit_schur_factorization_type","full",
      "-fieldsplit_velocity_pc_type","lu","-fieldsplit_pressure_ksp_rtol","1e-10","-fieldsplit_pressure_pc_type","jacobi",
      "-fieldsplit_pressure_inner_ksp_type","preonly","-fieldsplit_pressure_inner_pc_type","jacobi",
      "-fieldsplit_pressure_upper_ksp_type","preonly","-fieldsplit_pressure_upper_pc_type","jacobi"; tol = 1e-5)
  end

  # ── Multigrid tests ────────────────────────────────────────────────────────

  # 2 ranks — GMG on velocity block (geometric multigrid via dm_refine hierarchy)
  @testset "2 ranks — Q2/Q1 2D GMG on velocity (L²≈7e-7)" begin
    @test run_ex62_mpi(2,
      "-sol","quadratic","-dm_plex_simplex","0","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-petscds_jac_pre","0",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
      "-pc_fieldsplit_schur_factorization_type","full",
      "-pc_fieldsplit_schur_precondition","full",
      "-fieldsplit_velocity_ksp_type","fgmres",
      "-fieldsplit_velocity_pc_type","mg",
      "-fieldsplit_velocity_mg_coarse_pc_type","svd",
      "-fieldsplit_pressure_ksp_rtol","1e-9",
      "-fieldsplit_pressure_pc_type","jacobi",
      "-ksp_type","fgmres","-ksp_rtol","1e-9"; tol = 1e-5)
  end

  # 2 ranks — GAMG on pressure block
  @testset "2 ranks — Q2/Q1 2D GAMG on pressure (L²≈0.11)" begin
    @test run_ex62_mpi(2,
      "-sol","trig","-dm_plex_simplex","0","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pc_use_amat",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
      "-pc_fieldsplit_schur_factorization_type","full",
      "-pc_fieldsplit_schur_precondition","a11",
      "-fieldsplit_velocity_pc_type","lu",
      "-fieldsplit_pressure_ksp_rtol","1e-9",
      "-fieldsplit_pressure_pc_type","gamg",
      "-ksp_rtol","1e-9"; tol = 0.25)
  end

end
