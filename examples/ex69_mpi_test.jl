# INCLUDE IN MPI TEST
# MPI tests for examples/ex69.jl.
# Each test launches ex69.jl under multiple MPI ranks and checks L² error.
#
# Run via:
#   julia --project examples/ex69_mpi_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test
using MPI

const EX69 = joinpath(@__DIR__, "ex69.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

const SOLVER = ["-pc_use_amat",
                "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
                "-pc_fieldsplit_schur_factorization_type","full",
                "-pc_fieldsplit_schur_precondition","a11",
                "-fieldsplit_velocity_pc_type","lu",
                "-fieldsplit_pressure_ksp_rtol","1e-9",
                "-fieldsplit_pressure_pc_type","lu",
                "-ksp_rtol","1e-9"]

"""
    run_ex69_mpi(nranks, args...; tol=nothing) -> Bool

Run ex69.jl under `nranks` MPI ranks.  Returns `true` when the subprocess
exits cleanly.  If `tol` is given, also checks `L2 error:` < `tol`.
"""
function run_ex69_mpi(nranks::Int, args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$(MPI.mpiexec()) -n $nranks $JL $EX69 $args`;
                   stdout = out, stderr = devnull)
    ok = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true
    parse(Float64, m[1]) < tol
end

@testset "ex69 MPI" begin

  # 2 ranks — Q2/Q1 SolKx, 1 uniform refinement (L²≈0.0029)
  @testset "2 ranks — Q2/Q1 SolKx dm_refine 1" begin
    @test run_ex69_mpi(2,
      "-dm_plex_simplex","0","-dm_plex_separate_marker","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.01)
  end

  # 2 ranks — Q1/P0 SolKx, 1 uniform refinement (L²≈0.032)
  @testset "2 ranks — Q1/P0 SolKx dm_refine 1" begin
    @test run_ex69_mpi(2,
      "-dm_plex_simplex","0","-dm_plex_separate_marker","-dm_refine","1",
      "-vel_petscspace_degree","1","-pres_petscspace_degree","0",
      SOLVER...; tol = 0.1)
  end

  # 2 ranks — Q2/Q1 SolCx etaB=1000, 1 refinement (L²≈0.020)
  @testset "2 ranks — Q2/Q1 SolCx etaB=1000 dm_refine 1" begin
    @test run_ex69_mpi(2,
      "-dm_plex_simplex","0","-dm_plex_separate_marker","-dm_refine","1",
      "-sol_type","solcx","-etaB","1e3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.08)
  end

  # 4 ranks — Q2/Q1 SolKx, 2 uniform refinements (L²≈0.00069)
  @testset "4 ranks — Q2/Q1 SolKx dm_refine 2" begin
    @test run_ex69_mpi(4,
      "-dm_plex_simplex","0","-dm_plex_separate_marker","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.003)
  end

  # 2 ranks — P2/P1 simplex (triangle), 1 refinement (L²≈0.006)
  @testset "2 ranks — P2/P1 simplex SolKx dm_refine 1" begin
    @test run_ex69_mpi(2,
      "-dm_plex_separate_marker","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.02)
  end

  # ── FETI-DP tests (equivalent to C-code q2p1fetidp / p2p1fetidp testsets) ───
  # All use 2 MPI ranks (C tests use 4-5 ranks; 2 is sufficient to exercise the
  # FETI-DP communication paths while staying within available memory).

  # Shared base options for quad-mesh FETI-DP (discontinuous P1 pressure).
  FETIDP_QUAD = ["-dm_plex_simplex","0","-dm_plex_separate_marker",
                 "-petscpartitioner_type","simple",
                 "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
                 "-pres_petscspace_poly_tensor","0",
                 "-pres_petscdualspace_lagrange_continuity","0",
                 "-pres_petscdualspace_lagrange_node_endpoints","0",
                 "-dm_mat_type","is",
                 "-ksp_type","fetidp","-ksp_fetidp_saddlepoint",
                 "-fetidp_ksp_type","cg","-fetidp_ksp_norm_type","natural",
                 "-fetidp_bddc_pc_bddc_detect_disconnected",
                 "-fetidp_bddc_pc_bddc_symmetric",
                 "-fetidp_bddc_pc_bddc_vertex_size","3",
                 "-fetidp_bddc_pc_bddc_graph_maxcount","2",
                 "-fetidp_bddc_pc_bddc_coarse_redundant_pc_type","svd"]

  # Shared base options for simplex FETI-DP (continuous P1 pressure).
  FETIDP_SIMP = ["-dm_plex_separate_marker",
                 "-petscpartitioner_type","simple",
                 "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
                 "-dm_mat_type","is",
                 "-ksp_type","fetidp","-ksp_fetidp_saddlepoint",
                 "-fetidp_ksp_type","cg","-fetidp_ksp_norm_type","natural",
                 "-fetidp_pc_fieldsplit_schur_fact_type","diag",
                 "-fetidp_fieldsplit_p_pc_type","jacobi",
                 "-fetidp_fieldsplit_p_ksp_type","preonly",
                 "-fetidp_fieldsplit_lag_ksp_type","preonly",
                 "-fetidp_bddc_pc_bddc_detect_disconnected",
                 "-fetidp_bddc_pc_bddc_symmetric",
                 "-fetidp_bddc_pc_bddc_vertex_size","3",
                 "-fetidp_bddc_pc_bddc_graph_maxcount","2",
                 "-fetidp_bddc_pc_bddc_coarse_redundant_pc_type","cholesky",
                 "-fetidp_bddc_pc_bddc_dirichlet_pc_type","svd",
                 "-fetidp_bddc_pc_bddc_neumann_pc_type","svd"]

  @testset "FETI-DP q2p1fetidp aij (L²≈0.097)" begin
    @test run_ex69_mpi(2, FETIDP_QUAD...,
      "-mat_is_localmat_type","aij"; tol = 0.15)
  end

  @testset "FETI-DP q2p1fetidp deluxe MUMPS (L²≈0.097)" begin
    @test run_ex69_mpi(2, FETIDP_QUAD...,
      "-mat_is_localmat_type","aij",
      "-fetidp_bddc_pc_bddc_use_deluxe_scaling",
      "-fetidp_bddc_pc_bddc_deluxe_singlemat",
      "-fetidp_bddc_sub_schurs_mat_solver_type","mumps",
      "-fetidp_bddc_sub_schurs_mat_mumps_icntl_14","500"; tol = 0.15)
  end

  @testset "FETI-DP p2p1fetidp simplex (L²≈0.029)" begin
    @test run_ex69_mpi(2, FETIDP_SIMP...; tol = 0.1)
  end

  @testset "FETI-DP p2p1fetidp_allp (L²≈0.029)" begin
    @test run_ex69_mpi(2, FETIDP_SIMP...,
      "-ksp_fetidp_pressure_all"; tol = 0.1)
  end

  @testset "FETI-DP p2p1fetidp_discharm (L²≈0.029)" begin
    @test run_ex69_mpi(2, FETIDP_SIMP...,
      "-fetidp_pc_discrete_harmonic","-fetidp_harmonic_pc_type","cholesky",
      "-fetidp_bddc_pc_bddc_dirichlet_pc_type","none"; tol = 0.1)
  end

  @testset "FETI-DP p2p1fetidp_lumped (L²≈0.029)" begin
    @test run_ex69_mpi(2, FETIDP_SIMP...,
      "-fetidp_pc_lumped",
      "-fetidp_bddc_pc_bddc_dirichlet_pc_type","none"; tol = 0.1)
  end

  @testset "FETI-DP p2p1fetidp deluxe MUMPS (L²≈0.029)" begin
    @test run_ex69_mpi(2, FETIDP_SIMP...,
      "-fetidp_bddc_pc_bddc_use_deluxe_scaling",
      "-fetidp_bddc_pc_bddc_deluxe_singlemat",
      "-fetidp_bddc_sub_schurs_mat_solver_type","mumps",
      "-fetidp_bddc_sub_schurs_mat_mumps_icntl_14","500",
      "-fetidp_bddc_sub_schurs_posdef","0",
      "-fetidp_bddc_pc_bddc_dirichlet_pc_type","none",
      "-fetidp_bddc_pc_bddc_neumann_pc_type","svd"; tol = 0.1)
  end

  # ── BDDC benign-trick tests (equivalent to C-code q2p1bddc testset) ──────────

  BDDC_BASE = ["-dm_plex_simplex","0","-dm_plex_separate_marker",
               "-petscpartitioner_type","simple",
               "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
               "-pres_petscspace_poly_tensor","0",
               "-pres_petscdualspace_lagrange_continuity","0",
               "-pres_petscdualspace_lagrange_node_endpoints","0",
               "-petscds_jac_pre","0",
               "-dm_mat_type","is",
               "-ksp_type","cg","-ksp_norm_type","natural",
               "-pc_type","bddc",
               "-pc_bddc_benign_trick","-pc_bddc_nonetflux",
               "-pc_bddc_detect_disconnected",
               "-pc_bddc_vertex_size","2",
               "-pc_bddc_coarse_redundant_pc_type","svd",
               "-pc_bddc_use_qr_single"]

  @testset "BDDC q2p1bddc benign card (L²≈0.097)" begin
    @test run_ex69_mpi(2, BDDC_BASE...,
      "-pc_bddc_dirichlet_pc_type","svd",
      "-pc_bddc_neumann_pc_type","svd"; tol = 0.15)
  end

  @testset "BDDC q2p1bddc benign deluxe MUMPS (L²≈0.097)" begin
    @test run_ex69_mpi(2, BDDC_BASE...,
      "-pc_bddc_use_deluxe_scaling","-pc_bddc_deluxe_zerorows",
      "-sub_schurs_mat_solver_type","mumps",
      "-sub_schurs_mat_mumps_icntl_14","1000"; tol = 0.15)
  end

  @testset "BDDC q2p1bddc benign deluxe adaptive MUMPS (L²≈0.097)" begin
    @test run_ex69_mpi(2, BDDC_BASE...,
      "-pc_bddc_adaptive_threshold","1.7",
      "-pc_bddc_use_deluxe_scaling","-pc_bddc_deluxe_zerorows",
      "-sub_schurs_mat_solver_type","mumps",
      "-sub_schurs_mat_mumps_icntl_14","1000"; tol = 0.15)
  end

end
