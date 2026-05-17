# INCLUDE IN MPI TEST
# MPI tests for examples/ex69.jl.
# Each test launches ex69.jl under multiple MPI ranks and checks L² error.
#
# Run via:
#   julia --project examples/ex69_mpi_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test
using MPI

let rank = parse(Int, get(ENV, "PMI_RANK", get(ENV, "OMPI_COMM_WORLD_RANK", "0")))
    rank == 0 || exit(0)
end

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

    # ── Basic fieldsplit / Schur tests ────────────────────────────────────────

    @testset "2 ranks — Q2/Q1 SolKx (L²≈0.014)" begin
      @test run_ex69_mpi(2,
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
        "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
        SOLVER...; tol = 0.05)
    end

    @testset "2 ranks — Q1/P0 SolKx (L²≈0.062)" begin
      @test run_ex69_mpi(2,
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
        "-vel_petscspace_degree","1","-pres_petscspace_degree","0",
        SOLVER...; tol = 0.2)
    end

    @testset "2 ranks — Q2/Q1 SolCx etaB=1000 (L²≈0.039)" begin
      @test run_ex69_mpi(2,
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
        "-sol_type","solcx","-etaB","1e3",
        "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
        SOLVER...; tol = 0.15)
    end

    @testset "2 ranks — P2/P1 simplex SolKx (L²≈0.029)" begin
      @test run_ex69_mpi(2,
        "-dm_plex_separate_marker",
        "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
        SOLVER...; tol = 0.1)
    end

    # ── FETI-DP (quad disc-P1 pressure) ──────────────────────────────────────

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

    @testset "FETI-DP q2p1fetidp aij (L²≈0.097)" begin
      @test run_ex69_mpi(2, FETIDP_QUAD...,
        "-mat_is_localmat_type","aij"; tol = 0.15)
    end

    # ── FETI-DP (simplex cont-P1 pressure) ────────────────────────────────────

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

    @testset "FETI-DP p2p1fetidp simplex (L²≈0.029)" begin
      @test run_ex69_mpi(2, FETIDP_SIMP...; tol = 0.1)
    end

    # ── BDDC benign-trick ─────────────────────────────────────────────────────

    @testset "BDDC q2p1bddc benign card (L²≈0.097)" begin
      @test run_ex69_mpi(2,
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
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
        "-pc_bddc_use_qr_single",
        "-pc_bddc_dirichlet_pc_type","svd",
        "-pc_bddc_neumann_pc_type","svd"; tol = 0.15)
    end

end
