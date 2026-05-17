# Serial tests for examples/ex69.jl.
# Covers Q2/Q1 and Q1/P0 quad elements, SolKx and SolCx solution types,
# and a multi-refinement convergence check.
#
# Run via:
#   julia --project examples/ex69_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test

const EX69 = joinpath(@__DIR__, "ex69.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

# Common fieldsplit/Schur options used by all tests.
const SOLVER = ["-pc_use_amat",
                "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
                "-pc_fieldsplit_schur_factorization_type","full",
                "-pc_fieldsplit_schur_precondition","a11",
                "-fieldsplit_velocity_pc_type","lu",
                "-fieldsplit_pressure_ksp_rtol","1e-9",
                "-fieldsplit_pressure_pc_type","lu",
                "-ksp_rtol","1e-9"]

"""
    run_ex69(args...; tol=nothing) -> Bool

Run ex69.jl with `args`.  Returns `true` when the subprocess exits cleanly.
If `tol` is given, also checks that the printed `L2 error:` is less than `tol`.
"""
function run_ex69(args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$JL $EX69 $args`; stdout = out, stderr = devnull)
    ok  = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true   # no L2 line — trust exit code
    parse(Float64, m[1]) < tol
end

@testset "ex69 serial" begin

  # ── Q2/Q1 quad — SolKx (default, exponential viscosity) ─────────────────────

  @testset "Q2/Q1 SolKx — default 2×2 mesh (L²≈0.014)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.05)
  end

  @testset "Q2/Q1 SolKx — 4×4 mesh (L²≈0.0029)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.01)
  end

  # ── Q1/P0 quad — SolKx ───────────────────────────────────────────────────────

  @testset "Q1/P0 SolKx — default 2×2 mesh (L²≈0.062)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-vel_petscspace_degree","1","-pres_petscspace_degree","0",
      SOLVER...; tol = 0.2)
  end

  @testset "Q1/P0 SolKx — 4×4 mesh (L²≈0.032)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-dm_refine","1",
      "-vel_petscspace_degree","1","-pres_petscspace_degree","0",
      SOLVER...; tol = 0.1)
  end

  # ── Q2/Q1 quad — SolCx (piecewise-constant viscosity) ───────────────────────

  @testset "Q2/Q1 SolCx etaB=1 — default mesh (L²≈0.014)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-sol_type","solcx",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.05)
  end

  @testset "Q2/Q1 SolCx etaB=1000 — default mesh (L²≈0.039)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-sol_type","solcx","-etaB","1e3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.15)
  end

  @testset "Q2/Q1 SolCx etaB=1000 — 4×4 mesh (L²≈0.020)" begin
    @test run_ex69(
      "-dm_plex_simplex","0","-dm_plex_separate_marker",
      "-dm_refine","1",
      "-sol_type","solcx","-etaB","1e3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.08)
  end

  # ── Convergence rates (equivalent to C-code q2q1_conv / q1p0_conv tests) ────
  # Verify that L² error decreases by ~4× (Q2/Q1) and ~2× (Q1/P0) per refinement.

  @testset "Q2/Q1 SolKx convergence — 3 mesh levels" begin
    tols = [0.05, 0.015, 0.005]
    results = map(enumerate(tols)) do (r, tol)
      run_ex69(
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
        "-dm_refine","$(r-1)",
        "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
        SOLVER...; tol)
    end
    @test all(results)
  end

  @testset "Q1/P0 SolKx convergence — 3 mesh levels" begin
    tols = [0.2, 0.1, 0.05]
    results = map(enumerate(tols)) do (r, tol)
      run_ex69(
        "-dm_plex_simplex","0","-dm_plex_separate_marker",
        "-dm_refine","$(r-1)",
        "-vel_petscspace_degree","1","-pres_petscspace_degree","0",
        SOLVER...; tol)
    end
    @test all(results)
  end

  # ── P2/P1 simplex (triangle) — equivalent to C-code p2p1 / p2p1_conv tests ──

  @testset "P2/P1 simplex SolKx — default mesh (L²≈0.029)" begin
    @test run_ex69(
      "-dm_plex_separate_marker",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.1)
  end

  @testset "P2/P1 simplex SolKx convergence — 3 mesh levels" begin
    # C-code p2p1_conv: dm_refine 2 gives rate ~[3.0, 2.2] for vel/pres.
    tols = [0.1, 0.025, 0.007]
    results = map(enumerate(tols)) do (r, tol)
      run_ex69(
        "-dm_plex_separate_marker",
        "-dm_refine","$(r-1)",
        "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
        SOLVER...; tol)
    end
    @test all(results)
  end

  @testset "P2/P1 simplex SolCx etaB=1000 — default mesh" begin
    @test run_ex69(
      "-dm_plex_separate_marker",
      "-sol_type","solcx","-etaB","1e3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.2)
  end

end
