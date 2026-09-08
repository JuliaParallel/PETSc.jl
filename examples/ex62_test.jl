# Serial tests for examples/ex62.jl.
# Covers 2D and 3D, Q2/Q1 and P2/P1, both quadratic and trig MMS,
# block preconditioners, Vanka, and multigrid.
#
# Run via:
#   julia --project examples/ex62_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test

const EX62 = joinpath(@__DIR__, "ex62.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

# Common fieldsplit/Schur solver options.
const SOLVER = ["-pc_use_amat",
                "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
                "-pc_fieldsplit_schur_factorization_type","full",
                "-pc_fieldsplit_schur_precondition","a11",
                "-fieldsplit_velocity_pc_type","lu",
                "-fieldsplit_pressure_ksp_rtol","1e-9",
                "-fieldsplit_pressure_pc_type","lu",
                "-ksp_rtol","1e-9"]

"""
    run_ex62(args...; tol=nothing) -> Bool

Run ex62.jl with `args`.  Returns `true` when the subprocess exits cleanly and,
when `tol` is given, the printed `L2 error:` is below `tol`.
"""
function run_ex62(args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$JL $EX62 $args`; stdout = out, stderr = out)
    ok  = success(cmd)
    s   = String(take!(out))
    if !ok
        println("  ex62 exited non-zero: ", join(args, " "))
        println(s)
        return false
    end
    tol === nothing && return true
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    if m === nothing
        # A run that stops printing the error must not pass by default.
        println("  ex62 printed no L2 error: ", join(args, " "))
        println(s)
        return false
    end
    err = parse(Float64, m[1])
    println("  L2 error = ", err, " (tol ", tol, ")")
    err < tol && return true
    println("  ex62 L2 error exceeds tol: ", join(args, " "))
    return false
end

@testset "ex62 serial" begin

  # ── Quadratic MMS: L² ≈ machine eps for Q2/Q1 (exact polynomial) ───────────

  @testset "Q2/Q1 quad 2D — quadratic MMS (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  @testset "P2/P1 simplex 2D — quadratic MMS (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  @testset "P2/P1-disc simplex 2D — quadratic MMS Crouzeix-Raviart (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 1e-10)
  end

  @testset "Q2/Q1-disc quad 2D — quadratic MMS discontinuous pressure (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 1e-10)
  end

  @testset "Q2/Q1 hex 3D — quadratic MMS (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic",
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","3,3,3",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # ── dmsnes_check: Jacobian consistency via finite differences ────────────────
  # What these assert is that the run completes with a consistent Jacobian. 
  # They use the default solver, so the reported L2 error is set by the default KSP
  # tolerance and the 1e-4 finite-difference step rather than by the discretization. 
  # The tolerances below sit about two decades above the observed values so 
  # they bound the result without tracking round-off. 
  # Measured: 1.3e-9, 1.5e-10, 2.6e-9.

  @testset "Q2/Q1 quad 2D — dmsnes_check" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-dmsnes_check","0.0001"; tol = 1e-7)
  end

  @testset "P2/P1 simplex 2D — dmsnes_check" begin
    @test run_ex62(
      "-sol","quadratic",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-dmsnes_check","0.0001"; tol = 1e-8)
  end

  @testset "P3/P2 simplex 2D — quadratic MMS dmsnes_check" begin
    @test run_ex62(
      "-sol","quadratic",
      "-vel_petscspace_degree","3","-pres_petscspace_degree","2",
      "-dmsnes_check","0.0001"; tol = 1e-7)
  end

  # ── Trig MMS: measure convergence rate (Q2 → ~4–8× per refinement) ──────────

  @testset "Q2/Q1 quad 2D — trig MMS default mesh (L²≈0.65)" begin
    @test run_ex62(
      "-sol","trig","-dm_plex_simplex","0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1.0)
  end

  @testset "Q2/Q1 quad 2D — trig MMS dm_refine=1 (L²≈0.11)" begin
    @test run_ex62(
      "-sol","trig","-dm_plex_simplex","0","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.25)
  end

  @testset "Q2/Q1 quad 2D — trig MMS dm_refine=2 (L²≈0.024)" begin
    @test run_ex62(
      "-sol","trig","-dm_plex_simplex","0","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.06)
  end

  @testset "P2/P1 simplex 2D — trig MMS dm_refine=1 (L²≈0.15)" begin
    @test run_ex62(
      "-sol","trig","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 0.4)
  end

  @testset "P2/P1-disc simplex 2D — trig MMS dm_refine=1 Crouzeix-Raviart (L²≈0.15)" begin
    @test run_ex62(
      "-sol","trig","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 0.4)
  end

  @testset "Q2/Q1-disc quad 2D — trig MMS dm_refine=1 discontinuous pressure (L²≈0.11)" begin
    @test run_ex62(
      "-sol","trig","-dm_plex_simplex","0","-dm_refine","1",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-pres_petscdualspace_lagrange_continuity","0",
      SOLVER...; tol = 0.25)
  end

  # ── μ ≠ 1 ────────────────────────────────────────────────────────────────────

  @testset "Q2/Q1 quad 2D — quadratic MMS μ=2 (machine eps)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0","-mu","2.0",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      SOLVER...; tol = 1e-10)
  end

  # ── Block preconditioners ─────────────────────────────────────────────────────

  @testset "P2/P1 simplex 2D — block diagonal additive (L²≈2e-5)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1","-petscds_jac_pre","0",
      "-snes_error_if_not_converged",
      "-ksp_type","fgmres","-ksp_gmres_restart","100","-ksp_rtol","1.0e-4","-ksp_error_if_not_converged",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","additive",
      "-fieldsplit_velocity_pc_type","lu","-fieldsplit_pressure_pc_type","jacobi"; tol = 1e-3)
  end

  # The quadratic MMS is exact for P2/P1, so a fully converged solve leaves only round-off. 
  # At -ksp_rtol 1e-9 the residual has not fallen far enough for that,
  # and the reported error tracks where fgmres happened to stop. 
  # That depends on which kernels OpenBLAS dispatches for the runner CPU: 
  # measured 1.5e-7, 1.7e-7, 5.6e-7 and 3.3e-6 across four kernel paths, straddling this tolerance.
  # Solving to 1e-12 measures the discretization instead and collapses the spread to 1.5e-10 - 4.9e-10.
  @testset "P2/P1 simplex 2D — block triangular multiplicative (L²≈4e-10)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1","-petscds_jac_pre","0",
      "-snes_error_if_not_converged",
      "-ksp_type","fgmres","-ksp_gmres_restart","100","-ksp_rtol","1.0e-12","-ksp_error_if_not_converged",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","multiplicative",
      "-fieldsplit_velocity_pc_type","lu","-fieldsplit_pressure_pc_type","jacobi"; tol = 1e-6)
  end

  @testset "P2/P1 simplex 2D — Schur full factorization (L²≈0)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-snes_error_if_not_converged",
      "-ksp_type","fgmres","-ksp_gmres_restart","100","-ksp_rtol","1.0e-9","-ksp_error_if_not_converged",
      "-pc_use_amat",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
      "-pc_fieldsplit_schur_factorization_type","full","-pc_fieldsplit_off_diag_use_amat",
      "-fieldsplit_velocity_pc_type","lu",
      "-fieldsplit_pressure_ksp_rtol","1e-10","-fieldsplit_pressure_pc_type","jacobi"; tol = 1e-7)
  end

  # ── Vanka patch smoother ──────────────────────────────────────────────────────

  @testset "Q1/P0 quad 2D — Vanka patch (L²≈0.095)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0","-dm_refine","2",
      "-vel_petscspace_degree","1","-pres_petscspace_degree","0","-petscds_jac_pre","0",
      "-snes_rtol","1.0e-4",
      "-ksp_type","fgmres","-ksp_atol","1e-5","-ksp_error_if_not_converged",
      "-pc_type","patch","-pc_patch_partition_of_unity","0",
      "-pc_patch_construct_codim","0","-pc_patch_construct_type","vanka",
      "-sub_ksp_type","preonly","-sub_pc_type","lu"; tol = 0.2)
  end

  @testset "Q1/P0 quad 2D — Vanka dense inverse (L²≈0.095)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0","-dm_refine","2",
      "-vel_petscspace_degree","1","-pres_petscspace_degree","0","-petscds_jac_pre","0",
      "-snes_rtol","1.0e-4",
      "-ksp_type","fgmres","-ksp_atol","1e-5","-ksp_error_if_not_converged",
      "-pc_type","patch","-pc_patch_partition_of_unity","0",
      "-pc_patch_construct_codim","0","-pc_patch_construct_type","vanka",
      "-pc_patch_dense_inverse","-pc_patch_sub_mat_type","seqdense"; tol = 0.2)
  end

  # ── Multigrid (geometric multigrid on velocity, dm_refine for hierarchy) ────

  @testset "Q2/Q1 quad 2D — trig MMS GAMG on pressure" begin
    @test run_ex62(
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

  # Same reasoning as the multiplicative case above: solved to 1e-12 so the
  # reported error is the discretization error, not the iteration count.
  @testset "Q2/Q1 quad 2D — quadratic MMS GMG on velocity (L²≈3e-10)" begin
    @test run_ex62(
      "-sol","quadratic","-dm_plex_simplex","0","-dm_refine","2",
      "-vel_petscspace_degree","2","-pres_petscspace_degree","1",
      "-petscds_jac_pre","0",
      "-pc_type","fieldsplit","-pc_fieldsplit_type","schur",
      "-pc_fieldsplit_schur_factorization_type","full",
      "-pc_fieldsplit_schur_precondition","full",
      "-fieldsplit_velocity_ksp_type","fgmres",
      "-fieldsplit_velocity_pc_type","mg",
      "-fieldsplit_velocity_mg_coarse_pc_type","svd",
      "-fieldsplit_pressure_ksp_rtol","1e-12",
      "-fieldsplit_pressure_pc_type","jacobi",
      "-ksp_type","fgmres","-ksp_rtol","1e-12"; tol = 1e-5)
  end

end
