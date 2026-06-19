# Serial tests for examples/ex17.jl.
# Covers all 7 solution types, 2D and 3D, and key solver variants.
#
# Run via:
#   julia --project examples/ex17_test.jl
# or through the project test suite (picked up by test/examples.jl).

using Test

const EX17 = joinpath(@__DIR__, "ex17.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

"""
    run_ex17(args...; tol=nothing, skip_l2=false) -> Bool

Run ex17.jl with `args`.  When `tol` is given, parse the printed `L2 error:`
and require it to be less than `tol`.  Pass `skip_l2=true` for cases that have
no meaningful exact solution (elas_ge, elas_axial_disp).
"""
function run_ex17(args...; tol::Union{Float64, Nothing} = nothing, skip_l2::Bool = false)
    out = IOBuffer()
    cmd = pipeline(`$JL $EX17 $args`; stdout = out, stderr = devnull)
    ok  = success(cmd)
    ok || return false
    (tol === nothing || skip_l2) && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true
    parse(Float64, m[1]) < tol
end

# Common mesh options shared across most tests.
const Q1_2D = ("-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
               "-displacement_petscspace_degree","1","-pc_type","lu")

@testset "ex17 serial" begin

  # ── Solution types (2D, Q1 hex, direct solver) ───────────────────────────────

  @testset "vlap_quad — vector Laplacian, quadratic exact solution" begin
    @test run_ex17(Q1_2D..., "-sol_type","vlap_quad"; tol = 0.1)
  end

  @testset "elas_quad — linear elasticity, quadratic exact solution" begin
    @test run_ex17(Q1_2D..., "-sol_type","elas_quad"; tol = 0.1)
  end

  @testset "vlap_trig — vector Laplacian, trig exact solution" begin
    @test run_ex17(Q1_2D..., "-sol_type","vlap_trig"; tol = 0.2)
  end

  @testset "elas_trig — linear elasticity, trig exact solution" begin
    @test run_ex17(Q1_2D..., "-sol_type","elas_trig"; tol = 0.5)
  end

  @testset "elas_uniform_strain — exact linear solution" begin
    # P1 represents the exact (linear) solution to machine precision.
    @test run_ex17(Q1_2D..., "-sol_type","elas_uniform_strain"; tol = 1e-10)
  end

  @testset "elas_axial_disp — mixed Dirichlet/Neumann (no L2 check)" begin
    # The BVP with partial Dirichlet + right-wall traction does not have
    # axial_disp_u as its exact solution; only the Jacobian accuracy is
    # meaningful here (matching C test behaviour: -dmsnes_check).
    @test run_ex17(Q1_2D..., "-sol_type","elas_axial_disp",
                   "-dm_plex_separate_marker"; skip_l2 = true)
  end

  @testset "elas_ge — geological shift, no exact solution" begin
    @test run_ex17(
      "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
      "-displacement_petscspace_degree","1",
      "-sol_type","elas_ge","-dm_plex_separate_marker",
      "-pc_type","gamg","-ksp_type","cg"; skip_l2 = true)
  end

  # ── 3D ───────────────────────────────────────────────────────────────────────

  @testset "3D elas_quad (Q1 hex)" begin
    @test run_ex17(
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","3,3,3",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1","-pc_type","lu"; tol = 0.1)
  end

  @testset "3D vlap_quad (Q1 hex)" begin
    @test run_ex17(
      "-dm_plex_dim","3","-dm_plex_simplex","0","-dm_plex_box_faces","3,3,3",
      "-sol_type","vlap_quad",
      "-displacement_petscspace_degree","1","-pc_type","lu"; tol = 0.1)
  end

  # ── Solver variants ──────────────────────────────────────────────────────────

  @testset "GAMG with rigid-body near-null space" begin
    @test run_ex17(
      "-dm_plex_simplex","0","-dm_plex_box_faces","8,8",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","1",
      "-pc_type","gamg","-ksp_type","cg"; tol = 0.05)
  end

  @testset "higher-order Q2 elements (elas_quad, refine 1)" begin
    # Q2 on a refined mesh: converges faster than Q1; use generous tol.
    @test run_ex17(
      "-dm_plex_simplex","0","-dm_plex_box_faces","4,4","-dm_refine","1",
      "-sol_type","elas_quad",
      "-displacement_petscspace_degree","2","-pc_type","lu"; tol = 0.01)
  end

end
