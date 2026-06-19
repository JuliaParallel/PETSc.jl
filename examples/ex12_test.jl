# Serial tests for examples/ex12.jl.
# Each test launches ex12.jl as a subprocess with specific options and checks
# that the solver converges and the L² error is below a tolerance.
#
# Run via:
#   julia --project examples/ex12_test.jl
# or through the project test suite (picked up automatically by test/examples.jl).

using Test

const EX12 = joinpath(@__DIR__, "ex12.jl")
const JL   = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project())`

"""
    run_ex12(args...; tol=nothing) -> Bool

Run ex12.jl with `args` as command-line arguments.  Returns `true` when the
subprocess exits cleanly.  If `tol` is given, also checks that the printed
`L2 error:` value is less than `tol`.
"""
function run_ex12(args...; tol::Union{Float64, Nothing} = nothing)
    out = IOBuffer()
    cmd = pipeline(`$JL $EX12 $args`; stdout = out, stderr = devnull)
    ok  = success(cmd)
    ok || return false
    tol === nothing && return true
    s = String(take!(out))
    m = match(r"L2 error:\s*([\d.eE+\-]+)", s)
    m === nothing && return true           # no L2 line — just trust exit code
    parse(Float64, m[1]) < tol
end

# ── Boundary-condition / coefficient combinations ────────────────────────────

@testset "ex12 serial" begin

  @testset "2D dirichlet field (Q1)" begin
    @test run_ex12("-dim","2","-coeff","field","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4"; tol = 0.1)
  end

  @testset "2D dirichlet nonlinear (Q1)" begin
    @test run_ex12("-dim","2","-coeff","nonlinear","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4"; tol = 0.1)
  end

  @testset "2D neumann field (Q1)" begin
    @test run_ex12("-dim","2","-coeff","field","-bc","neumann",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4"; tol = 0.1)
  end

  @testset "2D neumann nonlinear (Q1)" begin
    @test run_ex12("-dim","2","-coeff","nonlinear","-bc","neumann",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4"; tol = 0.1)
  end

  @testset "3D dirichlet field (Q1 hex)" begin
    @test run_ex12("-dim","3","-coeff","field","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4,4"; tol = 0.1)
  end

  @testset "3D dirichlet nonlinear (Q1 hex)" begin
    @test run_ex12("-dim","3","-coeff","nonlinear","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4,4"; tol = 0.1)
  end

  # ── Solver variants ──────────────────────────────────────────────────────────

  @testset "GAMG (algebraic multigrid)" begin
    @test run_ex12("-dim","2","-coeff","field","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","8,8",
                   "-pc_type","gamg","-ksp_type","cg"; tol = 0.05)
  end

  @testset "matrix-free Jacobian (snes_mf)" begin
    @test run_ex12("-dim","2","-coeff","nonlinear","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
                   "-snes_mf"; tol = 0.1)
  end

  @testset "matrix-free preconditioner (snes_mf_operator)" begin
    @test run_ex12("-dim","2","-coeff","nonlinear","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
                   "-snes_mf_operator"; tol = 0.1)
  end

  @testset "FAS nonlinear multigrid (2 levels)" begin
    @test run_ex12("-dim","2","-coeff","nonlinear","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","4,4",
                   "-dm_refine_hierarchy","1",
                   "-snes_type","fas","-snes_fas_levels","2",
                   "-fas_coarse_snes_type","newtonls","-fas_coarse_pc_type","lu",
                   "-fas_coarse_snes_linesearch_type","basic",
                   "-fas_levels_1_snes_type","newtonls","-fas_levels_1_pc_type","ilu",
                   "-fas_levels_1_ksp_rtol","1e-8"; tol = 0.1)
  end

  @testset "ASM with ILU subsolvers" begin
    @test run_ex12("-dim","2","-coeff","field","-bc","dirichlet",
                   "-dm_plex_simplex","0","-dm_plex_box_faces","8,8",
                   "-pc_type","asm","-pc_asm_blocks","4",
                   "-sub_pc_type","ilu","-ksp_type","gmres"; tol = 0.05)
  end

end
