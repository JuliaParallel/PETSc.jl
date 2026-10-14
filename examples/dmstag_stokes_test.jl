# Checks examples/dmstag_stokes.jl: the iterative solve matches a direct one at 16², and
# the iteration count stays flat under refinement.
#
#   julia --project examples/dmstag_stokes_test.jl

using Test

include(joinpath(@__DIR__, "dmstag_stokes.jl"))

# The iterative solve of the example: the fields and the iteration count
function iterative(n; levels, options...)
    dm, A, P, b = stokes_system(petsclib, comm, n, 1e4)
    ksp = stokes_solver(petsclib, comm, dm, A, P; levels, options...)
    x = PETSc.global_vec(dm)
    PETSc.solve!(x, ksp, b; check = true)
    sol = (; owned_fields(dm, x)..., iterations = PETSc.iteration_number(ksp))
    foreach(PETSc.destroy!, (ksp, x, b, P, A, dm))
    return sol
end

# LU on A. In a closed box A is singular, since the pressure is defined up to a constant,
# so the pressure in the first cell is fixed at zero.
function direct(n)
    dm, A, P, b = stokes_system(petsclib, comm, n, 1e4)
    PETSc.zero_rows_local!(A, dm, [PETSc.stencil(dm, PC, (1, 1))], 1.0)
    ksp = PETSc.KSP(A; ksp_type = "preonly", pc_type = "lu", pc_factor_mat_solver_type = "umfpack")
    x = PETSc.global_vec(dm)
    PETSc.solve!(x, ksp, b; check = true)
    sol = owned_fields(dm, x)
    foreach(PETSc.destroy!, (ksp, x, b, P, A, dm))
    return sol
end

demean(P) = P .- sum(P) / length(P)

@testset "dmstag_stokes" begin
    @testset "matches a direct solve" begin
        it = iterative(16; levels = 2, ksp_rtol = 1e-10)
        lu = direct(16)
        V = maximum(abs, lu.Vy)
        @test maximum(abs, it.Vx - lu.Vx) < 1e-7 * V
        @test maximum(abs, it.Vy - lu.Vy) < 1e-7 * V
        @test maximum(abs, demean(it.P) - demean(lu.P)) < 1e-7 * maximum(abs, lu.P)
        @test abs(sum(it.P)) / length(it.P) < 1e-10 * maximum(abs, it.P)  # the null space is removed
    end
    @testset "iterations do not grow with the grid" begin
        # From 64² on the sinker's outline is resolved and the count settles
        n64 = iterative(64; levels = 4).iterations
        n128 = iterative(128; levels = 5).iterations
        @test n64 < 40
        @test n128 <= 1.2 * n64
    end
end
