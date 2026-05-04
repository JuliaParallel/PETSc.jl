using Test
using PETSc
using SciMLBase
using DiffEqBase
using LinearAlgebra

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing

include(joinpath(dirname(dirname(@__DIR__)), "examples", "ex51_sciml.jl"))

@testset "ex51_sciml example" begin
    @testset "TSRK explicit Runge-Kutta" begin
        @testset "default algorithm (TSRK 5dp), fixed step, reaches final time and is accurate" begin
            result = solve_ex51(; verbose = false)
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 5.0e-8
        end

        @testset "TSRK 3bs is less accurate than 5dp (fixed step)" begin
            result_3bs = solve_ex51(; alg = PETSc.TSRK("3bs"), verbose = false)
            result_5dp = solve_ex51(; alg = PETSc.TSRK("5dp"), verbose = false)
            @test result_3bs.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test result_3bs.error < 5.0e-4
            @test result_5dp.error < result_3bs.error
        end

        @testset "repeated solves give the same result (fixed step)" begin
            result_1 = solve_ex51(; verbose = false)
            result_2 = solve_ex51(; verbose = false)
            @test result_1.error ≈ result_2.error rtol = 1e-12
        end

        @testset "smaller dt improves accuracy (fixed step)" begin
            result_coarse = solve_ex51(; dt = 0.25, verbose = false)
            result_fine   = solve_ex51(; dt = 0.05, verbose = false)
            @test result_fine.error < result_coarse.error
        end

        @testset "adaptive = true is accepted and reaches final time" begin
            result = solve_ex51(; adaptive = true, verbose = false)
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 5.0e-4
        end
    end

    @testset "TSRosW Rosenbrock-W" begin
        @testset "TSRosW(\"ra34pw2\") reaches final time and is accurate" begin
            result = solve_ex51(;
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                dt = 0.1, verbose = false,
            )
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 1.0e-3
        end

        @testset "TSRosW(\"rodas3\") reaches final time" begin
            result = solve_ex51(;
                alg = PETSc.TSRosW("rodas3", ["-snes_fd"]),
                dt = 0.1, verbose = false,
            )
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 1.0e-2
        end

        @testset "smaller dt improves accuracy for TSRosW" begin
            result_coarse = solve_ex51(;
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                dt = 0.2, verbose = false,
            )
            result_fine = solve_ex51(;
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            @test result_fine.error < result_coarse.error
        end

        @testset "repeated TSRosW solves give the same result" begin
            result_1 = solve_ex51(;
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                dt = 0.1, verbose = false,
            )
            result_2 = solve_ex51(;
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                dt = 0.1, verbose = false,
            )
            @test result_1.error ≈ result_2.error rtol = 1e-12
        end
    end

    @testset "TSIRK implicit Runge-Kutta (Gauss)" begin
        # TSIRK/Gauss requires an AIJ sparse Jacobian matrix so that PETSc can
        # form its internal MATKAIJ stage operator. The SciML extension currently
        # calls TSSetIFunction without TSSetIJacobian, so PETSc falls back to a
        # matrix-free (MFFD) Jacobian, which TSIRK/Gauss rejects — the same
        # limitation documented in examples/ex51_implicit.jl for -snes_mf.
        # These tests are marked @test_broken to track the gap. Once the
        # extension is extended to set up an AIJ Jacobian template and register
        # it via TSSetIJacobian (matching ex51_implicit.jl's
        # ex51_implicit_create_jacobian_template path), they should pass with the
        # error tolerances matching test/ts_ex51_implicit.jl.
        @testset "Gauss 1 stage (order 2): error < 5e-3" begin
            @test_broken begin
                result = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "1",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                isapprox(result.final_time, 1.0; atol = 100 * eps(Float64)) &&
                    length(result.solution) == 2 &&
                    result.error < 5.0e-3
            end
        end

        @testset "Gauss 2 stages (order 4): error < 5e-6" begin
            @test_broken begin
                result = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "2",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                isapprox(result.final_time, 1.0; atol = 100 * eps(Float64)) &&
                    length(result.solution) == 2 &&
                    result.error < 5.0e-6
            end
        end

        @testset "Gauss 3 stages (order 6): error < 5e-9" begin
            @test_broken begin
                result = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "3",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                isapprox(result.final_time, 1.0; atol = 100 * eps(Float64)) &&
                    length(result.solution) == 2 &&
                    result.error < 5.0e-9
            end
        end

        @testset "Gauss 4 stages (order 8): error < 5e-10" begin
            @test_broken begin
                result = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "4",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                isapprox(result.final_time, 1.0; atol = 100 * eps(Float64)) &&
                    length(result.solution) == 2 &&
                    result.error < 5.0e-10
            end
        end

        @testset "increasing stages strictly improve accuracy" begin
            @test_broken begin
                errors = map([1, 2, 3, 4]) do nstages
                    solve_ex51(;
                        alg = PETSc.TSGeneric("irk", [
                            "-ts_irk_type", "gauss",
                            "-ts_irk_nstages", string(nstages),
                            "-snes_fd_color",
                            "-ksp_type", "gmres",
                            "-pc_type", "none",
                        ]),
                        dt = 0.25, adaptive = false, verbose = false,
                    ).error
                end
                issorted(errors; rev = true)
            end
        end

        @testset "repeated solves give the same result (Gauss 2 stages)" begin
            @test_broken begin
                r1 = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "2",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                r2 = solve_ex51(;
                    alg = PETSc.TSGeneric("irk", [
                        "-ts_irk_type", "gauss",
                        "-ts_irk_nstages", "2",
                        "-snes_fd_color",
                        "-ksp_type", "gmres",
                        "-pc_type", "none",
                    ]),
                    dt = 0.25, adaptive = false, verbose = false,
                )
                isapprox(r1.error, r2.error; rtol = 1e-12)
            end
        end
    end

    @testset "TSImplicit fully implicit methods" begin
        @testset "TSImplicit(\"beuler\") reaches final time (1st-order)" begin
            result = solve_ex51(;
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                dt = 0.01, verbose = false,
            )
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 5.0e-2
        end

        @testset "TSImplicit(\"cn\") is more accurate than backward Euler at the same dt" begin
            result_beuler = solve_ex51(;
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            result_cn = solve_ex51(;
                alg = PETSc.TSImplicit("cn", ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            @test result_cn.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test result_cn.error < result_beuler.error
        end

        @testset "TSImplicit(\"theta\", 0.5) is second-order accurate" begin
            # theta=0.5 and "cn" are both second-order Crank-Nicolson variants but
            # differ internally (endpoint vs. midpoint staging), so their errors are
            # not numerically identical. Verify second-order accuracy directly by
            # checking that halving dt reduces the error by roughly a factor of 4.
            result_coarse = solve_ex51(;
                alg = PETSc.TSImplicit("theta", 0.5, ["-snes_fd"]),
                dt = 0.1, verbose = false,
            )
            result_fine = solve_ex51(;
                alg = PETSc.TSImplicit("theta", 0.5, ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            @test result_fine.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test result_fine.error < result_coarse.error / 2
        end

        @testset "TSImplicit(\"bdf\") reaches final time" begin
            result = solve_ex51(;
                alg = PETSc.TSImplicit("bdf", ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test result.error < 1.0e-2
        end

        @testset "smaller dt improves accuracy for TSImplicit beuler" begin
            result_coarse = solve_ex51(;
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                dt = 0.05, verbose = false,
            )
            result_fine = solve_ex51(;
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                dt = 0.01, verbose = false,
            )
            @test result_fine.error < result_coarse.error
        end
    end
end
