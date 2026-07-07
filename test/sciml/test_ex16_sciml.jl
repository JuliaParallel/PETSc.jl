using Test
using PETSc
using SciMLBase

include(joinpath(dirname(dirname(@__DIR__)), "examples", "ex16_sciml.jl"))

# Van der Pol ODE (ex16) — notes on the test design
#
# The problem is stiff for large mu; mu = 1000 is the upstream default.
# In the IMEX split (`imex = true`) the stiff part is isolated in f1 and
# handled implicitly, so TSARKIMEX converges reliably with -snes_fd even at
# mu = 1000. For the fully-implicit path (`imex = false`) the entire Jacobian
# is approximated by finite differences without the benefit of the analytic
# Jacobian registered in examples/ex16.jl, so the tests there use a more
# moderate mu = 100 to ensure robust SNES convergence.
#
# All solves use adaptive = true (the ex16_sciml default) so the step-size
# controller compensates for the stiffness automatically.

@testset "ex16_sciml example" begin

    @testset "TSARKIMEX IMEX (imex = true)" begin
        @testset "default solve reaches final time with a finite solution" begin
            result = solve_ex16(; verbose = false)
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test result.steps > 0
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
            @test result.imex == true
        end

        @testset "repeated solves give the same result" begin
            r1 = solve_ex16(; verbose = false)
            r2 = solve_ex16(; verbose = false)
            @test r1.solution ≈ r2.solution rtol = 1e-12
            @test r1.steps == r2.steps
        end

        @testset "TSARKIMEX subtypes all reach final time" begin
            for subtype in ("2e", "3", "4", "5", "ars122", "ars443")
                @testset "TSARKIMEX(\"$subtype\")" begin
                    result = solve_ex16(;
                        alg = PETSc.TSARKIMEX(subtype, ["-snes_fd"]),
                        verbose = false,
                    )
                    @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
                    @test result.steps > 0
                    @test length(result.solution) == 2
                    @test all(isfinite, result.solution)
                end
            end
        end
    end

    @testset "TSImplicit fully implicit (imex = false)" begin
        # Use mu = 100 (less stiff) so the finite-difference Jacobian converges
        # without the hand-coded analytical Jacobian from examples/ex16.jl.
        @testset "TSImplicit(\"beuler\") reaches final time" begin
            result = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                verbose = false,
            )
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test result.steps > 0
            @test !result.imex
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
        end

        @testset "TSImplicit(\"cn\") reaches final time" begin
            result = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSImplicit("cn", ["-snes_fd"]),
                verbose = false,
            )
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
        end

        @testset "TSImplicit(\"bdf\") reaches final time" begin
            result = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSImplicit("bdf", ["-snes_fd"]),
                verbose = false,
            )
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
        end

        @testset "repeated fully-implicit solves give the same result" begin
            r1 = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                verbose = false,
            )
            r2 = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                verbose = false,
            )
            @test r1.solution ≈ r2.solution rtol = 1e-12
        end
    end

    @testset "TSRosW Rosenbrock-W (imex = false)" begin
        @testset "TSRosW(\"ra34pw2\") reaches final time" begin
            result = solve_ex16(;
                imex = false, mu = 100.0,
                alg = PETSc.TSRosW("ra34pw2", ["-snes_fd"]),
                verbose = false,
            )
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test result.steps > 0
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
        end
    end

    @testset "mu parameter" begin
        @testset "mu = 100.0 returns the correct parameter in result" begin
            result = solve_ex16(; mu = 100.0, verbose = false)
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test result.mu ≈ 100.0
            @test length(result.solution) == 2
            @test all(isfinite, result.solution)
        end

        @testset "less stiff problem (mu = 100) needs no more steps than stiff (mu = 1000)" begin
            r_stiff = solve_ex16(; mu = 1000.0, verbose = false)
            r_mild  = solve_ex16(; mu = 100.0,  verbose = false)
            # A smaller stiffness parameter should require at most as many steps.
            @test r_mild.steps <= r_stiff.steps
        end

        @testset "mu = 1000.0 (default) yields a finite solution" begin
            result = solve_ex16(; mu = 1000.0, verbose = false)
            @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
            @test all(isfinite, result.solution)
        end
    end

    @testset "imex flag is reflected in the return value" begin
        r_imex    = solve_ex16(; imex = true,  verbose = false)
        r_noimex  = solve_ex16(; imex = false, mu = 100.0,
                                 alg = PETSc.TSImplicit("beuler", ["-snes_fd"]),
                                 verbose = false)
        @test r_imex.imex   == true
        @test r_noimex.imex == false
        # Both reach the requested final time.
        @test r_imex.final_time   ≈ 0.5 atol = 100 * eps(Float64)
        @test r_noimex.final_time ≈ 0.5 atol = 100 * eps(Float64)
    end
end
