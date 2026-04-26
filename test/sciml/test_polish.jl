using Test
using PETSc
using SciMLBase
using DiffEqBase
using DataStructures

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing

function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Step 9 — Polish" begin
    @testset "Algorithm types are exported from PETSc itself" begin
        # Top-level access: users should be able to write `PETSc.TSRK(...)`
        # without going through `Base.get_extension`.
        @test PETSc.TSRK isa Type
        @test PETSc.TSRosW isa Type
        @test PETSc.TSImplicit isa Type
        @test PETSc.TSARKIMEX isa Type
        @test PETSc.TSGeneric isa Type
        @test PETSc.PETScTSAlgorithm isa Type
        # And every concrete type subtypes the abstract one.
        for T in (PETSc.TSRK, PETSc.TSRosW, PETSc.TSImplicit, PETSc.TSARKIMEX, PETSc.TSGeneric)
            @test T <: PETSc.PETScTSAlgorithm
        end
    end

    @testset "Algorithm docstrings exist and mention the PETSc TS type" begin
        # These are user-facing public-API entry points; the docstrings should
        # mention the underlying PETSc TS type so users know what they get.
        @test occursin("rk", lowercase(string(@doc PETSc.TSRK)))
        @test occursin("rosw", lowercase(string(@doc PETSc.TSRosW)))
        @test occursin("beuler", lowercase(string(@doc PETSc.TSImplicit)))
        @test occursin("arkimex", lowercase(string(@doc PETSc.TSARKIMEX)))
        @test occursin("tssettype", lowercase(string(@doc PETSc.TSGeneric)))
    end

    @testset "TSGeneric pass-through with implicit \"beuler\"" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(decay!, u0, tspan)
        sol = solve(prob, PETSc.TSGeneric("beuler", ["-snes_fd"]); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end][1] ≈ exp(-1.0) atol = 5e-2
    end

    @testset "Float32 problem raises a clear ArgumentError" begin
        u0 = Float32[1.0]
        tspan = (0f0, 1f0)
        prob = ODEProblem(decay!, u0, tspan)
        @test_throws ArgumentError solve(prob, PETSc.TSRK("3bs"); dt = 0.1f0)
        # Exception text should call out the constraint so the user knows what
        # to pass.
        err = try
            solve(prob, PETSc.TSRK("3bs"); dt = 0.1f0)
            nothing
        catch e
            e
        end
        @test occursin("PetscReal", err.msg)
        @test occursin("Float64", err.msg)
    end
end
