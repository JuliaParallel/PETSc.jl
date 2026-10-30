# The API surface is fully accounted for by the register (naming.md §1.1).
#
# This is `scripts/api_surface.jl --check`, run from the test suite so that a new
# high-level binding cannot land without an entry in `scripts/renames.jl`. The
# fix for a failure is to add the name to RENAMES, INTERNAL or UNCHANGED_PUBLIC
# there, never to loosen the check.

using Test
using PETSc

include(joinpath(pkgdir(PETSc), "scripts", "api_surface.jl"))

@testset "api_surface --check" begin
    absent = APISurface.unregistered()
    if !isempty(absent)
        @info "bindings in no set of scripts/renames.jl" absent
    end
    @test isempty(absent)
    @test APISurface.check(devnull)
    # The surface is non-empty, so that a filter bug cannot make the check pass
    # by finding nothing at all.
    @test length(APISurface.surface()) > 100
end

@testset "handle types and type accessors are public" begin
    handles = (
        :AbstractPetscVec, :AbstractPetscMat, :AbstractPetscDM, :AbstractKSP,
        :AbstractPC, :AbstractSNES, :AbstractTS, :AbstractPetscOptions, :AbstractIS,
        :AbstractAO, :AbstractPF, :AbstractTao, :PetscDM, :PC, :IS, :AO, :PF, :Tao,
    )
    for name in handles
        @test Base.ispublic(PETSc, name)
        @test getfield(PETSc, name) === getfield(PETSc.LibPETSc, name)
        @test Docs.hasdoc(PETSc.LibPETSc, name)
    end
    # One function each, shared with LibPETSc
    for name in (:scalartype, :inttype, :realtype)
        @test Base.ispublic(PETSc, name)
        @test getfield(PETSc, name) === getfield(PETSc.LibPETSc, name)
    end
    petsclib = PETSc.petsclibs[1]
    @test PETSc.realtype(petsclib) === real(PETSc.scalartype(petsclib))
    @test PETSc.inttype(typeof(petsclib)) === petsclib.PetscInt
end

@testset "callback trampoline types are internal" begin
    trampolines = (
        :KSPComputeRHSFn, :KSPComputeOperatorsFn, :SNESSetFunctionFn, :SNESSetJacobianFn,
        :SNESSetConvergenceTestFn, :TSSetRHSFunctionFn, :TSSetRHSJacobianFn,
        :TSSetIFunctionFn, :TSSetIJacobianFn, :TSMonitorSetFn,
    )
    for name in trampolines
        @test isdefined(PETSc, name)
        @test !Base.ispublic(PETSc, name)
    end
    # the 0.4 spellings are gone
    for name in (:Fn_KSPComputeRHS, :Fn_KSPComputeOperators, :Fn_SNESSetFunction,
                 :Fn_SNESSetJacobian, :Fn_SNESSetConvergenceTest)
        @test !isdefined(PETSc, name)
    end
end
