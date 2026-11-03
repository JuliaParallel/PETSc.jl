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

@testset "wrap_local_array and unsafe_local_array are removed" begin
    for name in (:wrap_local_array, :wrap_localarray, :_unsafe_local_array,
                 :unsafe_local_array, :unsafe_localarray)
        @test !isdefined(PETSc, name)
        @test !Base.ispublic(PETSc, name)
    end
    # the scoped and the explicit checkouts that replace them
    for name in (:with_local_array!, :acquire_local_array, :release_local_array)
        @test Base.ispublic(PETSc, name)
    end
end

@testset "the 0.5 deprecations are removed" begin
    # No v0.4 name is a binding of PETSc's own; `setfield!` and the like resolve to Base.
    for (old, _) in APISurface.RENAMES
        owned = isdefined(PETSc, old) &&
                (startswith(String(old), "@") || parentmodule(getfield(PETSc, old)) === PETSc)
        @test !owned
    end
    for name in (:audit_file, :audit_petsc_file,
                 :local_arrays, :restore_local_arrays!, :local_coordinate_array,
                 :AbstractPETScMemBackend, :PetscKSP, :PetscSNES)
        @test !isdefined(PETSc, name)
        @test !Base.ispublic(PETSc, name)
    end
    @test !isdefined(PETSc.LibPETSc, :PetscKSP)

    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscScalar = PETSc.scalartype(petsclib)
    comm = PETSc.MPI.COMM_SELF
    v = PETSc.PetscVec(petsclib, 4)
    # The pre-0.5.4 DM keywords throw rather than become unread PETSc options
    for kw in (:setfromoptions, :dmsetfromoptions, :dmsetup)
        @test_throws "is called" PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (5,), 1;
                                            NamedTuple{(kw,)}((false,))...)
    end
    # Old argument forms have no method
    @test_throws MethodError PETSc.set_type!(v, "standard")
    @test_throws MethodError PETSc.ownership_range(v, false)
    @test_throws MethodError PETSc.ghost_update!(v, PETSc.INSERT_VALUES)
    @test_throws MethodError PETSc.PetscMat(petsclib, Int[0, 1], Int[0], PetscScalar[1]; comm)
    @test_throws MethodError macroexpand(@__MODULE__, :(PETSc.@simple_fn(f)))
    PETSc.destroy!(v)
    PETSc.finalize(petsclib)
end
