using Test
using PETSc, MPI
using PETSc: LibPETSc

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

# ============================================================================
#   Type stability of the generated wrappers
# ============================================================================
#
# 1. Static sweep: for every generated method of every LibPETSc function, replace the abstract
#    handle types in its signature by the concrete ones of the first library and ask inference
#    for the return type. It must be a single concrete type (no Union, no Any).
# 2. Dynamic checks on one representative wrapper per argument kind, for every library:
#    `@inferred` and zero allocations after warm-up (apart from arrays a wrapper has to return).

const HANDLE_CONCRETE = Dict(
    LibPETSc.AbstractPetscVec => LibPETSc.PetscVec, LibPETSc.AbstractPetscMat => LibPETSc.PetscMat,
    LibPETSc.AbstractPetscDM => LibPETSc.PetscDM, LibPETSc.AbstractPetscKSP => LibPETSc.PetscKSP,
    LibPETSc.AbstractPetscSNES => LibPETSc.PetscSNES, LibPETSc.AbstractPetscOptions => LibPETSc.PetscOptions,
    LibPETSc.AbstractIS => LibPETSc.IS, LibPETSc.AbstractPF => LibPETSc.PF, LibPETSc.AbstractTS => LibPETSc.TS,
    LibPETSc.AbstractAO => LibPETSc.AO, LibPETSc.AbstractTao => LibPETSc.Tao,
)

# concrete argument type for inference, given the library type
function concrete_arg(T, PetscLib)
    U = Base.unwrap_unionall(T)
    if U isa DataType
        for (A, C) in HANDLE_CONCRETE
            U.name === Base.unwrap_unionall(A).name && return C{PetscLib}
        end
        if U.name.name === :Union
            return T   # e.g. Union{Ptr, PetscVec}: leave as is, dispatch handles it
        end
        if U.name.name === :Array && !isempty(U.parameters) && U.parameters[1] isa TypeVar
            # Vector{<:AbstractPetscVec} and the like
            ub = Base.unwrap_unionall(U.parameters[1].ub)
            for (A, C) in HANDLE_CONCRETE
                ub isa DataType && ub.name === Base.unwrap_unionall(A).name && return Vector{C{PetscLib}}
            end
            return Vector{Any}
        end
    elseif T isa Union
        return T
    end
    return T
end

@testset "static return-type inference of every wrapper" begin
    petsclib = PETSc.petsclibs[1]
    PetscLib = typeof(petsclib)
    unstable = String[]
    checked = 0
    for name in names(LibPETSc; all = true)
        isdefined(LibPETSc, name) || continue
        f = getfield(LibPETSc, name)
        f isa Function || continue
        for m in methods(f)
            m.module === LibPETSc || continue
            sig = Base.unwrap_unionall(m.sig)
            sig isa DataType || continue
            params = sig.parameters[2:end]
            isempty(params) && continue
            # only the methods generated for the first library
            (params[1] isa Type && PetscLib <: params[1] && params[1] != Any) || continue
            argtypes = Any[PetscLib]
            for T in params[2:end]
                push!(argtypes, concrete_arg(T, PetscLib))
            end
            rts = try
                Base.return_types(f, Tuple{argtypes...})
            catch
                continue
            end
            checked += 1
            if length(rts) != 1 || !isconcretetype(rts[1]) && rts[1] !== Union{}
                push!(unstable, "$(m.name)$(Tuple{argtypes...}) -> $(rts)")
            end
        end
    end
    @info "static inference sweep" checked unstable = length(unstable)
    if !isempty(unstable)
        @info "first type-unstable wrappers" first(unstable, 20)
    end
    @test isempty(unstable)
end

@testset "representative wrappers: inferred and allocation free" begin
    comm = LibPETSc.PETSC_COMM_SELF
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscInt = petsclib.PetscInt
        PetscScalar = petsclib.PetscScalar
        n = PetscInt(8)

        # handle creation (Ref{CVec} output wrapped into a PetscVec)
        v = @inferred LibPETSc.VecCreateSeq(petsclib, comm, n)
        w = @inferred LibPETSc.VecDuplicate(petsclib, v)
        # scalar outputs
        @test (@inferred LibPETSc.VecGetSize(petsclib, v)) == n
        lo, hi = @inferred LibPETSc.VecGetOwnershipRange(petsclib, v)
        @test (lo, hi) == (0, n)
        # string output
        @test (@inferred LibPETSc.VecGetType(petsclib, v)) isa String
        # array input
        ix = PetscInt[0, 1]; vals = PetscScalar[1, 2]
        @inferred LibPETSc.VecSetValues(petsclib, v, PetscInt(2), ix, vals, LibPETSc.INSERT_VALUES)
        LibPETSc.VecAssemblyBegin(petsclib, v); LibPETSc.VecAssemblyEnd(petsclib, v)
        # PETSc-owned array output with a size rule
        a = @inferred LibPETSc.VecGetArray(petsclib, v)
        @test a isa Vector{PetscScalar} && a[1] == 1
        LibPETSc.VecRestoreArray(petsclib, v, a)
        # caller-allocated array output (sized by `ni`)
        y = @inferred LibPETSc.VecGetValues(petsclib, v, PetscInt(2), ix)
        @test y == vals
        # enum output and MPI_Comm output
        @test (@inferred LibPETSc.PetscObjectGetComm(petsclib, v)) isa MPI.Comm
        # opaque handle create/destroy by reference
        ns = @inferred LibPETSc.MatNullSpaceCreate(petsclib, comm, LibPETSc.PETSC_TRUE, PetscInt(0), LibPETSc.PetscVec{typeof(petsclib)}[])
        @test ns != C_NULL
        LibPETSc.MatNullSpaceDestroy(petsclib, ns)

        # allocations after warm-up
        LibPETSc.VecGetSize(petsclib, v)
        @test (@allocated LibPETSc.VecGetSize(petsclib, v)) == 0
        LibPETSc.VecGetOwnershipRange(petsclib, v)
        @test (@allocated LibPETSc.VecGetOwnershipRange(petsclib, v)) == 0
        LibPETSc.VecNorm(petsclib, v, LibPETSc.NORM_2)
        @test (@allocated LibPETSc.VecNorm(petsclib, v, LibPETSc.NORM_2)) == 0

        LibPETSc.VecDestroy(petsclib, w)
        LibPETSc.VecDestroy(petsclib, v)
        @test v.ptr == C_NULL
        PETSc.finalize(petsclib)
    end
end
