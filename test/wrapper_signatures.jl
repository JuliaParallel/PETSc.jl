using Test
using PETSc, MPI
using PETSc: LibPETSc

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

# ============================================================================
#   Wrapper input arguments take the abstract type
# ============================================================================
#
# A wrapper reads its inputs through `unsafe_convert`, so it has no reason to
# demand the plain wrapper struct. Annotating the abstract supertype lets a
# borrowed handle (`VecPtr`), a shell matrix, or a DM carrying its flavour in
# the type reach the same method.
#
# A regeneration that reverts this fails silently: the wrappers keep working
# for the plain types and only reject the subtypes, so the test checks the
# signatures rather than any one call.

# Concrete wrapper structs, which belong in a return position and not in an
# argument list.
const CONCRETE_NAMES = Set([
    :PetscVec,
    :PetscMat,
    :PetscDM,
    :PetscKSP,
    :PetscSNES,
    :PetscOptions,
    :IS,
    :PF,
    :TS,
    :AO,
    :Tao,
])

# Name of a type constructor, looking through `UnionAll` so both `PetscDM` and
# `PetscDM{PetscLib}` answer `:PetscDM`.
function type_name(@nospecialize(T))
    U = Base.unwrap_unionall(T)
    return U isa DataType ? U.name.name : nothing
end

# Argument types a method demands, including the element type of a vector
# argument, which is where the array-valued wrappers carry a wrapper type.
function argument_types(m::Method)
    sig = Base.unwrap_unionall(m.sig)
    sig isa DataType || return Any[]
    out = Any[]
    for T in sig.parameters[2:end]
        push!(out, T)
        U = Base.unwrap_unionall(T)
        if U isa DataType && U.name.name === :Array
            append!(out, U.parameters[1:1])
        end
    end
    return out
end

@testset "wrapper signatures take abstract types" begin
    offenders = String[]
    for name in names(LibPETSc; all = true)
        isdefined(LibPETSc, name) || continue
        f = getfield(LibPETSc, name)
        f isa Function || continue
        for m in methods(f)
            m.module in (LibPETSc, PETSc) || continue
            for T in argument_types(m)
                type_name(T) in CONCRETE_NAMES || continue
                push!(offenders, "$(m.name) at $(m.file):$(m.line) takes $(T)")
            end
        end
    end
    unique!(offenders)
    if !isempty(offenders)
        @info "wrapper arguments still annotated with a concrete type" count =
            length(offenders) first = first(offenders, 10)
    end
    @test isempty(offenders)
end

@testset "subtypes reach the wrappers" begin
    wrapper_comm = LibPETSc.PETSC_COMM_SELF
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscInt = petsclib.PetscInt

        # A borrowed vector handle, as PETScDiffEq.jl builds in every callback.
        v = LibPETSc.VecCreateSeq(petsclib, wrapper_comm, PetscInt(6))
        LibPETSc.VecZeroEntries(petsclib, v)
        borrowed = PETSc.VecPtr(petsclib, v.ptr, false)
        LibPETSc.VecSet(petsclib, borrowed, petsclib.PetscScalar(2))
        @test LibPETSc.VecNorm(petsclib, borrowed, LibPETSc.NORM_INFINITY) ≈ 2

        # The borrowed handle and its owner are one PETSc object.
        LibPETSc.VecScale(petsclib, borrowed, petsclib.PetscScalar(3))
        @test LibPETSc.VecNorm(petsclib, v, LibPETSc.NORM_INFINITY) ≈ 6

        # A shell matrix, whose type carries the wrapped Julia operator.
        shell = PETSc.MatShell(petsclib, (y, x) -> (y .= x), wrapper_comm, 4, 4)
        @test LibPETSc.MatGetLocalSize(petsclib, shell) == (4, 4)

        PETSc.destroy(v)
        PETSc.destroy(shell)
        PETSc.finalize(petsclib)
    end
end
