import .LibPETSc: AbstractKSP, CKSP, KSP, AbstractPetscDM   # KSP methods below are constructors of LibPETSc.KSP

# The Julia side of a KSP, kept with the PETSc object
mutable struct KSPState <: ObjectState
    computerhs!::Any
    computeops!::Any
    opts::Any
    alive::Bool
end
KSPState() = KSPState(
    x -> error("computerhs! not defined"),
    x -> error("computeops! not defined"),
    nothing,
    true,
)
state_type(::Type{<:KSP}) = KSPState

# `ksp.opts` and the callbacks reach the state
@inline Base.getproperty(ksp::KSP, name::Symbol) =
    (name === :ptr || name === :age || name === :own) ? getfield(ksp, name) :
    forward_getproperty(ksp, name)
@inline function Base.setproperty!(ksp::KSP, name::Symbol, value)
    (name === :ptr || name === :age || name === :own) &&
        return Base.setfield!(ksp, name, convert(fieldtype(typeof(ksp), name), value))
    return forward_setproperty!(ksp, name, value)
end
Base.propertynames(ksp::KSP, private::Bool = false) = forward_propertynames(typeof(ksp))


# positional constructor taking callbacks: the callbacks go to the state
function LibPETSc.KSP{PetscLib}(
    ptr::CKSP,
    age::Int,
    computerhs!::Function,
    computeops!::Function = x -> error("computeops! not defined"),
    opts = nothing,
) where {PetscLib}
    ptr == C_NULL && throw(ArgumentError("callbacks need a PETSc object; got a null pointer"))
    ksp = KSP{PetscLib}(ptr, age)
    ksp.computerhs! = computerhs!
    ksp.computeops! = computeops!
    ksp.opts = opts
    return ksp
end

# Custom display for REPL
function Base.show(io::IO, v::AbstractKSP{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc KSP (null pointer)")
        return
    else
        print(io, "PETSc KSP object")
    end
    return nothing
end


"""
    KSP(comm::MPI.Comm, A::AbstractPetscMat, P::AbstractPetscMat{PetscLib} = A; prefix="", options...)

Create a `KSP` using the matrix `A` and preconditioner construction matrix `P`
with optional `prefix` and `options`.

The communicator is obtained from `A` and if it has size `1` then the garbage
collector is set, otherwise the user is responsible for calling
[`destroy!`](@ref).

# External Links
$(doc_external("KSP/KSPCreate"))
$(doc_external("KSP/KSPSetOperators"))
$(doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(
    A::AbstractPetscMat{PetscLib},
    P::AbstractPetscMat{PetscLib} = A;
    prefix::String="",
    options...
) where {PetscLib}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    c = comm(A)
    ksp = LibPETSc.KSPCreate(petsclib, c)
    
    LibPETSc.KSPSetOperators(petsclib, ksp, A, P)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, prefix)
    end
    
    # Push options to PETSc options database
    if !isempty(options)
        opts = PetscOptions(petsclib; options...);
        push!(opts)
        LibPETSc.KSPSetFromOptions(petsclib, ksp)
        pop!(opts)
        ksp.opts = opts
    end

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(c) == 1
        finalizer(destroy!, ksp)
    end

    return ksp
end

"""
    KSP(dm::AbstractPetscDM; prefix="", options...)

Create a `KSP` associated with the `dm` with optional `prefix` and `options`.

The communicator is obtained from `dm`. The KSP can be used with geometric
multigrid when the DM provides grid hierarchy information.

# Arguments
- `dm::AbstractPetscDM`: The DM object to associate with the KSP
- `prefix::String`: Optional prefix for command-line options
- `options...`: Additional PETSc options as keyword arguments

# External Links
$(doc_external("KSP/KSPCreate"))
$(doc_external("KSP/KSPSetDM"))
$(doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(dm::AbstractPetscDM{PetscLib};
    prefix::String="",
    options...
) where {PetscLib}
    check_initialized(getlib(PetscLib))
    petsclib = getlib(PetscLib)
    c = comm(dm)
    ksp = LibPETSc.KSPCreate(petsclib, c)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, prefix)
    end
    
    LibPETSc.KSPSetDM(petsclib, ksp, dm)

    # Push options to PETSc options database
    if !isempty(options)
        opts = PetscOptions(petsclib; options...);
        push!(opts)
        LibPETSc.KSPSetFromOptions(petsclib, ksp)
        pop!(opts)
        ksp.opts = opts
    end

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(c) == 1
        finalizer(destroy!, ksp)
    end

    return ksp
end



"""
    KSP(petsclib, comm::MPI.Comm, A::SparseMatrixCSC; options...)

Create a [`KSP`](@ref) with the sparse matrix `A` using the `petsclib`. If
`petsclib` is not given, the default library will be used`.
"""
KSP(petsclib, comm, S::SparseMatrixCSC; kwargs...) 

function KSP(petsclib, comm, S::SparseMatrixCSC; kwargs...) 
    M = LibPETSc.PetscMat(petsclib, comm, S)
    ksp = KSP(M; kwargs...)
    destroy!(M)   # the KSP holds its own reference
    return ksp
end


function solve!(
    x::PetscVec{PetscLib},
    ksp::KSP{PetscLib},
    b::PetscVec{PetscLib},
) where {PetscLib}
    has_opts = !isnothing(ksp.opts)
    has_opts && push!(ksp.opts)
    try
        capture_callback_errors() do
            LibPETSc.KSPSolve(PetscLib, ksp, b, x)
        end
    finally
        has_opts && pop!(ksp.opts)
    end
    return x
end

function solve!(
    ksp::AbstractKSP{PetscLib},
) where {PetscLib}
    has_opts = hasproperty(ksp, :opts) && !isnothing(ksp.opts)
    has_opts && push!(ksp.opts)
    try
        capture_callback_errors() do
            LibPETSc.KSPSolve(PetscLib, ksp, C_NULL, C_NULL)
        end
    finally
        has_opts && pop!(ksp.opts)
    end
    return ksp
end

LinearAlgebra.ldiv!(x::PetscVec{PetscLib}, ksp::KSP{PetscLib}, b::PetscVec{PetscLib}) where {PetscLib} = solve!(x, ksp, b)

function Base.:\(ksp::KSP, b::PetscVec{PetscLib}) where {PetscLib}
    x = similar(b)
    ldiv!(x, ksp, b)
    return x
end

function Base.:\(
    ksp::KSP{PetscLib},
    b::Vector{PetscScalar},
) where {PetscLib, PetscScalar}
    PetscScalar === PetscLib.PetscScalar || throw(
        ArgumentError(
            "right-hand side has element type $PetscScalar, " *
            "but the library uses $(PetscLib.PetscScalar)",
        ),
    )
    c = comm(ksp)
    MPI.Comm_size(c) == 1 || throw(
        ArgumentError(
            "solving into a Julia Vector requires a sequential KSP, " *
            "but its communicator spans $(MPI.Comm_size(c)) ranks",
        ),
    )
    PetscInt = PetscLib.PetscInt

    # PETSc works on this copy of `b` in place, so it must outlive the solve
    b_copy = PetscScalar.(b)
    x = GC.@preserve b_copy begin
        petsc_b = LibPETSc.VecCreateSeqWithArray(getlib(PetscLib), c, PetscInt(1), PetscInt(length(b)), b_copy)
        petsc_x = ksp \ petsc_b
        x = petsc_x[:]
        destroy!(petsc_b)
        destroy!(petsc_x)
        x
    end

    return x
end


"""
    destroy!(ksp::KSP)

Destroy the solver `ksp` holds. Does nothing on a borrowed handle, such as the
one `ksp(ts)` hands back: see [`owns`](@ref).

# External Links
$(doc_external("KSP/KSPDestroy"))
"""
function destroy!(ksp::KSP{PetscLib}) where {PetscLib}
    owns(ksp) || return nothing
    if !isnothing(ksp.opts)
        destroy!(ksp.opts)
        ksp.opts = nothing
    end
    if isdestroyable(ksp, PetscLib)
        LibPETSc.KSPDestroy(PetscLib, ksp)
    end
    ksp.ptr = C_NULL
    return nothing
end



"""
    dm(ksp::AbstractKSP)

The DM attached to `ksp`, [`narrow`](@ref)ed to its flavour.

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetDM"))
"""
function dm(ksp::AbstractKSP{PetscLib}) where PetscLib
    return narrow(LibPETSc.KSPGetDM(getlib(PetscLib), ksp))
end

#
# Wrapper for calls to set_compute_rhs!
mutable struct KSPComputeRHSFn{PetscLib, PetscInt} end
function (w::KSPComputeRHSFn{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cb::CVec,
    ksp_ptr::Ptr{Cvoid},
) where {PetscLib, PetscInt}
    new_ksp = KSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age; own = false)
    b = PetscVec{PetscLib}(cb, getlib(PetscLib).age; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)::KSPState
    return run_callback("right-hand side rhs!") do
        ksp.computerhs!(b, new_ksp)
    end
end

"""
    set_compute_rhs!(rhs!::Function, ksp::AbstractKSP)

Define `rhs!` to be the right-hand side function of the `ksp`. A call to
`rhs!(b, new_ksp)` should set the elements of the PETSc vector `b` based on the
`new_ksp`.

!!! note

    The `new_ksp` passed to `rhs!` may not be the same as the `ksp` passed to
    `set_compute_rhs!`.

$(doc_callback())

# External Links
$(doc_external("KSP/KSPSetComputeRHS"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_compute_rhs! end
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function set_compute_rhs!(rhs!, ksp::AbstractKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        KSPComputeRHSFn{$PetscLib, $PetscInt}(),
        LibPETSc.PetscErrorCode,
        (CKSP, CVec, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computerhs! = rhs!
    LibPETSc.KSPSetComputeRHS($PetscLib, ksp, fptr, state_pointer(ksp))
    return ksp
end

# Wrapper for calls to set_compute_rhs!
mutable struct KSPComputeOperatorsFn{PetscLib, PetscInt} end
function (w::KSPComputeOperatorsFn{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cA::CMat,
    cP::CMat,
    ksp_ptr::Ptr{Cvoid},
) where {PetscLib, PetscInt}
    new_ksp = KSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age; own = false)
    A = PetscMat{PetscLib}(cA, getlib(PetscLib).age; own = false)
    P = PetscMat{PetscLib}(cP, getlib(PetscLib).age; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)::KSPState
    return run_callback("operators ops!") do
        ksp.computeops!(A, P, new_ksp)
    end
end

"""
    set_compute_operators!(ops!::Function, ksp::KSP)

Define `ops!` to be the compute operators function for the `ksp`. A call to
`ops!(A, P, new_ksp)` should set the elements of the PETSc matrix linear
operator `A` and preconditioning matrix `P` based on the `new_ksp`.

!!! note

    The `new_ksp` passed to `ops!` may not be the same as the `ksp` passed to
    `set_compute_operators!`.

$(doc_callback())

# External Links
$(doc_external("KSP/KSPSetComputeOperators"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_compute_operators! end
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function set_compute_operators!(ops!, ksp::AbstractKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        KSPComputeOperatorsFn{$PetscLib, $PetscInt}(),
        LibPETSc.PetscErrorCode,
        (CKSP, CMat, CMat, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computeops! = ops!
    LibPETSc.KSPSetComputeOperators($PetscLib, ksp, fptr, state_pointer(ksp))
    return ksp
end

"""
    solution(ksp::AbstractKSP)

The solution vector of `ksp`, as a `PetscVec`.

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetSolution"))
"""
solution(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetSolution(getlib(PetscLib), ksp)


"""
    type_name(ksp::AbstractKSP)

The name PETSc knows this solver by, as a `Symbol` (`:gmres`, `:cg`, …), or
`nothing` when no type has been set yet (docs/src/man/naming.md §3.1). v0.4
answered with a `String`; that is a break with no shim (§16).

# External Links
$(doc_external("KSP/KSPGetType"))
"""
type_name(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.KSPGetType(getlib(PetscLib), ksp))

"""
    set_type!(ksp::AbstractKSP, type::Symbol)

Set the Krylov method, for example `:gmres`, `:cg` or `:preonly`.

# External Links
$(doc_external("KSP/KSPSetType"))
"""
function set_type!(ksp::AbstractKSP{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.KSPSetType(getlib(PetscLib), ksp, String(type))
    return ksp
end

"""
    iteration_number(ksp::AbstractKSP)

The number of iterations the last [`solve!`](@ref) took.

# External Links
$(doc_external("KSP/KSPGetIterationNumber"))
"""
iteration_number(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetIterationNumber(getlib(PetscLib), ksp)

"""
    converged_reason(ksp::AbstractKSP)

Why the last [`solve!`](@ref) stopped, as a `LibPETSc.KSPConvergedReason`: positive
when it converged, negative when it diverged.

# External Links
$(doc_external("KSP/KSPGetConvergedReason"))
"""
converged_reason(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetConvergedReason(getlib(PetscLib), ksp)

"""
    set_operators!(ksp::AbstractKSP, A::AbstractPetscMat, P::AbstractPetscMat = A)

Solve with the operator `A`, building the preconditioner from `P`. Returns `ksp`.

`ksp` takes a reference to both matrices, so they stay valid inside it after the
caller destroys its own handles.

# External Links
$(doc_external("KSP/KSPSetOperators"))
"""
function set_operators!(
    ksp::AbstractKSP{PetscLib},
    A::AbstractPetscMat{PetscLib},
    P::AbstractPetscMat{PetscLib} = A,
) where {PetscLib}
    LibPETSc.KSPSetOperators(getlib(PetscLib), ksp, A, P)
    return ksp
end

const DM_ACTIVE_PARTS = (
    operator = LibPETSc.KSP_DMACTIVE_OPERATOR,
    rhs = LibPETSc.KSP_DMACTIVE_RHS,
    initial_guess = LibPETSc.KSP_DMACTIVE_INITIAL_GUESS,
    all = LibPETSc.KSP_DMACTIVE_ALL,
)

"""
    set_dm_active!(ksp::AbstractKSP, flag::Bool)
    set_dm_active!(ksp::AbstractKSP, part::Symbol, flag::Bool)

Whether the DM attached to `ksp` computes its operator, right-hand side and
initial guess. The first form sets all three; the second sets one, `part` being
`:operator`, `:rhs` or `:initial_guess` (`:all` is the first form). With `false`,
`ksp` keeps the DM for its geometry, for example for multigrid, while
[`set_operators!`](@ref) and [`solve!`](@ref) supply the rest. Returns `ksp`.

# External Links
$(doc_external("KSP/KSPSetDMActive"))
"""
function set_dm_active!(ksp::AbstractKSP{PetscLib}, part::Symbol, flag::Bool) where {PetscLib}
    haskey(DM_ACTIVE_PARTS, part) || throw(ArgumentError(
        "part must be one of $(join(repr.(keys(DM_ACTIVE_PARTS)), ", ")), got $(repr(part))",
    ))
    LibPETSc.KSPSetDMActive(getlib(PetscLib), ksp, DM_ACTIVE_PARTS[part], LibPETSc.PetscBool(flag))
    return ksp
end

set_dm_active!(ksp::AbstractKSP, flag::Bool) = set_dm_active!(ksp, :all, flag)
