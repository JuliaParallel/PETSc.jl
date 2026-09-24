import .LibPETSc: AbstractKSP, CKSP, KSP, AbstractPetscDM   # KSP methods below are constructors of LibPETSc.KSP

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
    return KSP(M; kwargs...)
end


function solve!(
    x::PetscVec{PetscLib},
    ksp::KSP{PetscLib},
    b::PetscVec{PetscLib},
) where {PetscLib}
    has_opts = !isnothing(ksp.opts)
    has_opts && push!(ksp.opts)
    try
        LibPETSc.KSPSolve(PetscLib, ksp, b, x)
    finally
        has_opts && pop!(ksp.opts)
    end
    return nothing
end

function solve!(
    ksp::AbstractKSP{PetscLib},
) where {PetscLib}
    has_opts = hasproperty(ksp, :opts) && !isnothing(ksp.opts)
    has_opts && push!(ksp.opts)
    try
        LibPETSc.KSPSolve(PetscLib, ksp, C_NULL, C_NULL)
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

    petsc_b = LibPETSc.VecCreateSeqWithArray(getlib(PetscLib), c, PetscInt(1), PetscInt(length(b)), PetscScalar.(b))
    petsc_x = ksp \ petsc_b
    x = petsc_x[:]
    destroy!(petsc_b)
    destroy!(petsc_x)

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
)::PetscInt where {PetscLib, PetscInt}
    PetscScalar = PetscLib.PetscScalar
    #new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)\
    #b = VecPtr(PetscLib, cb, false)
    new_ksp = KSP{PetscLib}(new_ksp_ptr, 0; own = false)
    b = PetscVec{PetscLib}(cb, 0; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)
    ierr = ksp.computerhs!(b, new_ksp)
    return PetscLib.PetscInt(ierr)
end

"""
    set_compute_rhs!(rhs!::Function, ksp::AbstractKSP)

Define `rhs!` to be the right-hand side function of the `ksp`. A call to
`rhs!(b, new_ksp)` should set the elements of the PETSc vector `b` based on the
`new_ksp`.

!!! note

    The `new_ksp` passed to `rhs!` may not be the same as the `ksp` passed to
    `set_compute_rhs!`.

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
        $PetscInt,
        (CKSP, CVec, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computerhs! = rhs!
    LibPETSc.KSPSetComputeRHS($PetscLib, ksp, fptr, pointer_from_objref(ksp))
    return ksp
end

# Wrapper for calls to set_compute_rhs!
mutable struct KSPComputeOperatorsFn{PetscLib, PetscInt} end
function (w::KSPComputeOperatorsFn{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cA::CMat,
    cP::CMat,
    ksp_ptr::Ptr{Cvoid},
)::PetscInt where {PetscLib, PetscInt}
    PetscScalar = PetscLib.PetscScalar
    #new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)
    new_ksp = KSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age; own = false)
    A = PetscMat{PetscLib}(cA, getlib(PetscLib).age; own = false)
    P = PetscMat{PetscLib}(cP, getlib(PetscLib).age; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)
    ierr = ksp.computeops!(A, P, new_ksp)
    return PetscLib.PetscInt(ierr)
end

"""
    set_compute_operators!(ops!::Function, ksp::KSP)

Define `ops!` to be the compute operators function for the `ksp`. A call to
`ops!(A, P, new_ksp)` should set the elements of the PETSc matrix linear
operator `A` and preconditioning matrix `P` based on the `new_ksp`.

!!! note

    The `new_ksp` passed to `ops!` may not be the same as the `ksp` passed to
    `set_compute_operators!`.

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
        $PetscInt,
        (CKSP, CMat, CMat, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computeops! = ops!
    LibPETSc.KSPSetComputeOperators($PetscLib, ksp, fptr, pointer_from_objref(ksp))
    return ksp
end

"""
    sol = solution(ksp::AbstractKSP)

Returns the solution vector associated with the KSP object.

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetSolution"))
"""
function solution(ksp::AbstractKSP{PetscLib}) where PetscLib
    petsclib = getlib(PetscLib)
    sol = LibPETSc.KSPGetSolution(petsclib, ksp)
    return VecPtr(petsclib, sol.ptr, false)   # owned by the KSP
end


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
    return nothing
end
