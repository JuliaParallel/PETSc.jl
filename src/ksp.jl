import .LibPETSc: AbstractPetscKSP, CKSP, PetscKSP, AbstractPetscDM

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscKSP{PetscLib}) where {PetscLib}
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
[`destroy`](@ref).

# External Links
$(_doc_external("KSP/KSPCreate"))
$(_doc_external("KSP/KSPSetOperators"))
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(
    A::AbstractPetscMat{PetscLib},
    P::AbstractPetscMat{PetscLib} = A;
    prefix::String="",
    options...
) where {PetscLib}
    @assert initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    comm = getcomm(A)
    ksp = LibPETSc.KSPCreate(petsclib,comm)
    
    LibPETSc.KSPSetOperators(petsclib, ksp, A, P)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, prefix)
    end
    
    # Push options to PETSc options database
    if !isempty(options)
        opts = PETSc.Options(petsclib; options...);
        push!(opts)
        LibPETSc.KSPSetFromOptions(petsclib, ksp)
        pop!(opts)
    end

    return ksp
end

"""
    KSP(comm::MPI.Comm, A::AbstractPetscMat, P::AbstractPetscMat{PetscLib} = A; prefix="", options...)

Create a `KSP` using the matrix `A` and preconditioner construction matrix `P`
with optional `prefix` and `options`.

The communicator is obtained from `A` and if it has size `1` then the garbage
collector is set, otherwise the user is responsible for calling
[`destroy`](@ref).

# External Links
$(_doc_external("KSP/KSPCreate"))
$(_doc_external("KSP/KSPSetOperators"))
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(dm::AbstractPetscDM{PetscLib};
    prefix::String="",
    options...
) where {PetscLib}
    @assert initialized(getlib(PetscLib))
    petsclib = getlib(PetscLib)
    comm = getcomm(dm)
    ksp = LibPETSc.KSPCreate(petsclib,comm)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, prefix)
    end
    
    LibPETSc.KSPSetDM(petsclib, ksp, dm)

    # Push options to PETSc options database
    if !isempty(options)
        opts = PETSc.Options(petsclib; options...);
        push!(opts)
        LibPETSc.KSPSetFromOptions(petsclib, ksp)
        pop!(opts)
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
    M = PETSc.MatCreateSeqAIJ(petsclib, comm, S)
    return KSP(M; kwargs...)
end


function solve!(
    x::PetscVec{PetscLib},
    ksp::PetscKSP{PetscLib},
    b::PetscVec{PetscLib},
) where {PetscLib}
    LibPETSc.KSPSolve(PetscLib, ksp, b, x)
    return nothing
end

function solve!(
    ksp::AbstractPetscKSP{PetscLib},
) where {PetscLib}
    #with(ksp.opts) do
    LibPETSc.KSPSolve(PetscLib, ksp, C_NULL, C_NULL)
    #end
    return ksp
end

LinearAlgebra.ldiv!(x::PetscVec{PetscLib}, ksp::PetscKSP{PetscLib}, b::PetscVec{PetscLib}) where {PetscLib} = solve!(x, ksp, b)

function Base.:\(ksp::PetscKSP, b::PetscVec{PetscLib}) where {PetscLib}
    x = similar(b)
    ldiv!(x, ksp, b)
    return x
end

function Base.:\(
    ksp::PetscKSP{PetscLib},
    b::Vector{PetscScalar},
) where {PetscLib, PetscScalar}
    @assert PetscScalar == PetscLib.PetscScalar
    comm = getcomm(ksp)
    @assert MPI.Comm_size(comm) == 1
    PetscInt = PetscLib.PetscInt

    petsc_b = LibPETSc.VecCreateSeqWithArray(getlib(PetscLib),comm, PetscInt(1), PetscInt(length(b)), PetscScalar.(b))
    petsc_x = ksp \ petsc_b
    x = petsc_x[:]
    destroy(petsc_b)
    destroy(petsc_x)

    return x
end


function destroy(ksp::PetscKSP{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) && ksp.ptr != C_NULL
        LibPETSc.KSPDestroy(PetscLib, ksp)
    end
    ksp.ptr = C_NULL
    return nothing
end



"""
    getDMDA(ksp::AbstractKSP)

Get `dmda` for `ksp`

The returned `dmda` is owned by the `ksp`

# External Links
$(_doc_external("KSP/KSPGetDM"))
"""
function getDMDA(ksp::AbstractPetscKSP{PetscLib}) where PetscLib
    dmda = LibPETSc.KSPGetDM(getlib(PetscLib),ksp)
    return dmda
end

#
# Wrapper for calls to setcomputerhs!
mutable struct Fn_KSPComputeRHS{PetscLib, PetscInt} end
function (w::Fn_KSPComputeRHS{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cb::CVec,
    ksp_ptr::Ptr{Cvoid},
)::PetscInt where {PetscLib, PetscInt}
    PetscScalar = PetscLib.PetscScalar
    #new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)\
    #b = VecPtr(PetscLib, cb, false)
    new_ksp = PetscKSP{PetscLib}(new_ksp_ptr, 0)
    b = PetscVec{PetscLib}(cb, 0)
    ksp = unsafe_pointer_to_objref(ksp_ptr)
    ierr = ksp.computerhs!(b, new_ksp)
    return PetscLib.PetscInt(ierr)
end

"""
    setcomputerhs!(ksp::AbstractKSP, rhs!::Function)
    setcomputerhs!(rhs!::Function, ksp::AbstractKSP)

Define `rhs!` to be the right-hand side function of the `ksp`. A call to
`rhs!(b, new_ksp)` should set the elements of the PETSc vector `b` based on the
`new_ksp`.

!!! note

    The `new_ksp` passed to `rhs!` may not be the same as the `ksp` passed to
    `setcomputerhs!`.

# External Links
$(_doc_external("KSP/KSPSetComputeRHS"))
"""
setcomputerhs!(ksp::AbstractPetscKSP, rhs!) = setcomputerhs!(rhs!, ksp)
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function setcomputerhs!(rhs!, ksp::AbstractPetscKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        Fn_KSPComputeRHS{$PetscLib, $PetscInt}(),
        $PetscInt,
        (CKSP, CVec, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computerhs! = rhs!
    LibPETSc.KSPSetComputeRHS($PetscLib, ksp, fptr, pointer_from_objref(ksp))
    return ksp
end

# Wrapper for calls to setcomputerhs!
mutable struct Fn_KSPComputeOperators{PetscLib, PetscInt} end
function (w::Fn_KSPComputeOperators{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cA::CMat,
    cP::CMat,
    ksp_ptr::Ptr{Cvoid},
)::PetscInt where {PetscLib, PetscInt}
    PetscScalar = PetscLib.PetscScalar
    #new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)
    new_ksp = PetscKSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age)
    A = PetscMat{PetscLib}(cA, getlib(PetscLib).age)
    P = PetscMat{PetscLib}(cP, getlib(PetscLib).age)
    ksp = unsafe_pointer_to_objref(ksp_ptr)
    ierr = ksp.computeops!(A, P, new_ksp)
    return PetscLib.PetscInt(ierr)
end

"""
    setcomputeoperators!(ksp::PetscKSP, ops!::Function)
    setcomputeoperators!(ops!::Function, ksp::PetscKSP)

Define `ops!` to be the compute operators function for the `ksp`. A call to
`ops!(A, P, new_ksp)` should set the elements of the PETSc matrix linear
operator `A` and preconditioning matrix `P` based on the `new_ksp`.

!!! note

    The `new_ksp` passed to `ops!` may not be the same as the `ksp` passed to
    `setcomputeoperators!`.

# External Links
$(_doc_external("KSP/KSPSetComputeOperators"))
"""
setcomputeoperators!(ksp::AbstractPetscKSP, ops!) = setcomputeoperators!(ops!, ksp)
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function setcomputeoperators!(ops!, ksp::AbstractPetscKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        Fn_KSPComputeOperators{$PetscLib, $PetscInt}(),
        $PetscInt,
        (CKSP, CMat, CMat, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computeops! = ops!
    LibPETSc.KSPSetComputeOperators($PetscLib, ksp, fptr, pointer_from_objref(ksp))
    return ksp
end

"""
    sol = get_solution(ksp::AbstractPetscKSP)
Returns the solution vector associated with the KSP object.
"""
function get_solution(ksp::AbstractPetscKSP{PetscLib}) where PetscLib
    sol = LibPETSc.KSPGetSolution(getlib(PetscLib),ksp)
    return sol
end