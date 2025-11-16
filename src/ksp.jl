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