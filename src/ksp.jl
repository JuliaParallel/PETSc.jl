const CKSP = Ptr{Cvoid}
const CKSPType = Cstring

abstract type AbstractKSP{PetscLib, PetscScalar} <: Factorization{PetscScalar} end

"""
    KSPPtr(petsclib, ksp::CKS, own)

Container type for a PETSc KSP that is just a raw pointer.
"""
mutable struct KSPPtr{PetscLib, PetscScalar} <:
               AbstractKSP{PetscLib, PetscScalar}
    ptr::CKSP
    age::Int
end

mutable struct KSP{PetscLib, PetscScalar} <: AbstractKSP{PetscLib, PetscScalar}
    ptr::CKSP
    opts::Options{PetscLib}
    age::Int
    computerhs! # Needed for KSPSetComputeRHS
    computeops! # Needed for KSPSetComputeOperators
    function KSP{PetscLib}(comm, opts) where {PetscLib}
        PetscScalar = PetscLib.PetscScalar
        ksp = new{PetscLib, PetscScalar}(
            C_NULL,
            opts,
            getlib(PetscLib).age,
            nothing,
            nothing,
        )
        with(ksp.opts) do
            LibPETSc.KSPCreate(PetscLib, comm, ksp)
        end

        # If there is only one rank we can finalize the KSP with GC
        if MPI.Comm_size(comm) == 1
            finalizer(destroy, ksp)
        end

        return ksp
    end
end

include("ksp_wrapped.jl")   


function setfromoptions!(ksp::AbstractKSP{PetscLib}) where {PetscLib}
    with(ksp.opts) do
        LibPETSc.KSPSetFromOptions(PetscLib, ksp)
    end
end

"""
    KSP(A::AbstractMat, P::AbstractMat{PetscLib} = A; options...)

Create a `KSP` using the matrix `A` and preconditioner construction matrix `P`
with the `options`.

The communicator is obtained from `A` and if it has size `1` then the garbage
collector is set, otherwise the user is responsible for calling
[`destroy`](@ref).

# External Links
$(_doc_external("KSP/KSPCreate"))
$(_doc_external("KSP/KSPSetOperators"))
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(
    A::AbstractMat{PetscLib},
    P::AbstractMat{PetscLib} = A;
    options...,
) where {PetscLib}
    @assert initialized(PetscLib)

    ksp = KSP{PetscLib}(getcomm(A), Options(PetscLib; options...))

    KSPSetOperators(ksp,A,P)

    setfromoptions!(ksp)

    return ksp
end


"""
    KSP([petsclib,] A::SparseMatrixCSC; options...)

Create a [`KSP`](@ref) with the sparse matrix `A` using the `petsclib`. If
`petsclib` is not given, the default library will be used`.
"""
KSP(petsclib, A::SparseMatrixCSC; kwargs...) =
    KSP(MatSeqAIJ(petsclib, A); kwargs...)
function KSP(A::SparseMatrixCSC{PetscScalar}; kwargs...) where {PetscScalar}
    KSP(MatSeqAIJ(getlib(; PetscScalar = PetscScalar), A); kwargs...)
end

"""
    KSP(da::AbstractDM; options...)

Construct a PETSc Krylov subspace solver from the distributed mesh

Any PETSc options prefixed with `ksp_` and `pc_` can be passed as keywords.

# External Links
$(_doc_external("KSP/KSPCreate"))
$(_doc_external("KSP/KSPSetDM"))
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(dm::AbstractDM{PetscLib}; options...) where {PetscLib}
    @assert initialized(PetscLib)

    comm = getcomm(dm)

    ksp = KSP{PetscLib}(comm, Options(PetscLib; options...))

    with(options) do
        KSPSetDM(ksp, dm)
    end

    setfromoptions!(ksp)

    return ksp
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
    new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)
    b = VecPtr(PetscLib, cb, false)
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
setcomputerhs!(ksp::AbstractKSP, rhs!) = setcomputerhs!(rhs!, ksp)
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function setcomputerhs!(rhs!, ksp::KSP{$PetscLib})
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
    new_ksp = KSPPtr{PetscLib, PetscScalar}(new_ksp_ptr, getlib(PetscLib).age)
    A = MatPtr{PetscLib, PetscScalar}(cA, getlib(PetscLib).age)
    P = MatPtr{PetscLib, PetscScalar}(cP, getlib(PetscLib).age)
    ksp = unsafe_pointer_to_objref(ksp_ptr)
    ierr = ksp.computeops!(A, P, new_ksp)
    return PetscLib.PetscInt(ierr)
end

"""
    setcomputeoperators!(ksp::AbstractKSP, ops!::Function)
    setcomputeoperators!(ops!::Function, ksp::AbstractKSP)

Define `ops!` to be the compute operators function for the `ksp`. A call to
`ops!(A, P, new_ksp)` should set the elements of the PETSc matrix linear
operator `A` and preconditioning matrix `P` based on the `new_ksp`.

!!! note

    The `new_ksp` passed to `ops!` may not be the same as the `ksp` passed to
    `setcomputeoperators!`.

# External Links
$(_doc_external("KSP/KSPSetComputeOperators"))
"""
setcomputeoperators!(ksp::AbstractKSP, ops!) = setcomputeoperators!(ops!, ksp)
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function setcomputeoperators!(ops!, ksp::KSP{$PetscLib})
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
    getDMDA(ksp::AbstractKSP)

Get `dmda` for `ksp`

The returned `dmda` is owned by the `ksp`

# External Links
$(_doc_external("KSP/KSPGetDM"))
"""
function getDMDA(ksp::AbstractKSP{PetscLib}) where PetscLib
    #t_dmda = Ref{CDM}()
    #LibPETSc.KSPGetDM(PetscLib, ksp, t_dmda)
    #dmda = DMDAPtr{PetscLib}(t_dmda[], getlib(PetscLib).age, false)
    dmda = KSPGetDM(ksp)
    return dmda
end


function destroy(ksp::AbstractKSP{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) &&
       ksp.age == getlib(PetscLib).age &&
       ksp.ptr != C_NULL
        LibPETSc.KSPDestroy(PetscLib, ksp)
    end
    ksp.ptr = C_NULL
    return nothing
end

function solve!(
    x::AbstractVec{PetscLib},
    ksp::AbstractKSP{PetscLib},
    b::AbstractVec{PetscLib},
) where {PetscLib}
    with(ksp.opts) do
        LibPETSc.KSPSolve(PetscLib, ksp, b, x)
    end
    return x
end

function solve!(
    ksp::AbstractKSP{PetscLib},
) where {PetscLib}
    with(ksp.opts) do
        LibPETSc.KSPSolve(PetscLib, ksp, C_NULL, C_NULL)
    end
    return ksp
end

"""
    createvecs(ksp::AbstractKSP; nright = 0, nleft = 0)

Create `nright` right and `nleft` left vectors compatible with the `ksp`.
Returned object `V` has `Tuple` members `V.right` and `V.left` containing the
vectors.

# External Links
$(_doc_external("KSP/KSPCreateVecs"))
"""
function createvecs(
    ksp::AbstractKSP{PetscLib};
    nright = 0,
    nleft = 0,
) where {PetscLib}
    # pointer of pointers to the base vectors
    r_right_vs = Ref{Ptr{CVec}}()
    r_left_vs = Ref{Ptr{CVec}}()

    # create 1 right and left vector
    LibPETSc.KSPCreateVecs(PetscLib, ksp, 1, r_right_vs, 1, r_left_vs)

    # create right vectors
    a_v = unsafe_wrap(Array, r_right_vs[], 1; own = false)
    v = VecPtr(PetscLib, a_v[1], false)
    right = ntuple(i -> similar(v), nright)

    # create left vectors
    a_v = unsafe_wrap(Array, r_left_vs[], 1; own = false)
    v = VecPtr(PetscLib, a_v[1], false)
    left = ntuple(i -> similar(v), nleft)

    LibPETSc.VecDestroyVecs(PetscLib, 1, r_right_vs)
    LibPETSc.VecDestroyVecs(PetscLib, 1, r_left_vs)

    (right = right, left = left)
end

function LinearAlgebra.ldiv!(x::AbstractVec, ksp::AbstractKSP, b::AbstractVec)
    solve!(x, ksp, b)
end

function Base.:\(ksp::AbstractKSP, b::AbstractVec)
    x = createvecs(ksp; nleft = 1).left[1]
    ldiv!(x, ksp, b)
    return x
end

function Base.:\(
    ksp::AbstractKSP{PetscLib},
    b::Vector{PetscScalar},
) where {PetscLib, PetscScalar}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert MPI.Comm_size(getcomm(ksp)) == 1

    petsc_b = VecSeqWithArray(PetscLib, b)
    petsc_x = createvecs(ksp; nleft = 1).left[1]

    ldiv!(petsc_x, ksp, petsc_b)

    x = similar(b, length(petsc_x))

    withlocalarray!(petsc_x; read = true, write = false) do y
        x .= y
    end

    destroy(petsc_b)
    destroy(petsc_x)

    return x
end


#=
#
# OLD WRAPPERS
#
struct WrappedKSP{T, PetscLib} <: AbstractKSP{T, PetscLib}
    ptr::CKSP
end

LinearAlgebra.transpose(ksp) = LinearAlgebra.Transpose(ksp)
LinearAlgebra.adjoint(ksp) = LinearAlgebra.Adjoint(ksp)

"""
    KSPSetComputeRHS!(
        ksp::KSP{Number},
        ComputeRHS!,
        ctx = C_NULL,
    )

Set the right-hand side function `ComputeRHS!` for the `ksp` using the user
`ctx`. `ComputeRHS!` should be callable with three arguments of type
`(::KSP{Number}, ::Vec, ::Ptr)`;
see [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetComputeRHS.html)
"""
function KSPSetComputeRHS! end

"""
    struct Fn_KSPComputeRHS{T} end

Type used to wrap `ComputeRHS!` functions in KSP
"""
struct Fn_KSPComputeRHS{T} end

"""
    KSPSetComputeOperators!(
        ksp::KSP{Number},
        ComputeOperators!,
        ctx = C_NULL,
    )

Set the linear operators function `ComputeOperators!` for the `ksp` using the
user `ctx`. `ComputeOperators!` should be callable with four arguments of type
`(::KSP{Number}, ::Mat, ::Mat, ::Ptr)`;
see [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetComputeOperators.html)
"""
function KSPSetComputeOperators! end

"""
    struct Fn_KSPComputeOperators{T} end

Type used to wrap `ComputeOperators!` functions in KSP
"""
struct Fn_KSPComputeOperators{T} end

@for_libpetsc begin

    function (::Fn_KSPComputeRHS{$PetscScalar})(
        new_ksp_ptr::CKSP,
        cb::CVec,
        ksp_ptr::Ptr{Cvoid}
       )::$PetscInt
        ksp = unsafe_pointer_to_objref(ksp_ptr)
        new_ksp = WrappedKSP{$PetscScalar, $PetscLib}(new_ksp_ptr)
        b = Vec{$PetscScalar}(cb)
        ierr = ksp.ComputeRHS!(new_ksp, b)
        return $PetscInt(ierr)
    end

    function KSPSetComputeRHS!(
        ksp::KSP{$PetscScalar},
        ComputeRHS!
    )
        ptr_ksp = pointer_from_objref(ksp)
        fptr = @cfunction(Fn_KSPComputeRHS{$PetscScalar}(),
                          $PetscInt,
                          (CKSP, CVec, Ptr{Cvoid}))
        with(ksp.opts) do
            @chk ccall((:KSPSetComputeRHS, $libpetsc), PetscErrorCode,
                (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
                ksp, fptr, ptr_ksp)
        end
        ksp.ComputeRHS! = ComputeRHS!
        return ksp
    end

    function (::Fn_KSPComputeOperators{$PetscScalar})(
        new_ksp_ptr::CKSP,
        cA::CMat,
        cP::CMat,
        ksp_ptr::Ptr{Cvoid}
       )::$PetscInt
        ksp = unsafe_pointer_to_objref(ksp_ptr)
        new_ksp = WrappedKSP{$PetscScalar, $PetscLib}(new_ksp_ptr)
        A = Mat{$PetscScalar}(cA)
        P = Mat{$PetscScalar}(cP)
        ierr = ksp.ComputeOperators!(new_ksp, A, P)
        return $PetscInt(ierr)
    end

    function KSPSetComputeOperators!(
        ksp::KSP{$PetscScalar},
        ComputeOperators!
    )
        ptr_ksp = pointer_from_objref(ksp)
        fptr = @cfunction(Fn_KSPComputeOperators{$PetscScalar}(),
                          $PetscInt,
                          (CKSP, CMat, CMat, Ptr{Cvoid}))
        with(ksp.opts) do
            @chk ccall((:KSPSetComputeOperators, $libpetsc), PetscErrorCode,
                (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
                ksp, fptr, ptr_ksp)
        end
        ksp.ComputeOperators! = ComputeOperators!
        return ksp
    end

    function DMDA(ksp::AbstractKSP{$PetscScalar})
        t_dm = Ref{CDM}()
        @chk ccall(
            (:KSPGetDM, $libpetsc),
            PetscErrorCode,
            (CKSP, Ptr{CDM}),
            ksp,
            t_dm,
        )
        dm = DMDA{$PetscLib}(t_dm[])
        return dm
    end

    function settolerances!(ksp::KSP{$PetscScalar}; rtol=PETSC_DEFAULT, atol=PETSC_DEFAULT, divtol=PETSC_DEFAULT, max_it=PETSC_DEFAULT)
        @chk ccall((:KSPSetTolerances, $libpetsc), PetscErrorCode, 
                    (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
                    ksp, rtol, atol, divtol, max_it)
        return nothing
    end


    function iters(ksp::KSP{$PetscScalar})
        r_its = Ref{$PetscInt}()
        @chk ccall((:KSPGetIterationNumber, $libpetsc), PetscErrorCode, 
        (KSP, Ptr{$PetscInt}), ksp, r_its)
        return r_its[]
    end

    function view(ksp::KSP{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, getcomm(ksp)))
        @chk ccall((:KSPView, $libpetsc), PetscErrorCode, 
                    (CKSP, CPetscViewer),
                ksp, viewer);
        return nothing
    end

    function resnorm(ksp::KSP{$PetscScalar})
        r_rnorm = Ref{$PetscReal}()
        @chk ccall((:KSPGetResidualNorm, $libpetsc), PetscErrorCode, 
        (KSP, Ptr{$PetscReal}), ksp, r_rnorm)
        return r_rnorm[]
    end

    function solve!(ksp::KSP{$PetscScalar})
        with(ksp.opts) do
            @chk ccall((:KSPSolve, $libpetsc), PetscErrorCode, 
            (CKSP, CVec, CVec), ksp, C_NULL, C_NULL)
        end
        return nothing
    end

    function solve!(x::AbstractVec{$PetscScalar}, tksp::Transpose{T,K}, b::AbstractVec{$PetscScalar}) where {T,K <: KSP{$PetscScalar}}
        ksp = parent(tksp)
        with(ksp.opts) do
            @chk ccall((:KSPSolveTranspose, $libpetsc), PetscErrorCode, 
            (CKSP, CVec, CVec), ksp, b, x)
        end
        return x
    end

end

# no generic Adjoint solve defined, but for Real we can use Adjoint
solve!(x::AbstractVec{T}, aksp::Adjoint{T,K}, b::AbstractVec{T}) where {K <: KSP{T}} where {T<:Real} =
    solve!(x, transpose(parent(aksp)), b)

const KSPAT{T, LT} = Union{KSP{T, LT}, Transpose{T, KSP{T, LT}}, Adjoint{T, KSP{T, LT}}}

LinearAlgebra.ldiv!(x::AbstractVec{T}, ksp::KSPAT{T, LT}, b::AbstractVec{T}) where {T, LT} = solve!(x, ksp, b)
function LinearAlgebra.ldiv!(x::AbstractVector{T}, ksp::KSPAT{T, LT}, b::AbstractVector{T}) where {T, LT}
    parent(solve!(AbstractVec(x), ksp, AbstractVec(b)))
end
Base.:\(ksp::KSPAT{T, LT}, b::AbstractVector{T}) where {T, LT} = ldiv!(similar(b), ksp, b)

Base.show(io::IO, ksp::KSP) = _show(io, ksp)

"""
    iters(ksp::KSP)

Gets the current iteration number; if the `solve!` is complete, returns the number of iterations used.

# External Links
$(_doc_external("KSP/KSPGetIterationNumber"))
"""
iters

"""
    resnorm(ksp::KSP)

Gets the last (approximate preconditioned) residual norm that has been computed.

# External Links
$(_doc_external("KSP/KSPGetResidualNorm"))
"""
resnorm
=#
