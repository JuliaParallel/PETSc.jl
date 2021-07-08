
const CKSP = Ptr{Cvoid}
const CKSPType = Cstring

abstract type AbstractKSP{T, PetscLib} <: Factorization{T} end

mutable struct KSP{T, PetscLib} <: AbstractKSP{T, PetscLib}
    ptr::CKSP
    __comm__::MPI.Comm # Do not access directly use `getcomm(ksp)`
    # keep around so that they don't get gc'ed
    A  # Operator
    P  # preconditioning operator
    opts::Options{PetscLib}
end

scalartype(::KSP{T}) where {T} = T

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CKSP}, obj::KSP) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CKSP}}, obj::KSP) =
    convert(Ptr{CKSP}, pointer_from_objref(obj))

Base.eltype(::KSP{T}) where {T} = T
LinearAlgebra.transpose(ksp) = LinearAlgebra.Transpose(ksp)
LinearAlgebra.adjoint(ksp) = LinearAlgebra.Adjoint(ksp)

@for_libpetsc begin

    function KSP{$PetscScalar}(comm::MPI.Comm; kwargs...)
        @assert initialized($petsclib)
        opts = Options($petsclib, kwargs...)
        ksp = KSP{$PetscScalar, $PetscLib}(C_NULL, comm, nothing, nothing, opts)
        @chk ccall((:KSPCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CKSP}), comm, ksp)
        if comm == MPI.COMM_SELF
            finalizer(destroy, ksp)
        end
        return ksp
    end

    function destroy(ksp::KSP{$PetscScalar})
        finalized($petsclib) ||
        @chk ccall((:KSPDestroy, $libpetsc), PetscErrorCode, (Ptr{CKSP},), ksp)
        return nothing
    end

    function setoperators!(ksp::KSP{$PetscScalar}, A::AbstractMat{$PetscScalar}, P::AbstractMat{$PetscScalar})
        @chk ccall((:KSPSetOperators, $libpetsc), PetscErrorCode, (CKSP, CMat, CMat), ksp, A, P)
        ksp.A = A
        ksp.P = P
        return nothing
    end

    function settolerances!(ksp::KSP{$PetscScalar}; rtol=PETSC_DEFAULT, atol=PETSC_DEFAULT, divtol=PETSC_DEFAULT, max_it=PETSC_DEFAULT)
        @chk ccall((:KSPSetTolerances, $libpetsc), PetscErrorCode, 
                    (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
                    ksp, rtol, atol, divtol, max_it)
        return nothing
    end

    function setfromoptions!(ksp::KSP{$PetscScalar})
        @chk ccall((:KSPSetFromOptions, $libpetsc), PetscErrorCode, (CKSP,), ksp)
    end

    function gettype(ksp::KSP{$PetscScalar})
        t_r = Ref{CKSPType}()
        @chk ccall((:KSPGetType, $libpetsc), PetscErrorCode, (CKSP, Ptr{CKSPType}), ksp, t_r)
        return unsafe_string(t_r[])
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

    function solve!(x::AbstractVec{$PetscScalar}, ksp::KSP{$PetscScalar}, b::AbstractVec{$PetscScalar})
        with(ksp.opts) do
            @chk ccall((:KSPSolve, $libpetsc), PetscErrorCode, 
            (CKSP, CVec, CVec), ksp, b, x)
        end
        return x
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


"""
    KSP(A, P; options...)

Construct a PETSc Krylov subspace solver.

Any PETSc options prefixed with `ksp_` and `pc_` can be passed as keywords.
"""
function KSP(A::AbstractMat{T}, P::AbstractMat{T}=A; kwargs...) where {T}
    ksp = KSP{T}(getcomm(A); kwargs...)
    setoperators!(ksp, A, P)
    with(ksp.opts) do
        setfromoptions!(ksp)
    end
    return ksp
end

Base.show(io::IO, ksp::KSP) = _show(io, ksp)


"""
    iters(ksp::KSP)

Gets the current iteration number; if the `solve!` is complete, returns the number of iterations used.

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetIterationNumber.html
"""
iters


"""
    resnorm(ksp::KSP)

Gets the last (approximate preconditioned) residual norm that has been computed.

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetResidualNorm.html
"""
resnorm

