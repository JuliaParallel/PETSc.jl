
const CKSP = Ptr{Cvoid}
const CKSPType = Cstring

abstract type AbstractKSP{T, PetscLib} <: Factorization{T} end

Base.@kwdef mutable struct KSP{T, PetscLib} <: AbstractKSP{T, PetscLib}
    ptr::CKSP = C_NULL
    opts::Options{PetscLib}
    # Stuff to keep around so that they don't get gc'ed
    _A = nothing
    _P = nothing
    _dm = nothing
    # Function pointers
    ComputeRHS! = nothing
    ComputeOperators! = nothing
end

struct WrappedKSP{T, PetscLib} <: AbstractKSP{T, PetscLib}
    ptr::CKSP
end

scalartype(::KSP{T}) where {T} = T
Base.eltype(::KSP{T}) where {T} = T

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

    function KSP{$PetscScalar}(comm::MPI.Comm; kwargs...)
        @assert initialized($petsclib)
        opts = Options($petsclib, kwargs...)
        ksp = KSP{$PetscScalar, $PetscLib}(opts=opts)
        with(ksp.opts) do
          @chk ccall((:KSPCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CKSP}), comm, ksp)
        end
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
        ksp._A = A
        ksp._P = P
        return nothing
    end

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

    function KSPSetDM!(ksp::KSP{$PetscScalar}, dm::AbstractDM{$PetscLib})
        with(ksp.opts) do
            @chk ccall((:KSPSetDM, $libpetsc), PetscErrorCode, (CKSP, CDM), ksp, dm)
        end
        ksp._dm = dm
        return nothing
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

    function setfromoptions!(ksp::KSP{$PetscScalar})
        with(ksp.opts) do
            @chk ccall((:KSPSetFromOptions, $libpetsc), PetscErrorCode, (CKSP,), ksp)
        end
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


"""
    KSP(A, P; options...)

Construct a PETSc Krylov subspace solver.

Any PETSc options prefixed with `ksp_` and `pc_` can be passed as keywords.
"""
function KSP(A::AbstractMat{T}, P::AbstractMat{T}=A; kwargs...) where {T}
    ksp = KSP{T}(getcomm(A); kwargs...)
    setoperators!(ksp, A, P)
    setfromoptions!(ksp)
    return ksp
end

"""
    KSP(da::AbstractDM; options...)

Construct a PETSc Krylov subspace solver from the distributed mesh

Any PETSc options prefixed with `ksp_` and `pc_` can be passed as keywords.

see [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetDM.html)
"""
function KSP(dm::AbstractDM{PetscLib}; kwargs...) where {PetscLib}
    T = scalartype(PetscLib)
    ksp = KSP{T}(getcomm(dm); kwargs...)
    KSPSetDM!(ksp, dm)
    setfromoptions!(ksp)
    return ksp
end

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

