
const CSNES = Ptr{Cvoid}
const CSNESType = Cstring


abstract type AbstractSNES{T, PetscLib} end

mutable struct SNES{T, PetscLib} <: AbstractSNES{T, PetscLib}
    ptr::CSNES
    opts::Options{PetscLib}
    fn!
    fn_vec
    update_jac!
    jac_A
    jac_P
    user_ctx    # Useful to transfer vectors and dms into the residual function
end

scalartype(::AbstractSNES{T}) where {T} = T

Base.eltype(::AbstractSNES{T}) where {T} = T

# How to handle Jacobians?
#  - https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESComputeJacobianDefault.html
#  - https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESComputeJacobianDefaultColor.html
#  -

struct SNESFn{T}
end

struct SNESJac{T}
end


#=
function _snesfn(csnes::CSNES, cx::CVec, cfx::CVec, ctx::Ptr{Cvoid})
    snes = unsafe_pointer_to_objref(ctx)
    snes.Feval(cfx, cx)
end

function _snesjac(csnes::CSNES, cx::CVec, cAmat::CMat, cPmat::CMat, ctx::Ptr{Cvoid})
    snes = unsafe_pointer_to_objref(ctx)
    snes.Jeval(cAmat, cPmat, cx)
end
=#
"""
    SNES{PetscScalar}(
        ::UnionPetscLib,
        comm::MPI.Comm; 
        snessetfromoptions = true,
        options...)

Initializes a `SNES` nonlinear solver object
"""
function SNES() end

@for_petsc  function SNES{$PetscScalar}(
        ::$UnionPetscLib,
        comm::MPI.Comm; 
        snessetfromoptions = true,
        options...)
        
        @assert initialized($petsclib)
        opts = Options($petsclib, options...)
        snes = SNES{$PetscScalar, $PetscLib}(C_NULL, opts, nothing, nothing, nothing, nothing, nothing, nothing)
        
        with(snes.opts) do
        
            @chk ccall(
                (:SNESCreate, $petsc_library), 
                PetscErrorCode, 
                (MPI.MPI_Comm, 
                Ptr{CSNES}), 
                comm,
                snes)

            snessetfromoptions && setfromoptions!(snes)

        end

        if comm == MPI.COMM_SELF
            finalizer(destroy, snes)
        end
        return snes
    end

@for_libpetsc begin

    function (::SNESFn{$PetscScalar})(csnes::CSNES, cx::CVec, cfx::CVec, ctx::Ptr{Cvoid})::$PetscInt
        snes = unsafe_pointer_to_objref(ctx)
        #x = unsafe_localarray($PetscScalar, cx; write=false)
        #fx = unsafe_localarray($PetscScalar, cfx; read=false)
        #snes.fn!(fx, x, snes.user_ctx)
        snes.fn!(cfx, cx, snes.user_ctx)
        #Base.finalize(x)
        #Base.finalize(fx)
        return $PetscInt(0)
    end

    function setfunction!(snes::AbstractSNES{$PetscScalar}, fn!, vec::AbstractVec{$PetscScalar})
        ctx = pointer_from_objref(snes)
        fptr = @cfunction(SNESFn{$PetscScalar}(), $PetscInt, (CSNES, CVec, CVec, Ptr{Cvoid}))
        with(snes.opts) do
            @chk ccall((:SNESSetFunction, $libpetsc), PetscErrorCode,
                (CSNES, CVec, Ptr{Cvoid}, Ptr{Cvoid}),
                snes, vec, fptr, ctx)
        end
        snes.fn_vec = vec
        snes.fn! = fn!
        return nothing
    end

    function destroy(snes::AbstractSNES{$PetscScalar})
        if snes.age == getlib(PetscLib).age && !(finalized(PetscLib)) && snes.ptr != C_NULL
            @chk ccall((:SNESDestroy, $libpetsc), PetscErrorCode, (Ptr{CSNES},), snes)
        end
        snes.ptr = C_NULL
        return nothing
    end

    function setfromoptions!(snes::AbstractSNES{$PetscScalar})
        @chk ccall((:SNESSetFromOptions, $libpetsc), PetscErrorCode, (CSNES,), snes)
    end

    function gettype(snes::AbstractSNES{$PetscScalar})
        t_r = Ref{CSNESType}()
        @chk ccall((:SNESGetType, $libpetsc), PetscErrorCode, (CSNES, Ptr{CSNESType}), snes, t_r)
        return unsafe_string(t_r[])
    end

    function view(snes::AbstractSNES{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, getcomm(snes)))
        @chk ccall((:SNESView, $libpetsc), PetscErrorCode,
                    (CSNES, CPetscViewer),
                snes, viewer);
        return nothing
    end



    function (::SNESJac{$PetscScalar})(csnes::CSNES, cx::CVec, cA::CMat, cP::CMat, ctx::Ptr{Cvoid})::$PetscInt
        snes = unsafe_pointer_to_objref(ctx)
        @assert snes.ptr == csnes
        @assert snes.jac_A.ptr == cA
        @assert snes.jac_P.ptr == cP
        #x = unsafe_localarray($PetscScalar, cx; write=false)
        #snes.update_jac!(x, snes.jac_A, snes.jac_P,snes.user_ctx)
        snes.update_jac!(cx, snes.jac_A, snes.jac_P,snes.user_ctx)
        #Base.finalize(x)
        return $PetscInt(0)
    end

    function setjacobian!(snes::AbstractSNES{$PetscScalar}, update_jac!, A::AbstractMat{$PetscScalar}, P::AbstractMat{$PetscScalar}=A)
        ctx = pointer_from_objref(snes)
        jacptr = @cfunction(SNESJac{$PetscScalar}(), $PetscInt, (CSNES, CVec, CMat, CMat, Ptr{Cvoid}))

        with(snes.opts) do
            @chk ccall((:SNESSetJacobian, $libpetsc), PetscErrorCode,
                (CSNES, CMat, CMat, Ptr{Cvoid}, Ptr{Cvoid}),
                snes, A, P, jacptr, ctx)
        end
        snes.update_jac! = update_jac!
        snes.jac_A = A
        snes.jac_P = P
        return nothing
    end


    function solve!(x::AbstractVec{$PetscScalar}, snes::AbstractSNES{$PetscScalar}, b::AbstractVec{$PetscScalar})
        with(snes.opts) do
            @chk ccall((:SNESSolve, $libpetsc), PetscErrorCode,
            (CSNES, CVec, CVec), snes, b, x)
        end
        return x
    end
    function solve!(x::AbstractVec{$PetscScalar}, snes::AbstractSNES{$PetscScalar})
        with(snes.opts) do
            @chk ccall((:SNESSolve, $libpetsc), PetscErrorCode,
            (CSNES, CVec, CVec), snes, C_NULL, x)
        end
        return x
    end

end

solve!(x::AbstractVector{T}, snes::AbstractSNES{T}) where {T} = parent(solve!(AbstractVec(x), snes))
