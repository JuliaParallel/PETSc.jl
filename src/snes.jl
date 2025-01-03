const CSNES = Ptr{Cvoid}
const CSNESType = Cstring

abstract type AbstractSNES{PetscLib} end

"""
    SNESPtr{PetscLib}(ptr, age)

Container type for a PETSc SNES that is just a raw pointer.
"""
mutable struct SNESPtr{PetscLib} <: AbstractSNES{PetscLib}
    ptr::CSNES
    age::Int
end

function destroy(snes::AbstractSNES{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) &&
       snes.age == getlib(PetscLib).age &&
       snes.ptr != C_NULL
        LibPETSc.SNESDestroy(PetscLib, snes)
    end
    snes.ptr = C_NULL
    return nothing
end

mutable struct SNES{PetscLib} <: AbstractSNES{PetscLib}
    ptr::CSNES
    opts::Options{PetscLib}
    age::Int
    # Garbage collection protections
    f!::Function
    updateJ!::Function
    function SNES{PetscLib}(ptr, opts, age) where {PetscLib}
        new{PetscLib}(
            ptr,
            opts,
            age,
            x -> error("function not defined"),
            x -> error("function not defined"),
        )
    end
end

"""
    SNES(
        petsclib::PetscLib
        comm::MPI.Comm;
        setfromoptions = true,
        options...
    )

Creates a PETSc nonlinear solve context on the communicator `comm`

If keyword argument `setfromoptions == true` then [`setfromoptions!`](@ref)
called on the created SNES context.

# External Links
$(_doc_external("SNES/SNESCreate"))
$(_doc_external("SNES/SNESSetFromOptions"))
"""
function SNES(
    petsclib::PetscLib,
    comm::MPI.Comm;
    setfromoptions = true,
    options...,
) where {PetscLib}
    opts = Options(petsclib; options...)

    snes = SNES{PetscLib}(C_NULL, opts, petsclib.age)

    with(snes.opts) do
        LibPETSc.SNESCreate(petsclib, comm, snes)
    end

    setfromoptions && setfromoptions!(snes)

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, snes)
    end

    return snes
end

function setfromoptions!(snes::AbstractSNES{PetscLib}) where {PetscLib}
    with(snes.opts) do
        LibPETSc.SNESSetFromOptions(PetscLib, snes)
    end
end

function gettype(snes::AbstractSNES{PetscLib}) where {PetscLib}
    r_type = Ref{CSNESType}()
    LibPETSc.SNESGetType(PetscLib, snes, r_type)
    return unsafe_string(r_type[])
end

function view(
    snes::AbstractSNES{PetscLib},
    viewer = LibPETSc.PETSC_VIEWER_STDOUT_(PetscLib, getcomm(snes)),
) where {PetscLib}
    LibPETSc.SNESView(PetscLib, snes, viewer)
    return nothing
end

"""
    setfunction!(snes::AbstractSNES, f!::Function, x::AbstractVec)
    setfunction!(f!::Function, snes::AbstractSNES, x::AbstractVec)

Define `f!` to be the function of the `snes`. A call to `f!(fx, snes, x)` should
set the elements of the PETSc vector `fx` based on the `x`.

# External Links
$(_doc_external("SNES/SNESSetFunction"))
"""
setfunction!(snes::AbstractSNES, rhs!, vec) = setfunction!(rhs!, snes, vec)

# Wrapper for calls to setfunction!
mutable struct Fn_SNESSetFunction{PetscLib, PetscInt} end
function (w::Fn_SNESSetFunction{PetscLib, PetscInt})(
    ::CSNES,
    r_x::CVec,
    r_fx::CVec,
    snes_ptr::Ptr{Cvoid},
)::PetscInt where {PetscLib, PetscInt}
    snes = unsafe_pointer_to_objref(snes_ptr)
    x = VecPtr(PetscLib, r_x, false)
    fx = VecPtr(PetscLib, r_fx, false)
    
    return snes.f!(fx, snes, x)
end

LibPETSc.@for_petsc function setfunction!(
    f!,
    snes::AbstractSNES{$PetscLib},
    vec::AbstractVec{$PetscLib},
) 
  
    ctx = pointer_from_objref(snes)
    PetscInt = $PetscLib.PetscInt
    fptr = @cfunction(
        Fn_SNESSetFunction{$PetscLib, $PetscInt}(),
        $PetscInt,
        (CSNES, CVec, CVec, Ptr{Cvoid})
    )
  
    with(snes.opts) do
        LibPETSc.SNESSetFunction($PetscLib, snes, vec, fptr, ctx)
    end
    snes.f! = f!
    return 0
end

"""
    setjacobian!(
        snes::AbstractSNES,
        updateJ!::Function,
        J::AbstractMat,
        P::AbstractMat = J
    )
    setjacobian!(
        updateJ!::Function,
        snes::AbstractSNES,
        J::AbstractMat,
        P::AbstractMat = J
    )

Define `updateJ!` to be the function that updates the Jacobian of the `snes`.

If `J == P` then a call to `updateJ!(J, snes, x)` should set the elements of the
PETSc Jacobian (approximation).

If `J ≠ P` then a call to `updateJ!(J, P, snes, x)` should set the elements of
the PETSc Jacobian (approximation) and preconditioning matrix `P`.

# External Links
$(_doc_external("SNES/SNESSetJacobian"))
"""
setjacobian!(snes::AbstractSNES, updateJ!, J, PJ = J) =
    setjacobian!(updateJ!, snes, J, PJ)

# Wrapper for calls to setjacobian!
mutable struct Fn_SNESSetJacobian{PetscLib, PetscInt} end
function (w::Fn_SNESSetJacobian{PetscLib, PetscInt})(
    ::CSNES,
    r_x::CVec,
    r_A::CMat,
    r_P::CMat,
    snes_ptr::Ptr{Cvoid},
)::PetscInt where {PetscLib, PetscInt}
    snes = unsafe_pointer_to_objref(snes_ptr)
    PetscScalar = PetscLib.PetscScalar
    petsclib = getlib(PetscLib)
    x = VecPtr(PetscLib, r_x, false)
    A = MatPtr{PetscLib, PetscScalar}(r_A, petsclib.age)
    P = MatPtr{PetscLib, PetscScalar}(r_P, petsclib.age)
    P == A
    return P == A ? snes.updateJ!(A, snes, x) : snes.updateJ!(A, P, snes, x)
end

LibPETSc.@for_petsc function setjacobian!(
    updateJ!,
    snes::AbstractSNES{$PetscLib},
    J::AbstractMat{$PetscLib},
    PJ::AbstractMat{$PetscLib} = J,
)
    ctx = pointer_from_objref(snes)
    fptr = @cfunction(
        Fn_SNESSetJacobian{$PetscLib, $PetscInt}(),
        $PetscInt,
        (CSNES, CVec, CMat, CMat, Ptr{Cvoid})
    )
    with(snes.opts) do
        LibPETSc.SNESSetJacobian($PetscLib, snes, J, PJ, fptr, ctx)
    end
    snes.updateJ! = updateJ!
    return nothing
end

function solve!(
    x::AbstractVec{PetscLib},
    snes::AbstractSNES{PetscLib},
    b::Union{Nothing, AbstractVec{PetscLib}} = nothing,
) where {PetscLib}
    with(snes.opts) do
        LibPETSc.SNESSolve(PetscLib, snes, isnothing(b) ? C_NULL : b, x)
    end
    return x
end

"""
    setDM!(snes::AbstractSNES, dm::AbstractDM)

Set `dm` for `snes`

# External Links
$(_doc_external("SNES/SNESSetDM"))
"""
function setDM!(
    snes::AbstractSNES{PetscLib},
    dm::AbstractDM{PetscLib},
) where {PetscLib}
    LibPETSc.SNESSetDM(PetscLib, snes, dm)
    return snes
end

"""
    getDMDA(snes::AbstractSNES)

Get `dmda` for `snes`

The returned `dmda` is owned by the `snes`

# External Links
$(_doc_external("SNES/SNESGetDM"))
"""
function getDMDA(
    snes::AbstractSNES{PetscLib},
) where {PetscLib}
    t_dmda = Ref{CDM}()
    LibPETSc.SNESGetDM(PetscLib, snes, t_dmda)
    dmda = DMDAPtr{PetscLib}(t_dmda[], getlib(PetscLib).age, false)
    return dmda
end


#=
@for_libpetsc begin

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

    function setjacobian!(snes::AbstractSNES{$PetscScalar}, update_jac!,
    A::AbstractMat{$PetscScalar}, P::AbstractMat{$PetscScalar}=A)
        ctx = pointer_from_objref(snes)
        jacptr = @cfunction(SNESJac{$PetscScalar}(), $PetscInt, (CSNES, CVec,
        CMat, CMat, Ptr{Cvoid}))

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
=#
