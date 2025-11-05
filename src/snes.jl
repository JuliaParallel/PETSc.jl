import .LibPETSc: AbstractPetscSNES, CSNES, PetscSNES

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscSNES{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc SNES (null pointer)")
        return
    else
        print(io, "PETSc SNES object")
    end
    return nothing
end

"""
    SNES(comm::MPI.Comm; prefix="", options...)

Creates a PETSc nonlinear solve context on the communicator `comm` with optional `prefix` and `options`.

The communicator is obtained from `A` and if it has size `1` then the garbage
collector is set, otherwise the user is responsible for calling
[`destroy`](@ref).

# External Links
$(_doc_external("SNES/SNESCreate"))
$(_doc_external("SNES/SNESSetFromOptions"))
"""
function SNES(
    petsclib::PetscLib,
    comm::MPI.Comm;
    prefix="",
    options...,
) where {PetscLib}
    @assert initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    snes = LibPETSc.SNESCreate(petsclib, comm)

    if !isempty(prefix)
        LibPETSc.SNESSetOptionsPrefix(petsclib, snes, prefix)
    end
    
    # Push options to PETSc options database
    if !isempty(options)
        opts = PETSc.Options(petsclib; options...);
        push!(opts)
        LibPETSc.SNESSetFromOptions(petsclib, snes)
        pop!(opts)
    end

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, snes)
    end

    return snes
end


function gettype(snes::AbstractPetscSNES{PetscLib}) where {PetscLib}
    return LibPETSc.SNESGetType(PetscLib, snes)
end

"""
    setfunction!(snes::AbstractPetscSNES, f!::Function, x::AbstractVec)
    setfunction!(f!::Function, snes::AbstractPetscSNES, x::AbstractVec)

Define `f!` to be the function of the `snes`. A call to `f!(fx, snes, x)` should
set the elements of the PETSc vector `fx` based on the `x`.

# External Links
$(_doc_external("SNES/SNESSetFunction"))
"""
setfunction!(snes::AbstractPetscSNES, rhs!, vec) = setfunction!(rhs!, snes, vec)

# Wrapper for calls to setfunction!
mutable struct Fn_SNESSetFunction{PetscLib} end
function (w::Fn_SNESSetFunction{PetscLib})(
    ::CSNES,
    r_x::CVec,
    r_fx::CVec,
    snes_ptr::Ptr{Cvoid},
) where {PetscLib}
    snes = unsafe_pointer_to_objref(snes_ptr)
    x  = PetscVec{PetscLib}(r_x)
    fx = PetscVec{PetscLib}(r_fx)

    return snes.f!(fx, snes, x)
end

LibPETSc.@for_petsc function setfunction!(
    f!,
    snes::AbstractPetscSNES{$PetscLib},
    vec::AbstractPetscVec{$PetscLib},
    ) 

    ctx = pointer_from_objref(snes)
    PetscInt = $PetscLib.PetscInt
    fptr = @cfunction(
        Fn_SNESSetFunction{$PetscLib}(),
        $PetscInt,
        (CSNES, CVec, CVec, Ptr{Cvoid})
    )
  
    #with(snes.opts) do
    LibPETSc.SNESSetFunction($PetscLib, snes, vec, fptr, ctx)
    #end
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
setjacobian!(snes::AbstractPetscSNES, updateJ!, J, PJ = J) =
    setjacobian!(updateJ!, snes, J, PJ)

# Wrapper for calls to setjacobian!
mutable struct Fn_SNESSetJacobian{PetscLib} end
function (w::Fn_SNESSetJacobian{PetscLib})(
    ::CSNES,
    r_x::CVec,
    r_A::CMat,
    r_P::CMat,
    snes_ptr::Ptr{Cvoid},
) where {PetscLib}
    snes = unsafe_pointer_to_objref(snes_ptr)
    x = PetscVec{PetscLib}(r_x)
    A = PetscMat{PetscLib}(r_A)
    P = PetscMat{PetscLib}(r_P)
    
    P == A

    return P == A ? snes.updateJ!(A, snes, x) : snes.updateJ!(A, P, snes, x)
end

LibPETSc.@for_petsc function setjacobian!(
    updateJ!,
    snes::AbstractPetscSNES{$PetscLib},
    J::AbstractPetscMat{$PetscLib},
    PJ::AbstractPetscMat{$PetscLib} = J,
)
    ctx = pointer_from_objref(snes)
    fptr = @cfunction(
        Fn_SNESSetJacobian{$PetscLib}(),
        $PetscInt,
        (CSNES, CVec, CMat, CMat, Ptr{Cvoid})
    )
    #with(snes.opts) do
        LibPETSc.SNESSetJacobian($PetscLib, snes, J, PJ, fptr, ctx)
    #end
    snes.updateJ! = updateJ!
    return nothing
end

function solve!(
    x::AbstractPetscVec{PetscLib},
    snes::AbstractPetscSNES{PetscLib},
    b::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
) where {PetscLib}
    #with(snes.opts) do
    LibPETSc.SNESSolve(PetscLib, snes, isnothing(b) ? C_NULL : b, x)
    #end
    return x
end


function destroy(snes::AbstractPetscSNES{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) && snes.ptr != C_NULL
        LibPETSc.SNESDestroy(PetscLib, snes)
    end
    snes.ptr = C_NULL
    return nothing
end