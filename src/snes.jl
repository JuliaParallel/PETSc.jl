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
    SNES(petsclib, comm::MPI.Comm; prefix="", options...)

Create a PETSc nonlinear solver (SNES) context on the communicator `comm`.

# Arguments
- `petsclib`: The PETSc library instance
- `comm::MPI.Comm`: MPI communicator
- `prefix::String`: Optional prefix for command-line options
- `options...`: Additional PETSc options as keyword arguments

If `comm` has size 1, the garbage collector will handle cleanup automatically.
Otherwise, the user is responsible for calling `destroy`.

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
    setfunction!(snes, f!, vec)
    setfunction!(f!, snes, vec)

Set the residual function `f!` for the nonlinear solver `snes`.

The function `f!` will be called as `f!(fx, snes, x)` where:
- `fx`: Output vector to store the residual F(x)
- `snes`: The SNES context
- `x`: Input vector with current solution

The `vec` argument is a template vector used for the residual.

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

    if Base.applicable(snes.f!, fx, snes, x, snes.user_ctx)
        return snes.f!(fx, snes, x, snes.user_ctx)
    else
        return snes.f!(fx, snes, x)
    end
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

If you set `snes.user_ctx`, then `updateJ!` may optionally accept that as an
additional last argument:

- `updateJ!(J, snes, x, user_ctx)` when `J == P`
- `updateJ!(J, P, snes, x, user_ctx)` when `J ≠ P`

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

    same_mat = (P.ptr == A.ptr)

    if same_mat
        if Base.applicable(snes.updateJ!, A, snes, x, snes.user_ctx)
            return snes.updateJ!(A, snes, x, snes.user_ctx)
        else
            return snes.updateJ!(A, snes, x)
        end
    else
        if Base.applicable(snes.updateJ!, A, P, snes, x, snes.user_ctx)
            return snes.updateJ!(A, P, snes, x, snes.user_ctx)
        else
            return snes.updateJ!(A, P, snes, x)
        end
    end
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


"""
    destroy(snes::AbstractPetscSNES)

Destroy a SNES (nonlinear solver) object and release associated resources.

This function is typically called automatically via finalizers when the object
is garbage collected, but can be called explicitly to free resources immediately.

# External Links
$(_doc_external("SNES/SNESDestroy"))
"""
function destroy(snes::AbstractPetscSNES{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) && snes.ptr != C_NULL
        LibPETSc.SNESDestroy(PetscLib, snes)
    end
    snes.ptr = C_NULL
    return nothing
end


"""
    dm = getDM(snes::AbstractPetscSNES)

Get `dmda` for `snes`

The returned `dmda` is owned by the `snes`

# External Links
$(_doc_external("SNES/SNESGetDM"))
"""
function getDM(
    snes::AbstractPetscSNES{PetscLib},
) where {PetscLib}
    dmda = LibPETSc.SNESGetDM(getlib(PetscLib), snes)
    return dmda
end


"""
    setDM!(snes::AbstractPetscSNES, dm::AbstractDM)

Set `dm` for `snes`

# External Links
$(_doc_external("SNES/SNESSetDM"))
"""
function setDM!(
    snes::AbstractPetscSNES{PetscLib},
    dm::AbstractPetscDM{PetscLib},
) where {PetscLib}
    LibPETSc.SNESSetDM(getlib(PetscLib), snes, dm)
    return nothing
end
