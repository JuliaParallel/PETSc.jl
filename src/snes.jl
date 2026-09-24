import .LibPETSc: AbstractSNES, CSNES, SNES

# The Julia side of a SNES, kept with the PETSc object (naming.md §18.3)
mutable struct SNESState <: ObjectState
    f!::Any
    updateJ!::Any
    convergence_test!::Any
    user_ctx::Any
    opts::Any
    alive::Bool
end
SNESState() = SNESState(
    x -> error("function not defined"),
    x -> error("function not defined"),
    nothing,
    nothing,
    nothing,
    true,
)
state_type(::Type{<:SNES}) = SNESState

# `snes.f!`, `snes.user_ctx`, `snes.opts` and the rest reach the state
@inline Base.getproperty(snes::SNES, name::Symbol) =
    (name === :ptr || name === :age || name === :own) ? getfield(snes, name) :
    forward_getproperty(snes, name)
@inline function Base.setproperty!(snes::SNES, name::Symbol, value)
    (name === :ptr || name === :age || name === :own) &&
        return Base.setfield!(snes, name, convert(fieldtype(typeof(snes), name), value))
    return forward_setproperty!(snes, name, value)
end
Base.propertynames(snes::SNES, private::Bool = false) = forward_propertynames(typeof(snes))


# positional constructors taking callbacks: the callbacks go to the state
function LibPETSc.SNES{PetscLib}(
    ptr::CSNES,
    age::Int,
    f!::Function,
    updateJ!::Function = x -> error("function not defined"),
    user_ctx = nothing,
    opts = nothing,
) where {PetscLib}
    ptr == C_NULL && throw(ArgumentError("callbacks need a PETSc object; got a null pointer"))
    snes = SNES{PetscLib}(ptr, age)
    snes.f! = f!
    snes.updateJ! = updateJ!
    snes.user_ctx = user_ctx
    snes.opts = opts
    return snes
end
LibPETSc.SNES(ptr::Ptr, lib::PetscLib, f!::Function, updateJ!::Function, user_ctx = nothing, age::Int = lib.age) where {PetscLib} =
    SNES{PetscLib}(ptr, age, f!, updateJ!, user_ctx)

# Custom display for REPL
function Base.show(io::IO, v::AbstractSNES{PetscLib}) where {PetscLib}
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
Otherwise, the user is responsible for calling `destroy!`.

# External Links
$(doc_external("SNES/SNESCreate"))
$(doc_external("SNES/SNESSetFromOptions"))
"""
function SNES(
    petsclib::PetscLib,
    comm::MPI.Comm;
    prefix="",
    options...,
) where {PetscLib <: LibPETSc.PetscLibType}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    snes = LibPETSc.SNESCreate(petsclib, comm)

    if !isempty(prefix)
        LibPETSc.SNESSetOptionsPrefix(petsclib, snes, prefix)
    end
    
    # Store options for deferred SNESSetFromOptions in solve!.
    # We do NOT call SetFromOptions here because the DM is typically
    # set after construction (via set_dm!), and PCs like MG need the
    # DM hierarchy to be available during SetFromOptions/SetUp.
    if !isempty(options)
        snes.opts = PetscOptions(petsclib; options...)
    end

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, snes)
    end

    return snes
end


"""
    type_name(snes::AbstractSNES)

The name PETSc knows this solver by, as a `Symbol` (`:newtonls`, `:fas`, …), or
`nothing` when no type has been set yet (docs/src/man/naming.md §3.1). v0.4
answered with a `String`; that is a break with no shim (§16).

# External Links
$(doc_external("SNES/SNESGetType"))
"""
type_name(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.SNESGetType(PetscLib, snes))

"""
    set_type!(snes::AbstractSNES, type::Symbol)

Set the nonlinear method, for example `:newtonls`, `:newtontr` or `:fas`.

# External Links
$(doc_external("SNES/SNESSetType"))
"""
function set_type!(snes::AbstractSNES{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.SNESSetType(getlib(PetscLib), snes, String(type))
    return snes
end

"""
    set_function!(f!, snes, vec)

Set the residual function `f!` for the nonlinear solver `snes`.

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).

The function `f!` will be called as `f!(fx, snes, x)` where:
- `fx`: Output vector to store the residual F(x)
- `snes`: The SNES context
- `x`: Input vector with current solution

The `vec` argument is a template vector used for the residual.

$(doc_callback())

# External Links
$(doc_external("SNES/SNESSetFunction"))
"""
function set_function! end

# Wrapper for calls to set_function!
mutable struct SNESSetFunctionFn{PetscLib} end
function (w::SNESSetFunctionFn{PetscLib})(
    actual_snes_ptr::CSNES,
    r_x::CVec,
    r_fx::CVec,
    snes_ptr::Ptr{Cvoid},
) where {PetscLib}
    snes = unsafe_pointer_to_objref(snes_ptr)::SNESState
    # Wrap the actual C SNES for the current MG level so that dm() inside
    # the callback returns the correct DM (matches the pattern in KSPComputeRHSFn).
    actual_snes = SNES{PetscLib}(actual_snes_ptr, getlib(PetscLib).age; own = false)
    x  = PetscVec{PetscLib}(r_x; own = false)
    fx = PetscVec{PetscLib}(r_fx; own = false)

    return run_callback("residual f!") do
        if Base.applicable(snes.f!, fx, actual_snes, x, snes.user_ctx)
            snes.f!(fx, actual_snes, x, snes.user_ctx)
        else
            snes.f!(fx, actual_snes, x)
        end
    end
end

LibPETSc.@for_petsc function set_function!(
    f!,
    snes::AbstractSNES{$PetscLib},
    vec::AbstractPetscVec{$PetscLib},
    ) 

    ctx = state_pointer(snes)
    PetscInt = $PetscLib.PetscInt
    fptr = @cfunction(
        SNESSetFunctionFn{$PetscLib}(),
        LibPETSc.PetscErrorCode,
        (CSNES, CVec, CVec, Ptr{Cvoid})
    )
  
    #with(snes.opts) do
    LibPETSc.SNESSetFunction($PetscLib, snes, vec, fptr, ctx)
    #end
    snes.f! = f!
    return snes
end

"""
    set_snes_jacobian!(
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

If a user context is set with [`set_user_ctx!`](@ref), then `updateJ!` may optionally accept that as an
additional last argument:

- `updateJ!(J, snes, x, user_ctx)` when `J == P`
- `updateJ!(J, P, snes, x, user_ctx)` when `J ≠ P`

$(doc_callback())

# External Links
$(doc_external("SNES/SNESSetJacobian"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_snes_jacobian! end

# Wrapper for calls to set_snes_jacobian!
mutable struct SNESSetJacobianFn{PetscLib} end
function (w::SNESSetJacobianFn{PetscLib})(
    actual_snes_ptr::CSNES,
    r_x::CVec,
    r_A::CMat,
    r_P::CMat,
    snes_ptr::Ptr{Cvoid},
) where {PetscLib}
    snes = unsafe_pointer_to_objref(snes_ptr)::SNESState
    actual_snes = SNES{PetscLib}(actual_snes_ptr, getlib(PetscLib).age; own = false)
    x = PetscVec{PetscLib}(r_x; own = false)
    A = PetscMat{PetscLib}(r_A; own = false)
    P = PetscMat{PetscLib}(r_P; own = false)

    same_mat = (P.ptr == A.ptr)

    return run_callback("Jacobian updateJ!") do
        if same_mat
            if Base.applicable(snes.updateJ!, A, actual_snes, x, snes.user_ctx)
                snes.updateJ!(A, actual_snes, x, snes.user_ctx)
            else
                snes.updateJ!(A, actual_snes, x)
            end
        else
            if Base.applicable(snes.updateJ!, A, P, actual_snes, x, snes.user_ctx)
                snes.updateJ!(A, P, actual_snes, x, snes.user_ctx)
            else
                snes.updateJ!(A, P, actual_snes, x)
            end
        end
    end
end

LibPETSc.@for_petsc function set_snes_jacobian!(
    updateJ!,
    snes::AbstractSNES{$PetscLib},
    J::AbstractPetscMat{$PetscLib},
    PJ::AbstractPetscMat{$PetscLib} = J,
)
    ctx = state_pointer(snes)
    fptr = @cfunction(
        SNESSetJacobianFn{$PetscLib}(),
        LibPETSc.PetscErrorCode,
        (CSNES, CVec, CMat, CMat, Ptr{Cvoid})
    )
    #with(snes.opts) do
        LibPETSc.SNESSetJacobian($PetscLib, snes, J, PJ, fptr, ctx)
    #end
    snes.updateJ! = updateJ!
    return snes
end

"""
    set_convergence_test!(test!::Function, snes::AbstractSNES)

Install a Julia closure as the `SNES` convergence test (`SNESSetConvergenceTest`).

`test!` is called as `test!(snes, it, xnorm, gnorm, fnorm)` at every iteration (`it`
starts at 0, before the first linear solve) and must return a `SNESConvergedReason`
(e.g. `LibPETSc.SNES_CONVERGED_ITERATING` to continue, a positive reason to report
convergence, or a negative reason to report divergence) — see `SNESConvergedReason`
in `LibPETSc`. `xnorm`/`gnorm`/`fnorm` are the current iterate/scaled-step/residual
2-norms, computed by PETSc exactly as for `SNESConvergedDefault`; a custom test
that instead needs its own residual (e.g. a normalised force residual computed
during `FormFunction`, not `‖F‖₂`) should ignore these and read whatever state it
cached during the residual evaluation via its own closure captures.

The closure is kept alive by `snes` until `snes` is destroyed or a new test is
installed. `snes.user_ctx` is left alone, so the residual and Jacobian callbacks
keep receiving it.

$(doc_callback("its return value is the `SNESConvergedReason` it decides"))

# External Links
$(doc_external("SNES/SNESSetConvergenceTest"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_convergence_test! end

# Wrapper for calls to set_convergence_test!. As for every SNES callback, the
# context pointer is the SNES's state, which holds the closure.
mutable struct SNESSetConvergenceTestFn{PetscLib} end
function (w::SNESSetConvergenceTestFn{PetscLib})(
    actual_snes_ptr::CSNES,
    it,
    xnorm,
    gnorm,
    fnorm,
    reason_ptr::Ptr{<:Integer},
    snes_ptr::Ptr{Cvoid},
) where {PetscLib}
    snes = unsafe_pointer_to_objref(snes_ptr)::SNESState
    actual_snes = SNES{PetscLib}(actual_snes_ptr, getlib(PetscLib).age; own = false)
    return run_callback("convergence test") do
        reason = snes.convergence_test!(actual_snes, Int(it), Float64(xnorm), Float64(gnorm), Float64(fnorm))
        unsafe_store!(reason_ptr, eltype(reason_ptr)(Int(reason)))
        return nothing
    end
end

LibPETSc.@for_petsc function set_convergence_test!(
    test!,
    snes::AbstractSNES{$PetscLib},
)
    snes.convergence_test! = test!
    ctx = state_pointer(snes)
    # PetscErrorCode is always Cint (see LibPETSc_const.jl), regardless of PetscInt's width;
    # SNESConvergedReason is a plain C enum, i.e. also Cint-sized.
    fptr = @cfunction(
        SNESSetConvergenceTestFn{$PetscLib}(),
        Cint,
        (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{Cint}, Ptr{Cvoid})
    )
    LibPETSc.SNESSetConvergenceTest($PetscLib, snes, fptr, ctx, C_NULL)
    return snes
end

"""
    user_ctx(snes::AbstractSNES)

Whatever was stored with [`set_user_ctx!`](@ref), or `nothing`.
"""
user_ctx(snes::AbstractSNES) = snes.user_ctx

"""
    set_user_ctx!(snes::AbstractSNES, ctx)

Attach `ctx` to `snes`, to be handed back as the last argument of the residual
and Jacobian callbacks that have a method accepting it. It is kept with the
PETSc object, so a borrowed handle such as `snes(ts)` sees it too. Returns `snes`.
"""
function set_user_ctx!(snes::AbstractSNES, ctx)
    snes.user_ctx = ctx
    return snes
end

function solve!(
    x::AbstractPetscVec{PetscLib},
    snes::AbstractSNES{PetscLib},
    b::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
) where {PetscLib}
    has_opts = !isnothing(snes.opts)
    if has_opts
        push!(snes.opts)
        # Call SetFromOptions here (rather than in the constructor) so that
        # the DM and callbacks are already attached. This is critical for
        # PCs like MG that need a DM hierarchy, and for FieldSplit sub-KSPs
        # that are created lazily during KSPSetUp/KSPSolve.
        LibPETSc.SNESSetFromOptions(PetscLib, snes)
    end
    try
        capture_callback_errors() do
            LibPETSc.SNESSolve(PetscLib, snes, isnothing(b) ? C_NULL : b, x)
        end
    finally
        has_opts && pop!(snes.opts)
    end
    return x
end


"""
    destroy!(snes::AbstractSNES)

Destroy a SNES (nonlinear solver) object and release associated resources.

This function is typically called automatically via finalizers when the object
is garbage collected, but can be called explicitly to free resources immediately.
Does nothing on a borrowed handle, such as the one [`snes(ts)`](@ref) hands back.

# External Links
$(doc_external("SNES/SNESDestroy"))
"""
function destroy!(snes::AbstractSNES{PetscLib}) where {PetscLib}
    owns(snes) || return nothing
    if !isnothing(snes.opts)
        destroy!(snes.opts)
        snes.opts = nothing
    end
    if isdestroyable(snes, PetscLib)
        LibPETSc.SNESDestroy(PetscLib, snes)
    end
    snes.ptr = C_NULL
    return nothing
end


"""
    d = dm(snes::AbstractSNES)

The DM attached to `snes`, [`narrow`](@ref)ed to its flavour.

$(doc_borrowed())

# External Links
$(doc_external("SNES/SNESGetDM"))
"""
function dm(
    snes::AbstractSNES{PetscLib},
) where {PetscLib}
    return narrow(LibPETSc.SNESGetDM(getlib(PetscLib), snes))
end


"""
    set_dm!(snes::AbstractSNES, dm::AbstractDM)

Set `dm` for `snes`

# External Links
$(doc_external("SNES/SNESSetDM"))
"""
function set_dm!(
    snes::AbstractSNES{PetscLib},
    dm::AbstractPetscDM{PetscLib},
) where {PetscLib}
    LibPETSc.SNESSetDM(getlib(PetscLib), snes, dm)
    return snes
end
