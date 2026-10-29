# The Julia side of a SNES, kept with the PETSc object (naming.md §18.3)
mutable struct SNESState <: ObjectState
    f!::Any
    updateJ!::Any
    convergence_test!::Any
    user_ctx::Any
    opts::Any
    options_at::Symbol
    options_applied::Bool
    alive::Bool
end
SNESState() = SNESState(
    x -> error("function not defined"),
    x -> error("function not defined"),
    nothing,
    nothing,
    nothing,
    :first_solve,
    false,
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
    SNES(petsclib, comm::MPI.Comm; prefix = "", options_at = :first_solve, options...)

Create a PETSc nonlinear solver (SNES) context on the communicator `comm`.

`options_at` decides when the options are applied, through `SNESSetFromOptions`:

- `:first_solve` (the default): once, at the first [`setup!`](@ref) or `solve!`,
  once the DM and the callbacks are attached. A setting made in code after that
  overrides them.
- `:each_solve`: at the start of every [`solve!`](@ref). They override a setting
  made in code for the same option, including a preconditioner configured
  between solves.
- `:construction`: here, once. A DM or callback attached afterwards is not
  visible to them.

The options stay pushed during `solve!`, and during `setup!` too unless
`options_at = :each_solve`, so the sub-solvers PETSc creates during setup read
them (see [`with_options`](@ref)).

# Arguments
- `petsclib`: The PETSc library instance
- `comm::MPI.Comm`: MPI communicator
- `prefix::AbstractString`: Optional prefix for command-line options
- `options_at::Symbol`: When the options are applied, see above
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
    prefix::AbstractString = "",
    options_at::Symbol = :first_solve,
    options...,
) where {PetscLib <: LibPETSc.PetscLibType}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    snes = LibPETSc.SNESCreate(petsclib, comm)

    if !isempty(prefix)
        LibPETSc.SNESSetOptionsPrefix(petsclib, snes, String(prefix))
    end
    
    store_options!(snes, petsclib, options_at, options)

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
    x  = PetscVec{PetscLib}(r_x, getlib(PetscLib).age; own = false)
    fx = PetscVec{PetscLib}(r_fx, getlib(PetscLib).age; own = false)

    return run_callback("residual f!") do
        if Base.applicable(snes.f!, fx, actual_snes, x, snes.user_ctx)
            warn_user_ctx(snes.user_ctx)
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

Deprecated: if a user context is set with [`set_user_ctx!`](@ref), `updateJ!` may
accept it as an additional last argument (`updateJ!(J, snes, x, user_ctx)` when
`J == P`, `updateJ!(J, P, snes, x, user_ctx)` when `J ≠ P`). Capture the context
in a closure instead; the trailing form goes in v0.6.

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
    x = PetscVec{PetscLib}(r_x, getlib(PetscLib).age; own = false)
    A = PetscMat{PetscLib}(r_A, getlib(PetscLib).age; own = false)
    P = PetscMat{PetscLib}(r_P, getlib(PetscLib).age; own = false)

    same_mat = (P.ptr == A.ptr)

    return run_callback("Jacobian updateJ!") do
        if same_mat
            if Base.applicable(snes.updateJ!, A, actual_snes, x, snes.user_ctx)
                warn_user_ctx(snes.user_ctx)
                snes.updateJ!(A, actual_snes, x, snes.user_ctx)
            else
                snes.updateJ!(A, actual_snes, x)
            end
        else
            if Base.applicable(snes.updateJ!, A, P, actual_snes, x, snes.user_ctx)
                warn_user_ctx(snes.user_ctx)
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
        PetscReal = getlib(PetscLib).PetscReal
        reason = snes.convergence_test!(actual_snes, Int(it), PetscReal(xnorm), PetscReal(gnorm), PetscReal(fnorm))
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

Deprecated: a closure carries the context without this extra argument, and the
callbacks then have one calling form. `user_ctx` and the trailing-context
callback forms go in v0.6.
"""
function user_ctx(snes::AbstractSNES)
    @warn "user_ctx is deprecated and goes in v0.6: capture the context in a closure" maxlog = 1
    return snes.user_ctx
end

"""
    set_user_ctx!(snes::AbstractSNES, ctx)

Attach `ctx` to `snes`, to be handed back as the last argument of the residual
and Jacobian callbacks that have a method accepting it. It is kept with the
PETSc object, so a borrowed handle such as `snes(ts)` sees it too. Returns `snes`.

Deprecated: a closure carries the context without this extra argument, and the
callbacks then have one calling form. `user_ctx` and the trailing-context
callback forms go in v0.6.
"""
function set_user_ctx!(snes::AbstractSNES, ctx)
    @warn "set_user_ctx! is deprecated and goes in v0.6: capture the context in a closure" maxlog = 1
    snes.user_ctx = ctx
    return snes
end

"""
    solve!(x::AbstractPetscVec, snes::SNES, b = nothing; check = true)

Solve the nonlinear system of `snes` starting from `x`, which holds the solution on
return, and return `x`. `b` is an optional right-hand side: the system solved is
F(x) = b. The options passed to the [`SNES`](@ref) constructor are applied at the
first `setup!` or `solve!` by default; its `options_at` keyword says when.

PETSc does not report a failed solve: `x` holds the last iterate and
[`converged_reason`](@ref) says why it stopped. `check` decides what happens then:

- `check = true` (the default) throws a [`ConvergenceError`](@ref).
- `check = false` returns normally, and the caller reads `converged_reason`.

# External Links
$(doc_external("SNES/SNESSolve"))
$(doc_external("SNES/SNESGetConvergedReason"))
"""
function solve!(
    x::AbstractPetscVec{PetscLib},
    snes::AbstractSNES{PetscLib},
    b::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing;
    check::Bool = true,
) where {PetscLib}
    # Under the default `options_at = :first_solve`, SNESSetFromOptions runs at the
    # first solve rather than in the constructor, so that it sees the DM hierarchy
    # that multigrid needs and the callbacks.
    with_solver_options(snes; solving = true) do
        with_convergence_check(snes, check) do
            capture_callback_errors() do
                LibPETSc.SNESSolve(PetscLib, snes, isnothing(b) ? C_NULL : b, x)
            end
        end
    end
    return x
end

"""
    setup!(snes::AbstractSNES)

Set up `snes`, as [`solve!`](@ref) does before its first solve, and return
`snes`. Under `options_at = :first_solve` or `:construction` the constructor
options are pushed meanwhile, and under `:first_solve` applied first (see
[`SNES`](@ref)).

# External Links
$(doc_external("SNES/SNESSetUp"))
"""
function setup!(snes::AbstractSNES{PetscLib}) where {PetscLib}
    with_solver_options(snes; solving = false) do
        capture_callback_errors(() -> LibPETSc.SNESSetUp(PetscLib, snes))
    end
    return snes
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
    owns(snes) && destroy_options!(snes)
    return destroy_handle!(LibPETSc.SNESDestroy, snes, PetscLib)
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

"""
    ksp(snes::AbstractSNES)

The linear solver `snes` uses for its Newton steps.

$(doc_borrowed())

# External Links
$(doc_external("SNES/SNESGetKSP"))
"""
ksp(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetKSP(getlib(PetscLib), snes)

"""
    solution(snes::AbstractSNES)

The vector `snes` solves for, as a `PetscVec`.

$(doc_borrowed())

# External Links
$(doc_external("SNES/SNESGetSolution"))
"""
solution(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetSolution(getlib(PetscLib), snes)

"""
    set_solution!(snes::AbstractSNES, x::AbstractPetscVec)

Set the vector `snes` solves for, which [`solution`](@ref) then returns.
[`solve!`](@ref)`(x, snes)` sets it too. Returns `snes`.

# External Links
$(doc_external("SNES/SNESSetSolution"))
"""
function set_solution!(
    snes::AbstractSNES{PetscLib},
    x::AbstractPetscVec{PetscLib},
) where {PetscLib}
    LibPETSc.SNESSetSolution(getlib(PetscLib), snes, x)
    return snes
end

"""
    iteration_number(snes::AbstractSNES)

The number of nonlinear iterations the last [`solve!`](@ref) took.

# External Links
$(doc_external("SNES/SNESGetIterationNumber"))
"""
iteration_number(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetIterationNumber(getlib(PetscLib), snes)

"""
    ksp_iterations(snes::AbstractSNES)

The number of linear iterations the last [`solve!`](@ref) took, summed over its
nonlinear iterations.

# External Links
$(doc_external("SNES/SNESGetLinearSolveIterations"))
"""
ksp_iterations(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetLinearSolveIterations(getlib(PetscLib), snes)

"""
    function_norm(snes::AbstractSNES)

The norm of the residual at the current iterate.

# External Links
$(doc_external("SNES/SNESGetFunctionNorm"))
"""
function_norm(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetFunctionNorm(getlib(PetscLib), snes)

"""
    converged_reason(snes::AbstractSNES)

Why the last [`solve!`](@ref) stopped, as a `LibPETSc.SNESConvergedReason`:
positive when it converged, negative when it diverged.

# External Links
$(doc_external("SNES/SNESGetConvergedReason"))
"""
converged_reason(snes::AbstractSNES{PetscLib}) where {PetscLib} =
    LibPETSc.SNESGetConvergedReason(getlib(PetscLib), snes)

"""
    tolerances(snes::AbstractSNES)

The convergence tolerances of `snes`, as `(; atol, rtol, stol, max_it, max_funcs)`:
the absolute and relative residual tolerances, the step tolerance (relative to the
iterate), the iteration limit and the limit on residual evaluations.

# External Links
$(doc_external("SNES/SNESGetTolerances"))
"""
function tolerances(snes::AbstractSNES{PetscLib}) where {PetscLib}
    atol, rtol, stol, max_it, max_funcs = LibPETSc.SNESGetTolerances(getlib(PetscLib), snes)
    return (; atol, rtol, stol, max_it = Int(max_it), max_funcs = Int(max_funcs))
end

"""
    set_tolerances!(snes::AbstractSNES; atol, rtol, stol, max_it, max_funcs)

Set the convergence tolerances of `snes`, named as [`tolerances`](@ref) returns them.
A keyword left at `nothing` keeps the value `snes` has. Returns `snes`.

An option for the same tolerance (`snes_rtol = 1e-8` to the constructor, or
`-snes_rtol` on the command line) is applied when the constructor's `options_at`
says, and overrides this setting if that comes later.

# External Links
$(doc_external("SNES/SNESSetTolerances"))
"""
function set_tolerances!(
    snes::AbstractSNES{PetscLib};
    atol = nothing,
    rtol = nothing,
    stol = nothing,
    max_it = nothing,
    max_funcs = nothing,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    cur = tolerances(snes)
    LibPETSc.SNESSetTolerances(
        petsclib,
        snes,
        PetscLib.PetscReal(something(atol, cur.atol)),
        PetscLib.PetscReal(something(rtol, cur.rtol)),
        PetscLib.PetscReal(something(stol, cur.stol)),
        PetscLib.PetscInt(something(max_it, cur.max_it)),
        PetscLib.PetscInt(something(max_funcs, cur.max_funcs)),
    )
    return snes
end

"""
    set_function_domain_error!(snes::AbstractSNES)

Tell `snes`, from inside the residual callback, that the iterate it was given
lies outside the function's domain (a negative pressure, say). A few solvers
then cut the step; the others stop, and [`converged_reason`](@ref) reports
`SNES_DIVERGED_FUNCTION_DOMAIN`. Returns `snes`.

Either way `solve!` returns normally. Throwing from the callback instead stops
the solve, and `solve!` rethrows the exception.

# External Links
$(doc_external("SNES/SNESSetFunctionDomainError"))
"""
function set_function_domain_error!(snes::AbstractSNES{PetscLib}) where {PetscLib}
    LibPETSc.SNESSetFunctionDomainError(getlib(PetscLib), snes)
    return snes
end


"""
    set_jacobian_nullspace!(snes, nullsp)

Retrieve the assembled Jacobian matrix from `snes` (after `SNESSetUp`) and
attach `nullsp` to it.  Must be called after `SNESSetUp` and before `SNESSolve`.
"""
function set_jacobian_nullspace! end

LibPETSc.@for_petsc function set_jacobian_nullspace!(
    snes::LibPETSc.SNES{$PetscLib},
    nullsp::LibPETSc.MatNullSpace,
)
    petsclib = getlib($PetscLib)
    J = LibPETSc.SNESGetJacobianMat(petsclib, snes)
    LibPETSc.MatSetNullSpace(petsclib, J, nullsp)
    return snes
end
