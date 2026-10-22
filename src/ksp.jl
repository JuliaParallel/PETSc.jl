import .LibPETSc: AbstractKSP, CKSP, KSP, AbstractPetscDM   # KSP methods below are constructors of LibPETSc.KSP

# The Julia side of a KSP, kept with the PETSc object
mutable struct KSPState <: ObjectState
    computerhs!::Any
    computeops!::Any
    pre_solve!::Any
    opts::Any
    options_at::Symbol
    options_applied::Bool
    alive::Bool
end
KSPState() = KSPState(
    x -> error("computerhs! not defined"),
    x -> error("computeops! not defined"),
    nothing,
    nothing,
    :construction,
    false,
    true,
)
state_type(::Type{<:KSP}) = KSPState

# `ksp.opts` and the callbacks reach the state
@inline Base.getproperty(ksp::KSP, name::Symbol) =
    (name === :ptr || name === :age || name === :own) ? getfield(ksp, name) :
    forward_getproperty(ksp, name)
@inline function Base.setproperty!(ksp::KSP, name::Symbol, value)
    (name === :ptr || name === :age || name === :own) &&
        return Base.setfield!(ksp, name, convert(fieldtype(typeof(ksp), name), value))
    return forward_setproperty!(ksp, name, value)
end
Base.propertynames(ksp::KSP, private::Bool = false) = forward_propertynames(typeof(ksp))


# positional constructor taking callbacks: the callbacks go to the state
function LibPETSc.KSP{PetscLib}(
    ptr::CKSP,
    age::Int,
    computerhs!::Function,
    computeops!::Function = x -> error("computeops! not defined"),
    opts = nothing,
) where {PetscLib}
    ptr == C_NULL && throw(ArgumentError("callbacks need a PETSc object; got a null pointer"))
    ksp = KSP{PetscLib}(ptr, age)
    ksp.computerhs! = computerhs!
    ksp.computeops! = computeops!
    ksp.opts = opts
    return ksp
end

# Custom display for REPL
function Base.show(io::IO, v::AbstractKSP{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc KSP (null pointer)")
        return
    else
        print(io, "PETSc KSP object")
    end
    return nothing
end


"""
    KSP(A::AbstractPetscMat, P::AbstractPetscMat = A; prefix = "", options_at = :construction, options...)

Create a `KSP` using the matrix `A` and preconditioner construction matrix `P`
with optional `prefix` and `options`.

`options_at` decides when the options are applied, through `KSPSetFromOptions`:

- `:construction` (the default): here, once. A setting made in code afterwards,
  such as [`set_tolerances!`](@ref), overrides them.
- `:each_solve`: at the start of every [`solve!`](@ref), so they override a
  setting made in code for the same option.
- `:first_solve`: once, at the first [`setup!`](@ref) or `solve!`, after the DM
  and callbacks are attached.

The options stay pushed during `solve!`, and during `setup!` too unless
`options_at = :each_solve`, so the sub-solvers PETSc creates during setup read
them (see [`with_options`](@ref)).

The communicator is obtained from `A` and if it has size `1` then the garbage
collector is set, otherwise the user is responsible for calling
[`destroy!`](@ref).

# External Links
$(doc_external("KSP/KSPCreate"))
$(doc_external("KSP/KSPSetOperators"))
$(doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(
    A::AbstractPetscMat{PetscLib},
    P::AbstractPetscMat{PetscLib} = A;
    prefix::AbstractString = "",
    options_at::Symbol = :construction,
    options...
) where {PetscLib}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    c = comm(A)
    ksp = LibPETSc.KSPCreate(petsclib, c)
    
    LibPETSc.KSPSetOperators(petsclib, ksp, A, P)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, String(prefix))
    end
    
    store_options!(ksp, petsclib, options_at, options)

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(c) == 1
        finalizer(destroy!, ksp)
    end

    return ksp
end

"""
    KSP(dm::AbstractPetscDM; prefix = "", options_at = :construction, options...)

Create a `KSP` associated with the `dm` with optional `prefix` and `options`.

`options_at` decides when the options are applied, through `KSPSetFromOptions`:

- `:construction` (the default): here, once. A setting made in code afterwards,
  such as [`set_tolerances!`](@ref), overrides them.
- `:each_solve`: at the start of every [`solve!`](@ref), so they override a
  setting made in code for the same option.
- `:first_solve`: once, at the first [`setup!`](@ref) or `solve!`, after the DM
  and callbacks are attached.

The options stay pushed during `solve!`, and during `setup!` too unless
`options_at = :each_solve`, so the sub-solvers PETSc creates during setup read
them (see [`with_options`](@ref)).

The communicator is obtained from `dm`. The KSP can be used with geometric
multigrid when the DM provides grid hierarchy information.

# Arguments
- `dm::AbstractPetscDM`: The DM object to associate with the KSP
- `prefix::AbstractString`: Optional prefix for command-line options
- `options_at::Symbol`: When the options are applied, see above
- `options...`: Additional PETSc options as keyword arguments

# External Links
$(doc_external("KSP/KSPCreate"))
$(doc_external("KSP/KSPSetDM"))
$(doc_external("KSP/KSPSetFromOptions"))
"""
function KSP(dm::AbstractPetscDM{PetscLib};
    prefix::AbstractString = "",
    options_at::Symbol = :construction,
    options...
) where {PetscLib}
    check_initialized(getlib(PetscLib))
    petsclib = getlib(PetscLib)
    c = comm(dm)
    ksp = LibPETSc.KSPCreate(petsclib, c)
    
    if !isempty(prefix)
        LibPETSc.KSPSetOptionsPrefix(petsclib, ksp, String(prefix))
    end
    
    LibPETSc.KSPSetDM(petsclib, ksp, dm)

    store_options!(ksp, petsclib, options_at, options)

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(c) == 1
        finalizer(destroy!, ksp)
    end

    return ksp
end



"""
    KSP(petsclib, comm::MPI.Comm, A::SparseMatrixCSC; options...)

Create a [`KSP`](@ref) for the sparse matrix `A`, on `comm`, using the library
`petsclib`. `A` is copied into a new PETSc matrix that the solver owns, so changing
`A` afterwards does not change the solver. The keyword arguments are PETSc options,
as for `KSP(A::AbstractPetscMat; options...)`.
"""
KSP(petsclib, comm, S::SparseMatrixCSC; kwargs...) 

function KSP(petsclib, comm, S::SparseMatrixCSC; kwargs...) 
    M = LibPETSc.PetscMat(petsclib, comm, S)
    ksp = KSP(M; kwargs...)
    destroy!(M)   # the KSP holds its own reference
    return ksp
end


"""
    solve!(x::AbstractPetscVec, ksp::AbstractKSP, b::AbstractPetscVec; check = nothing)
    solve!(ksp::KSP; check = nothing)

Solve the linear system of `ksp` for the right-hand side `b` into `x`, and return
`x`. The second form uses the vectors set on `ksp` (for example by
[`set_compute_rhs!`](@ref)) and returns `ksp`. 
The initial guess is zero unless `-ksp_initial_guess_nonzero` is set.

PETSc does not report a failed solve: `x` holds the last iterate and
[`converged_reason`](@ref) says why it stopped. `check` decides what happens then:

- `check = true` throws a [`ConvergenceError`](@ref).
- `check = false` returns normally, and the caller reads `converged_reason`.
- Not passed, it returns normally and warns once per call site that from v0.6 a
  failed solve throws.

# External Links
$(doc_external("KSP/KSPSolve"))
$(doc_external("KSP/KSPGetConvergedReason"))
"""
function solve!(
    x::AbstractPetscVec{PetscLib},
    ksp::AbstractKSP{PetscLib},
    b::AbstractPetscVec{PetscLib};
    check::Union{Nothing, Bool} = nothing,
) where {PetscLib}
    with_solver_options(ksp; solving = true) do
        with_convergence_check(ksp, check) do
            capture_callback_errors() do
                LibPETSc.KSPSolve(PetscLib, ksp, b, x)
            end
        end
    end
    return x
end

function solve!(
    ksp::AbstractKSP{PetscLib};
    check::Union{Nothing, Bool} = nothing,
) where {PetscLib}
    with_solver_options(ksp; solving = true) do
        with_convergence_check(ksp, check) do
            capture_callback_errors() do
                LibPETSc.KSPSolve(PetscLib, ksp, C_NULL, C_NULL)
            end
        end
    end
    return ksp
end

LinearAlgebra.ldiv!(x::AbstractPetscVec{PetscLib}, ksp::AbstractKSP{PetscLib}, b::AbstractPetscVec{PetscLib}) where {PetscLib} =
    solve!(x, ksp, b)

"""
    setup!(ksp::AbstractKSP)

Set up `ksp` and its preconditioner, as [`solve!`](@ref) does before its first
solve, and return `ksp`. Call it to build the preconditioner at a controlled
point, for example to time it, or to reach the sub-solvers of a fieldsplit
before the first solve. Unless `options_at = :each_solve`, the constructor options
are pushed meanwhile, and under `:first_solve` applied first (see [`KSP`](@ref)).

# External Links
$(doc_external("KSP/KSPSetUp"))
"""
function setup!(ksp::AbstractKSP{PetscLib}) where {PetscLib}
    with_solver_options(ksp; solving = false) do
        capture_callback_errors(() -> LibPETSc.KSPSetUp(PetscLib, ksp))
    end
    return ksp
end

"""
    set_pre_solve!(pre_solve!, ksp::AbstractKSP)

Call `pre_solve!(ksp, b, x)` at the start of every solve of `ksp`, before the
first iteration. `b` is the right-hand side and `x` the initial guess, both
borrowed. Returns `ksp`. A second call replaces the first callback.

It runs before the preconditioner is set up for that solve, so the parts of a
composite preconditioner aren't configured yet. To change those, use
[`set_post_setup!`](@ref) on the preconditioner.

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax works.

$(doc_callback())

# External Links
$(doc_external("KSP/KSPSetPreSolve"))
"""
function set_pre_solve! end

struct KSPSetPreSolveFn{PetscLib} end
function (::KSPSetPreSolveFn{PetscLib})(
    ksp_ptr::CKSP,
    b_ptr::CVec,
    x_ptr::CVec,
    state_ptr::Ptr{Cvoid},
) where {PetscLib}
    return run_callback("pre solve!") do
        age = getlib(PetscLib).age
        state = unsafe_pointer_to_objref(state_ptr)::KSPState
        state.pre_solve!(
            KSP{PetscLib}(ksp_ptr, age; own = false),
            PetscVec{PetscLib}(b_ptr, age; own = false),
            PetscVec{PetscLib}(x_ptr, age; own = false),
        )
    end
end

LibPETSc.@for_petsc function set_pre_solve!(pre_solve!, ksp::AbstractKSP{$PetscLib})
    object_state!(ksp).pre_solve! = pre_solve!
    fptr = @cfunction(
        KSPSetPreSolveFn{$PetscLib}(),
        LibPETSc.PetscErrorCode,
        (CKSP, CVec, CVec, Ptr{Cvoid})
    )
    LibPETSc.KSPSetPreSolve($PetscLib, ksp, fptr, state_pointer(ksp))
    return ksp
end

function Base.:\(ksp::AbstractKSP{PetscLib}, b::AbstractPetscVec{PetscLib}) where {PetscLib}
    x = similar(b)
    ldiv!(x, ksp, b)
    return x
end

# Solves into a copy of `b` in the library's scalar type, so any vector works, 
# and returns `x` with the axes of `b`
function Base.:\(ksp::KSP{PetscLib}, b::AbstractVector) where {PetscLib}
    c = comm(ksp)
    MPI.Comm_size(c) == 1 || throw(
        ArgumentError(
            "solving into a Julia Vector requires a sequential KSP, " *
            "but its communicator spans $(MPI.Comm_size(c)) ranks",
        ),
    )
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar

    # PETSc works on this copy of `b` in place, so it must outlive the solve
    b_copy = copyto!(Vector{PetscScalar}(undef, length(b)), b)
    x = GC.@preserve b_copy begin
        petsc_b = seq_vec_with_array(getlib(PetscLib), c, PetscInt(1), PetscInt(length(b)), b_copy)
        petsc_x = ksp \ petsc_b
        x = copyto!(similar(b, PetscScalar), petsc_x[:])
        destroy!(petsc_b)
        destroy!(petsc_x)
        x
    end

    return x
end


"""
    destroy!(ksp::KSP)

Destroy the solver `ksp` holds. Does nothing on a borrowed handle, such as the
one `ksp(ts)` hands back: see [`owns`](@ref).

# External Links
$(doc_external("KSP/KSPDestroy"))
"""
function destroy!(ksp::KSP{PetscLib}) where {PetscLib}
    owns(ksp) && destroy_options!(ksp)
    return destroy_handle!(LibPETSc.KSPDestroy, ksp, PetscLib)
end



"""
    dm(ksp::AbstractKSP)

The DM attached to `ksp`, [`narrow`](@ref)ed to its flavour.

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetDM"))
"""
function dm(ksp::AbstractKSP{PetscLib}) where PetscLib
    return narrow(LibPETSc.KSPGetDM(getlib(PetscLib), ksp))
end

#
# Wrapper for calls to set_compute_rhs!
mutable struct KSPComputeRHSFn{PetscLib, PetscInt} end
function (w::KSPComputeRHSFn{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cb::CVec,
    ksp_ptr::Ptr{Cvoid},
) where {PetscLib, PetscInt}
    new_ksp = KSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age; own = false)
    b = PetscVec{PetscLib}(cb, getlib(PetscLib).age; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)::KSPState
    return run_callback("right-hand side rhs!") do
        ksp.computerhs!(b, new_ksp)
    end
end

"""
    set_compute_rhs!(rhs!::Function, ksp::AbstractKSP)

Define `rhs!` to be the right-hand side function of the `ksp`. A call to
`rhs!(b, new_ksp)` should set the elements of the PETSc vector `b` based on the
`new_ksp`.

!!! note

    The `new_ksp` passed to `rhs!` may not be the same as the `ksp` passed to
    `set_compute_rhs!`.

$(doc_callback())

# External Links
$(doc_external("KSP/KSPSetComputeRHS"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_compute_rhs! end
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function set_compute_rhs!(rhs!, ksp::AbstractKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        KSPComputeRHSFn{$PetscLib, $PetscInt}(),
        LibPETSc.PetscErrorCode,
        (CKSP, CVec, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computerhs! = rhs!
    LibPETSc.KSPSetComputeRHS($PetscLib, ksp, fptr, state_pointer(ksp))
    return ksp
end

# Wrapper for calls to set_compute_rhs!
mutable struct KSPComputeOperatorsFn{PetscLib, PetscInt} end
function (w::KSPComputeOperatorsFn{PetscLib, PetscInt})(
    new_ksp_ptr::CKSP,
    cA::CMat,
    cP::CMat,
    ksp_ptr::Ptr{Cvoid},
) where {PetscLib, PetscInt}
    new_ksp = KSP{PetscLib}(new_ksp_ptr, getlib(PetscLib).age; own = false)
    A = PetscMat{PetscLib}(cA, getlib(PetscLib).age; own = false)
    P = PetscMat{PetscLib}(cP, getlib(PetscLib).age; own = false)
    ksp = unsafe_pointer_to_objref(ksp_ptr)::KSPState
    return run_callback("operators ops!") do
        ksp.computeops!(A, P, new_ksp)
    end
end

"""
    set_compute_operators!(ops!::Function, ksp::KSP)

Define `ops!` to be the compute operators function for the `ksp`. A call to
`ops!(A, P, new_ksp)` should set the elements of the PETSc matrix linear
operator `A` and preconditioning matrix `P` based on the `new_ksp`.

!!! note

    The `new_ksp` passed to `ops!` may not be the same as the `ksp` passed to
    `set_compute_operators!`.

$(doc_callback())

# External Links
$(doc_external("KSP/KSPSetComputeOperators"))

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax
works. v0.4 also accepted the subject-first order; that method is gone in
v0.5, and there is no shim for it (§16).
"""
function set_compute_operators! end
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function set_compute_operators!(ops!, ksp::AbstractKSP{$PetscLib})
    # We must wrap the user function in our own object
    fptr = @cfunction(
        KSPComputeOperatorsFn{$PetscLib, $PetscInt}(),
        LibPETSc.PetscErrorCode,
        (CKSP, CMat, CMat, Ptr{Cvoid})
    )
    # set the computerhs! in the ksp
    ksp.computeops! = ops!
    LibPETSc.KSPSetComputeOperators($PetscLib, ksp, fptr, state_pointer(ksp))
    return ksp
end

"""
    solution(ksp::AbstractKSP)

The solution vector of `ksp`, as a `PetscVec`.

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetSolution"))
"""
solution(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetSolution(getlib(PetscLib), ksp)


"""
    type_name(ksp::AbstractKSP)

The name PETSc knows this solver by, as a `Symbol` (`:gmres`, `:cg`, …), or
`nothing` when no type has been set yet (docs/src/man/naming.md §3.1). v0.4
answered with a `String`; that is a break with no shim (§16).

# External Links
$(doc_external("KSP/KSPGetType"))
"""
type_name(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.KSPGetType(getlib(PetscLib), ksp))

"""
    set_type!(ksp::AbstractKSP, type::Symbol)

Set the Krylov method, for example `:gmres`, `:cg` or `:preonly`.

# External Links
$(doc_external("KSP/KSPSetType"))
"""
function set_type!(ksp::AbstractKSP{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.KSPSetType(getlib(PetscLib), ksp, String(type))
    return ksp
end

"""
    iteration_number(ksp::AbstractKSP)

The number of iterations the last [`solve!`](@ref) took.

# External Links
$(doc_external("KSP/KSPGetIterationNumber"))
"""
iteration_number(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetIterationNumber(getlib(PetscLib), ksp)

"""
    converged_reason(ksp::AbstractKSP)

Why the last [`solve!`](@ref) stopped, as a `LibPETSc.KSPConvergedReason`: positive
when it converged, negative when it diverged.

# External Links
$(doc_external("KSP/KSPGetConvergedReason"))
"""
converged_reason(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetConvergedReason(getlib(PetscLib), ksp)

"""
    tolerances(ksp::AbstractKSP)

The convergence tolerances of `ksp`, as `(; rtol, atol, divtol, max_it)`: the relative
and absolute residual tolerances, the divergence tolerance (relative to the initial
residual) and the iteration limit.

# External Links
$(doc_external("KSP/KSPGetTolerances"))
"""
function tolerances(ksp::AbstractKSP{PetscLib}) where {PetscLib}
    rtol, atol, divtol, max_it = LibPETSc.KSPGetTolerances(getlib(PetscLib), ksp)
    return (; rtol, atol, divtol, max_it = Int(max_it))
end

"""
    set_tolerances!(ksp::AbstractKSP; rtol, atol, divtol, max_it)

Set the convergence tolerances of `ksp`, named as [`tolerances`](@ref) returns them.
A keyword left at `nothing` keeps the value `ksp` has. Returns `ksp`.

An option for the same tolerance (`ksp_rtol = 1e-8` to the constructor, or
`-ksp_rtol` on the command line) is applied when the constructor's `options_at`
says, and overrides this setting if that comes later.

# External Links
$(doc_external("KSP/KSPSetTolerances"))
"""
function set_tolerances!(
    ksp::AbstractKSP{PetscLib};
    rtol = nothing,
    atol = nothing,
    divtol = nothing,
    max_it = nothing,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    cur = tolerances(ksp)
    LibPETSc.KSPSetTolerances(
        petsclib,
        ksp,
        PetscLib.PetscReal(something(rtol, cur.rtol)),
        PetscLib.PetscReal(something(atol, cur.atol)),
        PetscLib.PetscReal(something(divtol, cur.divtol)),
        PetscLib.PetscInt(something(max_it, cur.max_it)),
    )
    return ksp
end

"""
    set_operators!(ksp::AbstractKSP, A::AbstractPetscMat, P::AbstractPetscMat = A)

Solve with the operator `A`, building the preconditioner from `P`. Returns `ksp`.

`ksp` takes a reference to both matrices, so they stay valid inside it after the
caller destroys its own handles.

# External Links
$(doc_external("KSP/KSPSetOperators"))
"""
function set_operators!(
    ksp::AbstractKSP{PetscLib},
    A::AbstractPetscMat{PetscLib},
    P::AbstractPetscMat{PetscLib} = A,
) where {PetscLib}
    LibPETSc.KSPSetOperators(getlib(PetscLib), ksp, A, P)
    return ksp
end

const DM_ACTIVE_PARTS = (
    operator = LibPETSc.KSP_DMACTIVE_OPERATOR,
    rhs = LibPETSc.KSP_DMACTIVE_RHS,
    initial_guess = LibPETSc.KSP_DMACTIVE_INITIAL_GUESS,
    all = LibPETSc.KSP_DMACTIVE_ALL,
)

"""
    set_dm_active!(ksp::AbstractKSP, flag::Bool)
    set_dm_active!(ksp::AbstractKSP, part::Symbol, flag::Bool)

Whether the DM attached to `ksp` computes its operator, right-hand side and
initial guess. The first form sets all three; the second sets one, `part` being
`:operator`, `:rhs` or `:initial_guess` (`:all` is the first form). With `false`,
`ksp` keeps the DM for its geometry, for example for multigrid, while
[`set_operators!`](@ref) and [`solve!`](@ref) supply the rest. Returns `ksp`.

# External Links
$(doc_external("KSP/KSPSetDMActive"))
"""
function set_dm_active!(ksp::AbstractKSP{PetscLib}, part::Symbol, flag::Bool) where {PetscLib}
    haskey(DM_ACTIVE_PARTS, part) || throw(ArgumentError(
        "part must be one of $(join(repr.(keys(DM_ACTIVE_PARTS)), ", ")), got $(repr(part))",
    ))
    LibPETSc.KSPSetDMActive(getlib(PetscLib), ksp, DM_ACTIVE_PARTS[part], LibPETSc.PetscBool(flag))
    return ksp
end

set_dm_active!(ksp::AbstractKSP, flag::Bool) = set_dm_active!(ksp, :all, flag)
