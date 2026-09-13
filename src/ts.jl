"""
    TSSetRHSFunction(petsclib, ts, r, fptr::Ptr{Cvoid}, ctx = C_NULL)

Convenience overload for low-level TS RHS callbacks created with `@cfunction`.

The generated bindings currently accept the PETSc function-wrapper type directly,
while Julia's `@cfunction` returns a raw pointer. This overload bridges that
gap so callback-based TS examples can use the low-level interface naturally.
"""
function LibPETSc.TSSetRHSFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    r::AbstractPetscVec,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetRHSFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    r::AbstractPetscVec{$PetscLib},
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSRHSFunctionFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetRHSFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{LibPETSc.TSRHSFunctionFn}, Ptr{Cvoid}),
        ts,
        r,
        typed_fptr,
        ctx,
    )
    return nothing
end

function LibPETSc.TSSetRHSFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetRHSFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSRHSFunctionFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetRHSFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{LibPETSc.TSRHSFunctionFn}, Ptr{Cvoid}),
        ts,
        C_NULL,
        typed_fptr,
        ctx,
    )
    return nothing
end

"""
    TSSetIFunction(petsclib, ts, r, fptr::Ptr{Cvoid}, ctx = C_NULL)

Convenience overload for low-level TS implicit-function callbacks created with
`@cfunction`.
"""
function LibPETSc.TSSetIFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    r::AbstractPetscVec,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetIFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    r::AbstractPetscVec{$PetscLib},
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSIFunctionFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetIFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{LibPETSc.TSIFunctionFn}, Ptr{Cvoid}),
        ts,
        r,
        typed_fptr,
        ctx,
    )
    return nothing
end

function LibPETSc.TSSetIFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetIFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSIFunctionFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetIFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{LibPETSc.TSIFunctionFn}, Ptr{Cvoid}),
        ts,
        C_NULL,
        typed_fptr,
        ctx,
    )
    return nothing
end

"""
    TSSetIJacobian(petsclib, ts, A, P, fptr::Ptr{Cvoid}, ctx = C_NULL)

Convenience overload for low-level TS implicit-Jacobian callbacks created with
`@cfunction`.
"""
function LibPETSc.TSSetIJacobian(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    A::AbstractPetscMat,
    P::AbstractPetscMat,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetIJacobian(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib},
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSIJacobianFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetIJacobian, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CMat, LibPETSc.CMat, Ptr{LibPETSc.TSIJacobianFn}, Ptr{Cvoid}),
        ts,
        A,
        P,
        typed_fptr,
        ctx,
    )
    return nothing
end

"""
    adapt = TSGetAdapt(petsclib, ts)

Return the adaptive time-step controller attached to `ts`.
"""
function LibPETSc.TSGetAdapt(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
) end

LibPETSc.@for_petsc function LibPETSc.TSGetAdapt(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
)
    adapt_ref = Ref{LibPETSc.TSAdapt}()
    LibPETSc.@chk ccall(
        (:TSGetAdapt, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{LibPETSc.TSAdapt}),
        ts,
        adapt_ref,
    )
    return adapt_ref[]
end

"""
    TSIRKGetNumStages(petsclib, ts)

Return the number of stages currently configured for a `TSIRK` method.
"""
function LibPETSc.TSIRKGetNumStages(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
) end

LibPETSc.@for_petsc function LibPETSc.TSIRKGetNumStages(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
)
    nstages_ref = Ref{$PetscInt}()
    LibPETSc.@chk ccall(
        (:TSIRKGetNumStages, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{$PetscInt}),
        ts,
        nstages_ref,
    )
    return nstages_ref[]
end

"""
    TSAdaptSetType(petsclib, adapt, type::String)

Convenience wrapper for setting the TS adaptivity controller using a Julia
string such as `"none"` or `"basic"`.
"""
function LibPETSc.TSAdaptSetType(
    petsclib::LibPETSc.PetscLibType,
    adapt::LibPETSc.TSAdapt,
    type::AbstractString,
)
    s = String(type)
    GC.@preserve s LibPETSc.TSAdaptSetType(petsclib, adapt, Base.unsafe_convert(Ptr{Cchar}, s))
    return nothing
end

"""
    TSMonitorSet(petsclib, ts, monitor::Ptr{Cvoid}, ctx = C_NULL, mdestroy = C_NULL)

Convenience overload for low-level TS monitor callbacks created with
`@cfunction`.
"""
function LibPETSc.TSMonitorSet(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    monitor::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
    mdestroy::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSMonitorSet(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    monitor::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
    mdestroy::Ptr{Cvoid} = C_NULL,
)
    typed_destroy = Ptr{LibPETSc.PetscCtxDestroyFn}(mdestroy)
    LibPETSc.@chk ccall(
        (:TSMonitorSet, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.external, Ptr{Cvoid}, Ptr{LibPETSc.PetscCtxDestroyFn}),
        ts,
        monitor,
        ctx,
        typed_destroy,
    )
    return nothing
end

"""
    TSARKIMEXRegister(
        petsclib,
        name::String,
        order,
        s,
        At,
        bt,
        ct,
        A,
        b,
        c,
        bembedt,
        bembed,
        pinterp,
        binterpt,
        binterp,
    )

Julia-friendly overload for registering a custom `TSARKIMEX` tableau. Optional
PETSc arrays may be passed as `nothing`, which is translated to `NULL`.

For the stage tables `At` and `A`, PETSc expects flat vectors in row-major
order, matching the layout used by C arrays. If you start from a Julia matrix,
do not pass `vec(A)` directly since Julia stores matrices column-major; flatten
row-by-row instead, for example with `vec(permutedims(A))`.
"""
function LibPETSc.TSARKIMEXRegister(
    petsclib::LibPETSc.PetscLibType,
    name::String,
    order::Integer,
    s::Integer,
    At::AbstractVector,
    bt::Union{Nothing, AbstractVector},
    ct::Union{Nothing, AbstractVector},
    A::AbstractVector,
    b::Union{Nothing, AbstractVector},
    c::Union{Nothing, AbstractVector},
    bembedt::Union{Nothing, AbstractVector},
    bembed::Union{Nothing, AbstractVector},
    pinterp::Integer,
    binterpt::Union{Nothing, AbstractVector},
    binterp::Union{Nothing, AbstractVector},
) end

LibPETSc.@for_petsc function LibPETSc.TSARKIMEXRegister(
    petsclib::$UnionPetscLib,
    name::String,
    order::Integer,
    s::Integer,
    At::AbstractVector,
    bt::Union{Nothing, AbstractVector},
    ct::Union{Nothing, AbstractVector},
    A::AbstractVector,
    b::Union{Nothing, AbstractVector},
    c::Union{Nothing, AbstractVector},
    bembedt::Union{Nothing, AbstractVector},
    bembed::Union{Nothing, AbstractVector},
    pinterp::Integer,
    binterpt::Union{Nothing, AbstractVector},
    binterp::Union{Nothing, AbstractVector},
)
    At_vals = $PetscReal.(At)
    A_vals = $PetscReal.(A)
    bt_vals = bt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bt)
    ct_vals = ct === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(ct)
    b_vals = b === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(b)
    c_vals = c === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(c)
    bembedt_vals =
        bembedt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bembedt)
    bembed_vals =
        bembed === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bembed)
    binterpt_vals =
        binterpt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(binterpt)
    binterp_vals =
        binterp === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(binterp)

    LibPETSc.@chk ccall(
        (:TSARKIMEXRegister, $petsc_library),
        LibPETSc.PetscErrorCode,
        (
            LibPETSc.TSARKIMEXType,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        name,
        $PetscInt(order),
        $PetscInt(s),
        At_vals,
        bt_vals,
        ct_vals,
        A_vals,
        b_vals,
        c_vals,
        bembedt_vals,
        bembed_vals,
        $PetscInt(pinterp),
        binterpt_vals,
        binterp_vals,
    )
    return nothing
end

# ============================================================================
#   High-level TS interface
# ============================================================================

import .LibPETSc:
    AbstractTS, CTS, TS, AbstractPetscDM, AbstractPetscVec, PetscVec, CVec

function Base.show(io::IO, ts::AbstractTS{PetscLib}) where {PetscLib}
    if ts.ptr == C_NULL
        print(io, "PETSc TS (null pointer)")
    else
        print(io, "PETSc TS object")
    end
    return nothing
end

"""
    TS(petsclib, comm::MPI.Comm; prefix = "", options...)

Create a PETSc time stepper on the communicator `comm`.

Options are stored and applied in [`solve!`](@ref) rather than here, so that a
DM and callbacks attached after construction are visible to `TSSetFromOptions`.

`exact_final_time` defaults to `TS_EXACTFINALTIME_MATCHSTEP`, so the last step
lands on the time set by [`set_max_time!`](@ref). 
PETSc's own default is `TS_EXACTFINALTIME_UNSPECIFIED`, which integrates 
to the wrong time without reporting an error.

If `comm` has size 1 the garbage collector calls [`destroy!`](@ref).
Otherwise destruction is the caller's responsibility.

# External Links
$(_doc_external("TS/TSCreate"))
$(_doc_external("TS/TSSetExactFinalTime"))
"""
function TS(
    petsclib::PetscLib,
    comm::MPI.Comm;
    prefix::String = "",
    exact_final_time::LibPETSc.TSExactFinalTimeOption = LibPETSc.TS_EXACTFINALTIME_MATCHSTEP,
    options...,
) where {PetscLib}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    ts = LibPETSc.TSCreate(petsclib, comm)

    if !isempty(prefix)
        LibPETSc.TSSetOptionsPrefix(petsclib, ts, prefix)
    end

    LibPETSc.TSSetExactFinalTime(petsclib, ts, exact_final_time)

    if !isempty(options)
        ts.opts = Options(petsclib; options...)
    end

    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, ts)
    end

    return ts
end

"""
    set_exact_final_time!(ts::AbstractTS, option)

Choose how the final step meets the time set by [`set_max_time!`](@ref):
`TS_EXACTFINALTIME_MATCHSTEP`, `TS_EXACTFINALTIME_INTERPOLATE` or `TS_EXACTFINALTIME_STEPOVER`.

# External Links
$(_doc_external("TS/TSSetExactFinalTime"))
"""
function set_exact_final_time!(
    ts::AbstractTS{PetscLib},
    option::LibPETSc.TSExactFinalTimeOption,
) where {PetscLib}
    LibPETSc.TSSetExactFinalTime(getlib(PetscLib), ts, option)
    return nothing
end

"""
    set_adapt_type!(ts::AbstractTS, type::Symbol)

Set the timestep adaptivity controller, for example `:none` to hold the step size fixed, 
or `:basic` for the default error-based controller.

# External Links
$(_doc_external("TS/TSAdaptSetType"))
"""
function set_adapt_type!(ts::AbstractTS{PetscLib}, type::Symbol) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.TSAdaptSetType(
        petsclib,
        LibPETSc.TSGetAdapt(petsclib, ts),
        String(type),
    )
    return nothing
end

"""
    destroy!(ts::AbstractTS)

Destroy `ts` and release the options database attached to it.

The call is a no-op when the library has been finalized or when `ts` predates
the current initialize/finalize cycle, so a stale handle never reaches `TSDestroy`.

# External Links
$(_doc_external("TS/TSDestroy"))
"""
function destroy!(ts::AbstractTS{PetscLib}) where {PetscLib}
    if !isnothing(ts.opts)
        destroy(ts.opts)
        ts.opts = nothing
    end
    if isdestroyable(ts, PetscLib)
        LibPETSc.TSDestroy(PetscLib, ts)
    end
    ts.ptr = C_NULL
    return nothing
end

"""
    destroy(ts::AbstractTS)

Destroy `ts`. Provided so that the spelling used by the rest of the package keeps working; 
[`destroy!`](@ref) is the name to prefer, since the call mutates `ts`.

# External Links
$(_doc_external("TS/TSDestroy"))
"""
destroy(ts::AbstractTS) = destroy!(ts)

"""
    comm(ts::AbstractTS)

The MPI communicator `ts` was built on.

# External Links
$(_doc_external("Sys/PetscObjectGetComm"))
"""
comm(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.PetscObjectGetComm(PetscLib, ts)

"""
    type(ts::AbstractTS)

The time-stepping method currently set on `ts`, as a `Symbol`, or `nothing` when none has been set yet.

PETSc reports an unset type as a null string, which the generated `TSGetType` cannot convert. 
Asking a freshly created stepper for its type is reasonable, so it is answered with `nothing` here.

# External Links
$(_doc_external("TS/TSGetType"))
"""
function type end

LibPETSc.@for_petsc function type(ts::AbstractTS{$PetscLib})
    r_type = Ref{Ptr{Cchar}}(C_NULL)
    LibPETSc.@chk ccall(
        (:TSGetType, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{Ptr{Cchar}}),
        ts,
        r_type,
    )
    return r_type[] == C_NULL ? nothing : Symbol(unsafe_string(r_type[]))
end

"""
    set_type!(ts::AbstractTS, type::Symbol)

Set the time-stepping method, for example `:bdf`, `:rk` or `:arkimex`.

# External Links
$(_doc_external("TS/TSSetType"))
"""
function set_type!(ts::AbstractTS{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.TSSetType(getlib(PetscLib), ts, String(type))
    return nothing
end

"""
    set_problem_type!(ts::AbstractTS, type)

Declare the problem as `LibPETSc.TS_LINEAR` or `LibPETSc.TS_NONLINEAR`.

# External Links
$(_doc_external("TS/TSSetProblemType"))
"""
function set_problem_type!(
    ts::AbstractTS{PetscLib},
    type::LibPETSc.TSProblemType,
) where {PetscLib}
    LibPETSc.TSSetProblemType(getlib(PetscLib), ts, type)
    return nothing
end

"""
    dm(ts::AbstractTS)

The DM attached to `ts`. The DM is owned by `ts`.

# External Links
$(_doc_external("TS/TSGetDM"))
"""
dm(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetDM(getlib(PetscLib), ts)

"""
    set_dm!(ts::AbstractTS, dm::AbstractPetscDM)

Attach `dm` to `ts`.

# External Links
$(_doc_external("TS/TSSetDM"))
"""
function set_dm!(
    ts::AbstractTS{PetscLib},
    dm::AbstractPetscDM{PetscLib},
) where {PetscLib}
    LibPETSc.TSSetDM(getlib(PetscLib), ts, dm)
    return nothing
end

"""
    TSGetSolution(petsclib, ts)

Return the solution vector held by `ts`.

The generated three-argument form takes the vector as an input and nulls 
the caller's handle, so it cannot be used to read the solution back. 
The vector is owned by `ts` and must not be destroyed.

# External Links
$(_doc_external("TS/TSGetSolution"))
"""
function LibPETSc.TSGetSolution(petsclib::LibPETSc.PetscLibType, ts::TS) end

LibPETSc.@for_petsc function LibPETSc.TSGetSolution(
    petsclib::$UnionPetscLib,
    ts::TS,
)
    v_ = Ref{CVec}(C_NULL)
    LibPETSc.@chk ccall(
        (:TSGetSolution, $petsc_library),
        LibPETSc.PetscErrorCode,
        (CTS, Ptr{CVec}),
        ts,
        v_,
    )
    return PetscVec(v_[], petsclib)
end

"""
    solution(ts::AbstractTS)

The solution vector held by `ts`. It is owned by `ts`, so do not destroy it.

# External Links
$(_doc_external("TS/TSGetSolution"))
"""
solution(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSolution(getlib(PetscLib), ts)

"""
    set_solution!(ts::AbstractTS, u::AbstractPetscVec)

Set the initial condition of `ts`.

# External Links
$(_doc_external("TS/TSSetSolution"))
"""
function set_solution!(
    ts::AbstractTS{PetscLib},
    u::AbstractPetscVec{PetscLib},
) where {PetscLib}
    LibPETSc.TSSetSolution(getlib(PetscLib), ts, u)
    return nothing
end

# Time and step controls
# ----------------------------------------------------------------------------

"""
    current_time(ts::AbstractTS)

The time `ts` has reached.

# External Links
$(_doc_external("TS/TSGetTime"))
"""
current_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetTime(getlib(PetscLib), ts)

"""
    set_time!(ts::AbstractTS, t)

Set the current time of `ts`.

# External Links
$(_doc_external("TS/TSSetTime"))
"""
function set_time!(ts::AbstractTS{PetscLib}, t) where {PetscLib}
    LibPETSc.TSSetTime(getlib(PetscLib), ts, PetscLib.PetscReal(t))
    return nothing
end

"""
    timestep(ts::AbstractTS)

The current step size.

# External Links
$(_doc_external("TS/TSGetTimeStep"))
"""
timestep(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetTimeStep(getlib(PetscLib), ts)

"""
    set_timestep!(ts::AbstractTS, dt)

Set the step size.

# External Links
$(_doc_external("TS/TSSetTimeStep"))
"""
function set_timestep!(ts::AbstractTS{PetscLib}, dt) where {PetscLib}
    LibPETSc.TSSetTimeStep(getlib(PetscLib), ts, PetscLib.PetscReal(dt))
    return nothing
end

"""
    max_time(ts::AbstractTS)

The time at which integration stops.

# External Links
$(_doc_external("TS/TSGetMaxTime"))
"""
max_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetMaxTime(getlib(PetscLib), ts)

"""
    set_max_time!(ts::AbstractTS, t)

Set the time at which integration stops.

# External Links
$(_doc_external("TS/TSSetMaxTime"))
"""
function set_max_time!(ts::AbstractTS{PetscLib}, t) where {PetscLib}
    LibPETSc.TSSetMaxTime(getlib(PetscLib), ts, PetscLib.PetscReal(t))
    return nothing
end

"""
    max_steps(ts::AbstractTS)

The step count at which integration stops.

# External Links
$(_doc_external("TS/TSGetMaxSteps"))
"""
max_steps(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetMaxSteps(getlib(PetscLib), ts)

"""
    set_max_steps!(ts::AbstractTS, n)

Set the step count at which integration stops.

# External Links
$(_doc_external("TS/TSSetMaxSteps"))
"""
function set_max_steps!(ts::AbstractTS{PetscLib}, n) where {PetscLib}
    LibPETSc.TSSetMaxSteps(getlib(PetscLib), ts, PetscLib.PetscInt(n))
    return nothing
end

"""
    step_number(ts::AbstractTS)

The number of steps taken so far.

# External Links
$(_doc_external("TS/TSGetStepNumber"))
"""
step_number(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetStepNumber(getlib(PetscLib), ts)

"""
    tolerances(ts::AbstractTS)

Local truncation error tolerances, as `(; atol, rtol, vatol, vrtol)`.

`vatol` and `vrtol` hold per-component tolerances and carry a null pointer when
only the scalar tolerances are set. Both are owned by `ts`.

# External Links
$(_doc_external("TS/TSGetTolerances"))
"""
function tolerances(ts::AbstractTS{PetscLib}) where {PetscLib}
    atol, vatol, rtol, vrtol = LibPETSc.TSGetTolerances(getlib(PetscLib), ts)
    return (; atol, rtol, vatol, vrtol)
end

"""
    set_tolerances!(ts::AbstractTS; atol, rtol, vatol, vrtol)

Set the local truncation error tolerances.

Pass `vatol` or `vrtol` to give per-component tolerances; 
leaving them at `nothing` selects the scalar tolerance.

# External Links
$(_doc_external("TS/TSSetTolerances"))
"""
function set_tolerances!(
    ts::AbstractTS{PetscLib};
    atol = 1e-8,
    rtol = 1e-6,
    vatol::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
    vrtol::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscReal = PetscLib.PetscReal
    null_vec = LibPETSc.PetscVec(petsclib)
    LibPETSc.TSSetTolerances(
        petsclib,
        ts,
        PetscReal(atol),
        isnothing(vatol) ? null_vec : vatol,
        PetscReal(rtol),
        isnothing(vrtol) ? null_vec : vrtol,
    )
    return nothing
end

"""
    converged_reason(ts::AbstractTS)

Why the integration stopped, as a `LibPETSc.TSConvergedReason`.

# External Links
$(_doc_external("TS/TSGetConvergedReason"))
"""
converged_reason(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetConvergedReason(getlib(PetscLib), ts)

"""
    solve_time(ts::AbstractTS)

The time reached by the last [`solve!`](@ref).

This is the time the integration actually stopped at, which is not the time
asked for by [`set_max_time!`](@ref) unless the final step was made to land on it; 
see [`set_exact_final_time!`](@ref).

# External Links
$(_doc_external("TS/TSGetSolveTime"))
"""
solve_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSolveTime(getlib(PetscLib), ts)

"""
    snes_iterations(ts::AbstractTS)

Total number of nonlinear iterations taken so far, summed over the steps.

# External Links
$(_doc_external("TS/TSGetSNESIterations"))
"""
snes_iterations(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNESIterations(getlib(PetscLib), ts)

"""
    ksp_iterations(ts::AbstractTS)

Total number of linear iterations taken so far, summed over the steps.

# External Links
$(_doc_external("TS/TSGetKSPIterations"))
"""
ksp_iterations(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetKSPIterations(getlib(PetscLib), ts)

"""
    step_rejections(ts::AbstractTS)

Number of steps the adaptivity controller has rejected.

# External Links
$(_doc_external("TS/TSGetStepRejections"))
"""
step_rejections(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetStepRejections(getlib(PetscLib), ts)

"""
    snes_failures(ts::AbstractTS)

Number of failed nonlinear solves.

# External Links
$(_doc_external("TS/TSGetSNESFailures"))
"""
snes_failures(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNESFailures(getlib(PetscLib), ts)

"""
    TSGetSNES(petsclib, ts)

Return the nonlinear solver held by `ts`.

The generated three-argument form takes the solver as an input and nulls the
caller's handle, so it cannot be used to read the solver back. 
The `SNES` is owned by `ts` and must not be destroyed.

# External Links
$(_doc_external("TS/TSGetSNES"))
"""
function LibPETSc.TSGetSNES(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
) end

LibPETSc.@for_petsc function LibPETSc.TSGetSNES(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
)
    snes_ = Ref{LibPETSc.CSNES}(C_NULL)
    LibPETSc.@chk ccall(
        (:TSGetSNES, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{LibPETSc.CSNES}),
        ts,
        snes_,
    )
    return LibPETSc.PetscSNES(snes_[], petsclib)
end

"""
    TSGetKSP(petsclib, ts)

Return the linear solver held by `ts`.

The generated three-argument form takes the solver as an input and nulls the
caller's handle, so it cannot be used to read the solver back. 
The `KSP` is owned by `ts` and must not be destroyed.

# External Links
$(_doc_external("TS/TSGetKSP"))
"""
function LibPETSc.TSGetKSP(petsclib::LibPETSc.PetscLibType, ts::LibPETSc.TS) end

LibPETSc.@for_petsc function LibPETSc.TSGetKSP(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
)
    ksp_ = Ref{LibPETSc.CKSP}(C_NULL)
    LibPETSc.@chk ccall(
        (:TSGetKSP, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{LibPETSc.CKSP}),
        ts,
        ksp_,
    )
    return LibPETSc.PetscKSP(ksp_[], petsclib)
end

"""
    snes(ts::AbstractTS)

The nonlinear solver `ts` steps with. It is owned by `ts`, so do not destroy it.

Only the implicit methods build one. Asking an explicit method for its `SNES`
creates an unused solver rather than reporting an error.

# External Links
$(_doc_external("TS/TSGetSNES"))
"""
snes(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNES(getlib(PetscLib), ts)

"""
    ksp(ts::AbstractTS)

The linear solver `ts` steps with. It is owned by `ts`, so do not destroy it.

PETSc only offers this for a problem declared `TS_LINEAR` with [`set_problem_type!`](@ref), 
and raises `PETSC_ERR_ARG_WRONG` otherwise. The linear solver of a nonlinear problem 
belongs to its `SNES`, so reach it through [`snes`](@ref).

# External Links
$(_doc_external("TS/TSGetKSP"))
"""
ksp(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetKSP(getlib(PetscLib), ts)

"""
    setup!(ts::AbstractTS)

Complete the setup of `ts`. [`solve!`](@ref) calls this, so it is only needed
when the setup must happen at a controlled point.

# External Links
$(_doc_external("TS/TSSetUp"))
"""
function setup!(ts::AbstractTS{PetscLib}) where {PetscLib}
    LibPETSc.TSSetUp(getlib(PetscLib), ts)
    return nothing
end

"""
    set_from_options!(ts::AbstractTS)

Apply the PETSc options database to `ts`.

# External Links
$(_doc_external("TS/TSSetFromOptions"))
"""
function set_from_options!(ts::AbstractTS{PetscLib}) where {PetscLib}
    LibPETSc.TSSetFromOptions(getlib(PetscLib), ts)
    return nothing
end

"""
    TSSolve(petsclib, ts, ::Nothing)

Integrate `ts` using the solution vector already set on it.

The generated binding requires a vector. PETSc reads the solution set by
`TSSetSolution` when it is handed `NULL` instead, which is what this overload passes.

# External Links
$(_doc_external("TS/TSSolve"))
"""
function LibPETSc.TSSolve(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    ::Nothing,
) end

LibPETSc.@for_petsc function LibPETSc.TSSolve(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    ::Nothing,
)
    LibPETSc.@chk ccall(
        (:TSSolve, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec),
        ts,
        C_NULL,
    )
    return nothing
end

# Options given to the constructor are applied here rather than there, so that
# a DM and callbacks attached in between are visible to `TSSetFromOptions`.
function _with_options(f, ts::AbstractTS{PetscLib}) where {PetscLib}
    isnothing(ts.opts) && return f()
    push!(ts.opts)
    try
        LibPETSc.TSSetFromOptions(PetscLib, ts)
        return f()
    finally
        pop!(ts.opts)
    end
end

"""
    solve!(u::AbstractPetscVec, ts::AbstractTS)
    solve!(ts::AbstractTS)

Integrate `ts`, starting from `u` and returning it. 
The second form uses the solution vector already set on `ts` and returns `ts`.

`u` is registered with `TSSetSolution` before the solve. 
Passing it to `TSSolve` alone leaves the stepper reading an uninitialized solution vector,
which the implicit methods see as an initial condition of zero.

Options passed to the [`TS`](@ref) constructor are applied here, 
once the DM and the callbacks are attached.

# External Links
$(_doc_external("TS/TSSolve"))
$(_doc_external("TS/TSSetSolution"))
"""
function solve!(
    u::AbstractPetscVec{PetscLib},
    ts::AbstractTS{PetscLib},
) where {PetscLib}
    LibPETSc.TSSetSolution(PetscLib, ts, u)
    _with_options(ts) do
        LibPETSc.TSSolve(PetscLib, ts, u)
    end
    return u
end

function solve!(ts::AbstractTS{PetscLib}) where {PetscLib}
    _with_options(ts) do
        LibPETSc.TSSolve(PetscLib, ts, nothing)
    end
    return ts
end

"""
    step!(ts::AbstractTS)

Take a single step. Unlike [`solve!`](@ref) this does not apply the options
given to the constructor, and it ignores the time set by [`set_max_time!`](@ref).

# External Links
$(_doc_external("TS/TSStep"))
"""
function step!(ts::AbstractTS{PetscLib}) where {PetscLib}
    LibPETSc.TSStep(getlib(PetscLib), ts)
    return nothing
end

"""
    reset!(ts::AbstractTS)

Release the work vectors and matrices `ts` allocated, keeping the callbacks
and the options.

The clock is not part of that state: the time and the step count survive, so
set them with [`set_time!`](@ref) and [`set_max_time!`](@ref) before integrating again.

# External Links
$(_doc_external("TS/TSReset"))
"""
function reset!(ts::AbstractTS{PetscLib}) where {PetscLib}
    LibPETSc.TSReset(getlib(PetscLib), ts)
    return nothing
end

"""
    interpolate!(u::AbstractPetscVec, ts::AbstractTS, t)

Fill `u` with the solution interpolated to time `t`, and return it.

Only the methods that keep a dense output can do this, 
and `t` must lie inside the step just taken.

# External Links
$(_doc_external("TS/TSInterpolate"))
"""
function interpolate!(
    u::AbstractPetscVec{PetscLib},
    ts::AbstractTS{PetscLib},
    t,
) where {PetscLib}
    LibPETSc.TSInterpolate(getlib(PetscLib), ts, PetscLib.PetscReal(t), u)
    return u
end

"""
    user_ctx(ts::AbstractTS)

Whatever was stored with [`set_user_ctx!`](@ref), or `nothing`.
"""
user_ctx(ts::AbstractTS) = ts.user_ctx

"""
    set_user_ctx!(ts::AbstractTS, ctx)

Attach `ctx` to `ts`, to be handed back as the last argument of every callback
that has a method accepting it.

The object is held by `ts` on the Julia side, so it is kept alive and needs no pinning.
"""
function set_user_ctx!(ts::AbstractTS, ctx)
    ts.user_ctx = ctx
    return nothing
end

# Callbacks
# ----------------------------------------------------------------------------
#
# Each setter stores the Julia function on the `ts` and hands PETSc a
# `@cfunction` trampoline plus a pointer to the `ts` itself as the context, 
# so the closure stays rooted for as long as the object lives. 

# A callback may return an error code; anything else is treated as success.

"""
    TSSetRHSJacobian(petsclib, ts, A, P, fptr::Ptr{Cvoid}, ctx = C_NULL)

Convenience overload for low-level TS RHS-Jacobian callbacks created with `@cfunction`.

# External Links
$(_doc_external("TS/TSSetRHSJacobian"))
"""
function LibPETSc.TSSetRHSJacobian(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.TS,
    A::AbstractPetscMat,
    P::AbstractPetscMat,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetRHSJacobian(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.TS,
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib},
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    typed_fptr = Base.unsafe_convert(Ptr{LibPETSc.TSRHSJacobianFn}, fptr)
    LibPETSc.@chk ccall(
        (:TSSetRHSJacobian, $petsc_library),
        LibPETSc.PetscErrorCode,
        (
            LibPETSc.CTS,
            LibPETSc.CMat,
            LibPETSc.CMat,
            Ptr{LibPETSc.TSRHSJacobianFn},
            Ptr{Cvoid},
        ),
        ts,
        A,
        P,
        typed_fptr,
        ctx,
    )
    return nothing
end

# A callback may return a PETSc error code; anything else counts as success.
_errorcode(r) =
    r isa Integer ? LibPETSc.PetscErrorCode(r) : LibPETSc.PetscErrorCode(0)

# "error in library called by PETSc", from `petscsystypes.h`.
const _PETSC_ERR_LIB = LibPETSc.PetscErrorCode(76)

# Call `f` and turn its result into a PETSc error code.
#
# A Julia exception must not cross the `@cfunction` boundary: PETSc is C and
# cannot unwind a Julia frame, so an escaping error takes the process down. 
# Log it here instead and report failure to PETSc, which unwinds its own stack 
# and leaves `@chk` to raise a `PetscError` from the enclosing `solve!`.
function _run_callback(f, name)
    try
        return _errorcode(f())
    catch e
        bt = catch_backtrace()
        # Reporting is itself Julia code, and nothing here may throw either.
        try
            @error "PETSc.jl: the $name callback failed" exception = (e, bt)
        catch
            Core.println("PETSc.jl: the ", name, " callback failed")
        end
        return _PETSC_ERR_LIB
    end
end

"""
    set_rhs_function!(f!, ts::AbstractTS, r = nothing)
    set_rhs_function!(ts::AbstractTS, f!, r = nothing)

Set the right-hand side ``G`` of an explicit problem ``du/dt = G(t, u)``.

`f!` is called as `f!(F, ts, t, u)`, filling the vector `F`. If `ts.user_ctx`
is set, `f!(F, ts, t, u, user_ctx)` is used instead when that method exists.

`r` is an optional template vector for the residual.

# External Links
$(_doc_external("TS/TSSetRHSFunction"))
"""
set_rhs_function!(ts::AbstractTS, f!, r = nothing) =
    set_rhs_function!(f!, ts, r)

mutable struct TSSetRHSFunctionFn{PetscLib, PetscReal} end
function (::TSSetRHSFunctionFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    F_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    ts = unsafe_pointer_to_objref(ctx)
    actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age)
    u = PetscVec{PetscLib}(u_ptr)
    F = PetscVec{PetscLib}(F_ptr)

    _run_callback("rhs_function!") do
        if Base.applicable(ts.rhs_function!, F, actual_ts, t, u, ts.user_ctx)
            ts.rhs_function!(F, actual_ts, t, u, ts.user_ctx)
        else
            ts.rhs_function!(F, actual_ts, t, u)
        end
    end
end

LibPETSc.@for_petsc function set_rhs_function!(
    f!,
    ts::AbstractTS{$PetscLib},
    r::Union{Nothing, AbstractPetscVec{$PetscLib}} = nothing,
)
    fptr = @cfunction(
        TSSetRHSFunctionFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, Ptr{Cvoid})
    )
    ts.rhs_function! = f!
    LibPETSc.TSSetRHSFunction($PetscLib, ts, r, fptr, pointer_from_objref(ts))
    return nothing
end

"""
    set_rhs_jacobian!(updateJ!, ts::AbstractTS, A, P = A)
    set_rhs_jacobian!(ts::AbstractTS, updateJ!, A, P = A)

Set the Jacobian of the right-hand side ``G``.

`updateJ!` is called as `updateJ!(A, P, ts, t, u)`, filling the Jacobian `A`
and the preconditioning matrix `P`. If `ts.user_ctx` is set,
`updateJ!(A, P, ts, t, u, user_ctx)` is used instead when that method exists.

# External Links
$(_doc_external("TS/TSSetRHSJacobian"))
"""
set_rhs_jacobian!(ts::AbstractTS, updateJ!, A, P = A) =
    set_rhs_jacobian!(updateJ!, ts, A, P)

mutable struct TSSetRHSJacobianFn{PetscLib, PetscReal} end
function (::TSSetRHSJacobianFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    A_ptr::CMat,
    P_ptr::CMat,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    ts = unsafe_pointer_to_objref(ctx)
    actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age)
    u = PetscVec{PetscLib}(u_ptr)
    A = PetscMat{PetscLib}(A_ptr)
    P = PetscMat{PetscLib}(P_ptr)

    _run_callback("rhs_jacobian!") do
        if Base.applicable(ts.rhs_jacobian!, A, P, actual_ts, t, u, ts.user_ctx)
            ts.rhs_jacobian!(A, P, actual_ts, t, u, ts.user_ctx)
        else
            ts.rhs_jacobian!(A, P, actual_ts, t, u)
        end
    end
end

LibPETSc.@for_petsc function set_rhs_jacobian!(
    updateJ!,
    ts::AbstractTS{$PetscLib},
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib} = A,
)
    fptr = @cfunction(
        TSSetRHSJacobianFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CMat, CMat, Ptr{Cvoid})
    )
    ts.rhs_jacobian! = updateJ!
    LibPETSc.TSSetRHSJacobian(
        $PetscLib,
        ts,
        A,
        P,
        fptr,
        pointer_from_objref(ts),
    )
    return nothing
end

"""
    set_ifunction!(f!, ts::AbstractTS, r = nothing)
    set_ifunction!(ts::AbstractTS, f!, r = nothing)

Set the residual ``F`` of an implicit problem ``F(t, u, du/dt) = 0``.

`f!` is called as `f!(F, ts, t, u, u_t)`, filling the vector `F`. 
If `ts.user_ctx` is set, `f!(F, ts, t, u, u_t, user_ctx)` is used instead 
when that method exists.

# External Links
$(_doc_external("TS/TSSetIFunction"))
"""
set_ifunction!(ts::AbstractTS, f!, r = nothing) = set_ifunction!(f!, ts, r)

mutable struct TSSetIFunctionFn{PetscLib, PetscReal} end
function (::TSSetIFunctionFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    udot_ptr::CVec,
    F_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    ts = unsafe_pointer_to_objref(ctx)
    actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age)
    u = PetscVec{PetscLib}(u_ptr)
    u_t = PetscVec{PetscLib}(udot_ptr)
    F = PetscVec{PetscLib}(F_ptr)

    _run_callback("ifunction!") do
        if Base.applicable(ts.ifunction!, F, actual_ts, t, u, u_t, ts.user_ctx)
            ts.ifunction!(F, actual_ts, t, u, u_t, ts.user_ctx)
        else
            ts.ifunction!(F, actual_ts, t, u, u_t)
        end
    end
end

LibPETSc.@for_petsc function set_ifunction!(
    f!,
    ts::AbstractTS{$PetscLib},
    r::Union{Nothing, AbstractPetscVec{$PetscLib}} = nothing,
)
    fptr = @cfunction(
        TSSetIFunctionFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, CVec, Ptr{Cvoid})
    )
    ts.ifunction! = f!
    LibPETSc.TSSetIFunction($PetscLib, ts, r, fptr, pointer_from_objref(ts))
    return nothing
end

"""
    set_ijacobian!(updateJ!, ts::AbstractTS, A, P = A)
    set_ijacobian!(ts::AbstractTS, updateJ!, A, P = A)

Set the Jacobian of the implicit residual ``F``.

`updateJ!` is called as `updateJ!(A, P, ts, t, u, u_t, shift)` and should fill
`A` with ``dF/du + shift * dF/du_t``. If `ts.user_ctx` is set, the method
taking a trailing `user_ctx` is used instead when it exists.

# External Links
$(_doc_external("TS/TSSetIJacobian"))
"""
set_ijacobian!(ts::AbstractTS, updateJ!, A, P = A) =
    set_ijacobian!(updateJ!, ts, A, P)

mutable struct TSSetIJacobianFn{PetscLib, PetscReal} end
function (::TSSetIJacobianFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    udot_ptr::CVec,
    shift::PetscReal,
    A_ptr::CMat,
    P_ptr::CMat,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    ts = unsafe_pointer_to_objref(ctx)
    actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age)
    u = PetscVec{PetscLib}(u_ptr)
    u_t = PetscVec{PetscLib}(udot_ptr)
    A = PetscMat{PetscLib}(A_ptr)
    P = PetscMat{PetscLib}(P_ptr)

    _run_callback("ijacobian!") do
        if Base.applicable(
            ts.ijacobian!,
            A,
            P,
            actual_ts,
            t,
            u,
            u_t,
            shift,
            ts.user_ctx,
        )
            ts.ijacobian!(A, P, actual_ts, t, u, u_t, shift, ts.user_ctx)
        else
            ts.ijacobian!(A, P, actual_ts, t, u, u_t, shift)
        end
    end
end

LibPETSc.@for_petsc function set_ijacobian!(
    updateJ!,
    ts::AbstractTS{$PetscLib},
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib} = A,
)
    fptr = @cfunction(
        TSSetIJacobianFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, CMat, Ptr{Cvoid})
    )
    ts.ijacobian! = updateJ!
    LibPETSc.TSSetIJacobian(
        $PetscLib,
        ts,
        A,
        P,
        fptr,
        pointer_from_objref(ts),
    )
    return nothing
end

"""
    set_monitor!(f, ts::AbstractTS)
    set_monitor!(ts::AbstractTS, f)

Call `f` once after every accepted step.

`f` is called as `f(ts, step, t, u)`, where `step` counts the steps taken, `t`
is the time reached and `u` holds the solution there. If `ts.user_ctx` is set,
`f(ts, step, t, u, user_ctx)` is used instead when that method exists. 
`u` is owned by `ts` and must not be destroyed.

Only one monitor can be set this way; a second call replaces the first. 
The monitors PETSc installs from the options database, such as `-ts_monitor`, 
are unaffected.

# External Links
$(_doc_external("TS/TSMonitorSet"))
"""
set_monitor!(ts::AbstractTS, f) = set_monitor!(f, ts)

mutable struct TSMonitorSetFn{PetscLib, PetscInt, PetscReal} end
function (::TSMonitorSetFn{PetscLib, PetscInt, PetscReal})(
    ts_ptr::CTS,
    step::PetscInt,
    t::PetscReal,
    u_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscInt, PetscReal}
    ts = unsafe_pointer_to_objref(ctx)
    actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age)
    u = PetscVec{PetscLib}(u_ptr)

    _run_callback("monitor") do
        if Base.applicable(ts.monitor, actual_ts, step, t, u, ts.user_ctx)
            ts.monitor(actual_ts, step, t, u, ts.user_ctx)
        else
            ts.monitor(actual_ts, step, t, u)
        end
    end
end

LibPETSc.@for_petsc function set_monitor!(f, ts::AbstractTS{$PetscLib})
    fptr = @cfunction(
        TSMonitorSetFn{$PetscLib, $PetscInt, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid})
    )
    ts.monitor = f
    LibPETSc.TSMonitorSet($PetscLib, ts, fptr, pointer_from_objref(ts))
    return nothing
end
