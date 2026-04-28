mutable struct DEOptions{SavT, TstopsT, CType, reltolType, abstolType}
    saveat::SavT
    tstops::TstopsT
    save_everystep::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
    save_discretes::Bool
    callback::CType
    reltol::reltolType
    abstol::abstolType
    maxiters::Int
end

_as_time_iter(::Nothing, ::Type) = ()
_as_time_iter(x::Number, ::Type{T}) where {T} = (T(x),)
_as_time_iter(x, ::Type{T}) where {T} = (T(t) for t in x)

# Materialize iterable `saveat` / `tstops` input into a concrete `Vector`
# *exactly once* up front. The wrapper then validates and filters that
# materialized data, so one-shot iterators (e.g. `Iterators.Stateful`)
# survive validation instead of being consumed before
# `_expand_saveat` / `_as_time_iter` can read them again.
_materialize_times(::Nothing) = nothing
_materialize_times(x::Number) = x
_materialize_times(x::AbstractVector) = x
_materialize_times(x::Tuple) = x
_materialize_times(x) = collect(x)

# Empty-check that works on both `Nothing` and any materialized iterable
# without consuming a stateful iterator. Use `_materialize_times` upstream
# so callers always pass a concrete container here.
_is_empty_times(::Nothing) = true
_is_empty_times(x::Number) = false
_is_empty_times(x) = isempty(x)

# Reject saveat input that PETSc would otherwise either silently turn into
# "save nothing" (scalar `0` / `Inf` / `NaN`) or convert to a non-finite
# timestamp (iterable element `Inf` / `NaN`). Validating here keeps the
# rest of the wrapper consistent with the `dt` / `dtmin` / `dtmax` /
# tolerance validators that fail loudly at the Julia boundary.
_validate_saveat(::Nothing) = nothing
_validate_saveat(::Tuple{}) = nothing
function _validate_saveat(saveat::Number)
    saveat isa Real || throw(ArgumentError(
        "PETSc.jl SciML extension: scalar `saveat = $(saveat)` must be a real " *
        "number, got type $(typeof(saveat)).",
    ))
    (isfinite(saveat) && saveat > 0) || throw(ArgumentError(
        "PETSc.jl SciML extension: scalar `saveat = $(saveat)` must be finite " *
        "and strictly positive (it is the save *spacing*).",
    ))
    return nothing
end
function _validate_saveat(saveat)
    for t in saveat
        (t isa Real && isfinite(t)) || throw(ArgumentError(
            "PETSc.jl SciML extension: every `saveat` entry must be a finite " *
            "real number, got $(repr(t)) of type $(typeof(t)).",
        ))
    end
    return nothing
end

# SciML semantics for scalar `saveat`: it is a *spacing*, not a single time.
# `saveat = 0.5` on tspan = (0.0, 2.0) saves at 0.5, 1.0, 1.5, 2.0 (subject to
# save_start / save_end handling). Expand a scalar value into the appropriate
# vector of timestamps in the integration direction.
function _expand_saveat(saveat::Number, tdir, tspan, ::Type{T}) where {T}
    spacing = abs(T(saveat))
    t0 = T(tspan[1])
    tf = T(tspan[2])
    return T[t0 + tdir * spacing * k
             for k in 1:floor(Int, abs(tf - t0) / spacing)]
end
_expand_saveat(saveat, _tdir, _tspan, ::Type{T}) where {T} =
    T[T(t) for t in saveat]
_expand_saveat(::Nothing, _tdir, _tspan, ::Type{T}) where {T} = T[]

function _build_opts(
    ::Type{tType},
    saveat,
    tstops,
    tdir,
    tspan;
    save_everystep::Bool,
    save_on::Bool,
    save_start::Bool,
    save_end::Bool,
    save_discretes::Bool,
    callback,
    reltol,
    abstol,
    maxiters::Integer,
) where {tType}
    t0 = tdir * tType(tspan[1])
    tf = tdir * tType(tspan[2])

    saveat_materialized = _materialize_times(saveat)
    tstops_materialized = _materialize_times(tstops)
    _validate_saveat(saveat_materialized)
    saveat_expanded = _expand_saveat(saveat_materialized, tdir, tspan, tType)
    saveat_data = tType[tdir * t for t in saveat_expanded if t0 < tdir * t <= tf]
    tstops_data = tType[
        tdir * t for t in _as_time_iter(tstops_materialized, tType) if t0 < tdir * t <= tf
    ]
    push!(tstops_data, tf)

    saveat_heap = BinaryMinHeap(saveat_data)
    tstops_heap = BinaryMinHeap(tstops_data)

    return DEOptions(
        saveat_heap,
        tstops_heap,
        save_everystep,
        save_on,
        save_start,
        save_end,
        save_discretes,
        callback,
        reltol,
        abstol,
        Int(maxiters),
    )
end
