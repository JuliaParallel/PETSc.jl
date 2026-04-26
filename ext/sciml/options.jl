mutable struct DEOptions{SavT, TstopsT, CType, reltolType, abstolType}
    saveat::SavT
    tstops::TstopsT
    save_everystep::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
    callback::CType
    reltol::reltolType
    abstol::abstolType
    maxiters::Int
    verbose::Bool
end

_as_time_iter(::Nothing, ::Type) = ()
_as_time_iter(x::Number, ::Type{T}) where {T} = (T(x),)
_as_time_iter(x, ::Type{T}) where {T} = (T(t) for t in x)

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
    callback,
    reltol,
    abstol,
    maxiters::Integer,
    verbose::Bool,
) where {tType}
    t0 = tdir * tType(tspan[1])
    tf = tdir * tType(tspan[2])

    saveat_data = tType[
        tdir * t for t in _as_time_iter(saveat, tType) if t0 < tdir * t <= tf
    ]
    tstops_data = tType[
        tdir * t for t in _as_time_iter(tstops, tType) if t0 < tdir * t <= tf
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
        callback,
        reltol,
        abstol,
        Int(maxiters),
        verbose,
    )
end
