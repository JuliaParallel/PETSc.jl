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
    type::String,
)
    c_str = Vector{UInt8}(type * "\0")
    ptr = Base.unsafe_convert(Ptr{Cchar}, pointer(c_str))
    LibPETSc.TSAdaptSetType(petsclib, adapt, ptr)
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
