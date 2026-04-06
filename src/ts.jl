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
