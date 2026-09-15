# override for DMProjectFunction; C signature: DMProjectFunction(<not in API snapshot>)
"""
    DMProjectFunction(petsclib, dm, time, funcs, ctxs, mode, X)

Project the pointwise functions `funcs` (one `Ptr{Cvoid}` per field, matching
`PetscSimplePointFn` signature) into the global vector `X`.  `ctxs` is a
matching `Vector{Ptr{Cvoid}}` of context pointers (use `C_NULL` entries for
no context).
"""
function DMProjectFunction(petsclib::PetscLibType, dm::AbstractPetscDM, time::Real,
                           funcs::Vector{Ptr{Cvoid}}, ctxs::Vector{Ptr{Cvoid}},
                           mode::InsertMode, X::AbstractPetscVec) end

@for_petsc function DMProjectFunction(petsclib::$UnionPetscLib, dm::AbstractPetscDM,
                                      time::Real,
                                      funcs::Vector{Ptr{Cvoid}},
                                      ctxs::Vector{Ptr{Cvoid}},
                                      mode::InsertMode, X::AbstractPetscVec)
    GC.@preserve funcs ctxs @chk ccall(
        (:DMProjectFunction, $petsc_library), PetscErrorCode,
        (CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, InsertMode, CVec),
        dm, $PetscReal(time), funcs, ctxs, mode, X)
    return nothing
end

