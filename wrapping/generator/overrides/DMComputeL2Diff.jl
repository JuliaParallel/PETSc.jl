# override for DMComputeL2Diff; C signature: DMComputeL2Diff(<not in API snapshot>)
"""
    DMComputeL2Diff(petsclib, dm, time, funcs, ctxs, X) -> PetscReal

Compute the L² norm of the difference between the global vector `X` and the
pointwise exact functions `funcs`.  `ctxs` is a `Vector{Ptr{Cvoid}}` of
context pointers (use `C_NULL` entries for no context).
"""
function DMComputeL2Diff(petsclib::PetscLibType, dm::AbstractPetscDM, time::Real,
                         funcs::Vector{Ptr{Cvoid}}, ctxs::Vector{Ptr{Cvoid}},
                         X::AbstractPetscVec) end

@for_petsc function DMComputeL2Diff(petsclib::$UnionPetscLib, dm::AbstractPetscDM,
                                    time::Real,
                                    funcs::Vector{Ptr{Cvoid}},
                                    ctxs::Vector{Ptr{Cvoid}},
                                    X::AbstractPetscVec)
    diff_ref = Ref{$PetscReal}(0)
    GC.@preserve funcs ctxs @chk ccall(
        (:DMComputeL2Diff, $petsc_library), PetscErrorCode,
        (CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, CVec, Ptr{$PetscReal}),
        dm, $PetscReal(time), funcs, ctxs, X, diff_ref)
    return diff_ref[]
end

