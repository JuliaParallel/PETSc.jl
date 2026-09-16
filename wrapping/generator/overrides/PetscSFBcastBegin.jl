# override for PetscSFBcastBegin; C signature: PetscSFBcastBegin(PetscSF sf, MPI_Datatype unit, const void *rootdata, void *leafdata, MPI_Op op)
"""
    PetscSFBcastBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, op::MPI.Op)

Begin broadcasting root values to leaves: for each leaf, `leafdata[leaf] = op(leafdata[leaf], rootdata[root])`.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`rootdata` and `leafdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.
`op` is usually `MPI.REPLACE`; use `MPI.SUM` etc. to combine with the existing leaf values. Both arrays must stay alive until `PetscSFBcastEnd` returns.

See also: `PetscSF`, `PetscSFBcastEnd()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFBcastBegin"))
"""
function PetscSFBcastBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFBcastBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata begin
        @chk ccall(
            (:PetscSFBcastBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata), op)
    end
    return nothing
end
