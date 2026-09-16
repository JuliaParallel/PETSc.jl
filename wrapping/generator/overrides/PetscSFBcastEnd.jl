# override for PetscSFBcastEnd; C signature: PetscSFBcastEnd(PetscSF sf, MPI_Datatype unit, const void *rootdata, void *leafdata, MPI_Op op)
"""
    PetscSFBcastEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, op::MPI.Op)

End a broadcast started with `PetscSFBcastBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`rootdata` and `leafdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.


See also: `PetscSF`, `PetscSFBcastBegin()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFBcastEnd"))
"""
function PetscSFBcastEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFBcastEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata begin
        @chk ccall(
            (:PetscSFBcastEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata), op)
    end
    return nothing
end
