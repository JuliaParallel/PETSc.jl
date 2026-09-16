# override for PetscSFReduceEnd; C signature: PetscSFReduceEnd(PetscSF sf, MPI_Datatype unit, const void *leafdata, void *rootdata, MPI_Op op)
"""
    PetscSFReduceEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata, rootdata, op::MPI.Op)

End a reduction started with `PetscSFReduceBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`leafdata` and `rootdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.


See also: `PetscSF`, `PetscSFReduceBegin()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFReduceEnd"))
"""
function PetscSFReduceEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFReduceEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve leafdata rootdata begin
        @chk ccall(
            (:PetscSFReduceEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, leafdata isa Ptr ? leafdata : pointer(leafdata), rootdata isa Ptr ? rootdata : pointer(rootdata), op)
    end
    return nothing
end
