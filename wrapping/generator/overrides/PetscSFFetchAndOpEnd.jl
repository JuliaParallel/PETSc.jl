# override for PetscSFFetchAndOpEnd; C signature: PetscSFFetchAndOpEnd(PetscSF sf, MPI_Datatype unit, void *rootdata, const void *leafdata, void *leafupdate, MPI_Op op)
"""
    PetscSFFetchAndOpEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, leafupdate, op::MPI.Op)

End a fetch-and-op started with `PetscSFFetchAndOpBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)`); `rootdata` has the SF's root count,
`leafdata` and `leafupdate` its leaf count. The arrays must stay alive until `PetscSFFetchAndOpEnd` returns.

See also: `PetscSF`, `PetscSFFetchAndOpBegin()`, `PetscSFReduceBegin()`

# External Links
$(_doc_external("Vec/PetscSFFetchAndOpEnd"))
"""
function PetscSFFetchAndOpEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFFetchAndOpEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata leafupdate begin
        @chk ccall(
            (:PetscSFFetchAndOpEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata),
            leafupdate isa Ptr ? leafupdate : pointer(leafupdate), op)
    end
    return nothing
end
