# override for PetscSFFetchAndOpBegin; C signature: PetscSFFetchAndOpBegin(PetscSF sf, MPI_Datatype unit, void *rootdata, const void *leafdata, void *leafupdate, MPI_Op op)
"""
    PetscSFFetchAndOpBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, leafupdate, op::MPI.Op)

Begin a fetch-and-op: every leaf fetches the current root value into `leafupdate` and then applies `rootdata[root] = op(rootdata[root], leafdata[leaf])`, atomically per root.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)`); `rootdata` has the SF's root count,
`leafdata` and `leafupdate` its leaf count. The arrays must stay alive until `PetscSFFetchAndOpEnd` returns.

See also: `PetscSF`, `PetscSFFetchAndOpEnd()`, `PetscSFReduceBegin()`

# External Links
$(_doc_external("Vec/PetscSFFetchAndOpBegin"))
"""
function PetscSFFetchAndOpBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFFetchAndOpBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata leafupdate begin
        @chk ccall(
            (:PetscSFFetchAndOpBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata),
            leafupdate isa Ptr ? leafupdate : pointer(leafupdate), op)
    end
    return nothing
end
