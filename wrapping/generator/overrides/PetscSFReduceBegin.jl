# override for PetscSFReduceBegin; C signature: PetscSFReduceBegin(PetscSF sf, MPI_Datatype unit, const void *leafdata, void *rootdata, MPI_Op op)
"""
    PetscSFReduceBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata, rootdata, op::MPI.Op)

Begin reducing leaf values into their roots: `rootdata[root] = op(rootdata[root], leafdata[leaf])` over all leaves of a root.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`leafdata` and `rootdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.
Typical `op`s are `MPI.SUM`, `MPI.MAX`, `MPI.MIN` and `MPI.REPLACE`. Both arrays must stay alive until `PetscSFReduceEnd` returns.

See also: `PetscSF`, `PetscSFReduceEnd()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFReduceBegin"))
"""
function PetscSFReduceBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFReduceBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve leafdata rootdata begin
        @chk ccall(
            (:PetscSFReduceBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, leafdata isa Ptr ? leafdata : pointer(leafdata), rootdata isa Ptr ? rootdata : pointer(rootdata), op)
    end
    return nothing
end
