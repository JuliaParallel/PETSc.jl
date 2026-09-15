# override for DMStagRestoreProductCoordinateArrays; C signature: DMStagRestoreProductCoordinateArrays(DM dm, void* arrX, void* arrY, void* arrZ)
"""
	DMStagRestoreProductCoordinateArrays(petsclib::PetscLibType,dm::AbstractPetscDM, arrX::PetscArray, arrY::PetscArray, arrZ::PetscArray) 
restore local array access

Logically Collective

Input Parameter:
- `dm` - the `DMSTAG` object
- `arrX` - local 1D coordinate arrays for x direction
- `arrY` - local 1D coordinate arrays for y direction
- `arrZ` - local 1D coordinate arrays for z direction

Level: intermediate

Notes:
This function does not automatically perform a local->global scatter to populate global coordinates from the local coordinates.
Thus, it may be required to explicitly perform these operations in some situations, as in the following partial example:
-vb
PetscCall(DMGetCoordinateDM(dm, &cdm));
for (PetscInt d = 0; d < 3; ++d) {
DM  subdm;
Vec coor, coor_local;

PetscCall(DMProductGetDM(cdm, d, &subdm));
PetscCall(DMGetCoordinates(subdm, &coor));
PetscCall(DMGetCoordinatesLocal(subdm, &coor_local));
PetscCall(DMLocalToGlobal(subdm, coor_local, INSERT_VALUES, coor));
PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Coordinates dim %" PetscInt_FMT ":\n", d));
PetscCall(VecView(coor, PETSC_VIEWER_STDOUT_WORLD));
}
-ve

See also: 
=== 
`DMSTAG`, `DMStagGetProductCoordinateArrays()`, `DMStagGetProductCoordinateArraysRead()`

# External Links
$(_doc_external("DMStag/DMStagRestoreProductCoordinateArrays"))
"""
function DMStagRestoreProductCoordinateArrays(petsclib::PetscLibType, dm::AbstractPetscDM, arrX::Union{PetscArray, Nothing}, arrY::Union{PetscArray, Nothing}, arrZ::Union{PetscArray, Nothing}) end

@for_petsc function DMStagRestoreProductCoordinateArrays(petsclib::$UnionPetscLib, dm::AbstractPetscDM, arrX::Union{PetscArray, Nothing}, arrY::Union{PetscArray, Nothing}, arrZ::Union{PetscArray, Nothing})

    refX = arrX === nothing ? C_NULL : Ref(arrX.ptr)
    refY = arrY === nothing ? C_NULL : Ref(arrY.ptr)
    refZ = arrZ === nothing ? C_NULL : Ref(arrZ.ptr)

    @chk ccall(
               (:DMStagRestoreProductCoordinateArrays, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{Ptr{Ptr{$PetscScalar}}}, Ptr{Ptr{Ptr{$PetscScalar}}}, Ptr{Ptr{Ptr{$PetscScalar}}}),
               dm, refX, refY, refZ
              )


	return nothing
end

