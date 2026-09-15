# override for DMStagRestoreProductCoordinateArraysRead; C signature: DMStagRestoreProductCoordinateArraysRead(DM dm, void* arrX, void* arrY, void* arrZ)
"""
	DMStagRestoreProductCoordinateArraysRead(petsclib::PetscLibType,dm::AbstractPetscDM, arrX::PetscArray, arrY::PetscArray, arrZ::PetscArray) 
restore local product array access, read

Logically Collective

Input Parameters:
- `dm`   - the `DMSTAG` object
- `arrX` - local 1D coordinate arrays for x direction
- `arrY` - local 1D coordinate arrays for y direction
- `arrZ` - local 1D coordinate arrays for z direction

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMStagGetProductCoordinateArrays()`, `DMStagGetProductCoordinateArraysRead()`

# External Links
$(_doc_external("DMStag/DMStagRestoreProductCoordinateArraysRead"))
"""
function DMStagRestoreProductCoordinateArraysRead(petsclib::PetscLibType, dm::AbstractPetscDM, arrX::Union{PetscArray, Nothing}, arrY::Union{PetscArray, Nothing}, arrZ::Union{PetscArray, Nothing}) end

@for_petsc function DMStagRestoreProductCoordinateArraysRead(petsclib::$UnionPetscLib, dm::AbstractPetscDM, arrX::Union{PetscArray, Nothing}, arrY::Union{PetscArray, Nothing}, arrZ::Union{PetscArray, Nothing} )

    refX = arrX === nothing ? C_NULL : Ref(arrX.ptr)
    refY = arrY === nothing ? C_NULL : Ref(arrY.ptr)
    refZ = arrZ === nothing ? C_NULL : Ref(arrZ.ptr)

    @chk ccall(
               (:DMStagRestoreProductCoordinateArrays, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{Ptr{Ptr{$PetscScalar}}}, Ptr{Ptr{Ptr{$PetscScalar}}}, Ptr{Ptr{Ptr{$PetscScalar}}}),
               dm, refX, refY, refZ,
              )

	return nothing
end

