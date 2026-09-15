# override for DMStagGetProductCoordinateArraysRead; C signature: DMStagGetProductCoordinateArraysRead(DM dm, void* arrX, void* arrY, void* arrZ)
"""
	arrX::PetscArray, arrY::PetscArray, arrZ::PetscArray = DMStagGetProductCoordinateArraysRead(petsclib::PetscLibType,dm::AbstractPetscDM) 
extract product coordinate arrays, read

Logically Collective

Input Parameter:
- `dm` - the `DMSTAG` object

Output Parameters:
- `arrX` - local 1D coordinate arrays for `x` direction
- `arrY` - local 1D coordinate arrays for `y` direction, not set for one dimensional problems
- `arrZ` - local 1D coordinate arrays for `z` direction, not set for one and two dimensional problems

Level: intermediate

Note:
See `DMStagGetProductCoordinateArrays()` for more information.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagGetProductCoordinateArrays()`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesProduct()`, `DMStagGetProductCoordinateLocationSlot()`

# External Links
$(_doc_external("DMStag/DMStagGetProductCoordinateArraysRead"))
"""
function DMStagGetProductCoordinateArraysRead(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMStagGetProductCoordinateArraysRead(petsclib::$UnionPetscLib, dm::AbstractPetscDM)
    arrX_ = Ref{Ptr{Ptr{$PetscScalar}}}()
    arrY_ = Ref{Ptr{Ptr{$PetscScalar}}}()
    arrZ_ = Ref{Ptr{Ptr{$PetscScalar}}}()

    xs,ys,zs,nx,ny,nz = DMStagGetGhostCorners(petsclib, dm)

    @chk ccall(
               (:DMStagGetProductCoordinateArraysRead, $petsc_library),
               PetscErrorCode,
               (CDM, Ref{Ptr{Ptr{$PetscScalar}}}, Ref{Ptr{Ptr{$PetscScalar}}}, Ref{Ptr{Ptr{$PetscScalar}}}),
               dm, arrX_, arrY_, arrZ_,
              )

    mat = unsafe_wrap(Array, unsafe_load(arrX_[], xs+1), (2,nx))
    mat = OffsetArray(PermutedDimsArray(mat, (2,1)), xs, 0)
    arrX = PetscArray(mat,arrX_[]) 

    if ny>0
        mat = unsafe_wrap(Array, unsafe_load(arrY_[], ys+1), (2,ny))
        mat = OffsetArray(PermutedDimsArray(mat, (2,1)), ys, 0)
        arrY = PetscArray(mat,arrY_[]) 
    else
        arrY = nothing
    end
    if nz>0
        mat = unsafe_wrap(Array, unsafe_load(arrZ_[], zs+1), (2,nz))
        mat = OffsetArray(PermutedDimsArray(mat, (2,1)), zs, 0)
        arrZ = PetscArray(mat,arrZ_[]) 
    else
        arrZ = nothing
    end

	return arrX,arrY,arrZ
end

