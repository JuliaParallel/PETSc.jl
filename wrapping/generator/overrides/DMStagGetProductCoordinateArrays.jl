# override for DMStagGetProductCoordinateArrays; C signature: DMStagGetProductCoordinateArrays(DM dm, void* arrX, void* arrY, void* arrZ)
"""
	arrX,arrY,arrZ = DMStagGetProductCoordinateArrays(petsclib::PetscLibType,dm::AbstractPetscDM) 
extract local product coordinate arrays, one per dimension

Logically Collective

Input Parameter:
- `dm` - the `DMSTAG` object

Output Parameters:
- `arrX` - local 1D coordinate arrays for x direction
- `arrY` - local 1D coordinate arrays for y direction, not set for one dimensional problems
- `arrZ` - local 1D coordinate arrays for z direction, not set for one and two dimensional problems

Level: intermediate

Notes:
A high-level helper function to quickly extract local coordinate arrays.

Note that 2-dimensional arrays are returned. See
`DMStagVecGetArray()`, which is called internally to produce these arrays
representing coordinates on elements and vertices (element boundaries)
for a 1-dimensional `DMSTAG` in each coordinate direction.

One should use `DMStagGetProductCoordinateLocationSlot()` to determine appropriate
indices for the second dimension in these returned arrays. This function
checks that the coordinate array is a suitable product of 1-dimensional
`DMSTAG` objects.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagGetProductCoordinateArraysRead()`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesProduct()`, `DMStagGetProductCoordinateLocationSlot()`

# External Links
$(_doc_external("DMStag/DMStagGetProductCoordinateArrays"))
"""
function DMStagGetProductCoordinateArrays(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMStagGetProductCoordinateArrays(petsclib::$UnionPetscLib, dm::AbstractPetscDM)
    arrX_ = Ref{Ptr{Ptr{$PetscScalar}}}()
    arrY_ = Ref{Ptr{Ptr{$PetscScalar}}}()
    arrZ_ = Ref{Ptr{Ptr{$PetscScalar}}}()

    xs,ys,zs,nx,ny,nz = DMStagGetGhostCorners(petsclib, dm)

    @chk ccall(
               (:DMStagGetProductCoordinateArrays, $petsc_library),
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

