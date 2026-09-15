# override for DMPlexGetFullMeet; C signature: DMPlexGetFullMeet(DM dm, PetscInt numPoints, PetscInt points[], PetscInt* numCoveredPoints, PetscInt* coveredPoints[])
"""
	numCoveredPoints::PetscInt,coveredPoints::Vector{PetscInt} = DMPlexGetFullMeet(petsclib::PetscLibType,dm::AbstractPetscDM, numPoints::PetscInt, points::Vector{PetscInt}) 
Get an array for the meet of the set of points

Not Collective

Input Parameters:
- `dm`        - The `DMPLEX` object
- `numPoints` - The number of input points for the meet
- `points`    - The input points, of length  `numPoints`

Output Parameters:
- `numCoveredPoints` - The number of points in the meet
- `coveredPoints`    - The points in the meet, of length  `numCoveredPoints`

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexGetMeet()`, `DMPlexRestoreMeet()`, `DMPlexGetJoin()`

# External Links
$(_doc_external("DMPlex/DMPlexGetFullMeet"))
"""
function DMPlexGetFullMeet(petsclib::PetscLibType, dm::AbstractPetscDM, numPoints::PetscInt, points::Vector{PetscInt}) end

@for_petsc function DMPlexGetFullMeet(petsclib::$UnionPetscLib, dm::AbstractPetscDM, numPoints::$PetscInt, points::Vector{$PetscInt} )
	numCoveringPoints_ = Ref{$PetscInt}()
	coveringPoints_ = Ref{Ptr{$PetscInt}}(C_NULL)
    @chk ccall(
               (:DMPlexGetFullMeet, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               dm, numPoints, points, numCoveringPoints_, coveringPoints_,
              )
	numCoveringPoints = numCoveringPoints_[]
	coveringPoints = unsafe_wrap(Array, coveringPoints_[], numCoveringPoints; own = false)
	return numCoveringPoints,coveringPoints
end

