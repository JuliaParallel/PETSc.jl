# autodefined type arguments for class ------
mutable struct _n_PetscGridHash end
const PetscGridHash = Ptr{_n_PetscGridHash}

# -------------------------------------------------------
"""
	box::PetscGridHash = PetscGridHashCreate(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, point::Vector{PetscScalar}) 

# External Links
$(_doc_external("Dm/PetscGridHashCreate"))
"""
function PetscGridHashCreate(petsclib::PetscLibType, comm::MPI_Comm, dim::PetscInt, point::Vector{PetscScalar}) end

@for_petsc function PetscGridHashCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, dim::$PetscInt, point::Vector{$PetscScalar} )
	box_ = Ref{PetscGridHash}()

    @chk ccall(
               (:PetscGridHashCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscScalar}, Ptr{PetscGridHash}),
               comm, dim, point, box_,
              )

	box = box_[]

	return box
end 

"""
	PetscGridHashEnlarge(petsclib::PetscLibType,box::PetscGridHash, point::Vector{PetscScalar}) 

# External Links
$(_doc_external("Dm/PetscGridHashEnlarge"))
"""
function PetscGridHashEnlarge(petsclib::PetscLibType, box::PetscGridHash, point::Vector{PetscScalar}) end

@for_petsc function PetscGridHashEnlarge(petsclib::$UnionPetscLib, box::PetscGridHash, point::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscGridHashEnlarge, $petsc_library),
               PetscErrorCode,
               (PetscGridHash, Ptr{$PetscScalar}),
               box, point,
              )


	return nothing
end 

"""
	PetscGridHashSetGrid(petsclib::PetscLibType,box::PetscGridHash, n::Vector{PetscInt}, h::Vector{PetscReal}) 
Divide the grid into boxes

Not Collective

Input Parameters:
- `box` - The grid hash object
- `n`   - The number of boxes in each dimension, may use `PETSC_DETERMINE` for the entries
- `h`   - The box size in each dimension, only used if n[d] == `PETSC_DETERMINE`, if not needed you can pass in `NULL`

Level: developer

-seealso: `DMPLEX`, `PetscGridHashCreate()`

# External Links
$(_doc_external("Dm/PetscGridHashSetGrid"))
"""
function PetscGridHashSetGrid(petsclib::PetscLibType, box::PetscGridHash, n::Vector{PetscInt}, h::Vector{PetscReal}) end

@for_petsc function PetscGridHashSetGrid(petsclib::$UnionPetscLib, box::PetscGridHash, n::Vector{$PetscInt}, h::Vector{$PetscReal} )

    @chk ccall(
               (:PetscGridHashSetGrid, $petsc_library),
               PetscErrorCode,
               (PetscGridHash, Ptr{$PetscInt}, Ptr{$PetscReal}),
               box, n, h,
              )


	return nothing
end 

"""
	dboxes::Vector{PetscInt},boxes::Vector{PetscInt} = PetscGridHashGetEnclosingBox(petsclib::PetscLibType,box::PetscGridHash, numPoints::PetscInt, points::Vector{PetscScalar}) 
Find the grid boxes containing each input point

Not Collective

Input Parameters:
- `box`       - The grid hash object
- `numPoints` - The number of input points
- `points`    - The input point coordinates

Output Parameters:
- `dboxes` - An array of `numPoints` x `dim` integers expressing the enclosing box as (i_0, i_1, ..., i_dim)
- `boxes`  - An array of `numPoints` integers expressing the enclosing box as single number, or `NULL`

Level: developer

-seealso: `DMPLEX`, `PetscGridHashCreate()`

# External Links
$(_doc_external("Dm/PetscGridHashGetEnclosingBox"))
"""
function PetscGridHashGetEnclosingBox(petsclib::PetscLibType, box::PetscGridHash, numPoints::PetscInt, points::Vector{PetscScalar}) end

@for_petsc function PetscGridHashGetEnclosingBox(petsclib::$UnionPetscLib, box::PetscGridHash, numPoints::$PetscInt, points::Vector{$PetscScalar} )
	dboxes = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!
	boxes = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscGridHashGetEnclosingBox, $petsc_library),
               PetscErrorCode,
               (PetscGridHash, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               box, numPoints, points, dboxes, boxes,
              )


	return dboxes,boxes
end 

"""
	PetscGridHashDestroy(petsclib::PetscLibType,box::PetscGridHash) 

# External Links
$(_doc_external("Dm/PetscGridHashDestroy"))
"""
function PetscGridHashDestroy(petsclib::PetscLibType, box::PetscGridHash) end

@for_petsc function PetscGridHashDestroy(petsclib::$UnionPetscLib, box::PetscGridHash )

    @chk ccall(
               (:PetscGridHashDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscGridHash},),
               box,
              )


	return nothing
end 

