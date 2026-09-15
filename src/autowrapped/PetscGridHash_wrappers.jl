"""
	box::PetscGridHash = PetscGridHashCreate(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, point::Vector{PetscScalar}) 
Create a `PetscGridHash` for spatially locating points in a mesh.

Collective

Input Parameters:
- `comm`  - the MPI communicator
- `dim`   - the spatial dimension
- `point` - an initial point used to seed the bounding box, or `NULL` for a zero-initialized box

Output Parameter:
- `box` - the newly created `PetscGridHash`

Level: developer

-seealso: `DMPLEX`, `PetscGridHash`, `PetscGridHashEnlarge()`, `PetscGridHashDestroy()`

# External Links
$(_doc_external("DMPlex/PetscGridHashCreate"))
"""
function PetscGridHashCreate(petsclib::PetscLibType, comm::MPI_Comm, dim::Integer, point::AbstractVector{<:Number})
    error("PetscGridHashCreate: no generated method for these argument types")
end

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
	PetscGridHashDestroy(petsclib::PetscLibType,box::Union{PetscGridHash, Ref{PetscGridHash}}) 
Destroy a `PetscGridHash` and free its resources.

Collective

Input Parameter:
- `box` - the `PetscGridHash` to destroy; set to `NULL` on return

Level: developer

-seealso: `DMPLEX`, `PetscGridHash`, `PetscGridHashCreate()`, `PetscGridHashEnlarge()`

# External Links
$(_doc_external("DMPlex/PetscGridHashDestroy"))
"""
function PetscGridHashDestroy(petsclib::PetscLibType, box::Union{PetscGridHash, Ref{PetscGridHash}})
    error("PetscGridHashDestroy: no generated method for these argument types")
end

@for_petsc function PetscGridHashDestroy(petsclib::$UnionPetscLib, box::Union{PetscGridHash, Ref{PetscGridHash}} )
	box_ = box isa Base.RefValue ? box : Ref{PetscGridHash}(box)

    @chk ccall(
               (:PetscGridHashDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscGridHash},),
               box_,
              )


	return nothing
end 

"""
	PetscGridHashEnlarge(petsclib::PetscLibType,box::PetscGridHash, point::Vector{PetscScalar}) 
Enlarge the bounding box of a `PetscGridHash` to include a new point.

Not Collective

Input Parameters:
- `box`   - the `PetscGridHash`
- `point` - the point whose coordinates extend the box's lower and upper bounds

Level: developer

-seealso: `DMPLEX`, `PetscGridHash`, `PetscGridHashCreate()`, `PetscGridHashDestroy()`

# External Links
$(_doc_external("DMPlex/PetscGridHashEnlarge"))
"""
function PetscGridHashEnlarge(petsclib::PetscLibType, box::PetscGridHash, point::AbstractVector{<:Number})
    error("PetscGridHashEnlarge: no generated method for these argument types")
end

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
	PetscGridHashGetEnclosingBox(petsclib::PetscLibType,box::PetscGridHash, numPoints::PetscInt, points::Vector{PetscScalar}, dboxes::Vector{PetscInt}, boxes::Vector{PetscInt}) 
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
$(_doc_external("DMPlex/PetscGridHashGetEnclosingBox"))
"""
function PetscGridHashGetEnclosingBox(petsclib::PetscLibType, box::PetscGridHash, numPoints::Integer, points::AbstractVector{<:Number}, dboxes::AbstractVector{<:Number}, boxes::AbstractVector{<:Number})
    error("PetscGridHashGetEnclosingBox: no generated method for these argument types")
end

@for_petsc function PetscGridHashGetEnclosingBox(petsclib::$UnionPetscLib, box::PetscGridHash, numPoints::$PetscInt, points::Vector{$PetscScalar}, dboxes::Vector{$PetscInt}, boxes::Vector{$PetscInt} )

    @chk ccall(
               (:PetscGridHashGetEnclosingBox, $petsc_library),
               PetscErrorCode,
               (PetscGridHash, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               box, numPoints, points, dboxes, boxes,
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
$(_doc_external("DMPlex/PetscGridHashSetGrid"))
"""
function PetscGridHashSetGrid(petsclib::PetscLibType, box::PetscGridHash, n::AbstractVector{<:Number}, h::AbstractVector{<:Number})
    error("PetscGridHashSetGrid: no generated method for these argument types")
end

@for_petsc function PetscGridHashSetGrid(petsclib::$UnionPetscLib, box::PetscGridHash, n::Vector{$PetscInt}, h::Vector{$PetscReal} )

    @chk ccall(
               (:PetscGridHashSetGrid, $petsc_library),
               PetscErrorCode,
               (PetscGridHash, Ptr{$PetscInt}, Ptr{$PetscReal}),
               box, n, h,
              )


	return nothing
end 

