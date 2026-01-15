# autodefined type arguments for class ------
mutable struct _n_PetscKDTree end
const PetscKDTree = Ptr{_n_PetscKDTree}

# -------------------------------------------------------
"""
	PetscKDTreeDestroy(petsclib::PetscLibType,tree::PetscKDTree) 
destroy a `PetscKDTree`

Not Collective, No Fortran Support

Input Parameters:
- `tree` - tree to destroy

Level: advanced

-seealso: `PetscKDTree`, `PetscKDTreeCreate()`

# External Links
$(_doc_external("Vec/PetscKDTreeDestroy"))
"""
function PetscKDTreeDestroy(petsclib::PetscLibType, tree::PetscKDTree) end

@for_petsc function PetscKDTreeDestroy(petsclib::$UnionPetscLib, tree::PetscKDTree )

    @chk ccall(
               (:PetscKDTreeDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscKDTree},),
               tree,
              )


	return nothing
end 

"""
	new_tree::PetscKDTree = PetscKDTreeCreate(petsclib::PetscLibType,num_coords::PetscCount, dim::PetscInt, coords::Vector{PetscReal}, copy_mode::PetscCopyMode, max_bucket_size::PetscInt) 
create a `PetscKDTree`

Not Collective, No Fortran Support

Input Parameters:
- `num_coords`      - number of coordinate points to build the `PetscKDTree`
- `dim`             - the dimension of the coordinates
- `coords`          - array of the coordinates, in point-major order
- `copy_mode`       - behavior handling `coords`, `PETSC_COPY_VALUES` generally more performant
- `max_bucket_size` - maximum number of points stored at each leaf

Output Parameter:
- `new_tree` - the resulting `PetscKDTree`

Level: advanced

-seealso: `PetscKDTree`, `PetscKDTreeDestroy()`, `PetscKDTreeQueryPointsNearestNeighbor()`

# External Links
$(_doc_external("Vec/PetscKDTreeCreate"))
"""
function PetscKDTreeCreate(petsclib::PetscLibType, num_coords::PetscCount, dim::PetscInt, coords::Vector{PetscReal}, copy_mode::PetscCopyMode, max_bucket_size::PetscInt) end

@for_petsc function PetscKDTreeCreate(petsclib::$UnionPetscLib, num_coords::PetscCount, dim::$PetscInt, coords::Vector{$PetscReal}, copy_mode::PetscCopyMode, max_bucket_size::$PetscInt )
	new_tree_ = Ref{PetscKDTree}()

    @chk ccall(
               (:PetscKDTreeCreate, $petsc_library),
               PetscErrorCode,
               (PetscCount, $PetscInt, Ptr{$PetscReal}, PetscCopyMode, $PetscInt, Ptr{PetscKDTree}),
               num_coords, dim, coords, copy_mode, max_bucket_size, new_tree_,
              )

	new_tree = new_tree_[]

	return new_tree
end 

"""
	distances::Vector{PetscReal} = PetscKDTreeQueryPointsNearestNeighbor(petsclib::PetscLibType,tree::PetscKDTree, num_points::PetscCount, points::Vector{PetscReal}, tolerance::PetscReal, indices::Vector{PetscCount}) 
find the nearest neighbor in a `PetscKDTree`

Not Collective, No Fortran Support

Input Parameters:
- `tree`       - tree to query
- `num_points` - number of points to query
- `points`     - array of the coordinates, in point-major order
- `tolerance`  - tolerance for nearest neighbor

Output Parameters:
- `indices`   - indices of the nearest neighbor to the query point
- `distances` - distance between the queried point and the nearest neighbor

Level: advanced

-seealso: `PetscKDTree`, `PetscKDTreeCreate()`

# External Links
$(_doc_external("Vec/PetscKDTreeQueryPointsNearestNeighbor"))
"""
function PetscKDTreeQueryPointsNearestNeighbor(petsclib::PetscLibType, tree::PetscKDTree, num_points::PetscCount, points::Vector{PetscReal}, tolerance::PetscReal, indices::Vector{PetscCount}) end

@for_petsc function PetscKDTreeQueryPointsNearestNeighbor(petsclib::$UnionPetscLib, tree::PetscKDTree, num_points::PetscCount, points::Vector{$PetscReal}, tolerance::$PetscReal, indices::Vector{PetscCount} )
	distances = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscKDTreeQueryPointsNearestNeighbor, $petsc_library),
               PetscErrorCode,
               (PetscKDTree, PetscCount, Ptr{$PetscReal}, $PetscReal, Ptr{PetscCount}, Ptr{$PetscReal}),
               tree, num_points, points, tolerance, indices, distances,
              )


	return distances
end 

"""
	PetscKDTreeView(petsclib::PetscLibType,tree::PetscKDTree, viewer::PetscViewer) 
view a `PetscKDTree`

Not Collective, No Fortran Support

Input Parameters:
- `tree`   - tree to view
- `viewer` - visualization context

Level: advanced

-seealso: `PetscKDTree`, `PetscKDTreeCreate()`, `PetscViewer`

# External Links
$(_doc_external("Vec/PetscKDTreeView"))
"""
function PetscKDTreeView(petsclib::PetscLibType, tree::PetscKDTree, viewer::PetscViewer) end

@for_petsc function PetscKDTreeView(petsclib::$UnionPetscLib, tree::PetscKDTree, viewer::PetscViewer )

    @chk ccall(
               (:PetscKDTreeView, $petsc_library),
               PetscErrorCode,
               (PetscKDTree, PetscViewer),
               tree, viewer,
              )


	return nothing
end 

