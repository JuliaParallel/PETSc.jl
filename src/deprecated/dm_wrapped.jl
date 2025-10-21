"""
	 UNTESTED !!!
	newdm = DMClone(dm::AbstractDM{PetscLib})

Creates a `DM` object with the same topology as the original.

Collective

Input Parameter:
===
- `dm` - The original `DM` object

Output Parameter:
===
- `newdm` - The new `DM` object

Level: beginner

Notes:
For some `DM` implementations this is a shallow clone, the result of which may share (reference counted) information with its parent. For example,
`DMClone()` applied to a `DMPLEX` object will result in a new `DMPLEX` that shares the topology with the original `DMPLEX`. It does not
share the `PetscSection` of the original `DM`.

The clone is considered set up if the original has been set up.

Use `DMConvert()` for a general way to create new `DM` from a given `DM`

See also: 
=== 
`DM`, `DMDestroy()`, `DMCreate()`, `DMSetType()`, `DMSetLocalSection()`, `DMSetGlobalSection()`, `DMPLEX`, `DMConvert()`

# External Links
$(_doc_external("DM/DMClone"))
"""
function DMClone(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	newdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMClone(
		PetscLib,
		dm,
		newdm,
	)

	return newdm
end
 
 
"""
	 UNTESTED !!!
	 DMSetType(dm::AbstractDM{PetscLib},method::DMType)

Builds a `DM`, for a particular `DM` implementation.

Collective

Input Parameters:
===
- `dm`     - The `DM` object
- `method` - The name of the `DMType`, for example `DMDA`, `DMPLEX`

Options Database Key:
===
- `-dm_type <type>` - Sets the `DM` type; use -help for a list of available types

Level: intermediate

Note:
Of the `DM` is constructed by directly calling a function to construct a particular `DM`, for example, `DMDACreate2d()` or `DMPlexCreateBoxMesh()`

See also: 
=== 
`DM`, `DMType`, `DMDA`, `DMPLEX`, `DMGetType()`, `DMCreate()`, `DMDACreate2d()`

# External Links
$(_doc_external("DM/DMSetType"))
"""
function DMSetType(dm::AbstractDM{PetscLib},method::DMType) where {PetscLib}

	LibPETSc.DMSetType(
		PetscLib,
		dm,
		method,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	type = DMGetType(dm::AbstractDM{PetscLib})

Gets the `DM` type name (as a string) from the `DM`.

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `type` - The `DMType` name

Level: intermediate

See also: 
=== 
`DM`, `DMType`, `DMDA`, `DMPLEX`, `DMSetType()`, `DMCreate()`

# External Links
$(_doc_external("DM/DMGetType"))
"""
function DMGetType(dm::AbstractDM{PetscLib}) where {PetscLib}
	r_type = Ref{PETSc.DMType}()

	LibPETSc.DMGetType(
		PetscLib,
		dm,
		r_type,
	)


	type = unsafe_string(r_type[])
	return type
end
 
 
"""
	 UNTESTED !!!
	 DMView(dm::AbstractDM{PetscLib},v::PetscViewer)

Views a `DM`. Depending on the `PetscViewer` and its `PetscViewerFormat` it may print some ASCII information about the `DM` to the screen or a file or
save the `DM` in a binary file to be loaded later or create a visualization of the `DM`

Collective

Input Parameters:
===
- `dm` - the `DM` object to view
- `v`  - the viewer

Level: beginner

Notes:

`PetscViewer` = `PETSCVIEWERHDF5` i.e. HDF5 format can be used with `PETSC_VIEWER_HDF5_PETSC` as the `PetscViewerFormat` to save multiple `DMPLEX`
meshes in a single HDF5 file. This in turn requires one to name the `DMPLEX` object with `PetscObjectSetName()`
before saving it with `DMView()` and before loading it with `DMLoad()` for identification of the mesh object.

`PetscViewer` = `PETSCVIEWEREXODUSII` i.e. ExodusII format assumes that element blocks (mapped to "Cell sets" labels)
consists of sequentially numbered cells.

If `dm` has been distributed, only the part of the `DM` on MPI rank 0 (including "ghost" cells and vertices) will be written.

Only TRI, TET, QUAD, and HEX cells are supported.

`DMPLEX` only represents geometry while most post-processing software expect that a mesh also provides information on the discretization space. This function assumes that the file represents Lagrange finite elements of order 1 or 2.
The order of the mesh shall be set using `PetscViewerExodusIISetOrder()`

Variable names can be set and querried using `PetscViewerExodusII[Set/Get][Nodal/Zonal]VariableNames[s]`.

See also: 
=== 
`DM`, `PetscViewer`, `PetscViewerFormat`, `PetscViewerSetFormat()`, `DMDestroy()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMLoad()`, `PetscObjectSetName()`

# External Links
$(_doc_external("DM/DMView"))
"""
function DMView(dm::AbstractDM{PetscLib},v::PetscViewer) where {PetscLib}

	LibPETSc.DMView(
		PetscLib,
		dm,
		v,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLoad(newdm::AbstractDM{PetscLib},viewer::PetscViewer)

Loads a DM that has been stored in binary  with `DMView()`.

Collective

Input Parameters:
===
- `newdm`  - the newly loaded `DM`, this needs to have been created with `DMCreate()` or
some related function before a call to `DMLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()` or
`PETSCVIEWERHDF5` file viewer, obtained from `PetscViewerHDF5Open()`

Level: intermediate

Notes:
The type is determined by the data in the file, any type set into the DM before this call is ignored.

Using `PETSCVIEWERHDF5` type with `PETSC_VIEWER_HDF5_PETSC` format, one can save multiple `DMPLEX`
meshes in a single HDF5 file. This in turn requires one to name the `DMPLEX` object with `PetscObjectSetName()`
before saving it with `DMView()` and before loading it with `DMLoad()` for identification of the mesh object.

See also: 
=== 
`DM`, `PetscViewerBinaryOpen()`, `DMView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("DM/DMLoad"))
"""
function DMLoad(newdm::AbstractDM{PetscLib},viewer::PetscViewer) where {PetscLib}

	LibPETSc.DMLoad(
		PetscLib,
		newdm,
		viewer,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dm = DMDestroy()

Destroys a `DM`.

Collective

Input Parameter:
===
- `dm` - the `DM` object to destroy

Level: developer

See also: 
=== 
`DM`, `DMCreate()`, `DMType`, `DMSetType()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`

# External Links
$(_doc_external("DM/DMDestroy"))
"""
function DMDestroy(dm::AbstractDM{PetscLib}) where {PetscLib}
	LibPETSc.DMDestroy(
		PetscLib,
		dm,
	)

	return dm
end
 
 
"""
	 UNTESTED !!!
	vec = DMCreateGlobalVector(dm::AbstractDM{PetscLib})

Creates a global vector from a `DM` object. A global vector is a parallel vector that has no duplicate values shared between MPI ranks,
that is it has no ghost locations.

Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `vec` - the global vector

Level: beginner

See also: 
=== 
`DM`, `Vec`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMDestroy()`, `DMView()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`,
`DMGlobalToLocalBegin()`, `DMGlobalToLocalEnd()`

# External Links
$(_doc_external("DM/DMCreateGlobalVector"))
"""
function DMCreateGlobalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	vec = CVec()

	LibPETSc.DMCreateGlobalVector(
		PetscLib,
		dm,
		vec,
	)

	return vec
end
 
 
"""
	 UNTESTED !!!
	vec = DMCreateLocalVector(dm::AbstractDM{PetscLib})

Creates a local vector from a `DM` object.

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `vec` - the local vector

Level: beginner

Note:
A local vector usually has ghost locations that contain values that are owned by different MPI ranks. A global vector has no ghost locations.

See also: 
=== 
`DM`, `Vec`, `DMCreateGlobalVector()`, `DMGetLocalVector()`, `DMDestroy()`, `DMView()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`
`DMGlobalToLocalBegin()`, `DMGlobalToLocalEnd()`

# External Links
$(_doc_external("DM/DMCreateLocalVector"))
"""
function DMCreateLocalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	vec = CVec()

	LibPETSc.DMCreateLocalVector(
		PetscLib,
		dm,
		vec,
	)

	return vec
end
 
 
"""
	 UNTESTED !!!
	g = DMGetLocalVector(dm::AbstractDM{PetscLib})

Gets a PETSc vector that may be used with the `DM` local routines. This vector has spaces for the ghost values.

Not Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `g` - the local vector

Level: beginner

Note:
The vector values are NOT initialized and may have garbage in them, so you may need
to zero them.

The output parameter, `g`, is a regular PETSc vector that should be returned with
`DMRestoreLocalVector()` DO NOT call `VecDestroy()` on it.

This is intended to be used for vectors you need for a short time, like within a single function call.
For vectors that you intend to keep around (for example in a C struct) or pass around large parts of your
code you should use `DMCreateLocalVector()`.

VecStride*() operations can be useful when using `DM` with dof > 1

-seealso: `DM`, `DMCreateGlobalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMGlobalToLocalBegin()`,
`DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`, `DMCreateLocalVector()`, `DMRestoreLocalVector()`,
`VecStrideMax()`, `VecStrideMin()`, `VecStrideNorm()`, `DMClearLocalVectors()`, `DMGetNamedGlobalVector()`, `DMGetNamedLocalVector()`

# External Links
$(_doc_external("DM/DMGetLocalVector"))
"""
function DMGetLocalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	g = CVec()

	LibPETSc.DMGetLocalVector(
		PetscLib,
		dm,
		g,
	)

	return g
end
 
 
"""
	 UNTESTED !!!
	g = DMRestoreLocalVector(dm::AbstractDM{PetscLib})

Returns a PETSc vector that was
obtained from `DMGetLocalVector()`. Do not use with vector obtained via
`DMCreateLocalVector()`.

Not Collective

Input Parameters:
===
- `dm` - the `DM`
- `g`  - the local vector

Level: beginner

-seealso: `DM`, `DMCreateGlobalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMGlobalToLocalBegin()`,
`DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`, `DMCreateLocalVector()`, `DMGetLocalVector()`, `DMClearLocalVectors()`

# External Links
$(_doc_external("DM/DMRestoreLocalVector"))
"""
function DMRestoreLocalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	g = CVec()

	LibPETSc.DMRestoreLocalVector(
		PetscLib,
		dm,
		g,
	)

	return g
end
 
 
"""
	 UNTESTED !!!
	g = DMGetGlobalVector(dm::AbstractDM{PetscLib})

Gets a PETSc vector that may be used with the `DM` global routines.

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `g` - the global vector

Level: beginner

Note:
The vector values are NOT initialized and may have garbage in them, so you may need
to zero them.

The output parameter, `g`, is a regular PETSc vector that should be returned with
`DMRestoreGlobalVector()` DO NOT call `VecDestroy()` on it.

This is intended to be used for vectors you need for a short time, like within a single function call.
For vectors that you intend to keep around (for example in a C struct) or pass around large parts of your
code you should use `DMCreateGlobalVector()`.

VecStride*() operations can be useful when using `DM` with dof > 1

-seealso: `DM`, `DMCreateGlobalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMGlobalToLocalBegin()`,
`DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`, `DMCreateLocalVector()`, `DMRestoreLocalVector()`
`VecStrideMax()`, `VecStrideMin()`, `VecStrideNorm()`, `DMClearGlobalVectors()`, `DMGetNamedGlobalVector()`, `DMGetNamedLocalVector()`

# External Links
$(_doc_external("DM/DMGetGlobalVector"))
"""
function DMGetGlobalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	g = CVec()

	LibPETSc.DMGetGlobalVector(
		PetscLib,
		dm,
		g,
	)

	return g
end
 
 
"""
	 UNTESTED !!!
	g = DMRestoreGlobalVector(dm::AbstractDM{PetscLib})

Returns a PETSc vector that
obtained from `DMGetGlobalVector()`. Do not use with vector obtained via
`DMCreateGlobalVector()`.

Not Collective

Input Parameters:
===
- `dm` - the `DM`
- `g`  - the global vector

Level: beginner

-seealso: `DM`, `DMCreateGlobalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMGlobalToGlobalBegin()`,
`DMGlobalToGlobalEnd()`, `DMGlobalToGlobal()`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMClearGlobalVectors()`

# External Links
$(_doc_external("DM/DMRestoreGlobalVector"))
"""
function DMRestoreGlobalVector(dm::AbstractDM{PetscLib}) where {PetscLib}
	g = CVec()

	LibPETSc.DMRestoreGlobalVector(
		PetscLib,
		dm,
		g,
	)

	return g
end
 
 
"""
	 UNTESTED !!!
	 DMClearGlobalVectors(dm::AbstractDM{PetscLib})

Destroys all the global vectors that have been created for `DMGetGlobalVector()` calls in this `DM`

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: developer

-seealso: `DM`, `DMCreateGlobalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMGlobalToLocalBegin()`,
`DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`, `DMCreateLocalVector()`, `DMRestoreLocalVector()`
`VecStrideMax()`, `VecStrideMin()`, `VecStrideNorm()`, `DMClearLocalVectors()`

# External Links
$(_doc_external("DM/DMClearGlobalVectors"))
"""
function DMClearGlobalVectors(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearGlobalVectors(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearLocalVectors(dm::AbstractDM{PetscLib})

Destroys all the local vectors that have been created for `DMGetLocalVector()` calls in this `DM`

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: developer

-seealso: `DM`, `DMCreateLocalVector()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMLocalToLocalBegin()`,
`DMLocalToLocalEnd()`, `DMRestoreLocalVector()`
`VecStrideMax()`, `VecStrideMin()`, `VecStrideNorm()`, `DMClearGlobalVectors()`

# External Links
$(_doc_external("DM/DMClearLocalVectors"))
"""
function DMClearLocalVectors(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearLocalVectors(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearNamedGlobalVectors(dm::AbstractDM{PetscLib})

Destroys all the named global vectors that have been created with `DMGetNamedGlobalVector()` in this `DM`

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: developer

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMGetNamedLocalVector()`, `DMClearNamedLocalVectors()`

# External Links
$(_doc_external("DM/DMClearNamedGlobalVectors"))
"""
function DMClearNamedGlobalVectors(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearNamedGlobalVectors(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearNamedLocalVectors(dm::AbstractDM{PetscLib})

Destroys all the named local vectors that have been created with `DMGetNamedLocalVector()` in this `DM`

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: developer

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMGetNamedLocalVector()`, `DMClearNamedGlobalVectors()`

# External Links
$(_doc_external("DM/DMClearNamedLocalVectors"))
"""
function DMClearNamedLocalVectors(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearNamedLocalVectors(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	exists = DMHasNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

check for a named, persistent global vector created with `DMGetNamedGlobalVector()`

Not Collective

Input Parameters:
===
- `dm`   - `DM` to hold named vectors
- `name` - unique name for `Vec`

Output Parameter:
===
- `exists` - true if the vector was previously created

Level: developer

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMRestoreNamedLocalVector()`, `DMClearNamedGlobalVectors()`

# External Links
$(_doc_external("DM/DMHasNamedGlobalVector"))
"""
function DMHasNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	exists = Ref{PetscBool}()

	LibPETSc.DMHasNamedGlobalVector(
		PetscLib,
		dm,
		name,
		exists,
	)

	return exists[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	X = DMGetNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

get access to a named, persistent global vector

Collective

Input Parameters:
===
- `dm`   - `DM` to hold named vectors
- `name` - unique name for `X`

Output Parameter:
===
- `X` - named `Vec`

Level: developer

Note:
If a `Vec` with the given name does not exist, it is created.

-seealso: `DM`, `DMRestoreNamedGlobalVector()`, `DMHasNamedGlobalVector()`, `DMClearNamedGlobalVectors()`, `DMGetGlobalVector()`, `DMGetLocalVector()`

# External Links
$(_doc_external("DM/DMGetNamedGlobalVector"))
"""
function DMGetNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	X = CVec()

	LibPETSc.DMGetNamedGlobalVector(
		PetscLib,
		dm,
		name,
		X,
	)

	return X
end
 
 
"""
	 UNTESTED !!!
	X = DMRestoreNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

restore access to a named, persistent global vector

Collective

Input Parameters:
===
- `dm`   - `DM` on which `X` was gotten
- `name` - name under which `X` was gotten
- `X`    - `Vec` to restore

Level: developer

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMClearNamedGlobalVectors()`

# External Links
$(_doc_external("DM/DMRestoreNamedGlobalVector"))
"""
function DMRestoreNamedGlobalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	X = CVec()

	LibPETSc.DMRestoreNamedGlobalVector(
		PetscLib,
		dm,
		name,
		X,
	)

	return X
end
 
 
"""
	 UNTESTED !!!
	exists = DMHasNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

check for a named, persistent local vector created with `DMGetNamedLocalVector()`

Not Collective

Input Parameters:
===
- `dm`   - `DM` to hold named vectors
- `name` - unique name for `Vec`

Output Parameter:
===
- `exists` - true if the vector was previously created

Level: developer

Note:
If a `Vec` with the given name does not exist, it is created.

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMRestoreNamedLocalVector()`, `DMClearNamedLocalVectors()`

# External Links
$(_doc_external("DM/DMHasNamedLocalVector"))
"""
function DMHasNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	exists = Ref{PetscBool}()

	LibPETSc.DMHasNamedLocalVector(
		PetscLib,
		dm,
		name,
		exists,
	)

	return exists[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	X = DMGetNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

get access to a named, persistent local vector

Not Collective

Input Parameters:
===
- `dm`   - `DM` to hold named vectors
- `name` - unique name for `X`

Output Parameter:
===
- `X` - named `Vec`

Level: developer

Note:
If a `Vec` with the given name does not exist, it is created.

-seealso: `DM`, `DMGetNamedGlobalVector()`, `DMRestoreNamedLocalVector()`, `DMHasNamedLocalVector()`, `DMClearNamedLocalVectors()`, `DMGetGlobalVector()`, `DMGetLocalVector()`

# External Links
$(_doc_external("DM/DMGetNamedLocalVector"))
"""
function DMGetNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	X = CVec()

	LibPETSc.DMGetNamedLocalVector(
		PetscLib,
		dm,
		name,
		X,
	)

	return X
end
 
 
"""
	 UNTESTED !!!
	X = DMRestoreNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char})

restore access to a named, persistent local vector obtained with `DMGetNamedLocalVector()`

Not Collective

Input Parameters:
===
- `dm`   - `DM` on which `X` was gotten
- `name` - name under which `X` was gotten
- `X`    - `Vec` to restore

Level: developer

-seealso: `DM`, `DMRestoreNamedGlobalVector()`, `DMGetNamedLocalVector()`, `DMClearNamedLocalVectors()`

# External Links
$(_doc_external("DM/DMRestoreNamedLocalVector"))
"""
function DMRestoreNamedLocalVector(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	X = CVec()

	LibPETSc.DMRestoreNamedLocalVector(
		PetscLib,
		dm,
		name,
		X,
	)

	return X
end
 
 
"""
	 UNTESTED !!!
	ltog = DMGetLocalToGlobalMapping(dm::AbstractDM{PetscLib})

Accesses the local

Collective

Input Parameter:
===
- `dm` - the `DM` that provides the mapping

Output Parameter:
===
- `ltog` - the mapping

Level: advanced

Notes:
The global to local mapping allows one to set values into the global vector or matrix using `VecSetValuesLocal()` and `MatSetValuesLocal()`

Vectors obtained with  `DMCreateGlobalVector()` and matrices obtained with `DMCreateMatrix()` already contain the global mapping so you do
need to use this function with those objects.

This mapping can then be used by `VecSetLocalToGlobalMapping()` or `MatSetLocalToGlobalMapping()`.

See also: 
=== 
`DM`, `DMCreateLocalVector()`, `DMCreateGlobalVector()`, `VecSetLocalToGlobalMapping()`, `MatSetLocalToGlobalMapping()`,
`DMCreateMatrix()`

# External Links
$(_doc_external("DM/DMGetLocalToGlobalMapping"))
"""
function DMGetLocalToGlobalMapping(dm::AbstractDM{PetscLib}) where {PetscLib}
	ltog = LibPETSc.ISLocalToGlobalMapping()

	LibPETSc.DMGetLocalToGlobalMapping(
		PetscLib,
		dm,
		ltog,
	)

	return ltog
end
 
 
"""
	 UNTESTED !!!
	bs = DMGetBlockSize(dm::AbstractDM{PetscLib})

Gets the inherent block size associated with a `DM`

Not Collective

Input Parameter:
===
- `dm` - the `DM` with block structure

Output Parameter:
===
- `bs` - the block size, 1 implies no exploitable block structure

Level: intermediate

Notes:
This might be the number of degrees of freedom at each grid point for a structured grid.

Complex `DM` that represent multiphysics or staggered grids or mixed-methods do not generally have a single inherent block size, but
rather different locations in the vectors may have a different block size.

See also: 
=== 
`DM`, `ISCreateBlock()`, `VecSetBlockSize()`, `MatSetBlockSize()`, `DMGetLocalToGlobalMapping()`

# External Links
$(_doc_external("DM/DMGetBlockSize"))
"""
function DMGetBlockSize(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	bs = [PetscInt(1)]

	LibPETSc.DMGetBlockSize(
		PetscLib,
		dm,
		Ref(bs,1),
	)

	return bs[1]
end
 
 
"""
	 UNTESTED !!!
	coloring = DMCreateColoring(dm::AbstractDM{PetscLib},ctype::ISColoringType)

Gets coloring of a graph associated with the `DM`. Often the graph represents the operator matrix associated with the discretization
of a PDE on the `DM`.

Collective

Input Parameters:
===
- `dm`    - the `DM` object
- `ctype` - `IS_COLORING_LOCAL` or `IS_COLORING_GLOBAL`

Output Parameter:
===
- `coloring` - the coloring

Level: developer

Notes:
Coloring of matrices can also be computed directly from the sparse matrix nonzero structure via the `MatColoring` object or from the mesh from which the
matrix comes from (what this function provides). In general using the mesh produces a more optimal coloring (fewer colors).

This produces a coloring with the distance of 2, see `MatSetColoringDistance()` which can be used for efficiently computing Jacobians with `MatFDColoringCreate()`
For `DMDA` in three dimensions with periodic boundary conditions the number of grid points in each dimension must be divisible by 2*stencil_width + 1,
otherwise an error will be generated.

See also: 
=== 
`DM`, `ISColoring`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatType()`, `MatColoring`, `MatFDColoringCreate()`

# External Links
$(_doc_external("DM/DMCreateColoring"))
"""
function DMCreateColoring(dm::AbstractDM{PetscLib},ctype::ISColoringType) where {PetscLib}
	coloring = LibPETSc.ISColoring()

	LibPETSc.DMCreateColoring(
		PetscLib,
		dm,
		ctype,
		coloring,
	)

	return coloring
end
 
 
"""
	 UNTESTED !!!
	 DMSetMatrixPreallocateSkip(dm::AbstractDM{PetscLib},skip::PetscBool)

When `DMCreateMatrix()` is called the matrix sizes and
`ISLocalToGlobalMapping` will be properly set, but the data structures to store values in the
matrices will not be preallocated.

Logically Collective

Input Parameters:
===
- `dm`   - the `DM`
- `skip` - `PETSC_TRUE` to skip preallocation

Level: developer

Note:
This is most useful to reduce initialization costs when `MatSetPreallocationCOO()` and
`MatSetValuesCOO()` will be used.

See also: 
=== 
`DM`, `DMCreateMatrix()`, `DMSetMatrixStructureOnly()`, `DMSetMatrixPreallocateOnly()`

# External Links
$(_doc_external("DM/DMSetMatrixPreallocateSkip"))
"""
function DMSetMatrixPreallocateSkip(dm::AbstractDM{PetscLib},skip::PetscBool) where {PetscLib}

	LibPETSc.DMSetMatrixPreallocateSkip(
		PetscLib,
		dm,
		skip,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetMatrixPreallocateOnly(dm::AbstractDM{PetscLib},only::PetscBool)

When `DMCreateMatrix()` is called the matrix will be properly
preallocated but the nonzero structure and zero values will not be set.

Logically Collective

Input Parameters:
===
- `dm`   - the `DM`
- `only` - `PETSC_TRUE` if only want preallocation

Options Database Key:
===
- `-dm_preallocate_only` - Only preallocate the matrix for `DMCreateMatrix()`, `DMCreateMassMatrix()`, but do not fill it with zeros

Level: developer

See also: 
=== 
`DM`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatrixStructureOnly()`, `DMSetMatrixPreallocateSkip()`

# External Links
$(_doc_external("DM/DMSetMatrixPreallocateOnly"))
"""
function DMSetMatrixPreallocateOnly(dm::AbstractDM{PetscLib},only::PetscBool) where {PetscLib}

	LibPETSc.DMSetMatrixPreallocateOnly(
		PetscLib,
		dm,
		only,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetMatrixStructureOnly(dm::AbstractDM{PetscLib},only::PetscBool)

When `DMCreateMatrix()` is called, the matrix structure will be created
but the array for numerical values will not be allocated.

Logically Collective

Input Parameters:
===
- `dm`   - the `DM`
- `only` - `PETSC_TRUE` if you only want matrix structure

Level: developer

See also: 
=== 
`DM`, `DMCreateMatrix()`, `DMSetMatrixPreallocateOnly()`, `DMSetMatrixPreallocateSkip()`

# External Links
$(_doc_external("DM/DMSetMatrixStructureOnly"))
"""
function DMSetMatrixStructureOnly(dm::AbstractDM{PetscLib},only::PetscBool) where {PetscLib}

	LibPETSc.DMSetMatrixStructureOnly(
		PetscLib,
		dm,
		only,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetBlockingType(dm::AbstractDM{PetscLib},btype::DMBlockingType)

set the blocking granularity to be used for variable block size `DMCreateMatrix()` is called

Logically Collective

Input Parameters:
===
- `dm`    - the `DM`
- `btype` - block by topological point or field node

Options Database Key:
===
- `-dm_blocking_type [topological_point, field_node]` - use topological point blocking or field node blocking

Level: advanced

See also: 
=== 
`DM`, `DMCreateMatrix()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("DM/DMSetBlockingType"))
"""
function DMSetBlockingType(dm::AbstractDM{PetscLib},btype::DMBlockingType) where {PetscLib}

	LibPETSc.DMSetBlockingType(
		PetscLib,
		dm,
		btype,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	btype = DMGetBlockingType(dm::AbstractDM{PetscLib})

get the blocking granularity to be used for variable block size `DMCreateMatrix()` is called

Not Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `btype` - block by topological point or field node

Level: advanced

See also: 
=== 
`DM`, `DMCreateMatrix()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("DM/DMGetBlockingType"))
"""
function DMGetBlockingType(dm::AbstractDM{PetscLib}) where {PetscLib}
	btype = LibPETSc.DMBlockingType()

	LibPETSc.DMGetBlockingType(
		PetscLib,
		dm,
		btype,
	)

	return btype
end
 
 
"""
	 UNTESTED !!!
	mat,vec = DMCreateInterpolation(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib})

Gets the interpolation matrix between two `DM` objects. The resulting matrix map degrees of freedom in the vector obtained by
`DMCreateGlobalVector()` on the coarse `DM` to similar vectors on the fine grid `DM`.

Collective

Input Parameters:
===
- `dmc` - the `DM` object
- `dmf` - the second, finer `DM` object

Output Parameters:
===
- `mat` - the interpolation
- `vec` - the scaling (optional, pass `NULL` if not needed), see `DMCreateInterpolationScale()`

Level: developer

Notes:
For `DMDA` objects this only works for "uniform refinement", that is the refined mesh was obtained `DMRefine()` or the coarse mesh was obtained by
DMCoarsen(). The coordinates set into the `DMDA` are completely ignored in computing the interpolation.

For `DMDA` objects you can use this interpolation (more precisely the interpolation from the `DMGetCoordinateDM()`) to interpolate the mesh coordinate
vectors EXCEPT in the periodic case where it does not make sense since the coordinate vectors are not periodic.

See also: 
=== 
`DM`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMRefine()`, `DMCoarsen()`, `DMCreateRestriction()`, `DMCreateInterpolationScale()`

# External Links
$(_doc_external("DM/DMCreateInterpolation"))
"""
function DMCreateInterpolation(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib}) where {PetscLib}
	vec = CVec()

	LibPETSc.DMCreateInterpolation(
		PetscLib,
		dmc,
		dmf,
		mat,
		vec,
	)

	return mat,vec
end
 
 
"""
	 UNTESTED !!!
	mat = DMCreateRestriction(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib})

Gets restriction matrix between two `DM` objects. The resulting matrix map degrees of freedom in the vector obtained by
`DMCreateGlobalVector()` on the fine `DM` to similar vectors on the coarse grid `DM`.

Collective

Input Parameters:
===
- `dmc` - the `DM` object
- `dmf` - the second, finer `DM` object

Output Parameter:
===
- `mat` - the restriction

Level: developer

Note:
This only works for `DMSTAG`. For many situations either the transpose of the operator obtained with `DMCreateInterpolation()` or that
matrix multiplied by the vector obtained with `DMCreateInterpolationScale()` provides the desired object.

See also: 
=== 
`DM`, `DMRestrict()`, `DMInterpolate()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMRefine()`, `DMCoarsen()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMCreateRestriction"))
"""
function DMCreateRestriction(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCreateRestriction(
		PetscLib,
		dmc,
		dmf,
		mat,
	)

	return mat
end
 
 
"""
	 UNTESTED !!!
	mat = DMCreateInjection(dac::AbstractDM{PetscLib},daf::AbstractDM{PetscLib})

Gets injection matrix between two `DM` objects.

Collective

Input Parameters:
===
- `dac` - the `DM` object
- `daf` - the second, finer `DM` object

Output Parameter:
===
- `mat` - the injection

Level: developer

Notes:
This is an operator that applied to a vector obtained with `DMCreateGlobalVector()` on the
fine grid maps the values to a vector on the vector on the coarse `DM` by simply selecting
the values on the coarse grid points. This compares to the operator obtained by
`DMCreateRestriction()` or the transpose of the operator obtained by
`DMCreateInterpolation()` that uses a "local weighted average" of the values around the
coarse grid point as the coarse grid value.

For `DMDA` objects this only works for "uniform refinement", that is the refined mesh was obtained `DMRefine()` or the coarse mesh was obtained by
`DMCoarsen()`. The coordinates set into the `DMDA` are completely ignored in computing the injection.

See also: 
=== 
`DM`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMCreateInterpolation()`,
`DMCreateRestriction()`, `MatRestrict()`, `MatInterpolate()`

# External Links
$(_doc_external("DM/DMCreateInjection"))
"""
function DMCreateInjection(dac::AbstractDM{PetscLib},daf::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCreateInjection(
		PetscLib,
		dac,
		daf,
		mat,
	)

	return mat
end
 
 
"""
	 UNTESTED !!!
	mat = DMCreateMassMatrix(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib})

Gets the mass matrix between two `DM` objects, Mᵢⱼ = ∫ϕᵢϕⱼ	 where the ϕ are Galerkin basis functions for a
a Galerkin finite element model on the `DM`

Collective

Input Parameters:
===
- `dmc` - the target `DM` object
- `dmf` - the source `DM` object, can be `NULL`

Output Parameter:
===
- `mat` - the mass matrix

Level: developer

Notes:
For `DMPLEX` the finite element model for the `DM` must have been already provided.

if `dmc` is `dmf` or `NULL`, then x^t M x is an approximation to the L2 norm of the vector x which is obtained by `DMCreateGlobalVector()`

See also: 
=== 
`DM`, `DMCreateMassMatrixLumped()`, `DMCreateMatrix()`, `DMRefine()`, `DMCoarsen()`, `DMCreateRestriction()`, `DMCreateInterpolation()`, `DMCreateInjection()`

# External Links
$(_doc_external("DM/DMCreateMassMatrix"))
"""
function DMCreateMassMatrix(dmc::AbstractDM{PetscLib},dmf::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCreateMassMatrix(
		PetscLib,
		dmc,
		dmf,
		mat,
	)

	return mat
end
 
 
"""
	 UNTESTED !!!
	llm,lm = DMCreateMassMatrixLumped(dm::AbstractDM{PetscLib})

Gets the lumped mass matrix for a given `DM`

Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameters:
===
- `llm` - the local lumped mass matrix, which is a diagonal matrix, represented as a vector
- `lm`  - the global lumped mass matrix, which is a diagonal matrix, represented as a vector

Level: developer

Note:
See `DMCreateMassMatrix()` for how to create the non-lumped version of the mass matrix.

See also: 
=== 
`DM`, `DMCreateMassMatrix()`, `DMCreateMatrix()`, `DMRefine()`, `DMCoarsen()`, `DMCreateRestriction()`, `DMCreateInterpolation()`, `DMCreateInjection()`

# External Links
$(_doc_external("DM/DMCreateMassMatrixLumped"))
"""
function DMCreateMassMatrixLumped(dm::AbstractDM{PetscLib}) where {PetscLib}
	llm = CVec()
	lm = CVec()

	LibPETSc.DMCreateMassMatrixLumped(
		PetscLib,
		dm,
		llm,
		lm,
	)

	return llm,lm
end
 
 
"""
	 UNTESTED !!!
	mem = DMGetWorkArray(dm::AbstractDM{PetscLib},count::Int,dtype::MPI_Datatype)

Gets a work array guaranteed to be at least the input size, restore with `DMRestoreWorkArray()`

Not Collective

Input Parameters:
===
- `dm`    - the `DM` object
- `count` - The minimum size
- `dtype` - MPI data type, often `MPIU_REAL`, `MPIU_SCALAR`, or `MPIU_INT`)

Output Parameter:
===
- `mem` - the work array

Level: developer

Notes:
A `DM` may stash the array between instantiations so using this routine may be more efficient than calling `PetscMalloc()`

The array may contain nonzero values

See also: 
=== 
`DM`, `DMDestroy()`, `DMCreate()`, `DMRestoreWorkArray()`, `PetscMalloc()`

# External Links
$(_doc_external("DM/DMGetWorkArray"))
"""
function DMGetWorkArray(dm::AbstractDM{PetscLib},count::Int,dtype::MPI_Datatype) where {PetscLib}
	Float64 = PetscLib.Float64
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_mem = PETSc_RefPtr(dims, Float64)

	LibPETSc.DMGetWorkArray(
		PetscLib,
		dm,
		count,
		dtype,
		r_mem,
	)

	mem = PETSc_unsafe_wrap(r_mem, dims; own=false)

	return mem
end
 
 
"""
	 UNTESTED !!!
	mem = DMRestoreWorkArray(dm::AbstractDM{PetscLib},count::Int,dtype::MPI_Datatype)

Restores a work array obtained with `DMCreateWorkArray()`

Not Collective

Input Parameters:
===
- `dm`    - the `DM` object
- `count` - The minimum size
- `dtype` - MPI data type, often `MPIU_REAL`, `MPIU_SCALAR`, `MPIU_INT`

Output Parameter:
===
- `mem` - the work array

Level: developer

Developer Note:
count and dtype are ignored, they are only needed for `DMGetWorkArray()`

See also: 
=== 
`DM`, `DMDestroy()`, `DMCreate()`, `DMGetWorkArray()`

# External Links
$(_doc_external("DM/DMRestoreWorkArray"))
"""
function DMRestoreWorkArray(dm::AbstractDM{PetscLib},count::Int,dtype::MPI_Datatype) where {PetscLib}
	Float64 = PetscLib.Float64
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_mem = PETSc_RefPtr(dims, Float64)

	LibPETSc.DMRestoreWorkArray(
		PetscLib,
		dm,
		count,
		dtype,
		r_mem,
	)

	mem = PETSc_unsafe_wrap(r_mem, dims; own=false)

	return mem
end
 
 
"""
	 UNTESTED !!!
	cdm = DMGetCoarseDM(dm::AbstractDM{PetscLib})

Get the coarse `DM`from which this `DM` was obtained by refinement

Not Collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `cdm` - The coarse `DM`

Level: intermediate

See also: 
=== 
`DM`, `DMSetCoarseDM()`, `DMCoarsen()`

# External Links
$(_doc_external("DM/DMGetCoarseDM"))
"""
function DMGetCoarseDM(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	cdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMGetCoarseDM(
		PetscLib,
		dm,
		cdm,
	)

	return cdm
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoarseDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib})

Set the coarse `DM` from which this `DM` was obtained by refinement

Input Parameters:
===
- `dm`  - The `DM` object
- `cdm` - The coarse `DM`

Level: intermediate

Note:
Normally this is set automatically by `DMRefine()`

See also: 
=== 
`DM`, `DMGetCoarseDM()`, `DMCoarsen()`, `DMSetRefine()`, `DMSetFineDM()`

# External Links
$(_doc_external("DM/DMSetCoarseDM"))
"""
function DMSetCoarseDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetCoarseDM(
		PetscLib,
		dm,
		cdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	fdm = DMGetFineDM(dm::AbstractDM{PetscLib})

Get the fine mesh from which this `DM` was obtained by coarsening

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `fdm` - The fine `DM`

Level: intermediate

See also: 
=== 
`DM`, `DMSetFineDM()`, `DMCoarsen()`, `DMRefine()`

# External Links
$(_doc_external("DM/DMGetFineDM"))
"""
function DMGetFineDM(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	fdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMGetFineDM(
		PetscLib,
		dm,
		fdm,
	)

	return fdm
end
 
 
"""
	 UNTESTED !!!
	 DMSetFineDM(dm::AbstractDM{PetscLib},fdm::AbstractDM{PetscLib})

Set the fine mesh from which this was obtained by coarsening

Input Parameters:
===
- `dm`  - The `DM` object
- `fdm` - The fine `DM`

Level: developer

Note:
Normally this is set automatically by `DMCoarsen()`

See also: 
=== 
`DM`, `DMGetFineDM()`, `DMCoarsen()`, `DMRefine()`

# External Links
$(_doc_external("DM/DMSetFineDM"))
"""
function DMSetFineDM(dm::AbstractDM{PetscLib},fdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetFineDM(
		PetscLib,
		dm,
		fdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMRefineHierarchy(dm::AbstractDM{PetscLib},nlevels::Int,dmf::Vector{DM})

Refines a `DM` object, all levels at once

Collective

Input Parameters:
===
- `dm`      - the `DM` object
- `nlevels` - the number of levels of refinement

Output Parameter:
===
- `dmf` - the refined `DM` hierarchy

Level: developer

See also: 
=== 
`DM`, `DMCoarsen()`, `DMCoarsenHierarchy()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMRefineHierarchy"))
"""
function DMRefineHierarchy(dm::AbstractDM{PetscLib},nlevels::Int,dmf::Vector{DM}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMRefineHierarchy(
		PetscLib,
		dm,
		nlevels,
		dmf,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCoarsenHierarchy(dm::AbstractDM{PetscLib},nlevels::Int,dmc::Vector{DM})

Coarsens a `DM` object, all levels at once

Collective

Input Parameters:
===
- `dm`      - the `DM` object
- `nlevels` - the number of levels of coarsening

Output Parameter:
===
- `dmc` - the coarsened `DM` hierarchy

Level: developer

See also: 
=== 
`DM`, `DMCoarsen()`, `DMRefineHierarchy()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMCoarsenHierarchy"))
"""
function DMCoarsenHierarchy(dm::AbstractDM{PetscLib},nlevels::Int,dmc::Vector{DM}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMCoarsenHierarchy(
		PetscLib,
		dm,
		nlevels,
		dmc,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMRestrict(fine::AbstractDM{PetscLib},restrct::AbstractMatrix,rscale::AbstractVector,inject::AbstractMatrix,coarse::AbstractDM{PetscLib})

restricts user

Collective if any hooks are

Input Parameters:
===
- `fine`    - finer `DM` from which the data is obtained
- `restrct` - restriction matrix, apply using `MatRestrict()`, usually the transpose of the interpolation
- `rscale`  - scaling vector for restriction
- `inject`  - injection matrix, also use `MatRestrict()`
- `coarse`  - coarser `DM` to update

Level: developer

Developer Note:
Though this routine is called `DMRestrict()` the hooks are added with `DMCoarsenHookAdd()`, a consistent terminology would be better

See also: 
=== 
`DM`, `DMCoarsenHookAdd()`, `MatRestrict()`, `DMInterpolate()`, `DMRefineHookAdd()`

# External Links
$(_doc_external("DM/DMRestrict"))
"""
function DMRestrict(fine::AbstractDM{PetscLib},restrct::AbstractMatrix,rscale::AbstractVector,inject::AbstractMatrix,coarse::AbstractDM{PetscLib}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMRestrict(
		PetscLib,
		fine,
		restrct,
		rscale,
		inject,
		coarse,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMInterpolate(coarse::AbstractDM{PetscLib},interp::AbstractMatrix,fine::AbstractDM{PetscLib})

interpolates user

Collective if any hooks are

Input Parameters:
===
- `coarse` - coarser `DM` to use as a base
- `interp` - interpolation matrix, apply using `MatInterpolate()`
- `fine`   - finer `DM` to update

Level: developer

Developer Note:
This routine is called `DMInterpolate()` while the hook is called `DMRefineHookAdd()`. It would be better to have an
an API with consistent terminology.

See also: 
=== 
`DM`, `DMRefineHookAdd()`, `MatInterpolate()`

# External Links
$(_doc_external("DM/DMInterpolate"))
"""
function DMInterpolate(coarse::AbstractDM{PetscLib},interp::AbstractMatrix,fine::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMInterpolate(
		PetscLib,
		coarse,
		interp,
		fine,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMInterpolateSolution(coarse::AbstractDM{PetscLib},fine::AbstractDM{PetscLib},interp::AbstractMatrix,coarseSol::AbstractVector,fineSol::AbstractVector)

Interpolates a solution from a coarse mesh to a fine mesh.

Collective

Input Parameters:
===
- `coarse`    - coarse `DM`
- `fine`      - fine `DM`
- `interp`    - (optional) the matrix computed by `DMCreateInterpolation()`.  Implementations may not need this, but if it
is available it can avoid some recomputation.  If it is provided, `MatInterpolate()` will be used if
the coarse `DM` does not have a specialized implementation.
- `coarseSol` - solution on the coarse mesh

Output Parameter:
===
- `fineSol` - the interpolation of coarseSol to the fine mesh

Level: developer

Note:
This function exists because the interpolation of a solution vector between meshes is not always a linear
map.  For example, if a boundary value problem has an inhomogeneous Dirichlet boundary condition that is compressed
out of the solution vector.  Or if interpolation is inherently a nonlinear operation, such as a method using
slope-limiting reconstruction.

Developer Note:
This doesn't just interpolate "solutions" so its API name is questionable.

See also: 
=== 
`DM`, `DMInterpolate()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMInterpolateSolution"))
"""
function DMInterpolateSolution(coarse::AbstractDM{PetscLib},fine::AbstractDM{PetscLib},interp::AbstractMatrix,coarseSol::AbstractVector,fineSol::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMInterpolateSolution(
		PetscLib,
		coarse,
		fine,
		interp,
		coarseSol,
		fineSol,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dme = DMExtrude(dm::AbstractDM{PetscLib},layers::Int)

Extrude a `DM` object from a surface

Collective

Input Parameters:
===
- `dm`     - the `DM` object
- `layers` - the number of extruded cell layers

Output Parameter:
===
- `dme` - the extruded `DM`, or `NULL`

Level: developer

Note:
If no extrusion was done, the return value is `NULL`

See also: 
=== 
`DM`, `DMRefine()`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`

# External Links
$(_doc_external("DM/DMExtrude"))
"""
function DMExtrude(dm::AbstractDM{PetscLib},layers::Int) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dme = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMExtrude(
		PetscLib,
		dm,
		layers,
		dme,
	)

	return dme
end
 
 
"""
	 UNTESTED !!!
	 DMSetFromOptions(dm::AbstractDM{PetscLib})

sets parameters in a `DM` from the options database

Collective

Input Parameter:
===
- `dm` - the `DM` object to set options for

Options Database Keys:
===
- `-dm_preallocate_only`                               - Only preallocate the matrix for `DMCreateMatrix()` and `DMCreateMassMatrix()`, but do not fill it with zeros
- `-dm_vec_type <type>`                                - type of vector to create inside `DM`
- `-dm_mat_type <type>`                                - type of matrix to create inside `DM`
- `-dm_is_coloring_type`                               - <global or local>
- `-dm_bind_below <n>`                                 - bind (force execution on CPU) for `Vec` and `Mat` objects with local size (number of vector entries or matrix rows) below n; currently only supported for `DMDA`
- `-dm_plex_filename <str>`                            - File containing a mesh
- `-dm_plex_boundary_filename <str>`                   - File containing a mesh boundary
- `-dm_plex_name <str>`                                - Name of the mesh in the file
- `-dm_plex_shape <shape>`                             - The domain shape, such as `BOX`, `SPHERE`, etc.
- `-dm_plex_cell <ct>`                                 - Cell shape
- `-dm_plex_reference_cell_domain <bool>`              - Use a reference cell domain
- `-dm_plex_dim <dim>`                                 - Set the topological dimension
- `-dm_plex_simplex <bool>`                            - `PETSC_TRUE` for simplex elements, `PETSC_FALSE` for tensor elements
- `-dm_plex_interpolate <bool>`                        - `PETSC_TRUE` turns on topological interpolation (creating edges and faces)
- `-dm_plex_scale <sc>`                                - Scale factor for mesh coordinates
- `-dm_coord_remap <bool>`                             - Map coordinates using a function
- `-dm_coord_map <mapname>`                            - Select a builtin coordinate map
- `-dm_coord_map_params <p0,p1,p2,...>`                - Set coordinate mapping parameters
- `-dm_plex_box_faces <m,n,p>`                         - Number of faces along each dimension
- `-dm_plex_box_lower <x,y,z>`                         - Specify lower-left-bottom coordinates for the box
- `-dm_plex_box_upper <x,y,z>`                         - Specify upper-right-top coordinates for the box
- `-dm_plex_box_bd <bx,by,bz>`                         - Specify the `DMBoundaryType` for each direction
- `-dm_plex_sphere_radius <r>`                         - The sphere radius
- `-dm_plex_ball_radius <r>`                           - Radius of the ball
- `-dm_plex_cylinder_bd <bz>`                          - Boundary type in the z direction
- `-dm_plex_cylinder_num_wedges <n>`                   - Number of wedges around the cylinder
- `-dm_plex_reorder <order>`                           - Reorder the mesh using the specified algorithm
- `-dm_refine_pre <n>`                                 - The number of refinements before distribution
- `-dm_refine_uniform_pre <bool>`                      - Flag for uniform refinement before distribution
- `-dm_refine_volume_limit_pre <v>`                    - The maximum cell volume after refinement before distribution
- `-dm_refine <n>`                                     - The number of refinements after distribution
- `-dm_extrude <l>`                                    - Activate extrusion and specify the number of layers to extrude
- `-dm_plex_transform_extrude_thickness <t>`           - The total thickness of extruded layers
- `-dm_plex_transform_extrude_use_tensor <bool>`       - Use tensor cells when extruding
- `-dm_plex_transform_extrude_symmetric <bool>`        - Extrude layers symmetrically about the surface
- `-dm_plex_transform_extrude_normal <n0,...,nd>`      - Specify the extrusion direction
- `-dm_plex_transform_extrude_thicknesses <t0,...,tl>` - Specify thickness of each layer
- `-dm_plex_create_fv_ghost_cells`                     - Flag to create finite volume ghost cells on the boundary
- `-dm_plex_fv_ghost_cells_label <name>`               - Label name for ghost cells boundary
- `-dm_distribute <bool>`                              - Flag to redistribute a mesh among processes
- `-dm_distribute_overlap <n>`                         - The size of the overlap halo
- `-dm_plex_adj_cone <bool>`                           - Set adjacency direction
- `-dm_plex_adj_closure <bool>`                        - Set adjacency size
- `-dm_plex_use_ceed <bool>`                           - Use LibCEED as the FEM backend
- `-dm_plex_check_symmetry`                            - Check that the adjacency information in the mesh is symmetric - `DMPlexCheckSymmetry()`
- `-dm_plex_check_skeleton`                            - Check that each cell has the correct number of vertices (only for homogeneous simplex or tensor meshes) - `DMPlexCheckSkeleton()`
- `-dm_plex_check_faces`                               - Check that the faces of each cell give a vertex order this is consistent with what we expect from the cell type - `DMPlexCheckFaces()`
- `-dm_plex_check_geometry`                            - Check that cells have positive volume - `DMPlexCheckGeometry()`
- `-dm_plex_check_pointsf`                             - Check some necessary conditions for `PointSF` - `DMPlexCheckPointSF()`
- `-dm_plex_check_interface_cones`                     - Check points on inter-partition interfaces have conforming order of cone points - `DMPlexCheckInterfaceCones()`
- `-dm_plex_check_all`                                 - Perform all the checks above

Level: intermediate

Note:
For some `DMType` such as `DMDA` this cannot be called after `DMSetUp()` has been called.

See also: 
=== 
`DM`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`,
`DMPlexCheckSymmetry()`, `DMPlexCheckSkeleton()`, `DMPlexCheckFaces()`, `DMPlexCheckGeometry()`, `DMPlexCheckPointSF()`, `DMPlexCheckInterfaceCones()`,
`DMSetOptionsPrefix()`, `DMType`, `DMPLEX`, `DMDA`, `DMSetUp()`

# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
function DMSetFromOptions(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetFromOptions(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMViewFromOptions(dm::AbstractDM{PetscLib},obj::PetscObject,name::Vector{Char})

View a `DM` in a particular way based on a request in the options database

Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `obj`  - optional object that provides the prefix for the options database (if `NULL` then the prefix in obj is used)
- `name` - option string that is used to activate viewing

Level: intermediate

Note:
See `PetscObjectViewFromOptions()` for a list of values that can be provided in the options database to determine how the `DM` is viewed

See also: 
=== 
`DM`, `DMView()`, `PetscObjectViewFromOptions()`, `DMCreate()`

# External Links
$(_doc_external("DM/DMViewFromOptions"))
"""
function DMViewFromOptions(dm::AbstractDM{PetscLib},obj::PetscObject,name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMViewFromOptions(
		PetscLib,
		dm,
		obj,
		name,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dmAdapt = DMAdaptLabel(dm::AbstractDM{PetscLib},label::DMLabel)

Adapt a `DM` based on a `DMLabel` with values interpreted as coarsening and refining flags.  Specific implementations of `DM` maybe have
specialized flags, but all implementations should accept flag values `DM_ADAPT_DETERMINE`, `DM_ADAPT_KEEP`, `DM_ADAPT_REFINE`, and,
`DM_ADAPT_COARSEN`.

Collective

Input Parameters:
===
- `dm`    - the pre-adaptation `DM` object
- `label` - label with the flags

Output Parameters:
===
- `dmAdapt` - the adapted `DM` object: may be `NULL` if an adapted `DM` could not be produced.

Level: intermediate

-seealso: `DM`, `DMAdaptMetric()`, `DMCoarsen()`, `DMRefine()`

# External Links
$(_doc_external("DM/DMAdaptLabel"))
"""
function DMAdaptLabel(dm::AbstractDM{PetscLib},label::DMLabel) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dmAdapt = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMAdaptLabel(
		PetscLib,
		dm,
		label,
		dmAdapt,
	)

	return dmAdapt
end
 
 
"""
	 UNTESTED !!!
	dmAdapt = DMAdaptMetric(dm::AbstractDM{PetscLib},metric::AbstractVector,bdLabel::DMLabel,rgLabel::DMLabel)

Generates a mesh adapted to the specified metric field.

Input Parameters:
===
- `dm`      - The DM object
- `metric`  - The metric to which the mesh is adapted, defined vertex-wise.
- `bdLabel` - Label for boundary tags, which will be preserved in the output mesh. `bdLabel` should be `NULL` if there is no such label, and should be different from "_boundary_".
- `rgLabel` - Label for cell tags, which will be preserved in the output mesh. `rgLabel` should be `NULL` if there is no such label, and should be different from "_regions_".

Output Parameter:
===
- `dmAdapt` - Pointer to the `DM` object containing the adapted mesh

Note:
The label in the adapted mesh will be registered under the name of the input `DMLabel` object

Level: advanced

-seealso: `DMAdaptLabel()`, `DMCoarsen()`, `DMRefine()`

# External Links
$(_doc_external("DM/DMAdaptMetric"))
"""
function DMAdaptMetric(dm::AbstractDM{PetscLib},metric::AbstractVector,bdLabel::DMLabel,rgLabel::DMLabel) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dmAdapt = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMAdaptMetric(
		PetscLib,
		dm,
		metric,
		bdLabel,
		rgLabel,
		dmAdapt,
	)

	return dmAdapt
end
 
 
"""
	 UNTESTED !!!
	 DMSetUp(dm::AbstractDM{PetscLib})

sets up the data structures inside a `DM` object

Collective

Input Parameter:
===
- `dm` - the `DM` object to setup

Level: intermediate

Note:
This is usually called after various parameter setting operations and `DMSetFromOptions()` are called on the `DM`

See also: 
=== 
`DM`, `DMCreate()`, `DMSetType()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`

# External Links
$(_doc_external("DM/DMSetUp"))
"""
function DMSetUp(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetUp(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	scale = DMCreateInterpolationScale(dac::AbstractDM{PetscLib},daf::AbstractDM{PetscLib},mat::AbstractMatrix)

Forms L = 1/(R*1) where 1 is the vector of all ones, and R is
the transpose of the interpolation between the `DM`.

Input Parameters:
===
- `dac` - `DM` that defines a coarse mesh
- `daf` - `DM` that defines a fine mesh
- `mat` - the restriction (or interpolation operator) from fine to coarse

Output Parameter:
===
- `scale` - the scaled vector

Level: advanced

Note:
xcoarse = diag(L)*R*xfine preserves scale and is thus suitable for state (versus residual)
restriction. In other words xcoarse is the coarse representation of xfine.

Developer Note:
If the fine-scale `DMDA` has the -dm_bind_below option set to true, then `DMCreateInterpolationScale()` calls `MatSetBindingPropagates()`
on the restriction/interpolation operator to set the bindingpropagates flag to true.

See also: 
=== 
`DM`, `MatRestrict()`, `MatInterpolate()`, `DMCreateInterpolation()`, `DMCreateRestriction()`, `DMCreateGlobalVector()`

# External Links
$(_doc_external("DM/DMCreateInterpolationScale"))
"""
function DMCreateInterpolationScale(dac::AbstractDM{PetscLib},daf::AbstractDM{PetscLib},mat::AbstractMatrix) where {PetscLib}
	scale = CVec()

	LibPETSc.DMCreateInterpolationScale(
		PetscLib,
		dac,
		daf,
		mat,
		scale,
	)

	return scale
end
 
 
"""
	 UNTESTED !!!
	 DMGlobalToLocalBegin(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector)

Begins updating local vectors from global vector

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `g`    - the global vector
- `mode` - `INSERT_VALUES` or `ADD_VALUES`
- `l`    - the local vector

Level: intermediate

Notes:
The operation is completed with `DMGlobalToLocalEnd()`

One can perform local computations between the `DMGlobalToLocalBegin()` and  `DMGlobalToLocalEnd()` to overlap communication and computation

`DMGlobalToLocal()` is a short form of  `DMGlobalToLocalBegin()` and  `DMGlobalToLocalEnd()`

`DMGlobalToLocalHookAdd()` may be used to provide additional operations that are performed during the update process.

See also: 
=== 
`DM`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocal()`, `DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`, `DMLocalToGlobal()`, `DMLocalToGlobalEnd()`

# External Links
$(_doc_external("DM/DMGlobalToLocalBegin"))
"""
function DMGlobalToLocalBegin(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGlobalToLocalBegin(
		PetscLib,
		dm,
		g,
		mode,
		l,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGlobalToLocalEnd(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector)

Ends updating local vectors from global vector

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `g`    - the global vector
- `mode` - `INSERT_VALUES` or `ADD_VALUES`
- `l`    - the local vector

Level: intermediate

Note:
See `DMGlobalToLocalBegin()` for details.

See also: 
=== 
`DM`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocal()`, `DMLocalToGlobalBegin()`, `DMLocalToGlobal()`, `DMLocalToGlobalEnd()`

# External Links
$(_doc_external("DM/DMGlobalToLocalEnd"))
"""
function DMGlobalToLocalEnd(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGlobalToLocalEnd(
		PetscLib,
		dm,
		g,
		mode,
		l,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalToGlobalBegin(dm::AbstractDM{PetscLib},l::AbstractVector,mode::InsertMode,g::AbstractVector)

begins updating global vectors from local vectors

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `l`    - the local vector
- `mode` - if `INSERT_VALUES` then no parallel communication is used, if `ADD_VALUES` then all ghost points from the same base point accumulate into that base point.
- `g`    - the global vector

Level: intermediate

Notes:
In the `ADD_VALUES` case you normally would zero the receiving vector before beginning this operation.

`INSERT_VALUES is` not supported for `DMDA`, in that case simply compute the values directly into a global vector instead of a local one.

Use `DMLocalToGlobalEnd()` to complete the communication process.

`DMLocalToGlobal()` is a short form of  `DMLocalToGlobalBegin()` and  `DMLocalToGlobalEnd()`

`DMLocalToGlobalHookAdd()` may be used to provide additional operations that are performed during the update process.

See also: 
=== 
`DM`, `DMLocalToGlobal()`, `DMLocalToGlobalEnd()`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocal()`, `DMGlobalToLocalEnd()`, `DMGlobalToLocalBegin()`

# External Links
$(_doc_external("DM/DMLocalToGlobalBegin"))
"""
function DMLocalToGlobalBegin(dm::AbstractDM{PetscLib},l::AbstractVector,mode::InsertMode,g::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMLocalToGlobalBegin(
		PetscLib,
		dm,
		l,
		mode,
		g,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalToGlobalEnd(dm::AbstractDM{PetscLib},l::AbstractVector,mode::InsertMode,g::AbstractVector)

updates global vectors from local vectors

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `l`    - the local vector
- `mode` - `INSERT_VALUES` or `ADD_VALUES`
- `g`    - the global vector

Level: intermediate

Note:
See `DMLocalToGlobalBegin()` for full details

See also: 
=== 
`DM`, `DMLocalToGlobalBegin()`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocalEnd()`

# External Links
$(_doc_external("DM/DMLocalToGlobalEnd"))
"""
function DMLocalToGlobalEnd(dm::AbstractDM{PetscLib},l::AbstractVector,mode::InsertMode,g::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMLocalToGlobalEnd(
		PetscLib,
		dm,
		l,
		mode,
		g,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalToLocalBegin(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector)

Begins the process of mapping values from a local vector (that include
ghost points that contain irrelevant values) to another local vector where the ghost points
in the second are set correctly from values on other MPI ranks.

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `g`    - the original local vector
- `mode` - one of `INSERT_VALUES` or `ADD_VALUES`

Output Parameter:
===
- `l` - the local vector with correct ghost values

Level: intermediate

Note:
Must be followed by `DMLocalToLocalEnd()`.

See also: 
=== 
`DM`, `DMLocalToLocalEnd()`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateLocalVector()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`

# External Links
$(_doc_external("DM/DMLocalToLocalBegin"))
"""
function DMLocalToLocalBegin(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMLocalToLocalBegin(
		PetscLib,
		dm,
		g,
		mode,
		l,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalToLocalEnd(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector)

Maps from a local vector to another local vector where the ghost
points in the second are set correctly. Must be preceded by `DMLocalToLocalBegin()`.

Neighbor-wise Collective

Input Parameters:
===
- `dm`   - the `DM` object
- `g`    - the original local vector
- `mode` - one of `INSERT_VALUES` or `ADD_VALUES`

Output Parameter:
===
- `l` - the local vector with correct ghost values

Level: intermediate

See also: 
=== 
`DM`, `DMLocalToLocalBegin()`, `DMCoarsen()`, `DMDestroy()`, `DMView()`, `DMCreateLocalVector()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMGlobalToLocalEnd()`, `DMLocalToGlobalBegin()`

# External Links
$(_doc_external("DM/DMLocalToLocalEnd"))
"""
function DMLocalToLocalEnd(dm::AbstractDM{PetscLib},g::AbstractVector,mode::InsertMode,l::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMLocalToLocalEnd(
		PetscLib,
		dm,
		g,
		mode,
		l,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	M = DMConvert(dm::AbstractDM{PetscLib},newtype::DMType)

Converts a `DM` to another `DM`, either of the same or different type.

Collective

Input Parameters:
===
- `dm`      - the `DM`
- `newtype` - new `DM` type (use "same" for the same type)

Output Parameter:
===
- `M` - pointer to new `DM`

Level: intermediate

Note:
Cannot be used to convert a sequential `DM` to a parallel or a parallel to sequential,
the MPI communicator of the generated `DM` is always the same as the communicator
of the input `DM`.

See also: 
=== 
`DM`, `DMSetType()`, `DMCreate()`, `DMClone()`

# External Links
$(_doc_external("DM/DMConvert"))
"""
function DMConvert(dm::AbstractDM{PetscLib},newtype::DMType) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	M = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMConvert(
		PetscLib,
		dm,
		newtype,
		M,
	)

	return M
end
 
 
"""
	 UNTESTED !!!
	dim = DMGetDimension(dm::AbstractDM{PetscLib})

Return the topological dimension of the `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `dim` - The topological dimension

Level: beginner

See also: 
=== 
`DM`, `DMSetDimension()`, `DMCreate()`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
function DMGetDimension(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	dim = [PetscInt(1)]

	LibPETSc.DMGetDimension(
		PetscLib,
		dm,
		Ref(dim,1),
	)

	return dim[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetDimension(dm::AbstractDM{PetscLib},dim::Int)

Set the topological dimension of the `DM`

Collective

Input Parameters:
===
- `dm`  - The `DM`
- `dim` - The topological dimension

Level: beginner

See also: 
=== 
`DM`, `DMGetDimension()`, `DMCreate()`

# External Links
$(_doc_external("DM/DMSetDimension"))
"""
function DMSetDimension(dm::AbstractDM{PetscLib},dim::Int) where {PetscLib}

	LibPETSc.DMSetDimension(
		PetscLib,
		dm,
		dim,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	pStart,pEnd = DMGetDimPoints(dm::AbstractDM{PetscLib},dim::Int)

Get the half

Collective

Input Parameters:
===
- `dm`  - the `DM`
- `dim` - the dimension

Output Parameters:
===
- `pStart` - The first point of the given dimension
- `pEnd`   - The first point following points of the given dimension

Level: intermediate

Note:
The points are vertices in the Hasse diagram encoding the topology. This is explained in
https://arxiv.org/abs/0908.4427. If no points exist of this dimension in the storage scheme,
then the interval is empty.

See also: 
=== 
`DM`, `DMPLEX`, `DMPlexGetDepthStratum()`, `DMPlexGetHeightStratum()`

# External Links
$(_doc_external("DM/DMGetDimPoints"))
"""
function DMGetDimPoints(dm::AbstractDM{PetscLib},dim::Int) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	pStart = [PetscInt(1)]
	pEnd = [PetscInt(1)]

	LibPETSc.DMGetDimPoints(
		PetscLib,
		dm,
		dim,
		Ref(pStart,1),
		Ref(pEnd,1),
	)

	return pStart[1],pEnd[1]
end
 
 
"""
	 UNTESTED !!!
	useNatural = DMGetUseNatural(dm::AbstractDM{PetscLib})

Get the flag for creating a mapping to the natural order when a `DM` is (re)distributed in parallel

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `useNatural` - `PETSC_TRUE` to build the mapping to a natural order during distribution

Level: beginner

See also: 
=== 
`DM`, `DMSetUseNatural()`, `DMCreate()`

# External Links
$(_doc_external("DM/DMGetUseNatural"))
"""
function DMGetUseNatural(dm::AbstractDM{PetscLib}) where {PetscLib}
	useNatural = Ref{PetscBool}()

	LibPETSc.DMGetUseNatural(
		PetscLib,
		dm,
		useNatural,
	)

	return useNatural[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMSetUseNatural(dm::AbstractDM{PetscLib},useNatural::PetscBool)

Set the flag for creating a mapping to the natural order when a `DM` is (re)distributed in parallel

Collective

Input Parameters:
===
- `dm`         - The `DM`
- `useNatural` - `PETSC_TRUE` to build the mapping to a natural order during distribution

Level: beginner

Note:
This also causes the map to be build after `DMCreateSubDM()` and `DMCreateSuperDM()`

See also: 
=== 
`DM`, `DMGetUseNatural()`, `DMCreate()`, `DMPlexDistribute()`, `DMCreateSubDM()`, `DMCreateSuperDM()`

# External Links
$(_doc_external("DM/DMSetUseNatural"))
"""
function DMSetUseNatural(dm::AbstractDM{PetscLib},useNatural::PetscBool) where {PetscLib}

	LibPETSc.DMSetUseNatural(
		PetscLib,
		dm,
		useNatural,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	nranks = DMGetNeighbors(dm::AbstractDM{PetscLib},ranks::Vector{PetscMPIInt})

Gets an array containing the MPI ranks of all the processes neighbors

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameters:
===
- `nranks` - the number of neighbours
- `ranks`  - the neighbors ranks

Level: beginner

Note:
Do not free the array, it is freed when the `DM` is destroyed.

See also: 
=== 
`DM`, `DMDAGetNeighbors()`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("DM/DMGetNeighbors"))
"""
function DMGetNeighbors(dm::AbstractDM{PetscLib},ranks::Vector{PetscMPIInt}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	nranks = [PetscInt(1)]

	LibPETSc.DMGetNeighbors(
		PetscLib,
		dm,
		ranks,
		Ref(nranks,1),
	)

	return nranks[1]
end
 
 
"""
	 UNTESTED !!!
	cdm = DMGetCoordinateDM(dm::AbstractDM{PetscLib})

Gets the `DM` that prescribes coordinate layout and scatters between global and local coordinates

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `cdm` - coordinate `DM`

Level: intermediate

-seealso: `DM`, `DMSetCoordinateDM()`, `DMSetCoordinates()`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`, `DMGSetCellCoordinateDM()`,


# External Links
$(_doc_external("DM/DMGetCoordinateDM"))
"""
function DMGetCoordinateDM(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	cdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMGetCoordinateDM(
		PetscLib,
		dm,
		cdm,
	)

	return cdm
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinateDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib})

Sets the `DM` that prescribes coordinate layout and scatters between global and local coordinates

Logically Collective

Input Parameters:
===
- `dm`  - the `DM`
- `cdm` - coordinate `DM`

Level: intermediate

-seealso: `DM`, `DMGetCoordinateDM()`, `DMSetCoordinates()`, `DMGetCellCoordinateDM()`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`,
`DMGSetCellCoordinateDM()`

# External Links
$(_doc_external("DM/DMSetCoordinateDM"))
"""
function DMSetCoordinateDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetCoordinateDM(
		PetscLib,
		dm,
		cdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	cdm = DMGetCellCoordinateDM(dm::AbstractDM{PetscLib})

Gets the `DM` that prescribes cellwise coordinate layout and scatters between global and local cellwise coordinates

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `cdm` - cellwise coordinate `DM`, or `NULL` if they are not defined

Level: intermediate

Note:
Call `DMLocalizeCoordinates()` to automatically create cellwise coordinates for periodic geometries.

-seealso: `DM`, `DMSetCellCoordinateDM()`, `DMSetCellCoordinates()`, `DMSetCellCoordinatesLocal()`, `DMGetCellCoordinates()`, `DMGetCellCoordinatesLocal()`,
`DMLocalizeCoordinates()`, `DMSetCoordinateDM()`, `DMGetCoordinateDM()`

# External Links
$(_doc_external("DM/DMGetCellCoordinateDM"))
"""
function DMGetCellCoordinateDM(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	cdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMGetCellCoordinateDM(
		PetscLib,
		dm,
		cdm,
	)

	return cdm
end
 
 
"""
	 UNTESTED !!!
	 DMSetCellCoordinateDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib})

Sets the `DM` that prescribes cellwise coordinate layout and scatters between global and local cellwise coordinates

Logically Collective

Input Parameters:
===
- `dm`  - the `DM`
- `cdm` - cellwise coordinate `DM`

Level: intermediate

Note:
As opposed to `DMSetCoordinateDM()` these coordinates are useful for discontinuous Galerkin methods since they support coordinate fields that are discontinuous at cell boundaries.

-seealso: `DMGetCellCoordinateDM()`, `DMSetCellCoordinates()`, `DMSetCellCoordinatesLocal()`, `DMGetCellCoordinates()`, `DMGetCellCoordinatesLocal()`,
`DMSetCoordinateDM()`, `DMGetCoordinateDM()`

# External Links
$(_doc_external("DM/DMSetCellCoordinateDM"))
"""
function DMSetCellCoordinateDM(dm::AbstractDM{PetscLib},cdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSetCellCoordinateDM(
		PetscLib,
		dm,
		cdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dim = DMGetCoordinateDim(dm::AbstractDM{PetscLib})

Retrieve the dimension of the embedding space for coordinate values. For example a mesh on the surface of a sphere would have a 3 dimensional embedding space

Not Collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `dim` - The embedding dimension

Level: intermediate

-seealso: `DM`, `DMSetCoordinateDim()`, `DMGetCoordinateSection()`, `DMGetCoordinateDM()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMGetCoordinateDim"))
"""
function DMGetCoordinateDim(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	dim = [PetscInt(1)]

	LibPETSc.DMGetCoordinateDim(
		PetscLib,
		dm,
		Ref(dim,1),
	)

	return dim[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinateDim(dm::AbstractDM{PetscLib},dim::Int)

Set the dimension of the embedding space for coordinate values.

Not Collective

Input Parameters:
===
- `dm`  - The `DM` object
- `dim` - The embedding dimension

Level: intermediate

-seealso: `DM`, `DMGetCoordinateDim()`, `DMSetCoordinateSection()`, `DMGetCoordinateSection()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMSetCoordinateDim"))
"""
function DMSetCoordinateDim(dm::AbstractDM{PetscLib},dim::Int) where {PetscLib}

	LibPETSc.DMSetCoordinateDim(
		PetscLib,
		dm,
		dim,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section = DMGetCoordinateSection(dm::AbstractDM{PetscLib})

Retrieve the `PetscSection` of coordinate values over the mesh.

Collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `section` - The `PetscSection` object

Level: intermediate

Note:
This just retrieves the local section from the coordinate `DM`. In other words,
-vb
DMGetCoordinateDM(dm, &cdm);
DMGetLocalSection(cdm, &section);
-ve

-seealso: `DM`, `DMGetCoordinateDM()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMGetCoordinateSection"))
"""
function DMGetCoordinateSection(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()

	LibPETSc.DMGetCoordinateSection(
		PetscLib,
		dm,
		section,
	)

	return section
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinateSection(dm::AbstractDM{PetscLib},dim::Int,section::PetscSection)

Set the `PetscSection` of coordinate values over the mesh.

Not Collective

Input Parameters:
===
- `dm`      - The `DM` object
- `dim`     - The embedding dimension, or `PETSC_DETERMINE`
- `section` - The `PetscSection` object

Level: intermediate

-seealso: `DM`, `DMGetCoordinateDim()`, `DMGetCoordinateSection()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMSetCoordinateSection"))
"""
function DMSetCoordinateSection(dm::AbstractDM{PetscLib},dim::Int,section::PetscSection) where {PetscLib}

	LibPETSc.DMSetCoordinateSection(
		PetscLib,
		dm,
		dim,
		section,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section = DMGetCellCoordinateSection(dm::AbstractDM{PetscLib})

Retrieve the `PetscSection` of cellwise coordinate values over the mesh.

Collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `section` - The `PetscSection` object, or `NULL` if no cellwise coordinates are defined

Level: intermediate

Note:
This just retrieves the local section from the cell coordinate `DM`. In other words,
-vb
DMGetCellCoordinateDM(dm, &cdm);
DMGetLocalSection(cdm, &section);
-ve

-seealso: `DM`, `DMGetCoordinateSection()`, `DMSetCellCoordinateSection()`, `DMGetCellCoordinateDM()`, `DMGetCoordinateDM()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMGetCellCoordinateSection"))
"""
function DMGetCellCoordinateSection(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()

	LibPETSc.DMGetCellCoordinateSection(
		PetscLib,
		dm,
		section,
	)

	return section
end
 
 
"""
	 UNTESTED !!!
	 DMSetCellCoordinateSection(dm::AbstractDM{PetscLib},dim::Int,section::PetscSection)

Set the `PetscSection` of cellwise coordinate values over the mesh.

Not Collective

Input Parameters:
===
- `dm`      - The `DM` object
- `dim`     - The embedding dimension, or `PETSC_DETERMINE`
- `section` - The `PetscSection` object for a cellwise layout

Level: intermediate

-seealso: `DM`, `DMGetCoordinateDim()`, `DMSetCoordinateSection()`, `DMGetCellCoordinateSection()`, `DMGetCoordinateSection()`, `DMGetCellCoordinateDM()`, `DMGetLocalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMSetCellCoordinateSection"))
"""
function DMSetCellCoordinateSection(dm::AbstractDM{PetscLib},dim::Int,section::PetscSection) where {PetscLib}

	LibPETSc.DMSetCellCoordinateSection(
		PetscLib,
		dm,
		dim,
		section,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCoordinates(dm::AbstractDM{PetscLib})

Gets a global vector with the coordinates associated with the `DM`.

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - global coordinate vector

Level: intermediate

Notes:
This is a borrowed reference, so the user should NOT destroy this vector. When the `DM` is
destroyed `c` will no longer be valid.

Each process has only the locally-owned portion of the global coordinates (does NOT have the ghost coordinates).

For `DMDA`, in two and three dimensions coordinates are interlaced (x_0,y_0,x_1,y_1,...)
and (x_0,y_0,z_0,x_1,y_1,z_1...)

-seealso: `DM`, `DMDA`, `DMSetCoordinates()`, `DMGetCoordinatesLocal()`, `DMGetCoordinateDM()`, `DMDASetUniformCoordinates()`

# External Links
$(_doc_external("DM/DMGetCoordinates"))
"""
function DMGetCoordinates(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCoordinates(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinates(dm::AbstractDM{PetscLib},c::AbstractVector)

Sets into the `DM` a global vector that holds the coordinates

Collective

Input Parameters:
===
- `dm` - the `DM`
- `c`  - coordinate vector

Level: intermediate

Notes:
The coordinates do not include those for ghost points, which are in the local vector.

The vector `c` can be destroyed after the call

-seealso: `DM`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`, `DMGetCoordinateDM()`, `DMDASetUniformCoordinates()`

# External Links
$(_doc_external("DM/DMSetCoordinates"))
"""
function DMSetCoordinates(dm::AbstractDM{PetscLib},c::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetCoordinates(
		PetscLib,
		dm,
		c,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCellCoordinates(dm::AbstractDM{PetscLib})

Gets a global vector with the cellwise coordinates associated with the `DM`.

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - global coordinate vector

Level: intermediate

Notes:
This is a borrowed reference, so the user should NOT destroy this vector. When the `DM` is
destroyed `c` will no longer be valid.

Each process has only the locally-owned portion of the global coordinates (does NOT have the ghost coordinates).

-seealso: `DM`, `DMGetCoordinates()`, `DMSetCellCoordinates()`, `DMGetCellCoordinatesLocal()`, `DMGetCellCoordinateDM()`

# External Links
$(_doc_external("DM/DMGetCellCoordinates"))
"""
function DMGetCellCoordinates(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCellCoordinates(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	 DMSetCellCoordinates(dm::AbstractDM{PetscLib},c::AbstractVector)

Sets into the `DM` a global vector that holds the cellwise coordinates

Collective

Input Parameters:
===
- `dm` - the `DM`
- `c`  - cellwise coordinate vector

Level: intermediate

Notes:
The coordinates do not include those for ghost points, which are in the local vector.

The vector `c` should be destroyed by the caller.

-seealso: `DM`, `DMGetCoordinates()`, `DMSetCellCoordinatesLocal()`, `DMGetCellCoordinates()`, `DMGetCellCoordinatesLocal()`, `DMGetCellCoordinateDM()`

# External Links
$(_doc_external("DM/DMSetCellCoordinates"))
"""
function DMSetCellCoordinates(dm::AbstractDM{PetscLib},c::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetCellCoordinates(
		PetscLib,
		dm,
		c,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetCoordinatesLocalSetUp(dm::AbstractDM{PetscLib})

Prepares a local vector of coordinates, so that `DMGetCoordinatesLocalNoncollective()` can be used as non

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: advanced

-seealso: `DM`, `DMSetCoordinates()`, `DMGetCoordinatesLocalNoncollective()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocalSetUp"))
"""
function DMGetCoordinatesLocalSetUp(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMGetCoordinatesLocalSetUp(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCoordinatesLocal(dm::AbstractDM{PetscLib})

Gets a local vector with the coordinates associated with the `DM`.

Collective the first time it is called

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - coordinate vector

Level: intermediate

Notes:
This is a borrowed reference, so the user should NOT destroy `c`

Each process has the local and ghost coordinates

For `DMDA`, in two and three dimensions coordinates are interlaced (x_0,y_0,x_1,y_1,...)
and (x_0,y_0,z_0,x_1,y_1,z_1...)

-seealso: `DM`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMSetCoordinates()`, `DMGetCoordinateDM()`, `DMGetCoordinatesLocalNoncollective()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function DMGetCoordinatesLocal(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCoordinatesLocal(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCoordinatesLocalNoncollective(dm::AbstractDM{PetscLib})

Non

Not Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - coordinate vector

Level: advanced

Note:
A previous call to  `DMGetCoordinatesLocal()` or `DMGetCoordinatesLocalSetUp()` ensures that a call to this function will not error.

-seealso: `DM`, `DMGetCoordinatesLocalSetUp()`, `DMGetCoordinatesLocal()`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMSetCoordinates()`, `DMGetCoordinateDM()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocalNoncollective"))
"""
function DMGetCoordinatesLocalNoncollective(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCoordinatesLocalNoncollective(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	pCoordSection,pCoord = DMGetCoordinatesLocalTuple(dm::AbstractDM{PetscLib},p::IS)

Gets a local vector with the coordinates of specified points and the section describing its layout.

Not Collective

Input Parameters:
===
- `dm` - the `DM`
- `p`  - the `IS` of points whose coordinates will be returned

Output Parameters:
===
- `pCoordSection` - the `PetscSection` describing the layout of pCoord, i.e. each point corresponds to one point in `p`, and DOFs correspond to coordinates
- `pCoord`        - the `Vec` with coordinates of points in `p`

Level: advanced

Notes:
`DMGetCoordinatesLocalSetUp()` must be called first. This function employs `DMGetCoordinatesLocalNoncollective()` so it is not collective.

This creates a new vector, so the user SHOULD destroy this vector

Each process has the local and ghost coordinates

For `DMDA`, in two and three dimensions coordinates are interlaced (x_0,y_0,x_1,y_1,...)
and (x_0,y_0,z_0,x_1,y_1,z_1...)

-seealso: `DM`, `DMDA`, `DMSetCoordinatesLocal()`, `DMGetCoordinatesLocal()`, `DMGetCoordinatesLocalNoncollective()`, `DMGetCoordinatesLocalSetUp()`, `DMGetCoordinates()`, `DMSetCoordinates()`, `DMGetCoordinateDM()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocalTuple"))
"""
function DMGetCoordinatesLocalTuple(dm::AbstractDM{PetscLib},p::IS) where {PetscLib}
	pCoordSection = LibPETSc.PetscSection()
	pCoord = CVec()

	LibPETSc.DMGetCoordinatesLocalTuple(
		PetscLib,
		dm,
		p,
		pCoordSection,
		pCoord,
	)

	return pCoordSection,pCoord
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinatesLocal(dm::AbstractDM{PetscLib},c::AbstractVector)

Sets into the `DM` a local vector, including ghost points, that holds the coordinates

Not Collective

Input Parameters:
===
- `dm` - the `DM`
- `c`  - coordinate vector

Level: intermediate

Notes:
The coordinates of ghost points can be set using `DMSetCoordinates()`
followed by `DMGetCoordinatesLocal()`. This is intended to enable the
setting of ghost coordinates outside of the domain.

The vector `c` should be destroyed by the caller.

-seealso: `DM`, `DMGetCoordinatesLocal()`, `DMSetCoordinates()`, `DMGetCoordinates()`, `DMGetCoordinateDM()`

# External Links
$(_doc_external("DM/DMSetCoordinatesLocal"))
"""
function DMSetCoordinatesLocal(dm::AbstractDM{PetscLib},c::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetCoordinatesLocal(
		PetscLib,
		dm,
		c,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetCellCoordinatesLocalSetUp(dm::AbstractDM{PetscLib})

Prepares a local vector of cellwise coordinates, so that `DMGetCellCoordinatesLocalNoncollective()` can be used as non

Collective

Input Parameter:
===
- `dm` - the `DM`

Level: advanced

-seealso: `DM`, `DMGetCellCoordinatesLocalNoncollective()`

# External Links
$(_doc_external("DM/DMGetCellCoordinatesLocalSetUp"))
"""
function DMGetCellCoordinatesLocalSetUp(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMGetCellCoordinatesLocalSetUp(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCellCoordinatesLocal(dm::AbstractDM{PetscLib})

Gets a local vector with the cellwise coordinates associated with the `DM`.

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - coordinate vector

Level: intermediate

Notes:
This is a borrowed reference, so the user should NOT destroy this vector

Each process has the local and ghost coordinates

-seealso: `DM`, `DMSetCellCoordinatesLocal()`, `DMGetCellCoordinates()`, `DMSetCellCoordinates()`, `DMGetCellCoordinateDM()`, `DMGetCellCoordinatesLocalNoncollective()`

# External Links
$(_doc_external("DM/DMGetCellCoordinatesLocal"))
"""
function DMGetCellCoordinatesLocal(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCellCoordinatesLocal(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	c = DMGetCellCoordinatesLocalNoncollective(dm::AbstractDM{PetscLib})

Non

Not Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameter:
===
- `c` - cellwise coordinate vector

Level: advanced

-seealso: `DM`, `DMGetCellCoordinatesLocalSetUp()`, `DMGetCellCoordinatesLocal()`, `DMSetCellCoordinatesLocal()`, `DMGetCellCoordinates()`, `DMSetCellCoordinates()`, `DMGetCellCoordinateDM()`

# External Links
$(_doc_external("DM/DMGetCellCoordinatesLocalNoncollective"))
"""
function DMGetCellCoordinatesLocalNoncollective(dm::AbstractDM{PetscLib}) where {PetscLib}
	c = CVec()

	LibPETSc.DMGetCellCoordinatesLocalNoncollective(
		PetscLib,
		dm,
		c,
	)

	return c
end
 
 
"""
	 UNTESTED !!!
	 DMSetCellCoordinatesLocal(dm::AbstractDM{PetscLib},c::AbstractVector)

Sets into the `DM` a local vector including ghost points that holds the cellwise coordinates

Not Collective

Input Parameters:
===
- `dm` - the `DM`
- `c`  - cellwise coordinate vector

Level: intermediate

Notes:
The coordinates of ghost points can be set using `DMSetCoordinates()`
followed by `DMGetCoordinatesLocal()`. This is intended to enable the
setting of ghost coordinates outside of the domain.

The vector `c` should be destroyed by the caller.

-seealso: `DM`, `DMGetCellCoordinatesLocal()`, `DMSetCellCoordinates()`, `DMGetCellCoordinates()`, `DMGetCellCoordinateDM()`

# External Links
$(_doc_external("DM/DMSetCellCoordinatesLocal"))
"""
function DMSetCellCoordinatesLocal(dm::AbstractDM{PetscLib},c::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetCellCoordinatesLocal(
		PetscLib,
		dm,
		c,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	field = DMGetCoordinateField(dm::AbstractDM{PetscLib})


# External Links
$(_doc_external("DM/DMGetCoordinateField"))
"""
function DMGetCoordinateField(dm::AbstractDM{PetscLib}) where {PetscLib}
	field = LibPETSc.DMField()

	LibPETSc.DMGetCoordinateField(
		PetscLib,
		dm,
		field,
	)

	return field
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinateField(dm::AbstractDM{PetscLib},field::DMField)


# External Links
$(_doc_external("DM/DMSetCoordinateField"))
"""
function DMSetCoordinateField(dm::AbstractDM{PetscLib},field::DMField) where {PetscLib}

	LibPETSc.DMSetCoordinateField(
		PetscLib,
		dm,
		field,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetLocalBoundingBox(dm::AbstractDM{PetscLib},lmin::Vector{Float64},lmax::Vector{Float64})

Returns the bounding box for the piece of the `DM` on this process.

Not Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameters:
===
- `lmin` - local minimum coordinates (length coord dim, optional)
- `lmax` - local maximum coordinates (length coord dim, optional)

Level: beginner

Note:
If the `DM` is a `DMDA` and has no coordinates, the index bounds are returned instead.

-seealso: `DM`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`, `DMGetBoundingBox()`

# External Links
$(_doc_external("DM/DMGetLocalBoundingBox"))
"""
function DMGetLocalBoundingBox(dm::AbstractDM{PetscLib},lmin::Vector{Float64},lmax::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetLocalBoundingBox(
		PetscLib,
		dm,
		lmin,
		lmax,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetBoundingBox(dm::AbstractDM{PetscLib},gmin::Vector{Float64},gmax::Vector{Float64})

Returns the global bounding box for the `DM`.

Collective

Input Parameter:
===
- `dm` - the `DM`

Output Parameters:
===
- `gmin` - global minimum coordinates (length coord dim, optional)
- `gmax` - global maximum coordinates (length coord dim, optional)

Level: beginner

-seealso: `DM`, `DMGetLocalBoundingBox()`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`

# External Links
$(_doc_external("DM/DMGetBoundingBox"))
"""
function DMGetBoundingBox(dm::AbstractDM{PetscLib},gmin::Vector{Float64},gmax::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetBoundingBox(
		PetscLib,
		dm,
		gmin,
		gmax,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoordinateDisc(dm::AbstractDM{PetscLib},disc::PetscFE,project::PetscBool)

Set a coordinate space

Input Parameters:
===
- `dm`      - The `DM` object
- `disc`    - The new coordinate discretization or `NULL` to ensure a coordinate discretization exists
- `project` - Project coordinates to new discretization

Level: intermediate

Notes:
A `PetscFE` defines an approximation space using a `PetscSpace`, which represents the basis functions, and a `PetscDualSpace`, which defines the interpolation operation in the space.

This function takes the current mesh coordinates, which are discretized using some `PetscFE` space, and projects this function into a new `PetscFE` space.
The coordinate projection is done on the continuous coordinates, but the discontinuous coordinates are not updated.

Developer Note:
With more effort, we could directly project the discontinuous coordinates also.

-seealso: `DM`, `PetscFE`, `DMGetCoordinateField()`

# External Links
$(_doc_external("DM/DMSetCoordinateDisc"))
"""
function DMSetCoordinateDisc(dm::AbstractDM{PetscLib},disc::PetscFE,project::PetscBool) where {PetscLib}

	LibPETSc.DMSetCoordinateDisc(
		PetscLib,
		dm,
		disc,
		project,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	cellSF = DMLocatePoints(dm::AbstractDM{PetscLib},v::AbstractVector,ltype::DMPointLocationType)

Locate the points in `v` in the mesh and return a `PetscSF` of the containing cells

Collective

Input Parameters:
===
- `dm`    - The `DM`
- `ltype` - The type of point location, e.g. `DM_POINTLOCATION_NONE` or `DM_POINTLOCATION_NEAREST`

Input/Output Parameters:
- `v`      - The `Vec` of points, on output contains the nearest mesh points to the given points if `DM_POINTLOCATION_NEAREST` is used
- `cellSF` - Points to either `NULL`, or a `PetscSF` with guesses for which cells contain each point;
on output, the `PetscSF` containing the MPI ranks and local indices of the containing points

Level: developer

Notes:
To do a search of the local cells of the mesh, `v` should have `PETSC_COMM_SELF` as its communicator.
To do a search of all the cells in the distributed mesh, `v` should have the same MPI communicator as `dm`.

Points will only be located in owned cells, not overlap cells arising from `DMPlexDistribute()` or other overlapping distributions.

If *cellSF is `NULL` on input, a `PetscSF` will be created.
If *cellSF is not `NULL` on input, it should point to an existing `PetscSF`, whose graph will be used as initial guesses.

An array that maps each point to its containing cell can be obtained with
-vb
const PetscSFNode *cells;
PetscInt           nFound;
const PetscInt    *found;

PetscSFGetGraph(cellSF,NULL,&nFound,&found,&cells);
-ve

Where cells[i].rank is the MPI rank of the process owning the cell containing point found[i] (or i if found == NULL), and cells[i].index is
the index of the cell in its MPI process' local numbering. This rank is in the communicator for `v`, so if `v` is on `PETSC_COMM_SELF` then the rank will always be 0.

-seealso: `DM`, `DMSetCoordinates()`, `DMSetCoordinatesLocal()`, `DMGetCoordinates()`, `DMGetCoordinatesLocal()`, `DMPointLocationType`

# External Links
$(_doc_external("DM/DMLocatePoints"))
"""
function DMLocatePoints(dm::AbstractDM{PetscLib},v::AbstractVector,ltype::DMPointLocationType) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	cellSF = LibPETSc.PetscSF()

	LibPETSc.DMLocatePoints(
		PetscLib,
		dm,
		v,
		ltype,
		cellSF,
	)

	return cellSF
end
 
 
"""
	 UNTESTED !!!
	 DMGetPeriodicity(dm::AbstractDM{PetscLib},maxCell::Vector{Float64},Lstart::Vector{Float64},L::Vector{Float64})

Get the description of mesh periodicity

Input Parameter:
===
- `dm` - The `DM` object

Output Parameters:
===
- `maxCell` - Over distances greater than this, we can assume a point has crossed over to another sheet, when trying to localize cell coordinates
- `Lstart`  - If we assume the mesh is a torus, this is the start of each coordinate, or `NULL` for 0.0
- `L`       - If we assume the mesh is a torus, this is the length of each coordinate, otherwise it is < 0.0

Level: developer

-seealso: `DM`

# External Links
$(_doc_external("DM/DMGetPeriodicity"))
"""
function DMGetPeriodicity(dm::AbstractDM{PetscLib},maxCell::Vector{Float64},Lstart::Vector{Float64},L::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetPeriodicity(
		PetscLib,
		dm,
		maxCell,
		Lstart,
		L,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetPeriodicity(dm::AbstractDM{PetscLib},maxCell::Vector{Float64},Lstart::Vector{Float64},L::Vector{Float64})

Set the description of mesh periodicity

Input Parameters:
===
- `dm`      - The `DM` object
- `maxCell` - Over distances greater than this, we can assume a point has crossed over to another sheet, when trying to localize cell coordinates. Pass `NULL` to remove such information.
- `Lstart`  - If we assume the mesh is a torus, this is the start of each coordinate, or `NULL` for 0.0
- `L`       - If we assume the mesh is a torus, this is the length of each coordinate, otherwise it is < 0.0

Level: developer

-seealso: `DM`, `DMGetPeriodicity()`

# External Links
$(_doc_external("DM/DMSetPeriodicity"))
"""
function DMSetPeriodicity(dm::AbstractDM{PetscLib},maxCell::Vector{Float64},Lstart::Vector{Float64},L::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetPeriodicity(
		PetscLib,
		dm,
		maxCell,
		Lstart,
		L,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalizeCoordinate(dm::AbstractDM{PetscLib},in::Vector{Float64},endpoint::PetscBool,out::Vector{Float64})

If a mesh is periodic (a torus with lengths L_i, some of which can be infinite), project the coordinate onto [0, L_i) in each dimension.

Input Parameters:
===
- `dm`       - The `DM`
- `in`       - The input coordinate point (dim numbers)
- `endpoint` - Include the endpoint L_i

Output Parameter:
===
- `out` - The localized coordinate point (dim numbers)

Level: developer

-seealso: `DM`, `DMLocalizeCoordinates()`, `DMLocalizeAddCoordinate()`

# External Links
$(_doc_external("DM/DMLocalizeCoordinate"))
"""
function DMLocalizeCoordinate(dm::AbstractDM{PetscLib},in::Vector{Float64},endpoint::PetscBool,out::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMLocalizeCoordinate(
		PetscLib,
		dm,
		in,
		endpoint,
		out,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMLocalizeCoordinates(dm::AbstractDM{PetscLib})

If a mesh is periodic, create local coordinates for cells having periodic faces

Collective

Input Parameter:
===
- `dm` - The `DM`

Level: developer

-seealso: `DM`, `DMSetPeriodicity()`, `DMLocalizeCoordinate()`, `DMLocalizeAddCoordinate()`

# External Links
$(_doc_external("DM/DMLocalizeCoordinates"))
"""
function DMLocalizeCoordinates(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMLocalizeCoordinates(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	areLocalized = DMGetCoordinatesLocalized(dm::AbstractDM{PetscLib})

Check if the `DM` coordinates have been localized for cells

Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `areLocalized` - `PETSC_TRUE` if localized

Level: developer

-seealso: `DM`, `DMLocalizeCoordinates()`, `DMSetPeriodicity()`, `DMGetCoordinatesLocalizedLocal()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocalized"))
"""
function DMGetCoordinatesLocalized(dm::AbstractDM{PetscLib}) where {PetscLib}
	areLocalized = Ref{PetscBool}()

	LibPETSc.DMGetCoordinatesLocalized(
		PetscLib,
		dm,
		areLocalized,
	)

	return areLocalized[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	areLocalized = DMGetCoordinatesLocalizedLocal(dm::AbstractDM{PetscLib})

Check if the `DM` coordinates have been localized for cells on this process

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `areLocalized` - `PETSC_TRUE` if localized

Level: developer

-seealso: `DM`, `DMLocalizeCoordinates()`, `DMGetCoordinatesLocalized()`, `DMSetPeriodicity()`

# External Links
$(_doc_external("DM/DMGetCoordinatesLocalizedLocal"))
"""
function DMGetCoordinatesLocalizedLocal(dm::AbstractDM{PetscLib}) where {PetscLib}
	areLocalized = Ref{PetscBool}()

	LibPETSc.DMGetCoordinatesLocalizedLocal(
		PetscLib,
		dm,
		areLocalized,
	)

	return areLocalized[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMSubDomainRestrict(global::AbstractDM{PetscLib},oscatter::VecScatter,gscatter::VecScatter,subdm::AbstractDM{PetscLib})

restricts user

Collective if any hooks are

Input Parameters:
===
- `global`   - The global `DM` to use as a base
- `oscatter` - The scatter from domain global vector filling subdomain global vector with overlap
- `gscatter` - The scatter from domain global vector filling subdomain local vector with ghosts
- `subdm`    - The subdomain `DM` to update

Level: developer

See also: 
=== 
`DM`, `DMCoarsenHookAdd()`, `MatRestrict()`, `DMCreateDomainDecomposition()`

# External Links
$(_doc_external("DM/DMSubDomainRestrict"))
"""
function DMSubDomainRestrict(v_global::AbstractDM{PetscLib},oscatter::VecScatter,gscatter::VecScatter,subdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMSubDomainRestrict(
		PetscLib,
		v_global,
		oscatter,
		gscatter,
		subdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char})

Sets the prefix prepended to all option names when searching through the options database

Logically Collective

Input Parameters:
===
- `dm`     - the `DM` context
- `prefix` - the prefix to prepend

Level: advanced

Note:
A hyphen (-) must NOT be given at the beginning of the prefix name.
The first character of all runtime options is AUTOMATICALLY the hyphen.

See also: 
=== 
`DM`, `PetscObjectSetOptionsPrefix()`, `DMSetFromOptions()`

# External Links
$(_doc_external("DM/DMSetOptionsPrefix"))
"""
function DMSetOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetOptionsPrefix(
		PetscLib,
		dm,
		prefix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMAppendOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char})

Appends an additional string to an already existing prefix used for searching for
`DM` options in the options database.

Logically Collective

Input Parameters:
===
- `dm`     - the `DM` context
- `prefix` - the string to append to the current prefix

Level: advanced

Note:
If the `DM` does not currently have an options prefix then this value is used alone as the prefix as if `DMSetOptionsPrefix()` had been called.
A hyphen (-) must NOT be given at the beginning of the prefix name.
The first character of all runtime options is AUTOMATICALLY the hyphen.

See also: 
=== 
`DM`, `DMSetOptionsPrefix()`, `DMGetOptionsPrefix()`, `PetscObjectAppendOptionsPrefix()`, `DMSetFromOptions()`

# External Links
$(_doc_external("DM/DMAppendOptionsPrefix"))
"""
function DMAppendOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMAppendOptionsPrefix(
		PetscLib,
		dm,
		prefix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char})

Gets the prefix used for searching for all
DM options in the options database.

Not Collective

Input Parameter:
===
- `dm` - the `DM` context

Output Parameter:
===
- `prefix` - pointer to the prefix string used is returned

Level: advanced

Fortran Note:
Pass in a string 'prefix' of
sufficient length to hold the prefix.

See also: 
=== 
`DM`, `DMSetOptionsPrefix()`, `DMAppendOptionsPrefix()`, `DMSetFromOptions()`

# External Links
$(_doc_external("DM/DMGetOptionsPrefix"))
"""
function DMGetOptionsPrefix(dm::AbstractDM{PetscLib},prefix::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetOptionsPrefix(
		PetscLib,
		dm,
		prefix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetVecType(dm::AbstractDM{PetscLib},ctype::VecType)

Sets the type of vector to be created with `DMCreateLocalVector()` and `DMCreateGlobalVector()`

Logically Collective

Input Parameters:
===
- `dm`    - initial distributed array
- `ctype` - the vector type, for example `VECSTANDARD`, `VECCUDA`, or `VECVIENNACL`

Options Database Key:
===
- `-dm_vec_type ctype` - the type of vector to create

Level: intermediate

See also: 
=== 
`DM`, `DMCreate()`, `DMDestroy()`, `DMDAInterpolationType`, `VecType`, `DMGetVecType()`, `DMSetMatType()`, `DMGetMatType()`,
`VECSTANDARD`, `VECCUDA`, `VECVIENNACL`, `DMCreateLocalVector()`, `DMCreateGlobalVector()`

# External Links
$(_doc_external("DM/DMSetVecType"))
"""
function DMSetVecType(dm::AbstractDM{PetscLib},ctype::VecType) where {PetscLib}

	LibPETSc.DMSetVecType(
		PetscLib,
		dm,
		ctype,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctype = DMGetVecType(da::AbstractDM{PetscLib})

Gets the type of vector created with `DMCreateLocalVector()` and `DMCreateGlobalVector()`

Logically Collective

Input Parameter:
===
- `da` - initial distributed array

Output Parameter:
===
- `ctype` - the vector type

Level: intermediate

See also: 
=== 
`DM`, `DMCreate()`, `DMDestroy()`, `DMDAInterpolationType`, `VecType`, `DMSetMatType()`, `DMGetMatType()`, `DMSetVecType()`

# External Links
$(_doc_external("DM/DMGetVecType"))
"""
function DMGetVecType(da::AbstractDM{PetscLib}) where {PetscLib}
	ctype = Ref{VecType}()

	LibPETSc.DMGetVecType(
		PetscLib,
		da,
		r_ctype,
	)


	ctype = unsafe_string(r_ctype[])
	return ctype
end
 
 
"""
	 UNTESTED !!!
	 DMSetMatType(dm::AbstractDM{PetscLib},ctype::MatType)

Sets the type of matrix created with `DMCreateMatrix()`

Logically Collective

Input Parameters:
===
- `dm`    - the `DM` context
- `ctype` - the matrix type, for example `MATMPIAIJ`

Options Database Key:
===
- `-dm_mat_type ctype` - the type of the matrix to create, for example mpiaij

Level: intermediate

See also: 
=== 
`DM`, `MatType`, `DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatrixPreallocateOnly()`, `DMGetMatType()`, `DMCreateGlobalVector()`, `DMCreateLocalVector()`

# External Links
$(_doc_external("DM/DMSetMatType"))
"""
function DMSetMatType(dm::AbstractDM{PetscLib},ctype::MatType) where {PetscLib}

	LibPETSc.DMSetMatType(
		PetscLib,
		dm,
		ctype,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctype = DMGetMatType(dm::AbstractDM{PetscLib})

Gets the type of matrix that would be created with `DMCreateMatrix()`

Logically Collective

Input Parameter:
===
- `dm` - the `DM` context

Output Parameter:
===
- `ctype` - the matrix type

Level: intermediate

See also: 
=== 
`DM`, `DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatrixPreallocateOnly()`, `MatType`, `DMSetMatType()`

# External Links
$(_doc_external("DM/DMGetMatType"))
"""
function DMGetMatType(dm::AbstractDM{PetscLib}) where {PetscLib}
	ctype = Ref{MatType}()

	LibPETSc.DMGetMatType(
		PetscLib,
		dm,
		r_ctype,
	)


	return ctype
end
 
 
"""
	 UNTESTED !!!
	 DMSetISColoringType(dm::AbstractDM{PetscLib},ctype::ISColoringType)

Sets the type of coloring, `IS_COLORING_GLOBAL` or `IS_COLORING_LOCAL` that is created by the `DM`

Logically Collective

Input Parameters:
===
- `dm`    - the `DM` context
- `ctype` - the matrix type

Options Database Key:
===
- `-dm_is_coloring_type` - global or local

Level: intermediate

See also: 
=== 
`DM`, `DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatrixPreallocateOnly()`, `MatType`, `DMGetMatType()`,
`DMGetISColoringType()`, `ISColoringType`, `IS_COLORING_GLOBAL`, `IS_COLORING_LOCAL`

# External Links
$(_doc_external("DM/DMSetISColoringType"))
"""
function DMSetISColoringType(dm::AbstractDM{PetscLib},ctype::ISColoringType) where {PetscLib}

	LibPETSc.DMSetISColoringType(
		PetscLib,
		dm,
		ctype,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctype = DMGetISColoringType(dm::AbstractDM{PetscLib})

Gets the type of coloring, `IS_COLORING_GLOBAL` or `IS_COLORING_LOCAL` that is created by the `DM`

Logically Collective

Input Parameter:
===
- `dm` - the `DM` context

Output Parameter:
===
- `ctype` - the matrix type

Options Database Key:
===
- `-dm_is_coloring_type` - global or local

Level: intermediate

See also: 
=== 
`DM`, `DMDACreate1d()`, `DMDACreate2d()`, `DMDACreate3d()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMSetMatrixPreallocateOnly()`, `MatType`, `DMGetMatType()`,
`ISColoringType`, `IS_COLORING_GLOBAL`, `IS_COLORING_LOCAL`

# External Links
$(_doc_external("DM/DMGetISColoringType"))
"""
function DMGetISColoringType(dm::AbstractDM{PetscLib}) where {PetscLib}
	ctype = Ref{ISColoringType}()

	LibPETSc.DMGetISColoringType(
		PetscLib,
		dm,
		ctype,
	)

	return ctype
end
 
 
"""
	 UNTESTED !!!
	ctx = DMGetApplicationContext(dm::AbstractDM{PetscLib})

Gets a user context from a `DM` object

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `ctx` - the user context

Level: intermediate

Note:
A user context is a way to pass problem specific information that is accessible whenever the `DM` is available

See also: 
=== 
`DM`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`

# External Links
$(_doc_external("DM/DMGetApplicationContext"))
"""
function DMGetApplicationContext(dm::AbstractDM{PetscLib}) where {PetscLib}
	Float64 = PetscLib.Float64
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, Float64)

	LibPETSc.DMGetApplicationContext(
		PetscLib,
		dm,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 
"""
	 UNTESTED !!!
	flg = DMHasVariableBounds(dm::AbstractDM{PetscLib})

does the `DM` object have a variable bounds function?

Not Collective

Input Parameter:
===
- `dm` - the `DM` object to destroy

Output Parameter:
===
- `flg` - `PETSC_TRUE` if the variable bounds function exists

Level: developer

See also: 
=== 
`DM`, `DMComputeVariableBounds()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMGetApplicationContext()`

# External Links
$(_doc_external("DM/DMHasVariableBounds"))
"""
function DMHasVariableBounds(dm::AbstractDM{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.DMHasVariableBounds(
		PetscLib,
		dm,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	flg = DMHasColoring(dm::AbstractDM{PetscLib})

does the `DM` object have a method of providing a coloring?

Not Collective

Input Parameter:
===
- `dm` - the DM object

Output Parameter:
===
- `flg` - `PETSC_TRUE` if the `DM` has facilities for `DMCreateColoring()`.

Level: developer

See also: 
=== 
`DM`, `DMCreateColoring()`

# External Links
$(_doc_external("DM/DMHasColoring"))
"""
function DMHasColoring(dm::AbstractDM{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.DMHasColoring(
		PetscLib,
		dm,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	flg = DMHasCreateRestriction(dm::AbstractDM{PetscLib})

does the `DM` object have a method of providing a restriction?

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `flg` - `PETSC_TRUE` if the `DM` has facilities for `DMCreateRestriction()`.

Level: developer

See also: 
=== 
`DM`, `DMCreateRestriction()`, `DMHasCreateInterpolation()`, `DMHasCreateInjection()`

# External Links
$(_doc_external("DM/DMHasCreateRestriction"))
"""
function DMHasCreateRestriction(dm::AbstractDM{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.DMHasCreateRestriction(
		PetscLib,
		dm,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	flg = DMHasCreateInjection(dm::AbstractDM{PetscLib})

does the `DM` object have a method of providing an injection?

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `flg` - `PETSC_TRUE` if the `DM` has facilities for `DMCreateInjection()`.

Level: developer

See also: 
=== 
`DM`, `DMCreateInjection()`, `DMHasCreateRestriction()`, `DMHasCreateInterpolation()`

# External Links
$(_doc_external("DM/DMHasCreateInjection"))
"""
function DMHasCreateInjection(dm::AbstractDM{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.DMHasCreateInjection(
		PetscLib,
		dm,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMComputeVariableBounds(dm::AbstractDM{PetscLib},xl::AbstractVector,xu::AbstractVector)

compute variable bounds used by `SNESVI`.

Logically Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameters:
===
- `xl` - lower bound
- `xu` - upper bound

Level: advanced

Note:
This is generally not called by users. It calls the function provided by the user with DMSetVariableBounds()

See also: 
=== 
`DM`, `DMHasVariableBounds()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMGetApplicationContext()`

# External Links
$(_doc_external("DM/DMComputeVariableBounds"))
"""
function DMComputeVariableBounds(dm::AbstractDM{PetscLib},xl::AbstractVector,xu::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMComputeVariableBounds(
		PetscLib,
		dm,
		xl,
		xu,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	is,subdm = DMCreateSubDM(dm::AbstractDM{PetscLib},numFields::Int,fields::Vector{Int})

Returns an `IS` and `DM` encapsulating a subproblem defined by the fields passed in.
The fields are defined by `DMCreateFieldIS()`.

Not collective

Input Parameters:
===
- `dm`        - The `DM` object
- `numFields` - The number of fields to select
- `fields`    - The field numbers of the selected fields

Output Parameters:
===
- `is`    - The global indices for all the degrees of freedom in the new sub `DM`, use `NULL` if not needed
- `subdm` - The `DM` for the subproblem, use `NULL` if not needed

Level: intermediate

Note:
You need to call `DMPlexSetMigrationSF()` on the original `DM` if you want the Global-To-Natural map to be automatically constructed

See also: 
=== 
`DM`, `DMCreateFieldIS()`, `DMCreateFieldDecomposition()`, `DMAddField()`, `DMCreateSuperDM()`, `IS`, `DMPlexSetMigrationSF()`, `DMDestroy()`, `DMView()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`

# External Links
$(_doc_external("DM/DMCreateSubDM"))
"""
function DMCreateSubDM(dm::AbstractDM{PetscLib},numFields::Int,fields::Vector{Int}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	is = LibPETSc.IS()
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	subdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMCreateSubDM(
		PetscLib,
		dm,
		numFields,
		fields,
		is,
		subdm,
	)

	return is,subdm
end
 
 
"""
	 UNTESTED !!!
	is,subdm = DMCreateSectionSubDM(dm::AbstractDM{PetscLib},numFields::Int,fields::Vector{Int},numComps::Vector{Int},comps::Vector{Int})

Returns an `IS` and `subDM` containing a `PetscSection` that encapsulates a subproblem defined by a subset of the fields in a `PetscSection` in the `DM`.

Not Collective

Input Parameters:
===
- `dm`        - The `DM` object
- `numFields` - The number of fields to incorporate into `subdm`
- `fields`    - The field numbers of the selected fields
- `numComps`  - The number of components from each field to incorporate into `subdm`, or PETSC_DECIDE for all components
- `comps`     - The component numbers of the selected fields (omitted for PTESC_DECIDE fields)

Output Parameters:
===
- `is`    - The global indices for the subproblem or `NULL`
- `subdm` - The `DM` for the subproblem, which must already have be cloned from `dm` or `NULL`

Level: intermediate

Notes:
If `is` and `subdm` are both `NULL` this does nothing

-seealso: `DMCreateSubDM()`, `DMGetLocalSection()`, `DMPlexSetMigrationSF()`, `DMView()`

# External Links
$(_doc_external("DM/DMCreateSectionSubDM"))
"""
function DMCreateSectionSubDM(dm::AbstractDM{PetscLib},numFields::Int,fields::Vector{Int},numComps::Vector{Int},comps::Vector{Int}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	is = LibPETSc.IS()
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	subdm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMCreateSectionSubDM(
		PetscLib,
		dm,
		numFields,
		fields,
		numComps,
		comps,
		is,
		subdm,
	)

	return is,subdm
end
 
 
"""
	 UNTESTED !!!
	level = DMGetRefineLevel(dm::AbstractDM{PetscLib})

Gets the number of refinements that have generated this `DM` from some initial `DM`.

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `level` - number of refinements

Level: developer

Note:
This can be used, by example, to set the number of coarser levels associated with this `DM` for a multigrid solver.

See also: 
=== 
`DM`, `DMRefine()`, `DMCoarsen()`, `DMGetCoarsenLevel()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMGetRefineLevel"))
"""
function DMGetRefineLevel(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	level = [PetscInt(1)]

	LibPETSc.DMGetRefineLevel(
		PetscLib,
		dm,
		Ref(level,1),
	)

	return level[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetRefineLevel(dm::AbstractDM{PetscLib},level::Int)

Sets the number of refinements that have generated this `DM`.

Not Collective

Input Parameters:
===
- `dm`    - the `DM` object
- `level` - number of refinements

Level: advanced

Notes:
This value is used by `PCMG` to determine how many multigrid levels to use

The values are usually set automatically by the process that is causing the refinements of an initial `DM` by calling this routine.

See also: 
=== 
`DM`, `DMGetRefineLevel()`, `DMCoarsen()`, `DMGetCoarsenLevel()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMSetRefineLevel"))
"""
function DMSetRefineLevel(dm::AbstractDM{PetscLib},level::Int) where {PetscLib}

	LibPETSc.DMSetRefineLevel(
		PetscLib,
		dm,
		level,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	level = DMGetCoarsenLevel(dm::AbstractDM{PetscLib})

Gets the number of coarsenings that have generated this `DM`.

Not Collective

Input Parameter:
===
- `dm` - the `DM` object

Output Parameter:
===
- `level` - number of coarsenings

Level: developer

See also: 
=== 
`DM`, `DMCoarsen()`, `DMSetCoarsenLevel()`, `DMGetRefineLevel()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMGetCoarsenLevel"))
"""
function DMGetCoarsenLevel(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	level = [PetscInt(1)]

	LibPETSc.DMGetCoarsenLevel(
		PetscLib,
		dm,
		Ref(level,1),
	)

	return level[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetCoarsenLevel(dm::AbstractDM{PetscLib},level::Int)

Sets the number of coarsenings that have generated this `DM`.

Collective

Input Parameters:
===
- `dm`    - the `DM` object
- `level` - number of coarsenings

Level: developer

Note:
This is rarely used directly, the information is automatically set when a `DM` is created with `DMCoarsen()`

See also: 
=== 
`DM`, `DMCoarsen()`, `DMGetCoarsenLevel()`, `DMGetRefineLevel()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateInterpolation()`

# External Links
$(_doc_external("DM/DMSetCoarsenLevel"))
"""
function DMSetCoarsenLevel(dm::AbstractDM{PetscLib},level::Int) where {PetscLib}

	LibPETSc.DMSetCoarsenLevel(
		PetscLib,
		dm,
		level,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dm = VecGetDM(v::AbstractVector)

Gets the `DM` defining the data layout of the vector

Not Collective

Input Parameter:
===
- `v` - The `Vec`

Output Parameter:
===
- `dm` - The `DM`

Level: intermediate

Note:
A `Vec` may not have a `DM` associated with it.

See also: 
=== 
`DM`, `VecSetDM()`, `DMGetLocalVector()`, `DMGetGlobalVector()`, `DMSetVecType()`

# External Links
$(_doc_external("DM/VecGetDM"))
"""
function VecGetDM(v::AbstractVector{PetscLib}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.VecGetDM(
		PetscLib,
		v,
		dm,
	)

	return dm
end
 
 
"""
	 UNTESTED !!!
	 VecSetDM(v::AbstractVector,dm::AbstractDM{PetscLib})

Sets the `DM` defining the data layout of the vector.

Not Collective

Input Parameters:
===
- `v`  - The `Vec`
- `dm` - The `DM`

Level: developer

Notes:
This is rarely used, generally one uses `DMGetLocalVector()` or  `DMGetGlobalVector()` to create a vector associated with a given `DM`

This is NOT the same as `DMCreateGlobalVector()` since it does not change the view methods or perform other customization, but merely sets the `DM` member.

See also: 
=== 
`DM`, `VecGetDM()`, `DMGetLocalVector()`, `DMGetGlobalVector()`, `DMSetVecType()`

# External Links
$(_doc_external("DM/VecSetDM"))
"""
function VecSetDM(v::AbstractVector,dm::AbstractDM{PetscLib}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.VecSetDM(
		PetscLib,
		v,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	dm = MatGetDM(A::AbstractMatrix)

Gets the `DM` defining the data layout of the matrix

Not Collective

Input Parameter:
===
- `A` - The `Mat`

Output Parameter:
===
- `dm` - The `DM`

Level: intermediate

Note:
A matrix may not have a `DM` associated with it

Developer Note:
Since the `Mat` class doesn't know about the `DM` class the `DM` object is associated with the `Mat` through a `PetscObjectCompose()` operation

See also: 
=== 
`DM`, `MatSetDM()`, `DMCreateMatrix()`, `DMSetMatType()`

# External Links
$(_doc_external("DM/MatGetDM"))
"""
function MatGetDM(A::AbstractMatrix{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.MatGetDM(
		PetscLib,
		A,
		dm,
	)

	return dm
end
 
 
"""
	 UNTESTED !!!
	 MatSetDM(A::AbstractMatrix,dm::AbstractDM{PetscLib})

Sets the `DM` defining the data layout of the matrix

Not Collective

Input Parameters:
===
- `A`  - The `Mat`
- `dm` - The `DM`

Level: developer

Note:
This is rarely used in practice, rather `DMCreateMatrix()` is used to create a matrix associated with a particular `DM`

Developer Note:
Since the `Mat` class doesn't know about the `DM` class the `DM` object is associated with
the `Mat` through a `PetscObjectCompose()` operation

See also: 
=== 
`DM`, `MatGetDM()`, `DMCreateMatrix()`, `DMSetMatType()`

# External Links
$(_doc_external("DM/MatSetDM"))
"""
function MatSetDM(A::AbstractMatrix,dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.MatSetDM(
		PetscLib,
		A,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 MatFDColoringUseDM(coloring::AbstractMatrix,fdcoloring::MatFDColoring)

allows a `MatFDColoring` object to use the `DM` associated with the matrix to compute a `IS_COLORING_LOCAL` coloring

Input Parameters:
===
- `coloring`   - The matrix to get the `DM` from
- `fdcoloring` - the `MatFDColoring` object

Level: advanced

Developer Note:
This routine exists because the PETSc `Mat` library does not know about the `DM` objects

See also: 
=== 
`DM`, `MatFDColoring`, `MatFDColoringCreate()`, `ISColoringType`

# External Links
$(_doc_external("DM/MatFDColoringUseDM"))
"""
function MatFDColoringUseDM(coloring::AbstractMatrix{PetscLib},fdcoloring::MatFDColoring) where {PetscLib}

	LibPETSc.MatFDColoringUseDM(
		PetscLib,
		coloring,
		fdcoloring,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMPrintCellIndices(PetscLib, c::Int,name::Vector{Char},len::Int,x::Vector{Int})


# External Links
$(_doc_external("DM/DMPrintCellIndices"))
"""
function DMPrintCellIndices(PetscLib, c::Int,name::Vector{Char},len::Int,x::Vector{Int}) 
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMPrintCellIndices(
		PetscLib,
		c,
		name,
		len,
		x,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMPrintCellVector(c::Int,name::Vector{Char},len::Int,x::Vector{Float64})


# External Links
$(_doc_external("DM/DMPrintCellVector"))
"""
function DMPrintCellVector(PetscLib, c::Int,name::Vector{Char},len::Int,x::Vector{Float64})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMPrintCellVector(
		PetscLib,
		c,
		name,
		len,
		x,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMPrintCellVectorReal(c::Int,name::Vector{Char},len::Int,x::Vector{Float64})


# External Links
$(_doc_external("DM/DMPrintCellVectorReal"))
"""
function DMPrintCellVectorReal(PetscLib, c::Int,name::Vector{Char},len::Int,x::Vector{Float64}) 
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMPrintCellVectorReal(
		PetscLib,
		c,
		name,
		len,
		x,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMPrintCellMatrix(c::Int,name::Vector{Char},rows::Int,cols::Int,A::Vector{Float64})


# External Links
$(_doc_external("DM/DMPrintCellMatrix"))
"""
function DMPrintCellMatrix(PetscLib, c::Int,name::Vector{Char},rows::Int,cols::Int,A::Vector{Float64})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMPrintCellMatrix(
		PetscLib,
		c,
		name,
		rows,
		cols,
		A,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMPrintLocalVec(dm::AbstractDM{PetscLib},name::Vector{Char},tol<:AbstractFloat,X::AbstractVector)


# External Links
$(_doc_external("DM/DMPrintLocalVec"))
"""
function DMPrintLocalVec(dm::AbstractDM{PetscLib},name::Vector{Char},tol::Float64,X::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMPrintLocalVec(
		PetscLib,
		dm,
		name,
		tol,
		X,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section = DMGetSection(dm::AbstractDM{PetscLib})

Get the `PetscSection` encoding the local data layout for the `DM`.   This is equivalent to `DMGetLocalSection()`. Deprecated in v3.12

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `section` - The `PetscSection`

Options Database Key:
===
- `-dm_petscsection_view` - View the `PetscSection` created by the `DM`

Level: advanced

Notes:
Use `DMGetLocalSection()` in new code.

This gets a borrowed reference, so the user should not destroy this `PetscSection`.

See also: 
=== 
`DM`, `DMGetLocalSection()`, `DMSetLocalSection()`, `DMGetGlobalSection()`

# External Links
$(_doc_external("DM/DMGetSection"))
"""
function DMGetSection(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()

	LibPETSc.DMGetSection(
		PetscLib,
		dm,
		section,
	)

	return section
end
 
 
"""
	 UNTESTED !!!
	 DMSetSection(dm::AbstractDM{PetscLib},section::PetscSection)

Set the `PetscSection` encoding the local data layout for the `DM`.  This is equivalent to `DMSetLocalSection()`. Deprecated in v3.12

Input Parameters:
===
- `dm`      - The `DM`
- `section` - The `PetscSection`

Level: advanced

Notes:
Use `DMSetLocalSection()` in new code.

Any existing `PetscSection` will be destroyed

See also: 
=== 
`DM`, `DMSetLocalSection()`, `DMGetLocalSection()`, `DMSetGlobalSection()`

# External Links
$(_doc_external("DM/DMSetSection"))
"""
function DMSetSection(dm::AbstractDM{PetscLib},section::PetscSection) where {PetscLib}

	LibPETSc.DMSetSection(
		PetscLib,
		dm,
		section,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section = DMGetLocalSection(dm::AbstractDM{PetscLib})

Get the `PetscSection` encoding the local data layout for the `DM`.

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `section` - The `PetscSection`

Options Database Key:
===
- `-dm_petscsection_view` - View the section created by the `DM`

Level: intermediate

Note:
This gets a borrowed reference, so the user should not destroy this `PetscSection`.

See also: 
=== 
`DM`, `DMSetLocalSection()`, `DMGetGlobalSection()`

# External Links
$(_doc_external("DM/DMGetLocalSection"))
"""
function DMGetLocalSection(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()

	LibPETSc.DMGetLocalSection(
		PetscLib,
		dm,
		section,
	)

	return section
end
 
 
"""
	 UNTESTED !!!
	 DMSetLocalSection(dm::AbstractDM{PetscLib},section::PetscSection)

Set the `PetscSection` encoding the local data layout for the `DM`.

Input Parameters:
===
- `dm`      - The `DM`
- `section` - The `PetscSection`

Level: intermediate

Note:
Any existing Section will be destroyed

See also: 
=== 
`DM`, `PetscSection`, `DMGetLocalSection()`, `DMSetGlobalSection()`

# External Links
$(_doc_external("DM/DMSetLocalSection"))
"""
function DMSetLocalSection(dm::AbstractDM{PetscLib},section::PetscSection) where {PetscLib}

	LibPETSc.DMSetLocalSection(
		PetscLib,
		dm,
		section,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section = DMGetGlobalSection(dm::AbstractDM{PetscLib})

Get the `PetscSection` encoding the global data layout for the `DM`.

Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `section` - The `PetscSection`

Level: intermediate

Note:
This gets a borrowed reference, so the user should not destroy this `PetscSection`.

See also: 
=== 
`DM`, `DMSetLocalSection()`, `DMGetLocalSection()`

# External Links
$(_doc_external("DM/DMGetGlobalSection"))
"""
function DMGetGlobalSection(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()

	LibPETSc.DMGetGlobalSection(
		PetscLib,
		dm,
		section,
	)

	return section
end
 
 
"""
	 UNTESTED !!!
	 DMSetGlobalSection(dm::AbstractDM{PetscLib},section::PetscSection)

Set the `PetscSection` encoding the global data layout for the `DM`.

Input Parameters:
===
- `dm`      - The `DM`
- `section` - The PetscSection, or `NULL`

Level: intermediate

Note:
Any existing `PetscSection` will be destroyed

See also: 
=== 
`DM`, `DMGetGlobalSection()`, `DMSetLocalSection()`

# External Links
$(_doc_external("DM/DMSetGlobalSection"))
"""
function DMSetGlobalSection(dm::AbstractDM{PetscLib},section::PetscSection) where {PetscLib}

	LibPETSc.DMSetGlobalSection(
		PetscLib,
		dm,
		section,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	perm,blockStarts = DMCreateSectionPermutation(dm::AbstractDM{PetscLib})

Create a permutation of the `PetscSection` chart and optionally a block structure.

Input Parameter:
===
- `dm` - The `DM`

Output Parameters:
===
- `perm`        - A permutation of the mesh points in the chart
- `blockStarts` - A high bit is set for the point that begins every block, or `NULL` for default blocking

Level: developer

See also: 
=== 
`DM`, `PetscSection`, `DMGetLocalSection()`, `DMGetGlobalSection()`

# External Links
$(_doc_external("DM/DMCreateSectionPermutation"))
"""
function DMCreateSectionPermutation(dm::AbstractDM{PetscLib}) where {PetscLib}
	perm = LibPETSc.IS()
	blockStarts = LibPETSc.PetscBT()

	LibPETSc.DMCreateSectionPermutation(
		PetscLib,
		dm,
		perm,
		blockStarts,
	)

	return perm,blockStarts
end
 
 
"""
	 UNTESTED !!!
	 DMReorderSectionSetDefault(dm::AbstractDM{PetscLib},reorder::DMReorderDefaultFlag)

Set flag indicating whether the local section should be reordered by default

Logically collective

Input Parameters:
===
- `dm`      - The DM
- `reorder` - Flag for reordering

Level: intermediate

-seealso: `DMReorderSectionGetDefault()`

# External Links
$(_doc_external("DM/DMReorderSectionSetDefault"))
"""
function DMReorderSectionSetDefault(dm::AbstractDM{PetscLib},reorder::DMReorderDefaultFlag) where {PetscLib}

	LibPETSc.DMReorderSectionSetDefault(
		PetscLib,
		dm,
		reorder,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMReorderSectionSetType(dm::AbstractDM{PetscLib},reorder::MatOrderingType)

Set the type of local section reordering

Logically collective

Input Parameters:
===
- `dm`      - The DM
- `reorder` - The reordering method

Level: intermediate

-seealso: `DMReorderSectionGetType()`, `DMReorderSectionSetDefault()`

# External Links
$(_doc_external("DM/DMReorderSectionSetType"))
"""
function DMReorderSectionSetType(dm::AbstractDM{PetscLib},reorder::MatOrderingType) where {PetscLib}

	LibPETSc.DMReorderSectionSetType(
		PetscLib,
		dm,
		reorder,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMUseTensorOrder(dm::AbstractDM{PetscLib},tensor::PetscBool)

Use a tensor product closure ordering for the default section

Input Parameters:
===
- `dm`     - The DM
- `tensor` - Flag for tensor order

Level: developer

-seealso: `DMPlexSetClosurePermutationTensor()`, `PetscSectionResetClosurePermutation()`

# External Links
$(_doc_external("DM/DMUseTensorOrder"))
"""
function DMUseTensorOrder(dm::AbstractDM{PetscLib},tensor::PetscBool) where {PetscLib}

	LibPETSc.DMUseTensorOrder(
		PetscLib,
		dm,
		tensor,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	sf = DMGetSectionSF(dm::AbstractDM{PetscLib})

Get the `PetscSF` encoding the parallel dof overlap for the `DM`. If it has not been set,
it is created from the default `PetscSection` layouts in the `DM`.

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `sf` - The `PetscSF`

Level: intermediate

Note:
This gets a borrowed reference, so the user should not destroy this `PetscSF`.

See also: 
=== 
`DM`, `DMSetSectionSF()`, `DMCreateSectionSF()`

# External Links
$(_doc_external("DM/DMGetSectionSF"))
"""
function DMGetSectionSF(dm::AbstractDM{PetscLib}) where {PetscLib}
	sf = LibPETSc.PetscSF()

	LibPETSc.DMGetSectionSF(
		PetscLib,
		dm,
		sf,
	)

	return sf
end
 
 
"""
	 UNTESTED !!!
	 DMSetSectionSF(dm::AbstractDM{PetscLib},sf::PetscSF)

Set the `PetscSF` encoding the parallel dof overlap for the `DM`

Input Parameters:
===
- `dm` - The `DM`
- `sf` - The `PetscSF`

Level: intermediate

Note:
Any previous `PetscSF` is destroyed

See also: 
=== 
`DM`, `DMGetSectionSF()`, `DMCreateSectionSF()`

# External Links
$(_doc_external("DM/DMSetSectionSF"))
"""
function DMSetSectionSF(dm::AbstractDM{PetscLib},sf::PetscSF) where {PetscLib}

	LibPETSc.DMSetSectionSF(
		PetscLib,
		dm,
		sf,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCreateSectionSF(dm::AbstractDM{PetscLib},localSection::PetscSection,globalSection::PetscSection)

Create the `PetscSF` encoding the parallel dof overlap for the `DM` based upon the `PetscSection`s
describing the data layout.

Input Parameters:
===
- `dm`            - The `DM`
- `localSection`  - `PetscSection` describing the local data layout
- `globalSection` - `PetscSection` describing the global data layout

Level: developer

Note:
One usually uses `DMGetSectionSF()` to obtain the `PetscSF`

Developer Note:
Since this routine has for arguments the two sections from the `DM` and puts the resulting `PetscSF`
directly into the `DM`, perhaps this function should not take the local and global sections as
input and should just obtain them from the `DM`? Plus PETSc creation functions return the thing
they create, this returns nothing

See also: 
=== 
`DM`, `DMGetSectionSF()`, `DMSetSectionSF()`, `DMGetLocalSection()`, `DMGetGlobalSection()`

# External Links
$(_doc_external("DM/DMCreateSectionSF"))
"""
function DMCreateSectionSF(dm::AbstractDM{PetscLib},localSection::PetscSection,globalSection::PetscSection) where {PetscLib}

	LibPETSc.DMCreateSectionSF(
		PetscLib,
		dm,
		localSection,
		globalSection,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	sf = DMGetPointSF(dm::AbstractDM{PetscLib})

Get the `PetscSF` encoding the parallel section point overlap for the `DM`.

Not collective but the resulting `PetscSF` is collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `sf` - The `PetscSF`

Level: intermediate

Note:
This gets a borrowed reference, so the user should not destroy this `PetscSF`.

See also: 
=== 
`DM`, `DMSetPointSF()`, `DMGetSectionSF()`, `DMSetSectionSF()`, `DMCreateSectionSF()`

# External Links
$(_doc_external("DM/DMGetPointSF"))
"""
function DMGetPointSF(dm::AbstractDM{PetscLib}) where {PetscLib}
	sf = LibPETSc.PetscSF()

	LibPETSc.DMGetPointSF(
		PetscLib,
		dm,
		sf,
	)

	return sf
end
 
 
"""
	 UNTESTED !!!
	 DMSetPointSF(dm::AbstractDM{PetscLib},sf::PetscSF)

Set the `PetscSF` encoding the parallel section point overlap for the `DM`.

Collective

Input Parameters:
===
- `dm` - The `DM`
- `sf` - The `PetscSF`

Level: intermediate

See also: 
=== 
`DM`, `DMGetPointSF()`, `DMGetSectionSF()`, `DMSetSectionSF()`, `DMCreateSectionSF()`

# External Links
$(_doc_external("DM/DMSetPointSF"))
"""
function DMSetPointSF(dm::AbstractDM{PetscLib},sf::PetscSF) where {PetscLib}

	LibPETSc.DMSetPointSF(
		PetscLib,
		dm,
		sf,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	sf = DMGetNaturalSF(dm::AbstractDM{PetscLib})

Get the `PetscSF` encoding the map back to the original mesh ordering

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `sf` - The `PetscSF`

Level: intermediate

Note:
This gets a borrowed reference, so the user should not destroy this `PetscSF`.

See also: 
=== 
`DM`, `DMSetNaturalSF()`, `DMSetUseNatural()`, `DMGetUseNatural()`, `DMPlexCreateGlobalToNaturalSF()`, `DMPlexDistribute()`

# External Links
$(_doc_external("DM/DMGetNaturalSF"))
"""
function DMGetNaturalSF(dm::AbstractDM{PetscLib}) where {PetscLib}
	sf = LibPETSc.PetscSF()

	LibPETSc.DMGetNaturalSF(
		PetscLib,
		dm,
		sf,
	)

	return sf
end
 
 
"""
	 UNTESTED !!!
	 DMSetNaturalSF(dm::AbstractDM{PetscLib},sf::PetscSF)

Set the PetscSF encoding the map back to the original mesh ordering

Input Parameters:
===
- `dm` - The DM
- `sf` - The PetscSF

Level: intermediate

See also: 
=== 
`DM`, `DMGetNaturalSF()`, `DMSetUseNatural()`, `DMGetUseNatural()`, `DMPlexCreateGlobalToNaturalSF()`, `DMPlexDistribute()`

# External Links
$(_doc_external("DM/DMSetNaturalSF"))
"""
function DMSetNaturalSF(dm::AbstractDM{PetscLib},sf::PetscSF) where {PetscLib}

	LibPETSc.DMSetNaturalSF(
		PetscLib,
		dm,
		sf,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	section,mat,bias = DMGetDefaultConstraints(dm::AbstractDM{PetscLib})

Get the `PetscSection` and `Mat` that specify the local constraint interpolation. See `DMSetDefaultConstraints()` for a description of the purpose of constraint interpolation.

not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameters:
===
- `section` - The `PetscSection` describing the range of the constraint matrix: relates rows of the constraint matrix to dofs of the default section.  Returns `NULL` if there are no local constraints.
- `mat`     - The `Mat` that interpolates local constraints: its width should be the layout size of the default section.  Returns `NULL` if there are no local constraints.
- `bias`    - Vector containing bias to be added to constrained dofs

Level: advanced

Note:
This gets borrowed references, so the user should not destroy the `PetscSection`, `Mat`, or `Vec`.

See also: 
=== 
`DM`, `DMSetDefaultConstraints()`

# External Links
$(_doc_external("DM/DMGetDefaultConstraints"))
"""
function DMGetDefaultConstraints(dm::AbstractDM{PetscLib}) where {PetscLib}
	section = LibPETSc.PetscSection()
	bias = CVec()

	LibPETSc.DMGetDefaultConstraints(
		PetscLib,
		dm,
		section,
		mat,
		bias,
	)

	return section,mat,bias
end
 
 
"""
	 UNTESTED !!!
	 DMSetDefaultConstraints(dm::AbstractDM{PetscLib},section::PetscSection,mat::AbstractMatrix,bias::AbstractVector)

Set the `PetscSection` and `Mat` that specify the local constraint interpolation.

Collective

Input Parameters:
===
- `dm`      - The `DM`
- `section` - The `PetscSection` describing the range of the constraint matrix: relates rows of the constraint matrix to dofs of the default section.  Must have a local communicator (`PETSC_COMM_SELF` or derivative).
- `mat`     - The `Mat` that interpolates local constraints: its width should be the layout size of the default section:  `NULL` indicates no constraints.  Must have a local communicator (`PETSC_COMM_SELF` or derivative).
- `bias`    - A bias vector to be added to constrained values in the local vector.  `NULL` indicates no bias.  Must have a local communicator (`PETSC_COMM_SELF` or derivative).

Level: advanced

Notes:
If a constraint matrix is specified, then it is applied during `DMGlobalToLocalEnd()` when mode is `INSERT_VALUES`, `INSERT_BC_VALUES`, or `INSERT_ALL_VALUES`.  Without a constraint matrix, the local vector l returned by `DMGlobalToLocalEnd()` contains values that have been scattered from a global vector without modification; with a constraint matrix A, l is modified by computing c = A * l + bias, l[s[i]] = c[i], where the scatter s is defined by the `PetscSection` returned by `DMGetDefaultConstraints()`.

If a constraint matrix is specified, then its adjoint is applied during `DMLocalToGlobalBegin()` when mode is `ADD_VALUES`, `ADD_BC_VALUES`, or `ADD_ALL_VALUES`.  Without a constraint matrix, the local vector l is accumulated into a global vector without modification; with a constraint matrix A, l is first modified by computing c[i] = l[s[i]], l[s[i]] = 0, l = l + A'*c, which is the adjoint of the operation described above.  Any bias, if specified, is ignored when accumulating.

This increments the references of the `PetscSection`, `Mat`, and `Vec`, so they user can destroy them.

See also: 
=== 
`DM`, `DMGetDefaultConstraints()`

# External Links
$(_doc_external("DM/DMSetDefaultConstraints"))
"""
function DMSetDefaultConstraints(dm::AbstractDM{PetscLib},section::PetscSection,mat::AbstractMatrix,bias::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetDefaultConstraints(
		PetscLib,
		dm,
		section,
		mat,
		bias,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	odm = DMGetOutputDM(dm::AbstractDM{PetscLib})

Retrieve the `DM` associated with the layout for output

Collective

Input Parameter:
===
- `dm` - The original `DM`

Output Parameter:
===
- `odm` - The `DM` which provides the layout for output

Level: intermediate

Note:
In some situations the vector obtained with `DMCreateGlobalVector()` excludes points for degrees of freedom that are associated with fixed (Dirichelet) boundary
conditions since the algebraic solver does not solve for those variables. The output `DM` includes these excluded points and its global vector contains the
locations for those dof so that they can be output to a file or other viewer along with the unconstrained dof.

See also: 
=== 
`DM`, `VecView()`, `DMGetGlobalSection()`, `DMCreateGlobalVector()`, `PetscSectionHasConstraints()`, `DMSetGlobalSection()`

# External Links
$(_doc_external("DM/DMGetOutputDM"))
"""
function DMGetOutputDM(dm::AbstractDM{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	odm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.DMGetOutputDM(
		PetscLib,
		dm,
		odm,
	)

	return odm
end
 
 
"""
	 UNTESTED !!!
	 DMSetOutputSequenceNumber(dm::AbstractDM{PetscLib},num::Int,val<:AbstractFloat)

Set the sequence number/value for output

Input Parameters:
===
- `dm`  - The original `DM`
- `num` - The output sequence number
- `val` - The output sequence value

Level: intermediate

Note:
This is intended for output that should appear in sequence, for instance
a set of timesteps in an `PETSCVIEWERHDF5` file, or a set of realizations of a stochastic system.

See also: 
=== 
`DM`, `VecView()`

# External Links
$(_doc_external("DM/DMSetOutputSequenceNumber"))
"""
function DMSetOutputSequenceNumber(dm::AbstractDM{PetscLib},num::Int,val::Float64) where {PetscLib}

	LibPETSc.DMSetOutputSequenceNumber(
		PetscLib,
		dm,
		num,
		val,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	numFields = DMGetNumFields(dm::AbstractDM{PetscLib})

Get the number of fields in the `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `numFields` - The number of fields

Level: intermediate

See also: 
=== 
`DM`, `DMSetNumFields()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMGetNumFields"))
"""
function DMGetNumFields(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	numFields = [PetscInt(1)]

	LibPETSc.DMGetNumFields(
		PetscLib,
		dm,
		Ref(numFields,1),
	)

	return numFields[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetNumFields(dm::AbstractDM{PetscLib},numFields::Int)

Set the number of fields in the `DM`

Logically Collective

Input Parameters:
===
- `dm`        - The `DM`
- `numFields` - The number of fields

Level: intermediate

See also: 
=== 
`DM`, `DMGetNumFields()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMSetNumFields"))
"""
function DMSetNumFields(dm::AbstractDM{PetscLib},numFields::Int) where {PetscLib}

	LibPETSc.DMSetNumFields(
		PetscLib,
		dm,
		numFields,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	label,disc = DMGetField(dm::AbstractDM{PetscLib},f::Int)

Return the `DMLabel` and discretization object for a given `DM` field

Not Collective

Input Parameters:
===
- `dm` - The `DM`
- `f`  - The field number

Output Parameters:
===
- `label` - The label indicating the support of the field, or `NULL` for the entire mesh (pass in `NULL` if not needed)
- `disc`  - The discretization object (pass in `NULL` if not needed)

Level: intermediate

See also: 
=== 
`DM`, `DMAddField()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMGetField"))
"""
function DMGetField(dm::AbstractDM{PetscLib},f::Int) where {PetscLib}
	label = LibPETSc.DMLabel()
	disc = LibPETSc.PetscObject()

	LibPETSc.DMGetField(
		PetscLib,
		dm,
		f,
		label,
		disc,
	)

	return label,disc
end
 
 
"""
	 UNTESTED !!!
	 DMSetField(dm::AbstractDM{PetscLib},f::Int,label::DMLabel,disc::PetscObject)

Set the discretization object for a given `DM` field. Usually one would call `DMAddField()` which automatically handles
the field numbering.

Logically Collective

Input Parameters:
===
- `dm`    - The `DM`
- `f`     - The field number
- `label` - The label indicating the support of the field, or `NULL` for the entire mesh
- `disc`  - The discretization object

Level: intermediate

See also: 
=== 
`DM`, `DMAddField()`, `DMGetField()`

# External Links
$(_doc_external("DM/DMSetField"))
"""
function DMSetField(dm::AbstractDM{PetscLib},f::Int,label::DMLabel,disc::PetscObject) where {PetscLib}

	LibPETSc.DMSetField(
		PetscLib,
		dm,
		f,
		label,
		disc,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMAddField(dm::AbstractDM{PetscLib},label::DMLabel,disc::PetscObject)

Add a field to a `DM` object. A field is a function space defined by of a set of discretization points (geometric entities)
and a discretization object that defines the function space associated with those points.

Logically Collective

Input Parameters:
===
- `dm`    - The `DM`
- `label` - The label indicating the support of the field, or `NULL` for the entire mesh
- `disc`  - The discretization object

Level: intermediate

Notes:
The label already exists or will be added to the `DM` with `DMSetLabel()`.

For example, a piecewise continuous pressure field can be defined by coefficients at the cell centers of a mesh and piecewise constant functions
within each cell. Thus a specific function in the space is defined by the combination of a `Vec` containing the coefficients, a `DM` defining the
geometry entities, a `DMLabel` indicating a subset of those geometric entities, and a discretization object, such as a `PetscFE`.

See also: 
=== 
`DM`, `DMSetLabel()`, `DMSetField()`, `DMGetField()`, `PetscFE`

# External Links
$(_doc_external("DM/DMAddField"))
"""
function DMAddField(dm::AbstractDM{PetscLib},label::DMLabel,disc::PetscObject) where {PetscLib}

	LibPETSc.DMAddField(
		PetscLib,
		dm,
		label,
		disc,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMSetFieldAvoidTensor(dm::AbstractDM{PetscLib},f::Int,avoidTensor::PetscBool)

Set flag to avoid defining the field on tensor cells

Logically Collective

Input Parameters:
===
- `dm`          - The `DM`
- `f`           - The field index
- `avoidTensor` - `PETSC_TRUE` to skip defining the field on tensor cells

Level: intermediate

See also: 
=== 
`DM`, `DMGetFieldAvoidTensor()`, `DMSetField()`, `DMGetField()`

# External Links
$(_doc_external("DM/DMSetFieldAvoidTensor"))
"""
function DMSetFieldAvoidTensor(dm::AbstractDM{PetscLib},f::Int,avoidTensor::PetscBool) where {PetscLib}

	LibPETSc.DMSetFieldAvoidTensor(
		PetscLib,
		dm,
		f,
		avoidTensor,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	avoidTensor = DMGetFieldAvoidTensor(dm::AbstractDM{PetscLib},f::Int)

Get flag to avoid defining the field on tensor cells

Not Collective

Input Parameters:
===
- `dm` - The `DM`
- `f`  - The field index

Output Parameter:
===
- `avoidTensor` - The flag to avoid defining the field on tensor cells

Level: intermediate

See also: 
=== 
`DM`, `DMAddField()`, `DMSetField()`, `DMGetField()`, `DMSetFieldAvoidTensor()`

# External Links
$(_doc_external("DM/DMGetFieldAvoidTensor"))
"""
function DMGetFieldAvoidTensor(dm::AbstractDM{PetscLib},f::Int) where {PetscLib}
	avoidTensor = Ref{PetscBool}()

	LibPETSc.DMGetFieldAvoidTensor(
		PetscLib,
		dm,
		f,
		avoidTensor,
	)

	return avoidTensor[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMClearFields(dm::AbstractDM{PetscLib})

Remove all fields from the `DM`

Logically Collective

Input Parameter:
===
- `dm` - The `DM`

Level: intermediate

See also: 
=== 
`DM`, `DMGetNumFields()`, `DMSetNumFields()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMClearFields"))
"""
function DMClearFields(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearFields(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCopyFields(dm::AbstractDM{PetscLib},minDegree::Int,maxDegree::Int,newdm::AbstractDM{PetscLib})

Copy the discretizations for the `DM` into another `DM`

Collective

Input Parameters:
===
- `dm`        - The `DM`
- `minDegree` - Minimum degree for a discretization, or `PETSC_DETERMINE` for no limit
- `maxDegree` - Maximum degree for a discretization, or `PETSC_DETERMINE` for no limit

Output Parameter:
===
- `newdm` - The `DM`

Level: advanced

See also: 
=== 
`DM`, `DMGetField()`, `DMSetField()`, `DMAddField()`, `DMCopyDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMCopyFields"))
"""
function DMCopyFields(dm::AbstractDM{PetscLib},minDegree::Int,maxDegree::Int,newdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCopyFields(
		PetscLib,
		dm,
		minDegree,
		maxDegree,
		newdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	useCone,useClosure = DMGetAdjacency(dm::AbstractDM{PetscLib},f::Int)

Returns the flags for determining variable influence

Not Collective

Input Parameters:
===
- `dm` - The `DM` object
- `f`  - The field number, or `PETSC_DEFAULT` for the default adjacency

Output Parameters:
===
- `useCone`    - Flag for variable influence starting with the cone operation
- `useClosure` - Flag for variable influence using transitive closure

Level: developer

Notes:
-vb
FEM:   Two points p and q are adjacent if q in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
FVM:   Two points p and q are adjacent if q in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
FVM++: Two points p and q are adjacent if q in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
-ve
Further explanation can be found in the User's Manual Section on the Influence of Variables on One Another.

See also: 
=== 
`DM`, `DMSetAdjacency()`, `DMGetField()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMGetAdjacency"))
"""
function DMGetAdjacency(dm::AbstractDM{PetscLib},f::Int) where {PetscLib}
	useCone = Ref{PetscBool}()
	useClosure = Ref{PetscBool}()

	LibPETSc.DMGetAdjacency(
		PetscLib,
		dm,
		f,
		useCone,
		useClosure,
	)

	return useCone[] == PETSC_TRUE,useClosure[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMSetAdjacency(dm::AbstractDM{PetscLib},f::Int,useCone::PetscBool,useClosure::PetscBool)

Set the flags for determining variable influence

Not Collective

Input Parameters:
===
- `dm`         - The `DM` object
- `f`          - The field number
- `useCone`    - Flag for variable influence starting with the cone operation
- `useClosure` - Flag for variable influence using transitive closure

Level: developer

Notes:
-vb
FEM:   Two points p and q are adjacent if q in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
FVM:   Two points p and q are adjacent if q in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
FVM++: Two points p and q are adjacent if q in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
-ve
Further explanation can be found in the User's Manual Section on the Influence of Variables on One Another.

See also: 
=== 
`DM`, `DMGetAdjacency()`, `DMGetField()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMSetAdjacency"))
"""
function DMSetAdjacency(dm::AbstractDM{PetscLib},f::Int,useCone::PetscBool,useClosure::PetscBool) where {PetscLib}

	LibPETSc.DMSetAdjacency(
		PetscLib,
		dm,
		f,
		useCone,
		useClosure,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	useCone,useClosure = DMGetBasicAdjacency(dm::AbstractDM{PetscLib})

Returns the flags for determining variable influence, using either the default or field 0 if it is defined

Not collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameters:
===
- `useCone`    - Flag for variable influence starting with the cone operation
- `useClosure` - Flag for variable influence using transitive closure

Level: developer

Notes:
-vb
FEM:   Two points p and q are adjacent if q in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
FVM:   Two points p and q are adjacent if q in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
FVM++: Two points p and q are adjacent if q in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
-ve

See also: 
=== 
`DM`, `DMSetBasicAdjacency()`, `DMGetField()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMGetBasicAdjacency"))
"""
function DMGetBasicAdjacency(dm::AbstractDM{PetscLib}) where {PetscLib}
	useCone = Ref{PetscBool}()
	useClosure = Ref{PetscBool}()

	LibPETSc.DMGetBasicAdjacency(
		PetscLib,
		dm,
		useCone,
		useClosure,
	)

	return useCone[] == PETSC_TRUE,useClosure[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMSetBasicAdjacency(dm::AbstractDM{PetscLib},useCone::PetscBool,useClosure::PetscBool)

Set the flags for determining variable influence, using either the default or field 0 if it is defined

Not Collective

Input Parameters:
===
- `dm`         - The `DM` object
- `useCone`    - Flag for variable influence starting with the cone operation
- `useClosure` - Flag for variable influence using transitive closure

Level: developer

Notes:
-vb
FEM:   Two points p and q are adjacent if q in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
FVM:   Two points p and q are adjacent if q in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
FVM++: Two points p and q are adjacent if q in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
-ve

See also: 
=== 
`DM`, `DMGetBasicAdjacency()`, `DMGetField()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMSetBasicAdjacency"))
"""
function DMSetBasicAdjacency(dm::AbstractDM{PetscLib},useCone::PetscBool,useClosure::PetscBool) where {PetscLib}

	LibPETSc.DMSetBasicAdjacency(
		PetscLib,
		dm,
		useCone,
		useClosure,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	Nds = DMGetNumDS(dm::AbstractDM{PetscLib})

Get the number of discrete systems in the `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `Nds` - The number of `PetscDS` objects

Level: intermediate

See also: 
=== 
`DM`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMGetNumDS"))
"""
function DMGetNumDS(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	Nds = [PetscInt(1)]

	LibPETSc.DMGetNumDS(
		PetscLib,
		dm,
		Ref(Nds,1),
	)

	return Nds[1]
end
 
 
"""
	 UNTESTED !!!
	ds = DMGetDS(dm::AbstractDM{PetscLib})

Get the default `PetscDS`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `ds` - The default `PetscDS`

Level: intermediate

See also: 
=== 
`DM`, `DMGetCellDS()`, `DMGetRegionDS()`

# External Links
$(_doc_external("DM/DMGetDS"))
"""
function DMGetDS(dm::AbstractDM{PetscLib}) where {PetscLib}
	ds = LibPETSc.PetscDS()

	LibPETSc.DMGetDS(
		PetscLib,
		dm,
		ds,
	)

	return ds
end
 
 
"""
	 UNTESTED !!!
	ds,dsIn = DMGetCellDS(dm::AbstractDM{PetscLib},point::Int)

Get the `PetscDS` defined on a given cell

Not Collective

Input Parameters:
===
- `dm`    - The `DM`
- `point` - Cell for the `PetscDS`

Output Parameters:
===
- `ds`   - The `PetscDS` defined on the given cell
- `dsIn` - The `PetscDS` for input on the given cell, or NULL if the same ds

Level: developer

See also: 
=== 
`DM`, `DMGetDS()`, `DMSetRegionDS()`

# External Links
$(_doc_external("DM/DMGetCellDS"))
"""
function DMGetCellDS(dm::AbstractDM{PetscLib},point::Int) where {PetscLib}
	ds = LibPETSc.PetscDS()
	dsIn = LibPETSc.PetscDS()

	LibPETSc.DMGetCellDS(
		PetscLib,
		dm,
		point,
		ds,
		dsIn,
	)

	return ds,dsIn
end
 
 
"""
	 UNTESTED !!!
	fields,ds,dsIn = DMGetRegionDS(dm::AbstractDM{PetscLib},label::DMLabel)

Get the `PetscDS` for a given mesh region, defined by a `DMLabel`

Not Collective

Input Parameters:
===
- `dm`    - The `DM`
- `label` - The `DMLabel` defining the mesh region, or `NULL` for the entire mesh

Output Parameters:
===
- `fields` - The `IS` containing the `DM` field numbers for the fields in this `PetscDS`, or `NULL`
- `ds`     - The `PetscDS` defined on the given region, or `NULL`
- `dsIn`   - The `PetscDS` for input in the given region, or `NULL`

Level: advanced

Note:
If a non-`NULL` label is given, but there is no `PetscDS` on that specific label,
the `PetscDS` for the full domain (if present) is returned. Returns with
fields = `NULL` and ds = `NULL` if there is no `PetscDS` for the full domain.

See also: 
=== 
`DM`, `DMGetRegionNumDS()`, `DMSetRegionDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMGetRegionDS"))
"""
function DMGetRegionDS(dm::AbstractDM{PetscLib},label::DMLabel) where {PetscLib}
	fields = LibPETSc.IS()
	ds = LibPETSc.PetscDS()
	dsIn = LibPETSc.PetscDS()

	LibPETSc.DMGetRegionDS(
		PetscLib,
		dm,
		label,
		fields,
		ds,
		dsIn,
	)

	return fields,ds,dsIn
end
 
 
"""
	 UNTESTED !!!
	 DMSetRegionDS(dm::AbstractDM{PetscLib},label::DMLabel,fields::IS,ds::PetscDS,dsIn::PetscDS)

Set the `PetscDS` for a given mesh region, defined by a `DMLabel`

Collective

Input Parameters:
===
- `dm`     - The `DM`
- `label`  - The `DMLabel` defining the mesh region, or `NULL` for the entire mesh
- `fields` - The `IS` containing the `DM` field numbers for the fields in this `PetscDS`, or `NULL` for all fields
- `ds`     - The `PetscDS` defined on the given region
- `dsIn`   - The `PetscDS` for input on the given cell, or `NULL` if it is the same `PetscDS`

Level: advanced

Note:
If the label has a `PetscDS` defined, it will be replaced. Otherwise, it will be added to the `DM`. If the `PetscDS` is replaced,
the fields argument is ignored.

See also: 
=== 
`DM`, `DMGetRegionDS()`, `DMSetRegionNumDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMSetRegionDS"))
"""
function DMSetRegionDS(dm::AbstractDM{PetscLib},label::DMLabel,fields::IS,ds::PetscDS,dsIn::PetscDS) where {PetscLib}

	LibPETSc.DMSetRegionDS(
		PetscLib,
		dm,
		label,
		fields,
		ds,
		dsIn,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	label,fields,ds,dsIn = DMGetRegionNumDS(dm::AbstractDM{PetscLib},num::Int)

Get the `PetscDS` for a given mesh region, defined by the region number

Not Collective

Input Parameters:
===
- `dm`  - The `DM`
- `num` - The region number, in [0, Nds)

Output Parameters:
===
- `label`  - The region label, or `NULL`
- `fields` - The `IS` containing the `DM` field numbers for the fields in this `PetscDS`, or `NULL`
- `ds`     - The `PetscDS` defined on the given region, or `NULL`
- `dsIn`   - The `PetscDS` for input in the given region, or `NULL`

Level: advanced

See also: 
=== 
`DM`, `DMGetRegionDS()`, `DMSetRegionDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMGetRegionNumDS"))
"""
function DMGetRegionNumDS(dm::AbstractDM{PetscLib},num::Int) where {PetscLib}
	label = LibPETSc.DMLabel()
	fields = LibPETSc.IS()
	ds = LibPETSc.PetscDS()
	dsIn = LibPETSc.PetscDS()

	LibPETSc.DMGetRegionNumDS(
		PetscLib,
		dm,
		num,
		label,
		fields,
		ds,
		dsIn,
	)

	return label,fields,ds,dsIn
end
 
 
"""
	 UNTESTED !!!
	 DMSetRegionNumDS(dm::AbstractDM{PetscLib},num::Int,label::DMLabel,fields::IS,ds::PetscDS,dsIn::PetscDS)

Set the `PetscDS` for a given mesh region, defined by the region number

Not Collective

Input Parameters:
===
- `dm`     - The `DM`
- `num`    - The region number, in [0, Nds)
- `label`  - The region label, or `NULL`
- `fields` - The `IS` containing the `DM` field numbers for the fields in this `PetscDS`, or `NULL` to prevent setting
- `ds`     - The `PetscDS` defined on the given region, or `NULL` to prevent setting
- `dsIn`   - The `PetscDS` for input on the given cell, or `NULL` if it is the same `PetscDS`

Level: advanced

See also: 
=== 
`DM`, `DMGetRegionDS()`, `DMSetRegionDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMSetRegionNumDS"))
"""
function DMSetRegionNumDS(dm::AbstractDM{PetscLib},num::Int,label::DMLabel,fields::IS,ds::PetscDS,dsIn::PetscDS) where {PetscLib}

	LibPETSc.DMSetRegionNumDS(
		PetscLib,
		dm,
		num,
		label,
		fields,
		ds,
		dsIn,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	num = DMFindRegionNum(dm::AbstractDM{PetscLib},ds::PetscDS)

Find the region number for a given `PetscDS`, or

Not Collective

Input Parameters:
===
- `dm` - The `DM`
- `ds` - The `PetscDS` defined on the given region

Output Parameter:
===
- `num` - The region number, in [0, Nds), or -1 if not found

Level: advanced

See also: 
=== 
`DM`, `DMGetRegionNumDS()`, `DMGetRegionDS()`, `DMSetRegionDS()`, `DMGetDS()`, `DMGetCellDS()`

# External Links
$(_doc_external("DM/DMFindRegionNum"))
"""
function DMFindRegionNum(dm::AbstractDM{PetscLib},ds::PetscDS) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	num = [PetscInt(1)]

	LibPETSc.DMFindRegionNum(
		PetscLib,
		dm,
		ds,
		Ref(num,1),
	)

	return num[1]
end
 
 
"""
	 UNTESTED !!!
	 DMCreateDS(dm::AbstractDM{PetscLib})

Create the discrete systems for the `DM` based upon the fields added to the `DM`

Collective

Input Parameter:
===
- `dm` - The `DM`

Options Database Key:
===
- `-dm_petscds_view` - View all the `PetscDS` objects in this `DM`

Level: intermediate

See also: 
=== 
`DM`, `DMSetField`, `DMAddField()`, `DMGetDS()`, `DMGetCellDS()`, `DMGetRegionDS()`, `DMSetRegionDS()`

# External Links
$(_doc_external("DM/DMCreateDS"))
"""
function DMCreateDS(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCreateDS(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearDS(dm::AbstractDM{PetscLib})

Remove all discrete systems from the `DM`

Logically Collective

Input Parameter:
===
- `dm` - The `DM`

Level: intermediate

See also: 
=== 
`DM`, `DMGetNumDS()`, `DMGetDS()`, `DMSetField()`

# External Links
$(_doc_external("DM/DMClearDS"))
"""
function DMClearDS(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearDS(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCopyDS(dm::AbstractDM{PetscLib},minDegree::Int,maxDegree::Int,newdm::AbstractDM{PetscLib})

Copy the discrete systems for the `DM` into another `DM`

Collective

Input Parameters:
===
- `dm`        - The `DM`
- `minDegree` - Minimum degree for a discretization, or `PETSC_DETERMINE` for no limit
- `maxDegree` - Maximum degree for a discretization, or `PETSC_DETERMINE` for no limit

Output Parameter:
===
- `newdm` - The `DM`

Level: advanced

See also: 
=== 
`DM`, `DMCopyFields()`, `DMAddField()`, `DMGetDS()`, `DMGetCellDS()`, `DMGetRegionDS()`, `DMSetRegionDS()`

# External Links
$(_doc_external("DM/DMCopyDS"))
"""
function DMCopyDS(dm::AbstractDM{PetscLib},minDegree::Int,maxDegree::Int,newdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCopyDS(
		PetscLib,
		dm,
		minDegree,
		maxDegree,
		newdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCopyDisc(dm::AbstractDM{PetscLib},newdm::AbstractDM{PetscLib})

Copy the fields and discrete systems for the `DM` into another `DM`

Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `newdm` - The `DM`

Level: advanced

Developer Note:
Really ugly name, nothing in PETSc is called a `Disc` plus it is an ugly abbreviation

See also: 
=== 
`DM`, `DMCopyFields()`, `DMCopyDS()`

# External Links
$(_doc_external("DM/DMCopyDisc"))
"""
function DMCopyDisc(dm::AbstractDM{PetscLib},newdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCopyDisc(
		PetscLib,
		dm,
		newdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMComputeExactSolution(dm::AbstractDM{PetscLib},time<:AbstractFloat,u::AbstractVector,u_t::AbstractVector)

Compute the exact solution for a given `DM`, using the `PetscDS` information.

Collective

Input Parameters:
===
- `dm`   - The `DM`
- `time` - The time

Output Parameters:
===
- `u`   - The vector will be filled with exact solution values, or `NULL`
- `u_t` - The vector will be filled with the time derivative of exact solution values, or `NULL`

Level: developer

Note:
The user must call `PetscDSSetExactSolution()` before using this routine

See also: 
=== 
`DM`, `PetscDSSetExactSolution()`

# External Links
$(_doc_external("DM/DMComputeExactSolution"))
"""
function DMComputeExactSolution(dm::AbstractDM{PetscLib},time::Float64,u::AbstractVector,u_t::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMComputeExactSolution(
		PetscLib,
		dm,
		time,
		u,
		u_t,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	numAux = DMGetNumAuxiliaryVec(dm::AbstractDM{PetscLib})

Get the number of auxiliary vectors associated with this `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `numAux` - The number of auxiliary data vectors

Level: advanced

See also: 
=== 
`DM`, `DMClearAuxiliaryVec()`, `DMSetAuxiliaryVec()`, `DMGetAuxiliaryLabels()`, `DMGetAuxiliaryVec()`

# External Links
$(_doc_external("DM/DMGetNumAuxiliaryVec"))
"""
function DMGetNumAuxiliaryVec(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	numAux = [PetscInt(1)]

	LibPETSc.DMGetNumAuxiliaryVec(
		PetscLib,
		dm,
		Ref(numAux,1),
	)

	return numAux[1]
end
 
 
"""
	 UNTESTED !!!
	aux = DMGetAuxiliaryVec(dm::AbstractDM{PetscLib},label::DMLabel,value::Int,part::Int)

Get the auxiliary vector for region specified by the given label and value, and equation part

Not Collective

Input Parameters:
===
- `dm`    - The `DM`
- `label` - The `DMLabel`
- `value` - The label value indicating the region
- `part`  - The equation part, or 0 if unused

Output Parameter:
===
- `aux` - The `Vec` holding auxiliary field data

Level: advanced

Note:
If no auxiliary vector is found for this (label, value), (NULL, 0, 0) is checked as well.

See also: 
=== 
`DM`, `DMClearAuxiliaryVec()`, `DMSetAuxiliaryVec()`, `DMGetNumAuxiliaryVec()`, `DMGetAuxiliaryLabels()`

# External Links
$(_doc_external("DM/DMGetAuxiliaryVec"))
"""
function DMGetAuxiliaryVec(dm::AbstractDM{PetscLib},label::DMLabel,value::Int,part::Int) where {PetscLib}
	aux = CVec()

	LibPETSc.DMGetAuxiliaryVec(
		PetscLib,
		dm,
		label,
		value,
		part,
		aux,
	)

	return aux
end
 
 
"""
	 UNTESTED !!!
	 DMSetAuxiliaryVec(dm::AbstractDM{PetscLib},label::DMLabel,value::Int,part::Int,aux::AbstractVector)

Set an auxiliary vector for region specified by the given label and value, and equation part

Not Collective because auxiliary vectors are not parallel

Input Parameters:
===
- `dm`    - The `DM`
- `label` - The `DMLabel`
- `value` - The label value indicating the region
- `part`  - The equation part, or 0 if unused
- `aux`   - The `Vec` holding auxiliary field data

Level: advanced

See also: 
=== 
`DM`, `DMClearAuxiliaryVec()`, `DMGetAuxiliaryVec()`, `DMGetAuxiliaryLabels()`, `DMCopyAuxiliaryVec()`

# External Links
$(_doc_external("DM/DMSetAuxiliaryVec"))
"""
function DMSetAuxiliaryVec(dm::AbstractDM{PetscLib},label::DMLabel,value::Int,part::Int,aux::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetAuxiliaryVec(
		PetscLib,
		dm,
		label,
		value,
		part,
		aux,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMGetAuxiliaryLabels(dm::AbstractDM{PetscLib},labels::Vector{DMLabel},values::Vector{Int},parts::Vector{Int})

Get the labels, values, and parts for all auxiliary vectors in this `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameters:
===
- `labels` - The `DMLabel`s for each `Vec`
- `values` - The label values for each `Vec`
- `parts`  - The equation parts for each `Vec`

Level: advanced

Note:
The arrays passed in must be at least as large as `DMGetNumAuxiliaryVec()`.

See also: 
=== 
`DM`, `DMClearAuxiliaryVec()`, `DMGetNumAuxiliaryVec()`, `DMGetAuxiliaryVec()`, `DMSetAuxiliaryVec()`, `DMCopyAuxiliaryVec()`

# External Links
$(_doc_external("DM/DMGetAuxiliaryLabels"))
"""
function DMGetAuxiliaryLabels(dm::AbstractDM{PetscLib},labels::Vector{DMLabel},values::Vector{Int},parts::Vector{Int}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetAuxiliaryLabels(
		PetscLib,
		dm,
		labels,
		values,
		parts,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCopyAuxiliaryVec(dm::AbstractDM{PetscLib},dmNew::AbstractDM{PetscLib})

Copy the auxiliary vector data on a `DM` to a new `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `dmNew` - The new `DM`, now with the same auxiliary data

Level: advanced

Note:
This is a shallow copy of the auxiliary vectors

See also: 
=== 
`DM`, `DMClearAuxiliaryVec()`, `DMGetNumAuxiliaryVec()`, `DMGetAuxiliaryVec()`, `DMSetAuxiliaryVec()`

# External Links
$(_doc_external("DM/DMCopyAuxiliaryVec"))
"""
function DMCopyAuxiliaryVec(dm::AbstractDM{PetscLib},dmNew::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCopyAuxiliaryVec(
		PetscLib,
		dm,
		dmNew,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearAuxiliaryVec(dm::AbstractDM{PetscLib})

Destroys the auxiliary vector information and creates a new empty one

Not Collective

Input Parameter:
===
- `dm` - The `DM`

Level: advanced

See also: 
=== 
`DM`, `DMCopyAuxiliaryVec()`, `DMGetNumAuxiliaryVec()`, `DMGetAuxiliaryVec()`, `DMSetAuxiliaryVec()`

# External Links
$(_doc_external("DM/DMClearAuxiliaryVec"))
"""
function DMClearAuxiliaryVec(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMClearAuxiliaryVec(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCreateLabel(dm::AbstractDM{PetscLib},name::Vector{Char})

Create a label of the given name if it does not already exist in the `DM`

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Level: intermediate

See also: 
=== 
`DM`, `DMLabelCreate()`, `DMHasLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMCreateLabel"))
"""
function DMCreateLabel(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMCreateLabel(
		PetscLib,
		dm,
		name,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMCreateLabelAtIndex(dm::AbstractDM{PetscLib},l::Int,name::Vector{Char})

Create a label of the given name at the given index. If it already exists in the `DM`, move it to this index.

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `l`    - The index for the label
- `name` - The label name

Level: intermediate

See also: 
=== 
`DM`, `DMCreateLabel()`, `DMLabelCreate()`, `DMHasLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMCreateLabelAtIndex"))
"""
function DMCreateLabelAtIndex(dm::AbstractDM{PetscLib},l::Int,name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMCreateLabelAtIndex(
		PetscLib,
		dm,
		l,
		name,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	value = DMGetLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int)

Get the value in a `DMLabel` for the given point, with

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name
- `point` - The mesh point

Output Parameter:
===
- `value` - The label value for this point, or -1 if the point is not in the label

Level: beginner

See also: 
=== 
`DM`, `DMLabelGetValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMGetLabelValue"))
"""
function DMGetLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	value = [PetscInt(1)]

	LibPETSc.DMGetLabelValue(
		PetscLib,
		dm,
		name,
		point,
		Ref(value,1),
	)

	return value[1]
end
 
 
"""
	 UNTESTED !!!
	 DMSetLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int,value::Int)

Add a point to a `DMLabel` with given value

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name
- `point` - The mesh point
- `value` - The label value for this point

Output Parameter:
===

Level: beginner

See also: 
=== 
`DM`, `DMLabelSetValue()`, `DMGetStratumIS()`, `DMClearLabelValue()`

# External Links
$(_doc_external("DM/DMSetLabelValue"))
"""
function DMSetLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int,value::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetLabelValue(
		PetscLib,
		dm,
		name,
		point,
		value,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int,value::Int)

Remove a point from a `DMLabel` with given value

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name
- `point` - The mesh point
- `value` - The label value for this point

Level: beginner

See also: 
=== 
`DM`, `DMLabelClearValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMClearLabelValue"))
"""
function DMClearLabelValue(dm::AbstractDM{PetscLib},name::Vector{Char},point::Int,value::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMClearLabelValue(
		PetscLib,
		dm,
		name,
		point,
		value,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	size = DMGetLabelSize(dm::AbstractDM{PetscLib},name::Vector{Char})

Get the value of `DMLabelGetNumValues()` of a `DMLabel` in the `DM`

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `size` - The number of different integer ids, or 0 if the label does not exist

Level: beginner

Developer Note:
This should be renamed to something like `DMGetLabelNumValues()` or removed.

See also: 
=== 
`DM`, `DMLabelGetNumValues()`, `DMSetLabelValue()`, `DMGetLabel()`

# External Links
$(_doc_external("DM/DMGetLabelSize"))
"""
function DMGetLabelSize(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	size = [PetscInt(1)]

	LibPETSc.DMGetLabelSize(
		PetscLib,
		dm,
		name,
		Ref(size,1),
	)

	return size[1]
end
 
 
"""
	 UNTESTED !!!
	ids = DMGetLabelIdIS(dm::AbstractDM{PetscLib},name::Vector{Char})

Get the `DMLabelGetValueIS()` from a `DMLabel` in the `DM`

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `ids` - The integer ids, or `NULL` if the label does not exist

Level: beginner

See also: 
=== 
`DM`, `DMLabelGetValueIS()`, `DMGetLabelSize()`

# External Links
$(_doc_external("DM/DMGetLabelIdIS"))
"""
function DMGetLabelIdIS(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	ids = LibPETSc.IS()

	LibPETSc.DMGetLabelIdIS(
		PetscLib,
		dm,
		name,
		ids,
	)

	return ids
end
 
 
"""
	 UNTESTED !!!
	size = DMGetStratumSize(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int)

Get the number of points in a label stratum

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name of the stratum
- `value` - The stratum value

Output Parameter:
===
- `size` - The number of points, also called the stratum size

Level: beginner

See also: 
=== 
`DM`, `DMLabelGetStratumSize()`, `DMGetLabelSize()`, `DMGetLabelIds()`

# External Links
$(_doc_external("DM/DMGetStratumSize"))
"""
function DMGetStratumSize(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	size = [PetscInt(1)]

	LibPETSc.DMGetStratumSize(
		PetscLib,
		dm,
		name,
		value,
		Ref(size,1),
	)

	return size[1]
end
 
 
"""
	 UNTESTED !!!
	points = DMGetStratumIS(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int)

Get the points in a label stratum

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name
- `value` - The stratum value

Output Parameter:
===
- `points` - The stratum points, or `NULL` if the label does not exist or does not have that value

Level: beginner

See also: 
=== 
`DM`, `DMLabelGetStratumIS()`, `DMGetStratumSize()`

# External Links
$(_doc_external("DM/DMGetStratumIS"))
"""
function DMGetStratumIS(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	points = LibPETSc.IS()

	LibPETSc.DMGetStratumIS(
		PetscLib,
		dm,
		name,
		value,
		points,
	)

	return points
end
 
 
"""
	 UNTESTED !!!
	 DMSetStratumIS(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int,points::IS)

Set the points in a label stratum

Not Collective

Input Parameters:
===
- `dm`     - The `DM` object
- `name`   - The label name
- `value`  - The stratum value
- `points` - The stratum points

Level: beginner

See also: 
=== 
`DM`, `DMLabel`, `DMClearLabelStratum()`, `DMLabelClearStratum()`, `DMLabelSetStratumIS()`, `DMGetStratumSize()`

# External Links
$(_doc_external("DM/DMSetStratumIS"))
"""
function DMSetStratumIS(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int,points::IS) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetStratumIS(
		PetscLib,
		dm,
		name,
		value,
		points,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 DMClearLabelStratum(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int)

Remove all points from a stratum from a `DMLabel`

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `name`  - The label name
- `value` - The label value for this point

Output Parameter:
===

Level: beginner

See also: 
=== 
`DM`, `DMLabel`, `DMLabelClearStratum()`, `DMSetLabelValue()`, `DMGetStratumIS()`, `DMClearLabelValue()`

# External Links
$(_doc_external("DM/DMClearLabelStratum"))
"""
function DMClearLabelStratum(dm::AbstractDM{PetscLib},name::Vector{Char},value::Int) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMClearLabelStratum(
		PetscLib,
		dm,
		name,
		value,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	output = DMGetLabelOutput(dm::AbstractDM{PetscLib},name::Vector{Char})

Get the output flag for a given label

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `output` - The flag for output

Level: developer

See also: 
=== 
`DM`, `DMLabel`, `DMSetLabelOutput()`, `DMCreateLabel()`, `DMHasLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMGetLabelOutput"))
"""
function DMGetLabelOutput(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	output = Ref{PetscBool}()

	LibPETSc.DMGetLabelOutput(
		PetscLib,
		dm,
		name,
		output,
	)

	return output[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMSetLabelOutput(dm::AbstractDM{PetscLib},name::Vector{Char},output::PetscBool)

Set if a given label should be saved to a `PetscViewer` in calls to `DMView()`

Not Collective

Input Parameters:
===
- `dm`     - The `DM` object
- `name`   - The label name
- `output` - `PETSC_TRUE` to save the label to the viewer

Level: developer

See also: 
=== 
`DM`, `DMLabel`, `DMGetOutputFlag()`, `DMGetLabelOutput()`, `DMCreateLabel()`, `DMHasLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMSetLabelOutput"))
"""
function DMSetLabelOutput(dm::AbstractDM{PetscLib},name::Vector{Char},output::PetscBool) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMSetLabelOutput(
		PetscLib,
		dm,
		name,
		output,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	numLabels = DMGetNumLabels(dm::AbstractDM{PetscLib})

Return the number of labels defined by on the `DM`

Not Collective

Input Parameter:
===
- `dm` - The `DM` object

Output Parameter:
===
- `numLabels` - the number of Labels

Level: intermediate

See also: 
=== 
`DM`, `DMLabel`, `DMGetLabelByNum()`, `DMGetLabelName()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMGetNumLabels"))
"""
function DMGetNumLabels(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	numLabels = [PetscInt(1)]

	LibPETSc.DMGetNumLabels(
		PetscLib,
		dm,
		Ref(numLabels,1),
	)

	return numLabels[1]
end
 
 
"""
	 UNTESTED !!!
	 DMGetLabelName(dm::AbstractDM{PetscLib},n::Int,name::Vector{Char})

Return the name of nth label

Not Collective

Input Parameters:
===
- `dm` - The `DM` object
- `n`  - the label number

Output Parameter:
===
- `name` - the label name

Level: intermediate

Developer Note:
Some of the functions that appropriate on labels using their number have the suffix ByNum, others do not.

See also: 
=== 
`DM`, `DMLabel`, `DMGetLabelByNum()`, `DMGetLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMGetLabelName"))
"""
function DMGetLabelName(dm::AbstractDM{PetscLib},n::Int,name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.DMGetLabelName(
		PetscLib,
		dm,
		n,
		name,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	hasLabel = DMHasLabel(dm::AbstractDM{PetscLib},name::Vector{Char})

Determine whether the `DM` has a label of a given name

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `hasLabel` - `PETSC_TRUE` if the label is present

Level: intermediate

See also: 
=== 
`DM`, `DMLabel`, `DMGetLabel()`, `DMGetLabelByNum()`, `DMCreateLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMHasLabel"))
"""
function DMHasLabel(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	hasLabel = Ref{PetscBool}()

	LibPETSc.DMHasLabel(
		PetscLib,
		dm,
		name,
		hasLabel,
	)

	return hasLabel[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	label = DMGetLabel(dm::AbstractDM{PetscLib},name::Vector{Char})

Return the label of a given name, or `NULL`, from a `DM`

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `label` - The `DMLabel`, or `NULL` if the label is absent

Default labels in a `DMPLEX`:
- `"depth"`       - Holds the depth (co-dimension) of each mesh point
- `"celltype"`    - Holds the topological type of each cell
- `"ghost"`       - If the DM is distributed with overlap, this marks the cells and faces in the overlap
- `"Cell Sets"`   - Mirrors the cell sets defined by GMsh and ExodusII
- `"Face Sets"`   - Mirrors the face sets defined by GMsh and ExodusII
- `"Vertex Sets"` - Mirrors the vertex sets defined by GMsh

Level: intermediate

See also: 
=== 
`DM`, `DMLabel`, `DMHasLabel()`, `DMGetLabelByNum()`, `DMAddLabel()`, `DMCreateLabel()`, `DMPlexGetDepthLabel()`, `DMPlexGetCellType()`

# External Links
$(_doc_external("DM/DMGetLabel"))
"""
function DMGetLabel(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	label = LibPETSc.DMLabel()

	LibPETSc.DMGetLabel(
		PetscLib,
		dm,
		name,
		label,
	)

	return label
end
 
 
"""
	 UNTESTED !!!
	 DMSetLabel(dm::AbstractDM{PetscLib},label::DMLabel)

Replaces the label of a given name, or ignores it if the name is not present

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `label` - The `DMLabel`, having the same name, to substitute

Default labels in a `DMPLEX`:
- `"depth"`       - Holds the depth (co-dimension) of each mesh point
- `"celltype"`    - Holds the topological type of each cell
- `"ghost"`       - If the DM is distributed with overlap, this marks the cells and faces in the overlap
- `"Cell Sets"`   - Mirrors the cell sets defined by GMsh and ExodusII
- `"Face Sets"`   - Mirrors the face sets defined by GMsh and ExodusII
- `"Vertex Sets"` - Mirrors the vertex sets defined by GMsh

Level: intermediate

See also: 
=== 
`DM`, `DMLabel`, `DMCreateLabel()`, `DMHasLabel()`, `DMPlexGetDepthLabel()`, `DMPlexGetCellType()`

# External Links
$(_doc_external("DM/DMSetLabel"))
"""
function DMSetLabel(dm::AbstractDM{PetscLib},label::DMLabel) where {PetscLib}

	LibPETSc.DMSetLabel(
		PetscLib,
		dm,
		label,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	label = DMGetLabelByNum(dm::AbstractDM{PetscLib},n::Int)

Return the nth label on a `DM`

Not Collective

Input Parameters:
===
- `dm` - The `DM` object
- `n`  - the label number

Output Parameter:
===
- `label` - the label

Level: intermediate

See also: 
=== 
`DM`, `DMLabel`, `DMAddLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMGetLabelByNum"))
"""
function DMGetLabelByNum(dm::AbstractDM{PetscLib},n::Int) where {PetscLib}
	label = LibPETSc.DMLabel()

	LibPETSc.DMGetLabelByNum(
		PetscLib,
		dm,
		n,
		label,
	)

	return label
end
 
 
"""
	 UNTESTED !!!
	 DMAddLabel(dm::AbstractDM{PetscLib},label::DMLabel)

Add the label to this `DM`

Not Collective

Input Parameters:
===
- `dm`    - The `DM` object
- `label` - The `DMLabel`

Level: developer

See also: 
=== 
`DM`, `DMLabel`, `DMCreateLabel()`, `DMHasLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMGetStratumIS()`

# External Links
$(_doc_external("DM/DMAddLabel"))
"""
function DMAddLabel(dm::AbstractDM{PetscLib},label::DMLabel) where {PetscLib}

	LibPETSc.DMAddLabel(
		PetscLib,
		dm,
		label,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	label = DMRemoveLabel(dm::AbstractDM{PetscLib},name::Vector{Char})

Remove the label given by name from this `DM`

Not Collective

Input Parameters:
===
- `dm`   - The `DM` object
- `name` - The label name

Output Parameter:
===
- `label` - The `DMLabel`, or `NULL` if the label is absent. Pass in `NULL` to call `DMLabelDestroy()` on the label, otherwise the
caller is responsible for calling `DMLabelDestroy()`.

Level: developer

See also: 
=== 
`DM`, `DMLabel`, `DMCreateLabel()`, `DMHasLabel()`, `DMGetLabel()`, `DMGetLabelValue()`, `DMSetLabelValue()`, `DMLabelDestroy()`, `DMRemoveLabelBySelf()`

# External Links
$(_doc_external("DM/DMRemoveLabel"))
"""
function DMRemoveLabel(dm::AbstractDM{PetscLib},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	label = LibPETSc.DMLabel()

	LibPETSc.DMRemoveLabel(
		PetscLib,
		dm,
		name,
		label,
	)

	return label
end
 
 
"""
	 UNTESTED !!!
	label = DMRemoveLabelBySelf(dm::AbstractDM{PetscLib},failNotFound::PetscBool)

Remove the label from this `DM`

Not Collective

Input Parameters:
===
- `dm`           - The `DM` object
- `label`        - The `DMLabel` to be removed from the `DM`
- `failNotFound` - Should it fail if the label is not found in the `DM`?

Level: developer

Note:
Only exactly the same instance is removed if found, name match is ignored.
If the `DM` has an exclusive reference to the label, the label gets destroyed and
*label nullified.

See also: 
=== 
`DM`, `DMLabel`, `DMCreateLabel()`, `DMHasLabel()`, `DMGetLabel()` `DMGetLabelValue()`, `DMSetLabelValue()`, `DMLabelDestroy()`, `DMRemoveLabel()`

# External Links
$(_doc_external("DM/DMRemoveLabelBySelf"))
"""
function DMRemoveLabelBySelf(dm::AbstractDM{PetscLib},failNotFound::PetscBool) where {PetscLib}
	label = LibPETSc.DMLabel()

	LibPETSc.DMRemoveLabelBySelf(
		PetscLib,
		dm,
		failNotFound,
		label,
	)

	return label
end
 
 
"""
	 UNTESTED !!!
	 DMCopyLabels(dmA::AbstractDM{PetscLib},dmB::AbstractDM{PetscLib},mode::PetscCopyMode,all::PetscBool,emode::DMCopyLabelsMode)

Copy labels from one `DM` mesh to another `DM` with a superset of the points

Collective

Input Parameters:
===
- `dmA`   - The `DM` object with initial labels
- `dmB`   - The `DM` object to which labels are copied
- `mode`  - Copy labels by pointers (`PETSC_OWN_POINTER`) or duplicate them (`PETSC_COPY_VALUES`)
- `all`   - Copy all labels including "depth", "dim", and "celltype" (`PETSC_TRUE`) which are otherwise ignored (`PETSC_FALSE`)
- `emode` - How to behave when a `DMLabel` in the source and destination `DM`s with the same name is encountered (see `DMCopyLabelsMode`)

Level: intermediate

Note:
This is typically used when interpolating or otherwise adding to a mesh, or testing.

See also: 
=== 
`DM`, `DMLabel`, `DMAddLabel()`, `DMCopyLabelsMode`

# External Links
$(_doc_external("DM/DMCopyLabels"))
"""
function DMCopyLabels(dmA::AbstractDM{PetscLib},dmB::AbstractDM{PetscLib},mode::PetscCopyMode,all::PetscBool,emode::DMCopyLabelsMode) where {PetscLib}

	LibPETSc.DMCopyLabels(
		PetscLib,
		dmA,
		dmB,
		mode,
		all,
		emode,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	isBd = DMIsBoundaryPoint(dm::AbstractDM{PetscLib},point::Int)


# External Links
$(_doc_external("DM/DMIsBoundaryPoint"))
"""
function DMIsBoundaryPoint(dm::AbstractDM{PetscLib},point::Int) where {PetscLib}
	isBd = Ref{PetscBool}()

	LibPETSc.DMIsBoundaryPoint(
		PetscLib,
		dm,
		point,
		isBd,
	)

	return isBd[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	errorVec = DMComputeError(dm::AbstractDM{PetscLib},sol::AbstractVector,errors::Vector{Float64})

Computes the error assuming the user has provided the exact solution functions

Collective

Input Parameters:
===
- `dm`  - The `DM`
- `sol` - The solution vector

Input/Output Parameter:
- `errors` - An array of length Nf, the number of fields, or `NULL` for no output; on output
contains the error in each field

Output Parameter:
===
- `errorVec` - A vector to hold the cellwise error (may be `NULL`)

Level: developer

Note:
The exact solutions come from the `PetscDS` object, and the time comes from `DMGetOutputSequenceNumber()`.

See also: 
=== 
`DM`, `DMMonitorSet()`, `DMGetRegionNumDS()`, `PetscDSGetExactSolution()`, `DMGetOutputSequenceNumber()`

# External Links
$(_doc_external("DM/DMComputeError"))
"""
function DMComputeError(dm::AbstractDM{PetscLib},sol::AbstractVector,errors::Vector{Float64}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	errorVec = CVec()

	LibPETSc.DMComputeError(
		PetscLib,
		dm,
		sol,
		errors,
		errorVec,
	)

	return errorVec
end
 
 
"""
	 UNTESTED !!!
	flg = DMHasBasisTransform(dm::AbstractDM{PetscLib})

Whether the `DM` employs a basis transformation from functions in global vectors to functions in local vectors

Input Parameter:
===
- `dm` - The `DM`

Output Parameter:
===
- `flg` - `PETSC_TRUE` if a basis transformation should be done

Level: developer

See also: 
=== 
`DM`, `DMPlexGlobalToLocalBasis()`, `DMPlexLocalToGlobalBasis()`, `DMPlexCreateBasisRotation()`

# External Links
$(_doc_external("DM/DMHasBasisTransform"))
"""
function DMHasBasisTransform(dm::AbstractDM{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.DMHasBasisTransform(
		PetscLib,
		dm,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMCopyTransform(dm::AbstractDM{PetscLib},newdm::AbstractDM{PetscLib})


# External Links
$(_doc_external("DM/DMCopyTransform"))
"""
function DMCopyTransform(dm::AbstractDM{PetscLib},newdm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMCopyTransform(
		PetscLib,
		dm,
		newdm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	compatible,set = DMGetCompatibility(dm1::AbstractDM{PetscLib},dm2::AbstractDM{PetscLib})

determine if two `DM`s are compatible

Collective

Input Parameters:
===
- `dm1` - the first `DM`
- `dm2` - the second `DM`

Output Parameters:
===
- `compatible` - whether or not the two `DM`s are compatible
- `set`        - whether or not the compatible value was actually determined and set

Level: advanced

Notes:
Two `DM`s are deemed compatible if they represent the same parallel decomposition
of the same topology. This implies that the section (field data) on one
"makes sense" with respect to the topology and parallel decomposition of the other.
Loosely speaking, compatible `DM`s represent the same domain and parallel
decomposition, but hold different data.

Typically, one would confirm compatibility if intending to simultaneously iterate
over a pair of vectors obtained from different `DM`s.

For example, two `DMDA` objects are compatible if they have the same local
and global sizes and the same stencil width. They can have different numbers
of degrees of freedom per node. Thus, one could use the node numbering from
either `DM` in bounds for a loop over vectors derived from either `DM`.

Consider the operation of summing data living on a 2-dof `DMDA` to data living
on a 1-dof `DMDA`, which should be compatible, as in the following snippet.
-vb
-..
PetscCall(DMGetCompatibility(da1,da2,&compatible,&set));
if (set && compatible)  {
PetscCall(DMDAVecGetArrayDOF(da1,vec1,&arr1));
PetscCall(DMDAVecGetArrayDOF(da2,vec2,&arr2));
PetscCall(DMDAGetCorners(da1,&x,&y,NULL,&m,&n,NULL));
for (j=y; j<y+n; ++j) {
for (i=x; i<x+m, ++i) {
arr1[j][i][0] = arr2[j][i][0] + arr2[j][i][1];
}
}
PetscCall(DMDAVecRestoreArrayDOF(da1,vec1,&arr1));
PetscCall(DMDAVecRestoreArrayDOF(da2,vec2,&arr2));
} else {
SETERRQ(PetscObjectComm((PetscObject)da1,PETSC_ERR_ARG_INCOMP,"DMDA objects incompatible");
}
-..
-ve

Checking compatibility might be expensive for a given implementation of `DM`,
or might be impossible to unambiguously confirm or deny. For this reason,
this function may decline to determine compatibility, and hence users should
always check the "set" output parameter.

A `DM` is always compatible with itself.

In the current implementation, `DM`s which live on "unequal" communicators
(MPI_UNEQUAL in the terminology of MPI_Comm_compare()) are always deemed
incompatible.

This function is labeled "Collective," as information about all subdomains
is required on each rank. However, in `DM` implementations which store all this
information locally, this function may be merely "Logically Collective".

Developer Note:
Compatibility is assumed to be a symmetric concept; `DM` A is compatible with `DM` B
iff B is compatible with A. Thus, this function checks the implementations
of both dm and dmc (if they are of different types), attempting to determine
compatibility. It is left to `DM` implementers to ensure that symmetry is
preserved. The simplest way to do this is, when implementing type-specific
logic for this function, is to check for existing logic in the implementation
of other `DM` types and let *set = PETSC_FALSE if found.

See also: 
=== 
`DM`, `DMDACreateCompatibleDMDA()`, `DMStagCreateCompatibleDMStag()`

# External Links
$(_doc_external("DM/DMGetCompatibility"))
"""
function DMGetCompatibility(dm1::AbstractDM{PetscLib},dm2::AbstractDM{PetscLib}) where {PetscLib}
	compatible = Ref{PetscBool}()
	set = Ref{PetscBool}()

	LibPETSc.DMGetCompatibility(
		PetscLib,
		dm1,
		dm2,
		compatible,
		set,
	)

	return compatible[] == PETSC_TRUE,set[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 DMMonitorCancel(dm::AbstractDM{PetscLib})

Clears all the monitor functions for a `DM` object.

Logically Collective

Input Parameter:
===
- `dm` - the DM

Options Database Key:
===
- `-dm_monitor_cancel` - cancels all monitors that have been hardwired
into a code by calls to `DMonitorSet()`, but does not cancel those
set via the options database

Level: intermediate

Note:
There is no way to clear one specific monitor from a `DM` object.

See also: 
=== 
`DM`, `DMMonitorSet()`, `DMMonitorSetFromOptions()`, `DMMonitor()`

# External Links
$(_doc_external("DM/DMMonitorCancel"))
"""
function DMMonitorCancel(dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.DMMonitorCancel(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ornt,found = DMPolytopeMatchOrientation(ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})

Determine an orientation (transformation) that takes the source face arrangement to the target face arrangement

Not Collective

Input Parameters:
===
- `ct`         - The `DMPolytopeType`
- `sourceCone` - The source arrangement of faces
- `targetCone` - The target arrangement of faces

Output Parameters:
===
- `ornt`  - The orientation (transformation) which will take the source arrangement to the target arrangement
- `found` - Flag indicating that a suitable orientation was found

Level: advanced

Note:
An arrangement is a face order combined with an orientation for each face

Each orientation (transformation) is labeled with an integer from negative `DMPolytopeTypeGetNumArrangements(ct)`/2 to `DMPolytopeTypeGetNumArrangements(ct)`/2
that labels each arrangement (face ordering plus orientation for each face).

See `DMPolytopeMatchVertexOrientation()` to find a new vertex orientation that takes the source vertex arrangement to the target vertex arrangement

See also: 
=== 
`DM`, `DMPolytopeGetOrientation()`, `DMPolytopeMatchVertexOrientation()`, `DMPolytopeGetVertexOrientation()`

# External Links
$(_doc_external("DM/DMPolytopeMatchOrientation"))
"""
function DMPolytopeMatchOrientation(PetscLib, ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	ornt = [PetscInt(1)]
	found = Ref{PetscBool}()

	LibPETSc.DMPolytopeMatchOrientation(
		PetscLib,
		ct,
		sourceCone,
		targetCone,
		Ref(ornt,1),
		found,
	)

	return ornt[1],found[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	ornt,found = DMPolytopeMatchVertexOrientation(PetscLib,ct::DMPolytopeType,sourceVert::Vector{Int},targetVert::Vector{Int})

Determine an orientation (transformation) that takes the source vertex arrangement to the target vertex arrangement

Not Collective

Input Parameters:
===
- `ct`         - The `DMPolytopeType`
- `sourceVert` - The source arrangement of vertices
- `targetVert` - The target arrangement of vertices

Output Parameters:
===
- `ornt`  - The orientation (transformation) which will take the source arrangement to the target arrangement
- `found` - Flag indicating that a suitable orientation was found

Level: advanced

Notes:
An arrangement is a vertex order

Each orientation (transformation) is labeled with an integer from negative `DMPolytopeTypeGetNumArrangements(ct)`/2 to `DMPolytopeTypeGetNumArrangements(ct)`/2
that labels each arrangement (vertex ordering).

See `DMPolytopeMatchOrientation()` to find a new face orientation that takes the source face arrangement to the target face arrangement

See also: 
=== 
`DM`, `DMPolytopeType`, `DMPolytopeGetOrientation()`, `DMPolytopeMatchOrientation()`, `DMPolytopeTypeGetNumVertices()`, `DMPolytopeTypeGetVertexArrangement()`

# External Links
$(_doc_external("DM/DMPolytopeMatchVertexOrientation"))
"""
function DMPolytopeMatchVertexOrientation(PetscLib, ct::DMPolytopeType,sourceVert::Vector{Int},targetVert::Vector{Int})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	ornt = [PetscInt(1)]
	found = Ref{PetscBool}()

	LibPETSc.DMPolytopeMatchVertexOrientation(
		PetscLib,
		ct,
		sourceVert,
		targetVert,
		Ref(ornt,1),
		found,
	)

	return ornt[1],found[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	ornt = DMPolytopeGetOrientation(PetscLib, ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})

Determine an orientation (transformation) that takes the source face arrangement to the target face arrangement

Not Collective

Input Parameters:
===
- `ct`         - The `DMPolytopeType`
- `sourceCone` - The source arrangement of faces
- `targetCone` - The target arrangement of faces

Output Parameter:
===
- `ornt` - The orientation (transformation) which will take the source arrangement to the target arrangement

Level: advanced

Note:
This function is the same as `DMPolytopeMatchOrientation()` except it will generate an error if no suitable orientation can be found.

Developer Note:
It is unclear why this function needs to exist since one can simply call `DMPolytopeMatchOrientation()` and error if none is found

See also: 
=== 
`DM`, `DMPolytopeType`, `DMPolytopeMatchOrientation()`, `DMPolytopeGetVertexOrientation()`, `DMPolytopeMatchVertexOrientation()`

# External Links
$(_doc_external("DM/DMPolytopeGetOrientation"))
"""
function DMPolytopeGetOrientation(PetscLib, ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	ornt = [PetscInt(1)]

	LibPETSc.DMPolytopeGetOrientation(
		PetscLib,
		ct,
		sourceCone,
		targetCone,
		Ref(ornt,1),
	)

	return ornt[1]
end
 
 
"""
	 UNTESTED !!!
	ornt = DMPolytopeGetVertexOrientation(PetscLib, ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})

Determine an orientation (transformation) that takes the source vertex arrangement to the target vertex arrangement

Not Collective

Input Parameters:
===
- `ct`         - The `DMPolytopeType`
- `sourceCone` - The source arrangement of vertices
- `targetCone` - The target arrangement of vertices

Output Parameter:
===
- `ornt` - The orientation (transformation) which will take the source arrangement to the target arrangement

Level: advanced

Note:
This function is the same as `DMPolytopeMatchVertexOrientation()` except it errors if not orientation is possible.

Developer Note:
It is unclear why this function needs to exist since one can simply call `DMPolytopeMatchVertexOrientation()` and error if none is found

See also: 
=== 
`DM`, `DMPolytopeType`, `DMPolytopeMatchVertexOrientation()`, `DMPolytopeGetOrientation()`

# External Links
$(_doc_external("DM/DMPolytopeGetVertexOrientation"))
"""
function DMPolytopeGetVertexOrientation(PetscLib, ct::DMPolytopeType,sourceCone::Vector{Int},targetCone::Vector{Int})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	ornt = [PetscInt(1)]

	LibPETSc.DMPolytopeGetVertexOrientation(
		PetscLib,
		ct,
		sourceCone,
		targetCone,
		Ref(ornt,1),
	)

	return ornt[1]
end
 
 
"""
	 UNTESTED !!!
	inside = DMPolytopeInCellTest(PetscLib, ct::DMPolytopeType,point::Vector{Float64})

Check whether a point lies inside the reference cell of given type

Not Collective

Input Parameters:
===
- `ct`    - The `DMPolytopeType`
- `point` - Coordinates of the point

Output Parameter:
===
- `inside` - Flag indicating whether the point is inside the reference cell of given type

Level: advanced

See also: 
=== 
`DM`, `DMPolytopeType`, `DMLocatePoints()`

# External Links
$(_doc_external("DM/DMPolytopeInCellTest"))
"""
function DMPolytopeInCellTest(PetscLib, ct::DMPolytopeType,point::Vector{Float64})
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	inside = Ref{PetscBool}()

	LibPETSc.DMPolytopeInCellTest(
		PetscLib,
		ct,
		point,
		inside,
	)

	return inside[] == PETSC_TRUE
end
 
 
