"""
	dm = DMStagCreate1d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,M::Int,dof0::Int,dof1::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...)

Create an object to manage data living on the elements and vertices of a parallelized regular 1D grid.

Collective

Input Parameters:
===
- petsclib     - the PETSc library
- `comm`         - MPI communicator
- `bndx`         - boundary type: `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `M`            - global number of elements
- `dof0`         - number of degrees of freedom per vertex/0-cell
- `dof1`         - number of degrees of freedom per element/1-cell
- `stencilType`  - ghost/halo region type: `DMSTAG_STENCIL_BOX` or `DMSTAG_STENCIL_NONE`
- `stencilWidth` - width, in elements, of halo/ghost region
- `lx`           - array of local sizes, of length equal to the comm size, summing to `M` or `NULL`
- `dmsetfromoptions` - call set from options
- `dmsetup`          - call setup
- `options`          - additional options

Output Parameter:
===
- `dm` - the new `DMSTAG` object

Options Database Keys:
===
- `-dm_view`                                      - calls `DMViewFromOptions()` at the conclusion of `DMSetUp()`
- `-stag_grid_x <nx>`                             - number of elements in the x direction
- `-stag_ghost_stencil_width`                     - width of ghost region, in elements
- `-stag_boundary_type_x <none,ghosted,periodic>` - `DMBoundaryType` value

Level: beginner

Notes:
You must call `DMSetUp()` after this call before using the `DM`.
If you wish to use the options database (see the keys above) to change values in the `DMSTAG`, you must call
`DMSetFromOptions()` after this function but before `DMSetUp()`.

See also: 
=== 
`DMSTAG`, `DMStagCreate2d()`, `DMStagCreate3d()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateLocalVector()`, `DMLocalToGlobalBegin()`, `DMDACreate1d()`

# External Links
$(_doc_external("DMStag/DMStagCreate1d"))
"""
function DMStagCreate1d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,M::Int,dof0::Int,dof1::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...) where {PetscLib}
	opts = Options(petsclib; options...)
	dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

	with(dm.opts) do
		LibPETSc.DMStagCreate1d(
			petsclib,
			comm,
			bndx,
			M,
			dof0,
			dof1,
			stencilType,
			stencilWidth,
			lx,
			dm,
		)
	end
	dmsetfromoptions && setfromoptions!(dm)
	dmsetup && setup!(dm)

	return dm
end
 
 
"""
	dm = DMStagCreate2d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,bndy::DMBoundaryType,M::Int,N::Int,m::Int,n::Int,dof0::Int,dof1::Int,dof2::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},ly::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...)

Create an object to manage data living on the elements, faces, and vertices of a parallelized regular 2D grid.

Collective

Input Parameters:
===
- petsclib     - the PETSc library
- `comm`         - MPI communicator
- `bndx`         - x boundary type, `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `bndy`         - y boundary type, `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `M`            - global number of elements in x direction
- `N`            - global number of elements in y direction
- `m`            - number of ranks in the x direction (may be `PETSC_DECIDE`)
- `n`            - number of ranks in the y direction (may be `PETSC_DECIDE`)
- `dof0`         - number of degrees of freedom per vertex/0-cell
- `dof1`         - number of degrees of freedom per face/1-cell
- `dof2`         - number of degrees of freedom per element/2-cell
- `stencilType`  - ghost/halo region type: `DMSTAG_STENCIL_NONE`, `DMSTAG_STENCIL_BOX`, or `DMSTAG_STENCIL_STAR`
- `stencilWidth` - width, in elements, of halo/ghost region
- `lx`           - array of local x element counts, of length equal to `m`, summing to `M`, or `NULL`
- `ly`           - array of local y element counts, of length equal to `n`, summing to `N`, or `NULL`
- `dmsetfromoptions` - call set from options
- `dmsetup`          - call setup
- `options`          - additional options

Output Parameter:
===
- `dm` - the new `DMSTAG` object

Options Database Keys:
===
- `-dm_view`                                      - calls `DMViewFromOptions()` at the conclusion of `DMSetUp()`
- `-stag_grid_x <nx>`                             - number of elements in the x direction
- `-stag_grid_y <ny>`                             - number of elements in the y direction
- `-stag_ranks_x <rx>`                            - number of ranks in the x direction
- `-stag_ranks_y <ry>`                            - number of ranks in the y direction
- `-stag_ghost_stencil_width`                     - width of ghost region, in elements
- `-stag_boundary_type_x <none,ghosted,periodic>` - `DMBoundaryType` value
- `-stag_boundary_type_y <none,ghosted,periodic>` - `DMBoundaryType` value

Level: beginner

Notes:
You must call `DMSetUp()` after this call, before using the `DM`.
If you wish to use the options database (see the keys above) to change values in the `DMSTAG`, you must call
`DMSetFromOptions()` after this function but before `DMSetUp()`.

See also: 
=== 
`DMSTAG`, `DMStagCreate1d()`, `DMStagCreate3d()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateLocalVector()`, `DMLocalToGlobalBegin()`, `DMDACreate2d()`

# External Links
$(_doc_external("DMStag/DMStagCreate2d"))
"""
function DMStagCreate2d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,bndy::DMBoundaryType,M::Int,N::Int,m::Int,n::Int,dof0::Int,dof1::Int,dof2::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},ly::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...) where {PetscLib}
	opts = Options(petsclib; options...)
	dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

	with(dm.opts) do
		LibPETSc.DMStagCreate2d(
			petsclib,
			comm,
			bndx,
			bndy,
			M,
			N,
			m,
			n,
			dof0,
			dof1,
			dof2,
			stencilType,
			stencilWidth,
			lx,
			ly,
			dm,
		)
	end
	dmsetfromoptions && setfromoptions!(dm)
	dmsetup && setup!(dm)

	return dm
end
 
 
"""
	dm = DMStagCreate3d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,bndy::DMBoundaryType,bndz::DMBoundaryType,M::Int,N::Int,P::Int,m::Int,n::Int,p::Int,dof0::Int,dof1::Int,dof2::Int,dof3::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},ly::Vector{Int},lz::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...)

Create an object to manage data living on the elements, faces, edges, and vertices of a parallelized regular 3D grid.

Collective

Input Parameters:
===
- petsclib     - the PETSc library
- `comm`         - MPI communicator
- `bndx`         - x boundary type, `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `bndy`         - y boundary type, `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `bndz`         - z boundary type, `DM_BOUNDARY_NONE`, `DM_BOUNDARY_PERIODIC`, or `DM_BOUNDARY_GHOSTED`
- `M`            - global number of elements in x direction
- `N`            - global number of elements in y direction
- `P`            - global number of elements in z direction
- `m`            - number of ranks in the x direction (may be `PETSC_DECIDE`)
- `n`            - number of ranks in the y direction (may be `PETSC_DECIDE`)
- `p`            - number of ranks in the z direction (may be `PETSC_DECIDE`)
- `dof0`         - number of degrees of freedom per vertex/0-cell
- `dof1`         - number of degrees of freedom per edge/1-cell
- `dof2`         - number of degrees of freedom per face/2-cell
- `dof3`         - number of degrees of freedom per element/3-cell
- `stencilType`  - ghost/halo region type: `DMSTAG_STENCIL_NONE`, `DMSTAG_STENCIL_BOX`, or `DMSTAG_STENCIL_STAR`
- `stencilWidth` - width, in elements, of halo/ghost region
- `lx`           - array of local x  element counts, of length equal to `m`, summing to `M`, or `NULL`
- `ly`           - arrays of local y element counts, of length equal to `n`, summing to `N`, or `NULL`
- `lz`           - arrays of local z element counts, of length equal to `p`, summing to `P`, or `NULL`
- `dmsetfromoptions` - call set from options
- `dmsetup`          - call setup
- `options`          - additional options

Output Parameter:
===
- `dm` - the new `DMSTAG` object

Options Database Keys:
===
- `-dm_view`                                      - calls `DMViewFromOptions()` at the conclusion of `DMSetUp()`
- `-stag_grid_x <nx>`                             - number of elements in the x direction
- `-stag_grid_y <ny>`                             - number of elements in the y direction
- `-stag_grid_z <nz>`                             - number of elements in the z direction
- `-stag_ranks_x <rx>`                            - number of ranks in the x direction
- `-stag_ranks_y <ry>`                            - number of ranks in the y direction
- `-stag_ranks_z <rz>`                            - number of ranks in the z direction
- `-stag_ghost_stencil_width`                     - width of ghost region, in elements
- `-stag_boundary_type x <none,ghosted,periodic>` - `DMBoundaryType` value
- `-stag_boundary_type y <none,ghosted,periodic>` - `DMBoundaryType` value
- `-stag_boundary_type z <none,ghosted,periodic>` - `DMBoundaryType` value

Level: beginner

Notes:
You must call `DMSetUp()` after this call before using the `DM`.
If you wish to use the options database (see the keys above) to change values in the `DMSTAG`, you must call
`DMSetFromOptions()` after this function but before `DMSetUp()`.

See also: 
=== 
`DMSTAG`, `DMStagCreate1d()`, `DMStagCreate2d()`, `DMDestroy()`, `DMView()`, `DMCreateGlobalVector()`, `DMCreateLocalVector()`, `DMLocalToGlobalBegin()`, `DMDACreate3d()`

# External Links
$(_doc_external("DMStag/DMStagCreate3d"))
"""
function DMStagCreate3d(petsclib::PetscLib,comm::MPI_Comm,bndx::DMBoundaryType,bndy::DMBoundaryType,bndz::DMBoundaryType,M::Int,N::Int,P::Int,m::Int,n::Int,p::Int,dof0::Int,dof1::Int,dof2::Int,dof3::Int,stencilType::DMStagStencilType,stencilWidth::Int,lx::Vector{Int},ly::Vector{Int},lz::Vector{Int},dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...) where {PetscLib}
	opts = Options(petsclib; options...)
	dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

	with(dm.opts) do
		LibPETSc.DMStagCreate3d(
			petsclib,
			comm,
			bndx,
			bndy,
			bndz,
			M,
			N,
			P,
			m,
			n,
			p,
			dof0,
			dof1,
			dof2,
			dof3,
			stencilType,
			stencilWidth,
			lx,
			ly,
			lz,
			dm,
		)
	end
	dmsetfromoptions && setfromoptions!(dm)
	dmsetup && setup!(dm)

	return dm
end
 
 
"""
	newdm = DMStagCreateCompatibleDMStag(dm::AbstractDMStag{PetscLib},dof0::Int,dof1::Int,dof2::Int,dof3::Int,dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...)

create a compatible `DMSTAG` with different dof/stratum

Collective

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `dof0` - number of dof on the first stratum in the new `DMSTAG`
- `dof1` - number of dof on the second stratum in the new `DMSTAG`
- `dof2` - number of dof on the third stratum in the new `DMSTAG`
- `dof3` - number of dof on the fourth stratum in the new `DMSTAG`
- `dmsetfromoptions` - call set from options
- `dmsetup`          - call setup
- `options`          - additional options

Output Parameter:
===
- `newdm` - the new, compatible `DMSTAG`

Level: intermediate

Notes:
DOF supplied for strata too big for the dimension are ignored; these may be set to `0`.
For example, for a 2-dimensional `DMSTAG`, `dof2` sets the number of dof per element,
and `dof3` is unused. For a 3-dimensional `DMSTAG`, `dof3` sets the number of DOF per element.

In contrast to `DMDACreateCompatibleDMDA()`, coordinates are not reused.

See also: 
=== 
`DMSTAG`, `DMDACreateCompatibleDMDA()`, `DMGetCompatibility()`, `DMStagMigrateVec()`

# External Links
$(_doc_external("DMStag/DMStagCreateCompatibleDMStag"))
"""
function DMStagCreateCompatibleDMStag(dm::AbstractDMStag{PetscLib},dof0::Int,dof1::Int,dof2::Int,dof3::Int,dmsetfromoptions::Bool=true,dmsetup::Bool=true,options...) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib; options...)
	newdm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

	with(dm.opts) do
		LibPETSc.DMStagCreateCompatibleDMStag(
			PetscLib,
			dm,
			dof0,
			dof1,
			dof2,
			dof3,
			newdm,
		)
	end

	dmsetfromoptions && setfromoptions!(newdm)
	dmsetup && setup!(newdm)

	return newdm
end
 
 
"""
	boundaryTypeX,boundaryTypeY,boundaryTypeZ = DMStagGetBoundaryTypes(dm::AbstractDMStag{PetscLib})

get boundary types

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `boundaryTypeX` - boundary type for x direction
- `boundaryTypeY` - boundary type for y direction, not set for one dimensional problems
- `boundaryTypeZ` - boundary type for z direction, not set for one and two dimensional problems

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMBoundaryType`

# External Links
$(_doc_external("DMStag/DMStagGetBoundaryTypes"))
"""
function DMStagGetBoundaryTypes(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	boundaryTypeX = Ref{DMBoundaryType}(DM_BOUNDARY_NONE)
	boundaryTypeY = Ref{DMBoundaryType}(DM_BOUNDARY_NONE)
	boundaryTypeZ = Ref{DMBoundaryType}(DM_BOUNDARY_NONE)

	LibPETSc.DMStagGetBoundaryTypes(
		PetscLib,
		dm,
		boundaryTypeX,
		boundaryTypeY,
		boundaryTypeZ,
	)

	return boundaryTypeX[],boundaryTypeY[],boundaryTypeZ[]
end
 
 
"""
	x,y,z,m,n,p,nExtrax,nExtray,nExtraz = DMStagGetCorners(dm::AbstractDMStag{PetscLib})

return global element indices of the local region (excluding ghost points)

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `x`       - starting element index in first direction
- `y`       - starting element index in second direction
- `z`       - starting element index in third direction
- `m`       - element width in first direction
- `n`       - element width in second direction
- `p`       - element width in third direction
- `nExtrax` - number of extra partial elements in first direction
- `nExtray` - number of extra partial elements in second direction
- `nExtraz` - number of extra partial elements in third direction

Level: beginner

Notes:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

The number of extra partial elements is either 1 or 0.
The value is 1 on right, top, and front non-periodic domain ("physical") boundaries,
in the x, y, and z directions respectively, and otherwise 0.

See also: 
=== 
`DMSTAG`, `DMStagGetGhostCorners()`, `DMDAGetCorners()`

# External Links
$(_doc_external("DMStag/DMStagGetCorners"))
"""
function DMStagGetCorners(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	x = [PetscInt(1)]
	y = [PetscInt(1)]
	z = [PetscInt(1)]
	m = [PetscInt(1)]
	n = [PetscInt(1)]
	p = [PetscInt(1)]
	nExtrax = [PetscInt(1)]
	nExtray = [PetscInt(1)]
	nExtraz = [PetscInt(1)]

	LibPETSc.DMStagGetCorners(
		PetscLib,
		dm,
		Ref(x,1),
		Ref(y,1),
		Ref(z,1),
		Ref(m,1),
		Ref(n,1),
		Ref(p,1),
		Ref(nExtrax,1),
		Ref(nExtray,1),
		Ref(nExtraz,1),
	)

	return x[1],y[1],z[1],m[1],n[1],p[1],nExtrax[1],nExtray[1],nExtraz[1]
end
 
 
"""
	dof0,dof1,dof2,dof3 = DMStagGetDOF(dm::AbstractDMStag{PetscLib})

get number of DOF associated with each stratum of the grid

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `dof0` - the number of points per 0-cell (vertex/node)
- `dof1` - the number of points per 1-cell (element in 1D, edge in 2D and 3D)
- `dof2` - the number of points per 2-cell (element in 2D, face in 3D)
- `dof3` - the number of points per 3-cell (element in 3D)

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagGetCorners()`, `DMStagGetGhostCorners()`, `DMStagGetGlobalSizes()`, `DMStagGetStencilWidth()`, `DMStagGetBoundaryTypes()`, `DMStagGetLocationDOF()`, `DMDAGetDof()`

# External Links
$(_doc_external("DMStag/DMStagGetDOF"))
"""
function DMStagGetDOF(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	dof0 = [PetscInt(1)]
	dof1 = [PetscInt(1)]
	dof2 = [PetscInt(1)]
	dof3 = [PetscInt(1)]

	LibPETSc.DMStagGetDOF(
		PetscLib,
		dm,
		Ref(dof0,1),
		Ref(dof1,1),
		Ref(dof2,1),
		Ref(dof3,1),
	)

	return dof0[1],dof1[1],dof2[1],dof3[1]
end
 
 
"""
	entriesPerElement = DMStagGetEntries(dm::AbstractDMStag{PetscLib})

get number of native entries in the global representation

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameter:
===
- `entries` - number of rank-native entries in the global representation

Level: developer

Note:
This is the number of entries on this rank for a global vector associated with `dm`.
That is, it is value of `size` returned by `VecGetLocalSize(vec,&size)` when
`DMCreateGlobalVector(dm,&vec) is used to create a `Vec`. Users would typically
use these functions.

See also: 
=== 
`DMSTAG`, `DMStagGetDOF()`, `DMStagGetEntriesLocal()`, `DMStagGetEntriesPerElement()`, `DMCreateLocalVector()`

# External Links
$(_doc_external("DMStag/DMStagGetEntries"))
"""
function DMStagGetEntries(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	entriesPerElement = [PetscInt(1)]

	LibPETSc.DMStagGetEntries(
		PetscLib,
		dm,
		Ref(entriesPerElement,1),
	)

	return entriesPerElement[1]
end
 
 
"""
	entries = DMStagGetEntriesLocal(dm::AbstractDMStag{PetscLib})

get number of entries in the local representation

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameter:
===
- `entries` - number of entries in the local representation

Level: developer

Note:
This is the number of entries on this rank in the local representation.
That is, it is value of `size` returned by `VecGetSize(vec,&size)` or
`VecGetLocalSize(vec,&size)` when `DMCreateLocalVector(dm,&vec)` is used to
create a `Vec`. Users would typically use these functions.

See also: 
=== 
`DMSTAG`, `DMStagGetDOF()`, `DMStagGetEntries()`, `DMStagGetEntriesPerElement()`, `DMCreateLocalVector()`

# External Links
$(_doc_external("DMStag/DMStagGetEntriesLocal"))
"""
function DMStagGetEntriesLocal(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	entries = [PetscInt(1)]

	LibPETSc.DMStagGetEntriesLocal(
		PetscLib,
		dm,
		Ref(entries,1),
	)

	return entries[1]
end
 
 
"""
	entriesPerElement = DMStagGetEntriesPerElement(dm::AbstractDMStag{PetscLib})

get number of entries per element in the local representation

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameter:
===
- `entriesPerElement` - number of entries associated with each element in the local representation

Level: developer

Notes:
This is the natural block size for most local operations. In 1D it is equal to `dof0` + `dof1`,
in 2D it is equal to `dof0` + 2`dof1` + `dof2`, and in 3D it is equal to `dof0` + 3`dof1` + 3`dof2` + `dof3`

See also: 
=== 
`DMSTAG`, `DMStagGetDOF()`

# External Links
$(_doc_external("DMStag/DMStagGetEntriesPerElement"))
"""
function DMStagGetEntriesPerElement(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	entriesPerElement = [PetscInt(1)]

	LibPETSc.DMStagGetEntriesPerElement(
		PetscLib,
		dm,
		Ref(entriesPerElement,1),
	)

	return entriesPerElement[1]
end
 
 
"""
	x,y,z,m,n,p = DMStagGetGhostCorners(dm::AbstractDMStag{PetscLib})

return global element indices of the local region, including ghost points

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `x` - the starting element index in the first direction
- `y` - the starting element index in the second direction
- `z` - the starting element index in the third direction
- `m` - the element width in the first direction
- `n` - the element width in the second direction
- `p` - the element width in the third direction

Level: beginner

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagGetCorners()`, `DMDAGetGhostCorners()`

# External Links
$(_doc_external("DMStag/DMStagGetGhostCorners"))
"""
function DMStagGetGhostCorners(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	x = [PetscInt(1)]
	y = [PetscInt(1)]
	z = [PetscInt(1)]
	m = [PetscInt(1)]
	n = [PetscInt(1)]
	p = [PetscInt(1)]

	LibPETSc.DMStagGetGhostCorners(
		PetscLib,
		dm,
		Ref(x,1),
		Ref(y,1),
		Ref(z,1),
		Ref(m,1),
		Ref(n,1),
		Ref(p,1),
	)

	return x[1],y[1],z[1],m[1],n[1],p[1]
end
 
 
"""
	M,N,P = DMStagGetGlobalSizes(dm::AbstractDMStag{PetscLib})

get global element counts

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `M` - global element counts in the x direction
- `N` - global element counts in the y direction
- `P` - global element counts in the z direction

Level: beginner

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagGetLocalSizes()`, `DMDAGetInfo()`

# External Links
$(_doc_external("DMStag/DMStagGetGlobalSizes"))
"""
function DMStagGetGlobalSizes(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	M = [PetscInt(1)]
	N = [PetscInt(1)]
	P = [PetscInt(1)]

	LibPETSc.DMStagGetGlobalSizes(
		PetscLib,
		dm,
		Ref(M,1),
		Ref(N,1),
		Ref(P,1),
	)

	return M[1],N[1],P[1]
end
 
 
"""
	isFirstRank0,isFirstRank1,isFirstRank2 = DMStagGetIsFirstRank(dm::AbstractDMStag{PetscLib})

get boolean value for whether this rank is first in each direction in the rank grid

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `isFirstRank0` - whether this rank is first in the x direction
- `isFirstRank1` - whether this rank is first in the y direction
- `isFirstRank2` - whether this rank is first in the z direction

Level: intermediate

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagGetIsLastRank()`

# External Links
$(_doc_external("DMStag/DMStagGetIsFirstRank"))
"""
function DMStagGetIsFirstRank(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	isFirstRank0 = Ref{PetscBool}()
	isFirstRank1 = Ref{PetscBool}()
	isFirstRank2 = Ref{PetscBool}()

	LibPETSc.DMStagGetIsFirstRank(
		PetscLib,
		dm,
		isFirstRank0,
		isFirstRank1,
		isFirstRank2,
	)

	return isFirstRank0[] == PETSC_TRUE,isFirstRank1[] == PETSC_TRUE,isFirstRank2[] == PETSC_TRUE
end
 
 
"""
	isLastRank0,isLastRank1,isLastRank2 = DMStagGetIsLastRank(dm::AbstractDMStag{PetscLib})

get boolean value for whether this rank is last in each direction in the rank grid

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `isLastRank0` - whether this rank is last in the x direction
- `isLastRank1` - whether this rank is last in the y direction
- `isLastRank2` - whether this rank is last in the z direction

Level: intermediate

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagGetIsFirstRank()`

# External Links
$(_doc_external("DMStag/DMStagGetIsLastRank"))
"""
function DMStagGetIsLastRank(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	isLastRank0 = Ref{PetscBool}()
	isLastRank1 = Ref{PetscBool}()
	isLastRank2 = Ref{PetscBool}()

	LibPETSc.DMStagGetIsLastRank(
		PetscLib,
		dm,
		isLastRank0,
		isLastRank1,
		isLastRank2,
	)

	return isLastRank0[] == PETSC_TRUE,isLastRank1[] == PETSC_TRUE,isLastRank2[] == PETSC_TRUE
end
 
 
"""
	m,n,p = DMStagGetLocalSizes(dm::AbstractDMStag{PetscLib})

get local elementwise sizes

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `m` - local element counts (excluding ghosts) in the x direction
- `n` - local element counts (excluding ghosts) in the y direction
- `p` - local element counts (excluding ghosts) in the z direction

Level: beginner

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagGetGlobalSizes()`, `DMStagGetDOF()`, `DMStagGetNumRanks()`, `DMDAGetLocalInfo()`

# External Links
$(_doc_external("DMStag/DMStagGetLocalSizes"))
"""
function DMStagGetLocalSizes(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	m = [PetscInt(1)]
	n = [PetscInt(1)]
	p = [PetscInt(1)]

	LibPETSc.DMStagGetLocalSizes(
		PetscLib,
		dm,
		Ref(m,1),
		Ref(n,1),
		Ref(p,1),
	)

	return m[1],n[1],p[1]
end
 
 
"""
	dof = DMStagGetLocationDOF(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation)

Get number of DOF associated with a given point in a `DMSTAG` grid

Not Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `loc` - grid point (see `DMStagStencilLocation`)

Output Parameter:
===
- `dof` - the number of DOF (components) living at `loc` in `dm`

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMStagStencilLocation`, `DMStagStencil`, `DMDAGetDof()`

# External Links
$(_doc_external("DMStag/DMStagGetLocationDOF"))
"""
function DMStagGetLocationDOF(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	dof = [PetscInt(1)]

	LibPETSc.DMStagGetLocationDOF(
		PetscLib,
		dm,
		loc,
		Ref(dof,1),
	)

	return dof[1]
end
 
 
"""
	slot = DMStagGetLocationSlot(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation,c::Int)

get index to use in accessing raw local arrays

Not Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `loc` - location relative to an element
- `c`   - component

Output Parameter:
===
- `slot` - index to use

Level: beginner

Notes:
Provides an appropriate index to use with `DMStagVecGetArray()` and friends.
This is required so that the user doesn't need to know about the ordering of
dof associated with each local element.

See also: 
=== 
`DMSTAG`, `DMStagVecGetArray()`, `DMStagVecGetArrayRead()`, `DMStagGetDOF()`, `DMStagGetEntriesPerElement()`

# External Links
$(_doc_external("DMStag/DMStagGetLocationSlot"))
"""
function DMStagGetLocationSlot(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation,c::Int) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	slot = [PetscInt(1)]

	LibPETSc.DMStagGetLocationSlot(
		PetscLib,
		dm,
		loc,
		c,
		Ref(slot,1),
	)

	return slot[1]
end
 
 
"""
	nRanks0,nRanks1,nRanks2 = DMStagGetNumRanks(dm::AbstractDMStag{PetscLib})

get number of ranks in each direction in the global grid decomposition

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `nRanks0` - number of ranks in the x direction in the grid decomposition
- `nRanks1` - number of ranks in the y direction in the grid decomposition
- `nRanks2` - number of ranks in the z direction in the grid decomposition

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMStagGetGlobalSizes()`, `DMStagGetLocalSize()`, `DMStagSetNumRanks()`, `DMDAGetInfo()`

# External Links
$(_doc_external("DMStag/DMStagGetNumRanks"))
"""
function DMStagGetNumRanks(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	nRanks0 = [PetscInt(1)]
	nRanks1 = [PetscInt(1)]
	nRanks2 = [PetscInt(1)]

	LibPETSc.DMStagGetNumRanks(
		PetscLib,
		dm,
		Ref(nRanks0,1),
		Ref(nRanks1,1),
		Ref(nRanks2,1),
	)

	return nRanks0[1],nRanks1[1],nRanks2[1]
end
 
 
"""
	lx,ly,lz = DMStagGetOwnershipRanges(dm::AbstractDMStag{PetscLib})

get elements per rank in each direction

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `lx` - ownership along x direction (optional)
- `ly` - ownership along y direction (optional)
- `lz` - ownership along z direction (optional)

Level: intermediate

Notes:
These correspond to the optional final arguments passed to `DMStagCreate1d()`, `DMStagCreate2d()`, and `DMStagCreate3d()`.

Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

In C you should not free these arrays, nor change the values in them.
They will only have valid values while the `DMSTAG` they came from still exists (has not been destroyed).

See also: 
=== 
`DMSTAG`, `DMStagSetGlobalSizes()`, `DMStagSetOwnershipRanges()`, `DMStagCreate1d()`, `DMStagCreate2d()`, `DMStagCreate3d()`, `DMDAGetOwnershipRanges()`

# External Links
$(_doc_external("DMStag/DMStagGetOwnershipRanges"))
"""
function DMStagGetOwnershipRanges(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	n_lx,n_ly,n_lz = DMStagGetNumRanks(dm::AbstractDMStag{PetscLib})
	nX = (n_lx,)
	nY = (n_ly,)
	nZ = (n_lz,)
	
	r_lx = PETSc_RefPtr(nX, PetscInt)
	r_ly = PETSc_RefPtr(nY, PetscInt)
	r_lz = PETSc_RefPtr(nZ, PetscInt)

	LibPETSc.DMStagGetOwnershipRanges(
		PetscLib,
		dm,
		r_lx,
		r_ly,
		r_lz,
	)

	arrX = PETSc_unsafe_wrap(r_lx, nX; own=false)
	arrY = PETSc_unsafe_wrap(r_ly, nY; own=false)
	arrZ = PETSc_unsafe_wrap(r_lz, nZ; own=false)

	return arrX, arrY, arrZ
end
 
 
"""
	arrX,arrY,arrZ,p_arrx,p_arry,p_arrz = DMStagGetProductCoordinateArrays(dm::AbstractDMStag{PetscLib})

extract local product coordinate arrays, one per dimension

Logically Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `arrX` - local 1D coordinate arrays for x direction
- `arrY` - local 1D coordinate arrays for y direction, not set for one dimensional problems
- `arrZ` - local 1D coordinate arrays for z direction, not set for one and two dimensional problems
- `p_arrx` - pointer to local 1D coordinate arrays for x direction
- `p_arry` - pointer to local 1D coordinate arrays for y direction
- `p_arrz` - pointer to local 1D coordinate arrays for z direction

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
function DMStagGetProductCoordinateArrays(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	PetscInt    = PetscLib.PetscInt
	dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	nX = (dims[1],PetscInt(1))
	nY = (dims[2],PetscInt(1))
	nZ = (dims[3],PetscInt(1))
	r_arrX = PETSc_RefPtr(nX, PetscScalar)
	r_arrY = PETSc_RefPtr(nY, PetscScalar)
	r_arrZ = PETSc_RefPtr(nZ, PetscScalar)

	LibPETSc.DMStagGetProductCoordinateArrays(
		PetscLib,
		dm,
		r_arrX,
		r_arrY,
		r_arrZ,
	)

	arrX = PETSc_unsafe_wrap(r_arrX, nX; own=false)
	arrY = PETSc_unsafe_wrap(r_arrY, nY; own=false)
	arrZ = PETSc_unsafe_wrap(r_arrZ, nZ; own=false)

	return arrX[:],arrY[:],arrZ[:], r_arrX, r_arrY, r_arrZ
end
 
 
"""
	arrX,arrY,arrZ,p_arrx,p_arry,p_arrz = DMStagGetProductCoordinateArraysRead(dm::AbstractDMStag{PetscLib})

extract product coordinate arrays, read

Logically Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `arrX` - local 1D coordinate arrays for x direction
- `arrY` - local 1D coordinate arrays for y direction, not set for one dimensional problems
- `arrZ` - local 1D coordinate arrays for z direction, not set for one and two dimensional problems
- `p_arrx` - pointer to local 1D coordinate arrays for x direction
- `p_arry` - pointer to local 1D coordinate arrays for y direction
- `p_arrz` - pointer to local 1D coordinate arrays for z direction

Level: intermediate

Note:
See `DMStagGetProductCoordinateArrays()` for more information.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagGetProductCoordinateArrays()`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesProduct()`, `DMStagGetProductCoordinateLocationSlot()`

# External Links
$(_doc_external("DMStag/DMStagGetProductCoordinateArraysRead"))
"""
function DMStagGetProductCoordinateArraysRead(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	PetscInt    = PetscLib.PetscInt
	
	dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	nX = (dims[1],PetscInt(1))
	nY = (dims[2],PetscInt(1))
	nZ = (dims[3],PetscInt(1))
	r_arrX = PETSc_RefPtr(nX, PetscScalar)
	r_arrY = PETSc_RefPtr(nY, PetscScalar)
	r_arrZ = PETSc_RefPtr(nZ, PetscScalar)
	
	LibPETSc.DMStagGetProductCoordinateArraysRead(
		PetscLib,
		dm,
		r_arrX,
		r_arrY,
		r_arrZ,
	)

	arrX = PETSc_unsafe_wrap(r_arrX, nX; own=false)
	arrY = PETSc_unsafe_wrap(r_arrY, nY; own=false)
	arrZ = PETSc_unsafe_wrap(r_arrZ, nZ; own=false)

	return arrX[:],arrY[:],arrZ[:], r_arrX, r_arrY, r_arrZ
end
 
 
"""
	slot = DMStagGetProductCoordinateLocationSlot(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation)

get slot for use with local product coordinate arrays

Not Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `loc` - the grid location

Output Parameter:
===
- `slot` - the index to use in local arrays

Level: intermediate

Notes:
High-level helper function to get slot indices for 1D coordinate `DM`s,
for use with `DMStagGetProductCoordinateArrays()` and related functions.

For `loc`, one should use `DMSTAG_LEFT`, `DMSTAG_ELEMENT`, or `DMSTAG_RIGHT` for "previous", "center" and "next"
locations, respectively, in each dimension.
One can equivalently use `DMSTAG_DOWN` or `DMSTAG_BACK` in place of `DMSTAG_LEFT`,
and `DMSTAG_UP` or `DMSTACK_FRONT` in place of `DMSTAG_RIGHT`;

This function checks that the coordinates are actually set up so that using the
slots from any of the 1D coordinate sub-`DM`s are valid for all the 1D coordinate sub-`DM`s.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagGetProductCoordinateArrays()`, `DMStagGetProductCoordinateArraysRead()`, `DMStagSetUniformCoordinates()`

# External Links
$(_doc_external("DMStag/DMStagGetProductCoordinateLocationSlot"))
"""
function DMStagGetProductCoordinateLocationSlot(dm::AbstractDMStag{PetscLib},loc::DMStagStencilLocation) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	slot = [PetscInt(1)]

	LibPETSc.DMStagGetProductCoordinateLocationSlot(
		PetscLib,
		dm,
		loc,
		Ref(slot,1),
	)

	return slot[1]
end
 
 
"""
	stencilType = DMStagGetStencilType(dm::AbstractDMStag{PetscLib})

get elementwise ghost/halo stencil type

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameter:
===
- `stencilType` - the elementwise ghost stencil type: `DMSTAG_STENCIL_BOX`, `DMSTAG_STENCIL_STAR`, or `DMSTAG_STENCIL_NONE`

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagSetStencilType()`, `DMStagGetStencilWidth`, `DMStagStencilType`

# External Links
$(_doc_external("DMStag/DMStagGetStencilType"))
"""
function DMStagGetStencilType(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	stencilType = Ref{DMStagStencilType}()

	LibPETSc.DMStagGetStencilType(
		PetscLib,
		dm,
		stencilType,
	)

	return string(stencilType[])
end
 
 
"""
	stencilWidth = DMStagGetStencilWidth(dm::AbstractDMStag{PetscLib})

get elementwise stencil width

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameter:
===
- `stencilWidth` - stencil/halo/ghost width in elements

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagSetStencilWidth()`, `DMStagGetStencilType()`, `DMDAGetStencilType()`

# External Links
$(_doc_external("DMStag/DMStagGetStencilWidth"))
"""
function DMStagGetStencilWidth(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	stencilWidth = [PetscInt(1)]

	LibPETSc.DMStagGetStencilWidth(
		PetscLib,
		dm,
		Ref(stencilWidth,1),
	)

	return stencilWidth[1]
end
 
 
"""
	refine_x,refine_y,refine_z = DMStagGetRefinementFactor(dm::AbstractDMStag{PetscLib})

get refinement ratios in each direction

Not Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object

Output Parameters:
===
- `refine_x` - ratio of fine grid to coarse in x-direction (2 by default)
- `refine_y` - ratio of fine grid to coarse in y-direction (2 by default)
- `refine_z` - ratio of fine grid to coarse in z-direction (2 by default)

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMRefine()`, `DMCoarsen()`, `DMStagSetRefinementFactor()`, `DMDASetRefinementFactor()`

# External Links
$(_doc_external("DMStag/DMStagGetRefinementFactor"))
"""
function DMStagGetRefinementFactor(dm::AbstractDMStag{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	refine_x = [PetscInt(1)]
	refine_y = [PetscInt(1)]
	refine_z = [PetscInt(1)]

	LibPETSc.DMStagGetRefinementFactor(
		PetscLib,
		dm,
		Ref(refine_x,1),
		Ref(refine_y,1),
		Ref(refine_z,1),
	)

	return refine_x[1],refine_y[1],refine_z[1]
end
 
 
"""
	val = DMStagMatGetValuesStencil(dm::AbstractDMStag{PetscLib},mat::AbstractMatrix,nRow::Int,posRow::Vector{DMStagStencil{Int}},nCol::Int,posCol::Vector{DMStagStencil{Int}})

retrieve local matrix entries using grid indexing

Not Collective

Input Parameters:
===
- `dm`     - the `DMSTAG` object
- `mat`    - the matrix
- `nRow`   - number of rows
- `posRow` - grid locations (including components) of rows
- `nCol`   - number of columns
- `posCol` - grid locations (including components) of columns

Output Parameter:
===
- `val` - logically two-dimensional array of values

Level: advanced

See also: 
=== 
`DMSTAG`, `DMStagStencil`, `DMStagStencilLocation`, `DMStagVecGetValuesStencil()`, `DMStagVecSetValuesStencil()`, `DMStagMatSetValuesStencil()`, `MatSetValuesStencil()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `DMCreateMatrix()`

# External Links
$(_doc_external("DMStag/DMStagMatGetValuesStencil"))
"""
function DMStagMatGetValuesStencil(dm::AbstractDMStag{PetscLib},mat::AbstractMatrix,nRow::Int,posRow::Vector{DMStagStencil{Int}},nCol::Int,posCol::Vector{DMStagStencil{Int}}) where {PetscLib}
	@assert length(posRow) == nRow 
	@assert length(posCol) == nCol
	
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscScalar = PetscLib.PetscScalar
	val = zeros(PetscScalar, nRow*nCol)

	LibPETSc.DMStagMatGetValuesStencil(
		PetscLib,
		dm,
		mat,
		nRow,
		posRow,
		nCol,
		posCol,
		Ref(val,nRow*nCol),
	)

	return val
end

DMStagMatGetValuesStencil(dm::AbstractDMStag{PetscLib},mat::AbstractMatrix,posRow::Vector{DMStagStencil{Int}},posCol::Vector{DMStagStencil{Int}}) where {PetscLib} =DMStagMatGetValuesStencil(dm,mat,length(posRow),posRow,length(posCol),posCol)

 
"""
	 DMStagMatSetValuesStencil(dm::AbstractDMStag{PetscLib},mat::AbstractMatrix,nRow::Int,posRow::Vector{DMStagStencil{Int}},nCol::Int,posCol::Vector{DMStagStencil{Int}},val::Vector{<:AbstractFloat},insertMode::InsertMode)

insert or add matrix entries using grid indexing

Not Collective

Input Parameters:
===
- `dm`         - the `DMSTAG` object
- `mat`        - the matrix
- `nRow`       - number of rows
- `posRow`     - grid locations (including components) of rows
- `nCol`       - number of columns
- `posCol`     - grid locations (including components) of columns
- `val`        - logically two-dimensional array of values
- `insertMode` - `INSERT_VALUES` or `ADD_VALUES`

Notes:
See notes for `MatSetValuesStencil()`

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMStagStencil`, `DMStagStencilLocation`, `DMStagVecGetValuesStencil()`, `DMStagVecSetValuesStencil()`, `DMStagMatGetValuesStencil()`, `MatSetValuesStencil()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `DMCreateMatrix()`

# External Links
$(_doc_external("DMStag/DMStagMatSetValuesStencil"))
"""
function DMStagMatSetValuesStencil(dm::AbstractDMStag{PetscLib},mat::AbstractMatrix,nRow::Int,posRow::Vector{DMStagStencil{Int}},nCol::Int,posCol::Vector{DMStagStencil{Int}},val::Vector{<:AbstractFloat},insertMode::InsertMode) where {PetscLib}
	@assert length(val) == nRow*nCol 
	LibPETSc.DMStagMatSetValuesStencil(
		PetscLib,
		dm,
		mat,
		nRow,
		posRow,
		nCol,
		posCol,
		val,
		insertMode,
	)

	return nothing
end
 
 
"""
	 DMStagMigrateVec(dm::AbstractDMStag{PetscLib},vec::AbstractVector,dmTo::AbstractDMStag{PetscLib},vecTo::AbstractVector)

transfer a vector associated with a `DMSTAG` to a vector associated with a compatible `DMSTAG`

Collective

Input Parameters:
===
- `dm`    - the source `DMSTAG` object
- `vec`   - the source vector, compatible with `dm`
- `dmTo`  - the compatible destination `DMSTAG` object
- `vecTo` - the destination vector, compatible with `dmTo`

Level: advanced

Notes:
Extra dof are ignored, and unfilled dof are zeroed.
Currently only implemented to migrate global vectors to global vectors.
For the definition of compatibility of `DM`s, see `DMGetCompatibility()`.

See also: 
=== 
`DMSTAG`, `DMStagCreateCompatibleDMStag()`, `DMGetCompatibility()`, `DMStagVecSplitToDMDA()`

# External Links
$(_doc_external("DMStag/DMStagMigrateVec"))
"""
function DMStagMigrateVec(dm::AbstractDMStag{PetscLib},vec::AbstractVector,dmTo::AbstractDMStag{PetscLib},vecTo::AbstractVector) where {PetscLib}

	LibPETSc.DMStagMigrateVec(
		PetscLib,
		dm,
		vec,
		dmTo,
		vecTo,
	)

	return nothing
end
 
 
"""
	 DMStagPopulateLocalToGlobalInjective(dm::AbstractDMStag{PetscLib})

populate an internal 1

Collective

Creates an internal object which explicitly maps a single local degree of
freedom to each global degree of freedom. This is used, if populated,
instead of SCATTER_REVERSE_LOCAL with the (1-to-many, in general)
global-to-local map, when DMLocalToGlobal() is called with INSERT_VALUES.
This allows usage, for example, even in the periodic, 1-rank case, where
the inverse of the global-to-local map, even when restricted to on-rank
communication, is non-injective. This is at the cost of storing an additional
VecScatter object inside each `DMSTAG` object.

Input Parameter:
===
- `dm` - the `DMSTAG` object

Level: developer

Notes:
In normal usage, library users shouldn't be concerned with this function,
as it is called during `DMSetUp()`, when required.

Returns immediately if the internal map is already populated.

Developer Notes:
This could, if desired, be moved up to a general `DM` routine. It would allow,
for example, `DMDA` to support `DMLocalToGlobal()` with `INSERT_VALUES`,
even in the single-rank periodic case.

See also: 
=== 
`DMSTAG`, `DMLocalToGlobal()`, `VecScatter`

# External Links
$(_doc_external("DMStag/DMStagPopulateLocalToGlobalInjective"))
"""
function DMStagPopulateLocalToGlobalInjective(dm::AbstractDMStag{PetscLib}) where {PetscLib}

	LibPETSc.DMStagPopulateLocalToGlobalInjective(
		PetscLib,
		dm,
	)

	return nothing
end
 
 
"""
	arrX,arrY,arrZ = DMStagRestoreProductCoordinateArrays(dm::AbstractDMStag{PetscLib}, p_arrX, p_arrY, p_arrZ)

restore local array access

Logically Collective

Input Parameter:
===
- `dm` - the `DMSTAG` object
- `p_arrX` - pointer to local 1D coordinate arrays for x direction
- `p_arrY` - pointer to local 1D coordinate arrays for y direction
- `p_arrZ` - pointer to local 1D coordinate arrays for z direction

Level: intermediate

Notes:
This function does not automatically perform a local->global scatter to populate global coordinates from the local coordinates.
Thus, it may be required to explicitly perform these operations in some situations, as in the following partial example:
```
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
```

See also: 
=== 
`DMSTAG`, `DMStagGetProductCoordinateArrays()`, `DMStagGetProductCoordinateArraysRead()`

# External Links
$(_doc_external("DMStag/DMStagRestoreProductCoordinateArrays"))
"""
function DMStagRestoreProductCoordinateArrays(dm::AbstractDMStag{PetscLib}, r_arrX, r_arrY, r_arrZ) where {PetscLib}
	
	LibPETSc.DMStagRestoreProductCoordinateArrays(
		PetscLib,
		dm,
		r_arrX,
		r_arrY,
		r_arrZ,
	)

	return nothing
end
 
 
"""
	arrX,arrY,arrZ = DMStagRestoreProductCoordinateArraysRead(dm::AbstractDMStag{PetscLib}, p_arrX, p_arrY, p_arrZ)

restore local product array access, read

Logically Collective

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `p_arrX` - pointer to local 1D coordinate arrays for x direction
- `p_arrY` - pointer to local 1D coordinate arrays for y direction
- `p_arrZ` - pointer to local 1D coordinate arrays for z direction

Level: intermediate

See also: 
=== 
`DMSTAG`, `DMStagGetProductCoordinateArrays()`, `DMStagGetProductCoordinateArraysRead()`

# External Links
$(_doc_external("DMStag/DMStagRestoreProductCoordinateArraysRead"))
"""
function DMStagRestoreProductCoordinateArraysRead(dm::AbstractDMStag{PetscLib},  r_arrX, r_arrY, r_arrZ) where {PetscLib}
	LibPETSc.DMStagRestoreProductCoordinateArraysRead(
		PetscLib,
		dm,
		r_arrX,
		r_arrY,
		r_arrZ,
	)

	return nothing
end
 
 
"""
	 xc = DMStagRestrictSimple(dmf::AbstractDMStag{PetscLib},xf::AbstractVector,dmc::AbstractDMStag{PetscLib})

restricts data from a fine to a coarse `DMSTAG`, in the simplest way

Values on coarse cells are averages of all fine cells that they cover.
Thus, values on vertices are injected, values on edges are averages
of the underlying two fine edges, and values on elements in
d dimensions are averages of 2^d underlying elements.

Input Parameters:
===
- `dmf` - fine `DM`
- `xf`  - data on fine `DM`
- `dmc` - coarse `DM`

Output Parameter:
===
- `xc` - data on coarse `DM`

Level: advanced

See also: 
=== 
`DMSTAG`, `DM`, `DMRestrict()`, `DMCoarsen()`, `DMCreateInjection()`

# External Links
$(_doc_external("DMStag/DMStagRestrictSimple"))
"""
function DMStagRestrictSimple(dmf::AbstractDMStag{PetscLib},xf::AbstractVector,dmc::AbstractDMStag{PetscLib}) where {PetscLib}

	xc = DMLocalVec(dmc)      

	LibPETSc.DMStagRestrictSimple(
		PetscLib,
		dmf,
		xf,
		dmc,
		xc
	)

	return xc
end
 
 
"""
	 DMStagSetBoundaryTypes(dm::AbstractDMStag{PetscLib},boundaryType0::DMBoundaryType,boundaryType1::DMBoundaryType,boundaryType2::DMBoundaryType)

set `DMSTAG` boundary types

Logically Collective; boundaryType0, boundaryType1, and boundaryType2 must contain common values

Input Parameters:
===
- `dm`            - the `DMSTAG` object
- `boundaryType2` - boundary type for x direction
- `boundaryType1` - boundary type for y direction, not set for one dimensional problems
- `boundaryType0` - boundary type for z direction, not set for one and two dimensional problems

Level: advanced

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

See also: 
=== 
`DMSTAG`, `DMBoundaryType`, `DMStagCreate1d()`, `DMStagCreate2d()`, `DMStagCreate3d()`, `DMDASetBoundaryType()`

# External Links
$(_doc_external("DMStag/DMStagSetBoundaryTypes"))
"""
function DMStagSetBoundaryTypes(dm::AbstractDMStag{PetscLib},boundaryType0::DMBoundaryType,boundaryType1::DMBoundaryType,boundaryType2::DMBoundaryType) where {PetscLib}

	LibPETSc.DMStagSetBoundaryTypes(
		PetscLib,
		dm,
		boundaryType0,
		boundaryType1,
		boundaryType2,
	)

	return nothing
end
 
 
"""
	 DMStagSetCoordinateDMType(dm::AbstractDMStag{PetscLib},dmtype::DMType)

set DM type to store coordinates

Logically Collective; `dmtype` must contain common value

Input Parameters:
===
- `dm`     - the `DMSTAG` object
- `dmtype` - `DMtype` for coordinates, either `DMSTAG` or `DMPRODUCT`

Level: advanced

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMGetCoordinateDM()`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesExplicit()`, `DMStagSetUniformCoordinatesProduct()`, `DMType`

# External Links
$(_doc_external("DMStag/DMStagSetCoordinateDMType"))
"""
function DMStagSetCoordinateDMType(dm::AbstractDMStag{PetscLib},dmtype::DMType) where {PetscLib}

	LibPETSc.DMStagSetCoordinateDMType(
		PetscLib,
		dm,
		dmtype,
	)

	return nothing
end
 
 
"""
	 DMStagSetDOF(dm::AbstractDMStag{PetscLib},dof0::Int,dof1::Int,dof2::Int,dof3::Int)

set dof/stratum

Logically Collective; `dof0`, `dof1`, `dof2`, and `dof3` must contain common values

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `dof0` - the number of points per 0-cell (vertex/node)
- `dof1` - the number of points per 1-cell (element in 1D, edge in 2D and 3D)
- `dof2` - the number of points per 2-cell (element in 2D, face in 3D)
- `dof3` - the number of points per 3-cell (element in 3D)

Level: advanced

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

See also: 
=== 
`DMSTAG`, `DMDASetDof()`

# External Links
$(_doc_external("DMStag/DMStagSetDOF"))
"""
function DMStagSetDOF(dm::AbstractDMStag{PetscLib},dof0::Int,dof1::Int,dof2::Int,dof3::Int) where {PetscLib}

	LibPETSc.DMStagSetDOF(
		PetscLib,
		dm,
		dof0,
		dof1,
		dof2,
		dof3,
	)

	return nothing
end
 
 
"""
	 DMStagSetGlobalSizes(dm::AbstractDMStag{PetscLib},N0::Int,N1::Int,N2::Int)

set global element counts in each direction

Logically Collective; `N0`, `N1`, and `N2` must contain common values

Input Parameters:
===
- `dm` - the `DMSTAG` object
- `N0` - global elementwise size in the x direction
- `N1` - global elementwise size in the y direction
- `N2` - global elementwise size in the z direction

Level: advanced

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

See also: 
=== 
`DMSTAG`, `DMStagGetGlobalSizes()`, `DMDASetSizes()`

# External Links
$(_doc_external("DMStag/DMStagSetGlobalSizes"))
"""
function DMStagSetGlobalSizes(dm::AbstractDMStag{PetscLib},N0::Int,N1::Int,N2::Int) where {PetscLib}

	LibPETSc.DMStagSetGlobalSizes(
		PetscLib,
		dm,
		N0,
		N1,
		N2,
	)

	return nothing
end
 
 
"""
	 DMStagSetNumRanks(dm::AbstractDMStag{PetscLib},nRanks0::Int,nRanks1::Int,nRanks2::Int)

set ranks in each direction in the global rank grid

Logically Collective; `nRanks0`, `nRanks1`, and `nRanks2` must contain common values

Input Parameters:
===
- `dm`      - the `DMSTAG` object
- `nRanks0` - number of ranks in the x direction
- `nRanks1` - number of ranks in the y direction
- `nRanks2` - number of ranks in the z direction

Level: developer

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

See also: 
=== 
`DMSTAG`, `DMDASetNumProcs()`

# External Links
$(_doc_external("DMStag/DMStagSetNumRanks"))
"""
function DMStagSetNumRanks(dm::AbstractDMStag{PetscLib},nRanks0::Int,nRanks1::Int,nRanks2::Int) where {PetscLib}

	LibPETSc.DMStagSetNumRanks(
		PetscLib,
		dm,
		nRanks0,
		nRanks1,
		nRanks2,
	)

	return nothing
end
 
 
"""
	 DMStagSetOwnershipRanges(dm::AbstractDMStag{PetscLib},lx::Vector{Int},ly::Vector{Int},lz::Vector{Int})

set elements per rank in each direction

Logically Collective; `lx`, `ly`, and `lz` must contain common values

Input Parameters:
===
- `dm` - the `DMSTAG` object
- `lx` - element counts for each rank in the x direction, may be `NULL`
- `ly` - element counts for each rank in the y direction, may be `NULL`
- `lz` - element counts for each rank in the z direction, may be `NULL`

Level: developer

Note:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

See also: 
=== 
`DMSTAG`, `DMStagSetGlobalSizes()`, `DMStagGetOwnershipRanges()`, `DMDASetOwnershipRanges()`

# External Links
$(_doc_external("DMStag/DMStagSetOwnershipRanges"))
"""
function DMStagSetOwnershipRanges(dm::AbstractDMStag{PetscLib},lx::Vector{Int},ly::Vector{Int},lz::Vector{Int}) where {PetscLib}
	LibPETSc.DMStagSetOwnershipRanges(
		PetscLib,
		dm,
		lx,
		ly,
		lz,
	)

	return nothing
end
 
 
"""
	 DMStagSetStencilType(dm::AbstractDMStag{PetscLib},stencilType::DMStagStencilType)

set elementwise ghost/halo stencil type

Logically Collective; `stencilType` must contain common value

Input Parameters:
===
- `dm`          - the `DMSTAG` object
- `stencilType` - the elementwise ghost stencil type: `DMSTAG_STENCIL_BOX`, `DMSTAG_STENCIL_STAR`, or `DMSTAG_STENCIL_NONE`

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagGetStencilType()`, `DMStagSetStencilWidth()`, `DMStagStencilType`

# External Links
$(_doc_external("DMStag/DMStagSetStencilType"))
"""
function DMStagSetStencilType(dm::AbstractDMStag{PetscLib},stencilType::DMStagStencilType) where {PetscLib}

	LibPETSc.DMStagSetStencilType(
		PetscLib,
		dm,
		stencilType,
	)

	return nothing
end
 
 
"""
	 DMStagSetStencilWidth(dm::AbstractDMStag{PetscLib},stencilWidth::Int)

set elementwise stencil width

Logically Collective; `stencilWidth` must contain common value

Input Parameters:
===
- `dm`           - the `DMSTAG` object
- `stencilWidth` - stencil/halo/ghost width in elements

Level: beginner

Note:
The width value is not used when `DMSTAG_STENCIL_NONE` is specified.

See also: 
=== 
`DMSTAG`, `DMStagGetStencilWidth()`, `DMStagGetStencilType()`, `DMStagStencilType`

# External Links
$(_doc_external("DMStag/DMStagSetStencilWidth"))
"""
function DMStagSetStencilWidth(dm::AbstractDMStag{PetscLib},stencilWidth::Int) where {PetscLib}

	LibPETSc.DMStagSetStencilWidth(
		PetscLib,
		dm,
		stencilWidth,
	)

	return nothing
end
 
 
"""
	 DMStagSetRefinementFactor(dm::AbstractDMStag{PetscLib},refine_x::Int,refine_y::Int,refine_z::Int)

set refinement ratios in each direction

Logically Collective

Input Parameters:
===
- `dm`       - the `DMSTAG` object
- `refine_x` - ratio of fine grid to coarse in x-direction (2 by default)
- `refine_y` - ratio of fine grid to coarse in y-direction (2 by default)
- `refine_z` - ratio of fine grid to coarse in z-direction (2 by default)

Level: intermediate

Note:
Pass `PETSC_IGNORE` to leave a value unchanged

See also: 
=== 
`DMSTAG`, `DMRefine()`, `DMCoarsen()`, `DMStagGetRefinementFactor()`, `DMDAGetRefinementFactor()`

# External Links
$(_doc_external("DMStag/DMStagSetRefinementFactor"))
"""
function DMStagSetRefinementFactor(dm::AbstractDMStag{PetscLib},refine_x::Int,refine_y::Int,refine_z::Int) where {PetscLib}

	LibPETSc.DMStagSetRefinementFactor(
		PetscLib,
		dm,
		refine_x,
		refine_y,
		refine_z,
	)

	return nothing
end
 
 
"""
	 DMStagSetUniformCoordinates(dm::AbstractDMStag{PetscLib},xmin<:AbstractFloat,xmax<:AbstractFloat,ymin<:AbstractFloat,ymax<:AbstractFloat,zmin<:AbstractFloat,zmax<:AbstractFloat)

set `DMSTAG` coordinates to be a uniform grid

Collective

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `xmin` - minimum global coordinate value in the x direction
- `xmax` - maximum global coordinate values in the x direction
- `ymin` - minimum global coordinate value in the y direction
- `ymax` - maximum global coordinate value in the y direction
- `zmin` - minimum global coordinate value in the z direction
- `zmax` - maximum global coordinate value in the z direction

Level: advanced

Notes:
`DMSTAG` supports 2 different types of coordinate `DM`: `DMSTAG` and `DMPRODUCT`.
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

Local coordinates are populated (using `DMSetCoordinatesLocal()`), linearly
extrapolated to ghost cells, including those outside the physical domain.
This is also done in case of periodic boundaries, meaning that the same
global point may have different coordinates in different local representations,
which are equivalent assuming a periodicity implied by the arguments to this function,
i.e. two points are equivalent if their difference is a multiple of (`xmax` - `xmin` )
in the x direction, ( `ymax` - `ymin` ) in the y direction, and ( `zmax` - `zmin` ) in the z direction.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagSetUniformCoordinatesExplicit()`, `DMStagSetUniformCoordinatesProduct()`, `DMStagSetCoordinateDMType()`, `DMGetCoordinateDM()`, `DMGetCoordinates()`, `DMDASetUniformCoordinates()`, `DMBoundaryType`

# External Links
$(_doc_external("DMStag/DMStagSetUniformCoordinates"))
"""
function DMStagSetUniformCoordinates(dm::AbstractDMStag{PetscLib},xmin,xmax,ymin,ymax,zmin,zmax) where {PetscLib}

	LibPETSc.DMStagSetUniformCoordinates(
		PetscLib,
		dm,
		xmin,
		xmax,
		ymin,
		ymax,
		zmin,
		zmax
	)

	return nothing
end 
 
"""
	 DMStagSetUniformCoordinatesExplicit(dm::AbstractDMStag{PetscLib},xmin<:AbstractFloat,xmax<:AbstractFloat,ymin<:AbstractFloat,ymax<:AbstractFloat,zmin<:AbstractFloat,zmax<:AbstractFloat)

set `DMSTAG` coordinates to be a uniform grid, storing all values

Collective

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `xmin` - minimum global coordinate value in the x direction
- `xmax` - maximum global coordinate value in the x direction
- `ymin` - minimum global coordinate value in the y direction
- `ymax` - maximum global coordinate value in the y direction
- `zmin` - minimum global coordinate value in the z direction
- `zmax` - maximum global coordinate value in the z direction

Level: beginner

Notes:
`DMSTAG` supports 2 different types of coordinate `DM`: either another `DMSTAG`, or a `DMPRODUCT`.
If the grid is orthogonal, using `DMPRODUCT` should be more efficient.

Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

See the manual page for `DMStagSetUniformCoordinates()` for information on how
coordinates for dummy cells outside the physical domain boundary are populated.

See also: 
=== 
`DMSTAG`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesProduct()`, `DMStagSetCoordinateDMType()`

# External Links
$(_doc_external("DMStag/DMStagSetUniformCoordinatesExplicit"))
"""
function DMStagSetUniformCoordinatesExplicit(dm::AbstractDMStag{PetscLib},xmin,xmax,ymin,ymax,zmin,zmax) where {PetscLib}

	LibPETSc.DMStagSetUniformCoordinatesExplicit(
		PetscLib,
		dm,
		xmin,
		xmax,
		ymin,
		ymax,
		zmin,
		zmax
	)

	return nothing
end

 
"""
	 DMStagSetUniformCoordinatesProduct(dm::AbstractDMStag{PetscLib},xmin<:AbstractFloat,xmax<:AbstractFloat,ymin<:AbstractFloat,ymax<:AbstractFloat,zmin<:AbstractFloat,zmax<:AbstractFloat)

create uniform coordinates, as a product of 1D arrays

Set the coordinate `DM` to be a `DMPRODUCT` of 1D `DMSTAG` objects, each of which have a coordinate `DM` (also a 1d `DMSTAG`) holding uniform coordinates.

Collective

Input Parameters:
===
- `dm`   - the `DMSTAG` object
- `xmin` - minimum global coordinate value in the x direction
- `xmax` - maximum global coordinate value in the x direction
- `ymin` - minimum global coordinate value in the y direction
- `ymax` - maximum global coordinate value in the y direction
- `zmin` - minimum global coordinate value in the z direction
- `zmax` - maximum global coordinate value in the z direction

Level: intermediate

Notes:
Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.

The per-dimension 1-dimensional `DMSTAG` objects that comprise the product
always have active 0-cells (vertices, element boundaries) and 1-cells
(element centers).

See the manual page for `DMStagSetUniformCoordinates()` for information on how
coordinates for dummy cells outside the physical domain boundary are populated.

See also: 
=== 
`DMSTAG`, `DMPRODUCT`, `DMStagSetUniformCoordinates()`, `DMStagSetUniformCoordinatesExplicit()`, `DMStagSetCoordinateDMType()`

# External Links
$(_doc_external("DMStag/DMStagSetUniformCoordinatesProduct"))
"""
function DMStagSetUniformCoordinatesProduct(dm::AbstractDMStag{PetscLib},xmin,xmax,ymin,ymax,zmin,zmax) where {PetscLib}

	LibPETSc.DMStagSetUniformCoordinatesProduct(
		PetscLib,
		dm,
		xmin,
		xmax,
		ymin,
		ymax,
		zmin,
		zmax
	)

	return nothing
end
 
 
"""
	ix = DMStagStencilToIndexLocal(dm::AbstractDMStag{PetscLib},dim::Int,n::Int,pos::Vector{DMStagStencil{Int}})

Convert an array of `DMStagStenci`l objects to an array of indices into a local vector.

Not Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `dim` - the dimension of the `DMSTAG` object
- `n`   - the number of `DMStagStencil` objects
- `pos` - an array of `n` `DMStagStencil` objects

Output Parameter:
===
- `ix` - output array of `n` indices

Notes:
The `DMStagStencil` objects in `pos` use global element indices.

The `.c` fields in `pos` must always be set (even if to `0`).

Developer Notes:
This is a "hot" function, and accepts the dimension redundantly to avoid having to perform any error checking inside the function.

Level: developer

See also: 
=== 
`DMSTAG`, `DMStagStencilLocation`, `DMStagStencil`, `DMGetLocalVector`, `DMCreateLocalVector`

# External Links
$(_doc_external("DMStag/DMStagStencilToIndexLocal"))
"""
function DMStagStencilToIndexLocal(dm::AbstractDMStag{PetscLib},dim::Int,n::Int,pos::Vector{DMStagStencil{Int}}) where {PetscLib}
	@assert length(pos) == n 
	
	PetscInt = PetscLib.PetscInt
	ix = zeros(PetscInt,n)

	LibPETSc.DMStagStencilToIndexLocal(
		PetscLib,
		dm,
		dim,
		n,
		pos,
		Ref(ix,n),
	)

	return ix
end
 
 
"""
	array, p_array = DMStagVecGetArray(dm::AbstractDMStag{PetscLib},vec::AbstractVector)

get access to local array

Logically Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object

Output Parameter:
===
- `array` - the array
- `p_array` - pointer to the array

Level: beginner

Note:
This function returns a (dim+1)-dimensional array for a dim-dimensional
`DMSTAG`.

The first 1-3 dimensions indicate an element in the global
numbering, using the standard C ordering.

The final dimension in this array corresponds to a degree
of freedom with respect to this element, for example corresponding to
the element or one of its neighboring faces, edges, or vertices.

For example, for a 3D `DMSTAG`, indexing is `array[k][j][i][idx]`, where `k` is the
index in the z-direction, `j` is the index in the y-direction, and `i` is the
index in the x-direction.

`idx` is obtained with `DMStagGetLocationSlot()`, since the correct offset
into the (d+1)-dimensional C array for a d-dimensional `DMSTAG` depends on the grid size and the number
of DOF stored at each location.

`DMStagVecRestoreArray()` must be called, once finished with the array

See also: 
=== 
`DMSTAG`, `DMStagVecGetArrayRead()`, `DMStagGetLocationSlot()`, `DMGetLocalVector()`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMCreateGlobalVector()`, `DMDAVecGetArray()`, `DMDAVecGetArrayDOF()`

# External Links
$(_doc_external("DMStag/DMStagVecGetArray"))
"""
function DMStagVecGetArray(dm::AbstractDMStag{PetscLib},vec::AbstractVector) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (dims...,dmE)
	r_array = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMStagVecGetArray(
		PetscLib,
		dm,
		vec,
		r_array,
	)

	array = PETSc_unsafe_wrap(r_array, dims; own=false)

	return array, r_array
end
 
 
"""
	array, p_array = DMStagVecGetArrayRead(dm::AbstractDMStag{PetscLib},vec::AbstractVector)

get read

Logically Collective

See the man page for `DMStagVecGetArray()` for more information.

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object

Output Parameter:
===
- `array` - the read-only array
- `p_array` - pointer to the array

Level: beginner

Note:
`DMStagVecRestoreArrayRead()` must be called, once finished with the array

See also: 
=== 
`DMSTAG`, `DMStagGetLocationSlot()`, `DMGetLocalVector()`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMCreateGlobalVector()`, `DMDAVecGetArrayRead()`, `DMDAVecGetArrayDOFRead()`

# External Links
$(_doc_external("DMStag/DMStagVecGetArrayRead"))
"""
function DMStagVecGetArrayRead(dm::AbstractDMStag{PetscLib},vec::AbstractVector) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (dims...,dmE)
	r_array = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMStagVecGetArrayRead(
		PetscLib,
		dm,
		vec,
		r_array,
	)

	array = PETSc_unsafe_wrap(r_array, dims; own=false)

	return array, r_array
end
 
 
"""
	val = DMStagVecGetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,n::Int,pos::Vector{DMStagStencil{Int}})

get vector values using grid indexing

Not Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - the vector object
- `n`   - the number of values to obtain
- `pos` - locations to obtain values from (as an array of `DMStagStencil` values)

Output Parameter:
===
- `val` - value at the point

Notes:
Accepts stencils which refer to global element numbers, but
only allows access to entries in the local representation (including ghosts).

This approach is not as efficient as getting values directly with `DMStagVecGetArray()`,
which is recommended for matrix-free operators.

Level: advanced

See also: 
=== 
`DMSTAG`, `DMStagStencil`, `DMStagStencilLocation`, `DMStagVecSetValuesStencil()`, `DMStagMatSetValuesStencil()`, `DMStagVecGetArray()`

# External Links
$(_doc_external("DMStag/DMStagVecGetValuesStencil"))
"""
function DMStagVecGetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,n::Int,pos::Vector{DMStagStencil{Int}}) where {PetscLib}
	@assert length(pos) == n 
	PetscScalar = PetscLib.PetscScalar
	val = zeros(PetscScalar,n)

	LibPETSc.DMStagVecGetValuesStencil(
		PetscLib,
		dm,
		vec,
		n,
		pos,
		Ref(val,n),
	)

	return val
end
 
DMStagVecGetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,pos::Vector{DMStagStencil{Int}}) where {PetscLib} = DMStagVecGetValuesStencil(dm,vec,length(pos),pos)


"""
	array = DMStagVecRestoreArray(dm::AbstractDMStag{PetscLib},vec::AbstractVector, p_array)

restore access to a raw array

Logically Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object
- `p_array` - pointer to the array

Output Parameter:
===
- `array` - the array

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagVecGetArray()`, `DMDAVecRestoreArray()`, `DMDAVecRestoreArrayDOF()`

# External Links
$(_doc_external("DMStag/DMStagVecRestoreArray"))
"""
function DMStagVecRestoreArray(dm::AbstractDMStag{PetscLib},vec::AbstractVector, r_array) where {PetscLib}
	
	LibPETSc.DMStagVecRestoreArray(
		PetscLib,
		dm,
		vec,
		r_array
	)

	return nothing
end
 
 
"""
	array = DMStagVecRestoreArrayRead(dm::AbstractDMStag{PetscLib},vec::AbstractVector, p_array)

restore read

Logically Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - the Vec object
- `p_array` - Pointer to array

Output Parameter:
===
- `array` - the read-only array

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagVecGetArrayRead()`, `DMDAVecRestoreArrayRead()`, `DMDAVecRestoreArrayDOFRead()`

# External Links
$(_doc_external("DMStag/DMStagVecRestoreArrayRead"))
"""
function DMStagVecRestoreArrayRead(dm::AbstractDMStag{PetscLib},vec::AbstractVector) where {PetscLib}

	LibPETSc.DMStagVecRestoreArrayRead(
		PetscLib,
		dm,
		vec,
		r_array,
	)

	return nothing
end
 
 
"""
	 DMStagVecSetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,n::Int,pos::Vector{DMStagStencil{Int}},val::Vector{<:AbstractFloat},insertMode::InsertMode)

Set `Vec` values using global grid indexing

Not Collective

Input Parameters:
===
- `dm`         - the `DMSTAG` object
- `vec`        - the `Vec`
- `n`          - the number of values to set
- `pos`        - the locations to set values, as an array of `DMStagStencil` structs
- `val`        - the values to set
- `insertMode` - `INSERT_VALUES` or `ADD_VALUES`

Notes:
The vector is expected to be a global vector compatible with the DM (usually obtained by `DMGetGlobalVector()` or `DMCreateGlobalVector()`).

This approach is not as efficient as setting values directly with `DMStagVecGetArray()`, which is recommended for matrix-free operators.
For assembling systems, where overhead may be less important than convenience, this routine could be helpful in assembling a righthand side and a matrix (using `DMStagMatSetValuesStencil()`).

Level: advanced

See also: 
=== 
`DMSTAG`, `Vec`, `DMStagStencil`, `DMStagStencilLocation`, `DMStagVecGetValuesStencil()`, `DMStagMatSetValuesStencil()`, `DMCreateGlobalVector()`, `DMGetLocalVector()`, `DMStagVecGetArray()`

# External Links
$(_doc_external("DMStag/DMStagVecSetValuesStencil"))
"""
function DMStagVecSetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,n::Int,pos::Vector{DMStagStencil{Int}},val::Vector{<:AbstractFloat},insertMode::InsertMode) where {PetscLib}
	@assert length(val) == length(pos) == n 
	
	LibPETSc.DMStagVecSetValuesStencil(
		PetscLib,
		dm,
		vec,
		n,
		pos,
		val,
		insertMode,
	)

	return nothing
end

"""
Same but no need to specify the length of input vectors
"""
DMStagVecSetValuesStencil(dm::AbstractDMStag{PetscLib},vec::AbstractVector,pos::Vector{DMStagStencil{Int}},val::Vector{<:AbstractFloat},insertMode::InsertMode) where {PetscLib} = DMStagVecSetValuesStencil(dm,vec,length(pos),pos,val,insertMode)


 
"""
	pda, pdavec = DMStagVecSplitToDMDA(dm::AbstractDMStag{PetscLib},vec::AbstractVector,loc::DMStagStencilLocation,c::Int)

create a `DMDA` and `Vec` from a subgrid of a `DMSTAG` and its `Vec`

Collective

Input Parameters:
===
- `dm`  - the `DMSTAG` object
- `vec` - `Vec` object associated with `dm`
- `loc` - which subgrid to extract (see `DMStagStencilLocation`)
- `c`   - which component to extract (see note below)

Output Parameters:
===
- `pda`    - the `DMDA`
- `pdavec` - the new `Vec`

Level: advanced

Notes:
If a `c` value of `-k` is provided, the first `k` DOF for that position are extracted,
padding with zero values if needed. If a non-negative value is provided, a single
DOF is extracted.

The caller is responsible for destroying the created `DMDA` and `Vec`.

See also: 
=== 
`DMSTAG`, `DMDA`, `DMStagStencilLocation`, `DM`, `Vec`, `DMStagMigrateVec()`, `DMStagCreateCompatibleDMStag()`

# External Links
$(_doc_external("DMStag/DMStagVecSplitToDMDA"))
"""
function DMStagVecSplitToDMDA(dm::AbstractDMStag{PetscLib},vec::AbstractVector,loc::DMStagStencilLocation,c::Int) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	pda 	= DMDA{PetscLib}(C_NULL, opts, petsclib.age)
	pdavec 	= VecSeq{PetscLib, petsclib.PetscScalar}(C_NULL, petsclib.age)

	LibPETSc.DMStagVecSplitToDMDA(
		PetscLib,
		dm,
		vec,
		loc,
		c,
		pda,
		pdavec,
	)

	return pda,pdavec
end
 
 
