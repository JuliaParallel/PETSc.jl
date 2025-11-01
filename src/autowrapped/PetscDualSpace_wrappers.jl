"""
	PetscDualSpaceRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscDualSpaceType`

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: `PetscDualSpace`, `PetscDualSpaceType`, `PetscDualSpaceRegisterAll()`, `PetscDualSpaceRegisterDestroy()`

# External Links
$(_doc_external("Dm/PetscDualSpaceRegister"))
"""
function PetscDualSpaceRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscDualSpaceRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscDualSpaceRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscDualSpaceSetType(petsclib::PetscLibType,sp::PetscDualSpace, name::PetscDualSpaceType) 
Builds a particular `PetscDualSpace` based on its `PetscDualSpaceType`

Collective

Input Parameters:
- `sp`   - The `PetscDualSpace` object
- `name` - The kind of space

Options Database Key:
- `-petscdualspace_type <type>` - Sets the PetscDualSpace type; use -help for a list of available types

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceType`, `PetscDualSpaceGetType()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetType"))
"""
function PetscDualSpaceSetType(petsclib::PetscLibType, sp::PetscDualSpace, name::PetscDualSpaceType) end

@for_petsc function PetscDualSpaceSetType(petsclib::$UnionPetscLib, sp::PetscDualSpace, name::PetscDualSpaceType )

    @chk ccall(
               (:PetscDualSpaceSetType, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDualSpaceType),
               sp, name,
              )


	return nothing
end 

"""
	name::PetscDualSpaceType = PetscDualSpaceGetType(petsclib::PetscLibType,sp::PetscDualSpace) 
Gets the `PetscDualSpaceType` name (as a string) from the object.

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `name` - The `PetscDualSpaceType` name

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceType`, `PetscDualSpaceSetType()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetType"))
"""
function PetscDualSpaceGetType(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetType(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	name_ = Ref{PetscDualSpaceType}()

    @chk ccall(
               (:PetscDualSpaceGetType, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscDualSpaceType}),
               sp, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscDualSpaceViewFromOptions(petsclib::PetscLibType,A::PetscDualSpace, obj::PetscObject, name::String) 
View a `PetscDualSpace` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscDualSpace` object
- `obj`  - Optional object, provides the options prefix
- `name` - command line option name

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceView()`, `PetscObjectViewFromOptions()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceViewFromOptions"))
"""
function PetscDualSpaceViewFromOptions(petsclib::PetscLibType, A::PetscDualSpace, obj::PetscObject, name::String) end

@for_petsc function PetscDualSpaceViewFromOptions(petsclib::$UnionPetscLib, A::PetscDualSpace, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscDualSpaceViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscDualSpaceView(petsclib::PetscLibType,sp::PetscDualSpace, v::PetscViewer) 
Views a `PetscDualSpace`

Collective

Input Parameters:
- `sp` - the `PetscDualSpace` object to view
- `v`  - the viewer

Level: beginner

-seealso: `PetscViewer`, `PetscDualSpaceDestroy()`, `PetscDualSpace`

# External Links
$(_doc_external("Dm/PetscDualSpaceView"))
"""
function PetscDualSpaceView(petsclib::PetscLibType, sp::PetscDualSpace, v::PetscViewer) end

@for_petsc function PetscDualSpaceView(petsclib::$UnionPetscLib, sp::PetscDualSpace, v::PetscViewer )

    @chk ccall(
               (:PetscDualSpaceView, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscViewer),
               sp, v,
              )


	return nothing
end 

"""
	PetscDualSpaceSetFromOptions(petsclib::PetscLibType,sp::PetscDualSpace) 
sets parameters in a `PetscDualSpace` from the options database

Collective

Input Parameter:
- `sp` - the `PetscDualSpace` object to set options for

Options Database Keys:
- `-petscdualspace_order <order>`                 - the approximation order of the space
- `-petscdualspace_form_degree <deg>`             - the form degree, say 0 for point evaluations, or 2 for area integrals
- `-petscdualspace_components <c>`                - the number of components, say d for a vector field
- `-petscdualspace_refcell <celltype>`            - Reference cell type name
- `-petscdualspace_lagrange_continuity`           - Flag for continuous element
- `-petscdualspace_lagrange_tensor`               - Flag for tensor dual space
- `-petscdualspace_lagrange_trimmed`              - Flag for trimmed dual space
- `-petscdualspace_lagrange_node_type <nodetype>` - Lagrange node location type
- `-petscdualspace_lagrange_node_endpoints`       - Flag for nodes that include endpoints
- `-petscdualspace_lagrange_node_exponent`        - Gauss-Jacobi weight function exponent
- `-petscdualspace_lagrange_use_moments`          - Use moments (where appropriate) for functionals
- `-petscdualspace_lagrange_moment_order <order>` - Quadrature order for moment functionals

Level: intermediate

-seealso: `PetscDualSpaceView()`, `PetscDualSpace`, `PetscObjectSetFromOptions()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetFromOptions"))
"""
function PetscDualSpaceSetFromOptions(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSetFromOptions(petsclib::$UnionPetscLib, sp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace,),
               sp,
              )


	return nothing
end 

"""
	PetscDualSpaceSetUp(petsclib::PetscLibType,sp::PetscDualSpace) 
Construct a basis for a `PetscDualSpace`

Collective

Input Parameter:
- `sp` - the `PetscDualSpace` object to setup

Level: intermediate

-seealso: `PetscDualSpaceView()`, `PetscDualSpaceDestroy()`, `PetscDualSpace`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetUp"))
"""
function PetscDualSpaceSetUp(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSetUp(petsclib::$UnionPetscLib, sp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceSetUp, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace,),
               sp,
              )


	return nothing
end 

"""
	PetscDualSpaceDestroy(petsclib::PetscLibType,sp::PetscDualSpace) 
Destroys a `PetscDualSpace` object

Collective

Input Parameter:
- `sp` - the `PetscDualSpace` object to destroy

Level: beginner

-seealso: `PetscDualSpace`, `PetscDualSpaceView()`, `PetscDualSpace()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceDestroy"))
"""
function PetscDualSpaceDestroy(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceDestroy(petsclib::$UnionPetscLib, sp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDualSpace},),
               sp,
              )


	return nothing
end 

"""
	sp::PetscDualSpace = PetscDualSpaceCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscDualSpace` object. The type can then be set with `PetscDualSpaceSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscDualSpace` object

Output Parameter:
- `sp` - The `PetscDualSpace` object

Level: beginner

-seealso: `PetscDualSpace`, `PetscDualSpaceSetType()`, `PETSCDUALSPACELAGRANGE`

# External Links
$(_doc_external("Dm/PetscDualSpaceCreate"))
"""
function PetscDualSpaceCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscDualSpaceCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	sp_ = Ref{PetscDualSpace}()

    @chk ccall(
               (:PetscDualSpaceCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDualSpace}),
               comm, sp_,
              )

	sp = sp_[]

	return sp
end 

"""
	spNew::PetscDualSpace = PetscDualSpaceDuplicate(petsclib::PetscLibType,sp::PetscDualSpace) 
Creates a duplicate `PetscDualSpace` object that is not setup.

Collective

Input Parameter:
- `sp` - The original `PetscDualSpace`

Output Parameter:
- `spNew` - The duplicate `PetscDualSpace`

Level: beginner

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`, `PetscDualSpaceSetType()`

# External Links
$(_doc_external("Dm/PetscDualSpaceDuplicate"))
"""
function PetscDualSpaceDuplicate(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceDuplicate(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	spNew_ = Ref{PetscDualSpace}()

    @chk ccall(
               (:PetscDualSpaceDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscDualSpace}),
               sp, spNew_,
              )

	spNew = spNew_[]

	return spNew
end 

"""
	PetscDualSpaceGetDM(petsclib::PetscLibType,sp::PetscDualSpace, dm::PetscDM) 
Get the `DM` representing the reference cell of a `PetscDualSpace`

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `dm` - The reference cell, that is a `DM` that consists of a single cell

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceSetDM()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetDM"))
"""
function PetscDualSpaceGetDM(petsclib::PetscLibType, sp::PetscDualSpace, dm::PetscDM) end

@for_petsc function PetscDualSpaceGetDM(petsclib::$UnionPetscLib, sp::PetscDualSpace, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:PetscDualSpaceGetDM, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{CDM}),
               sp, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	PetscDualSpaceSetDM(petsclib::PetscLibType,sp::PetscDualSpace, dm::PetscDM) 
Get the `DM` representing the reference cell

Not Collective

Input Parameters:
- `sp` - The `PetscDual`Space
- `dm` - The reference cell

Level: intermediate

-seealso: `PetscDualSpace`, `DM`, `PetscDualSpaceGetDM()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetDM"))
"""
function PetscDualSpaceSetDM(petsclib::PetscLibType, sp::PetscDualSpace, dm::PetscDM) end

@for_petsc function PetscDualSpaceSetDM(petsclib::$UnionPetscLib, sp::PetscDualSpace, dm::PetscDM )

    @chk ccall(
               (:PetscDualSpaceSetDM, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, CDM),
               sp, dm,
              )


	return nothing
end 

"""
	order::PetscInt = PetscDualSpaceGetOrder(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the order of the dual space

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `order` - The order

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceSetOrder()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetOrder"))
"""
function PetscDualSpaceGetOrder(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetOrder(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetOrder, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, order_,
              )

	order = order_[]

	return order
end 

"""
	PetscDualSpaceSetOrder(petsclib::PetscLibType,sp::PetscDualSpace, order::PetscInt) 
Set the order of the dual space

Not Collective

Input Parameters:
- `sp`    - The `PetscDualSpace`
- `order` - The order

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceGetOrder()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetOrder"))
"""
function PetscDualSpaceSetOrder(petsclib::PetscLibType, sp::PetscDualSpace, order::PetscInt) end

@for_petsc function PetscDualSpaceSetOrder(petsclib::$UnionPetscLib, sp::PetscDualSpace, order::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceSetOrder, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               sp, order,
              )


	return nothing
end 

"""
	Nc::PetscInt = PetscDualSpaceGetNumComponents(petsclib::PetscLibType,sp::PetscDualSpace) 
Return the number of components for this space

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscDualSpaceSetNumComponents()`, `PetscDualSpaceGetDimension()`, `PetscDualSpaceCreate()`, `PetscDualSpace`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetNumComponents"))
"""
function PetscDualSpaceGetNumComponents(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetNumComponents(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	Nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, Nc_,
              )

	Nc = Nc_[]

	return Nc
end 

"""
	PetscDualSpaceSetNumComponents(petsclib::PetscLibType,sp::PetscDualSpace, Nc::PetscInt) 
Set the number of components for this space

Input Parameters:
- `sp` - The `PetscDualSpace`
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscDualSpaceGetNumComponents()`, `PetscDualSpaceCreate()`, `PetscDualSpace`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetNumComponents"))
"""
function PetscDualSpaceSetNumComponents(petsclib::PetscLibType, sp::PetscDualSpace, Nc::PetscInt) end

@for_petsc function PetscDualSpaceSetNumComponents(petsclib::$UnionPetscLib, sp::PetscDualSpace, Nc::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceSetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               sp, Nc,
              )


	return nothing
end 

"""
	PetscDualSpaceGetFunctional(petsclib::PetscLibType,sp::PetscDualSpace, i::PetscInt, fncal::PetscQuadrature) 
Get the i

Not Collective

Input Parameters:
- `sp` - The `PetscDualSpace`
- `i`  - The basis number

Output Parameter:
- `functional` - The basis functional

Level: intermediate

-seealso: `PetscDualSpace`, `PetscQuadrature`, `PetscDualSpaceGetDimension()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetFunctional"))
"""
function PetscDualSpaceGetFunctional(petsclib::PetscLibType, sp::PetscDualSpace, i::PetscInt, fncal::PetscQuadrature) end

@for_petsc function PetscDualSpaceGetFunctional(petsclib::$UnionPetscLib, sp::PetscDualSpace, i::$PetscInt, fncal::PetscQuadrature )

    @chk ccall(
               (:PetscDualSpaceGetFunctional, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, Ptr{PetscQuadrature}),
               sp, i, fncal,
              )


	return nothing
end 

"""
	dim::PetscInt = PetscDualSpaceGetDimension(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the dimension of the dual space, i.e. the number of basis functionals

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `dim` - The dimension

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceGetFunctional()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetDimension"))
"""
function PetscDualSpaceGetDimension(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetDimension(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	intdim::PetscInt = PetscDualSpaceGetInteriorDimension(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the interior dimension of the dual space, i.e. the number of basis functionals assigned to the interior of the reference domain

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `intdim` - The dimension

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceGetFunctional()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetInteriorDimension"))
"""
function PetscDualSpaceGetInteriorDimension(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetInteriorDimension(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	intdim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetInteriorDimension, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, intdim_,
              )

	intdim = intdim_[]

	return intdim
end 

"""
	uniform::PetscBool = PetscDualSpaceGetUniform(petsclib::PetscLibType,sp::PetscDualSpace) 
Whether this dual space is uniform

Not Collective

Input Parameter:
- `sp` - A dual space

Output Parameter:
- `uniform` - `PETSC_TRUE` if (a) the dual space is the same for each point in a stratum of the reference `DMPLEX`, and
(b) every symmetry of each point in the reference `DMPLEX` is also a symmetry of the point's dual space.

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceGetPointSubspace()`, `PetscDualSpaceGetSymmetries()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetUniform"))
"""
function PetscDualSpaceGetUniform(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetUniform(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	uniform_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceGetUniform, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, uniform_,
              )

	uniform = uniform_[]

	return uniform
end 

"""
	numDof::Vector{PetscInt} = PetscDualSpaceGetNumDof(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the number of degrees of freedom for each spatial (topological) dimension

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `numDof` - An array of length dim+1 which holds the number of dofs for each dimension

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceGetFunctional()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetNumDof"))
"""
function PetscDualSpaceGetNumDof(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetNumDof(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	numDof_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDualSpaceGetNumDof, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{Ptr{$PetscInt}}),
               sp, numDof_,
              )

	numDof = unsafe_wrap(Array, numDof_[], VecGetLocalSize(petsclib, x); own = false)

	return numDof
end 

"""
	PetscDualSpaceGetSection(petsclib::PetscLibType,sp::PetscDualSpace, section::PetscSection) 
Create a `PetscSection` over the reference cell with the layout from this space

Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `section` - The section

Level: advanced

-seealso: `PetscDualSpace`, `PetscSection`, `PetscDualSpaceCreate()`, `DMPLEX`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetSection"))
"""
function PetscDualSpaceGetSection(petsclib::PetscLibType, sp::PetscDualSpace, section::PetscSection) end

@for_petsc function PetscDualSpaceGetSection(petsclib::$UnionPetscLib, sp::PetscDualSpace, section::PetscSection )

    @chk ccall(
               (:PetscDualSpaceGetSection, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscSection}),
               sp, section,
              )


	return nothing
end 

"""
	PetscDualSpaceGetInteriorSection(petsclib::PetscLibType,sp::PetscDualSpace, section::PetscSection) 
Create a `PetscSection` over the reference cell with the layout from this space
for interior degrees of freedom

Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `section` - The interior section

Level: advanced

-seealso: `PetscDualSpace`, `PetscSection`, `PetscDualSpaceCreate()`, `DMPLEX`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetInteriorSection"))
"""
function PetscDualSpaceGetInteriorSection(petsclib::PetscLibType, sp::PetscDualSpace, section::PetscSection) end

@for_petsc function PetscDualSpaceGetInteriorSection(petsclib::$UnionPetscLib, sp::PetscDualSpace, section::PetscSection )

    @chk ccall(
               (:PetscDualSpaceGetInteriorSection, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscSection}),
               sp, section,
              )


	return nothing
end 

"""
	value::PetscScalar = PetscDualSpaceApply(petsclib::PetscLibType,sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFEGeom, numComp::PetscInt, func::external, ctx::Cvoid) 
Apply a functional from the dual space basis to an input function

Input Parameters:
- `sp`      - The `PetscDualSpace` object
- `f`       - The basis functional index
- `time`    - The time
- `cgeom`   - A context with geometric information for this cell, we use v0 (the initial vertex) and J (the Jacobian) (or evaluated at the coordinates of the functional)
- `numComp` - The number of components for the function
- `func`    - The input function
- `ctx`     - A context for the function

Output Parameter:
- `value` - numComp output values

Calling sequence:
-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApply"))
"""
function PetscDualSpaceApply(petsclib::PetscLibType, sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFEGeom, numComp::PetscInt, func::external, ctx::Cvoid) end

@for_petsc function PetscDualSpaceApply(petsclib::$UnionPetscLib, sp::PetscDualSpace, f::$PetscInt, time::$PetscReal, cgeom::PetscFEGeom, numComp::$PetscInt, func::external, ctx::Cvoid )
	value_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApply, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, $PetscReal, Ptr{PetscFEGeom}, $PetscInt, external, Ptr{Cvoid}, Ptr{$PetscScalar}),
               sp, f, time, cgeom, numComp, func, ctx, value_,
              )

	value = value_[]

	return value
end 

"""
	spValue::PetscScalar = PetscDualSpaceApplyAll(petsclib::PetscLibType,sp::PetscDualSpace, pointEval::PetscScalar) 
Apply all functionals from the dual space basis to the result of an evaluation at the points returned by `PetscDualSpaceGetAllData()`

Input Parameters:
- `sp`        - The `PetscDualSpace` object
- `pointEval` - Evaluation at the points returned by `PetscDualSpaceGetAllData()`

Output Parameter:
- `spValue` - The values of all dual space functionals

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyAll"))
"""
function PetscDualSpaceApplyAll(petsclib::PetscLibType, sp::PetscDualSpace, pointEval::PetscScalar) end

@for_petsc function PetscDualSpaceApplyAll(petsclib::$UnionPetscLib, sp::PetscDualSpace, pointEval::$PetscScalar )
	spValue_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyAll, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               sp, pointEval, spValue_,
              )

	spValue = spValue_[]

	return spValue
end 

"""
	spValue::PetscScalar = PetscDualSpaceApplyInterior(petsclib::PetscLibType,sp::PetscDualSpace, pointEval::PetscScalar) 
Apply interior functionals from the dual space basis to the result of an evaluation at the points returned by `PetscDualSpaceGetInteriorData()`

Input Parameters:
- `sp`        - The `PetscDualSpace` object
- `pointEval` - Evaluation at the points returned by `PetscDualSpaceGetInteriorData()`

Output Parameter:
- `spValue` - The values of interior dual space functionals

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyInterior"))
"""
function PetscDualSpaceApplyInterior(petsclib::PetscLibType, sp::PetscDualSpace, pointEval::PetscScalar) end

@for_petsc function PetscDualSpaceApplyInterior(petsclib::$UnionPetscLib, sp::PetscDualSpace, pointEval::$PetscScalar )
	spValue_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyInterior, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               sp, pointEval, spValue_,
              )

	spValue = spValue_[]

	return spValue
end 

"""
	value::PetscScalar = PetscDualSpaceApplyDefault(petsclib::PetscLibType,sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFEGeom, Nc::PetscInt, func::external, ctx::Cvoid) 
Apply a functional from the dual space basis to an input function by assuming a point evaluation functional.

Input Parameters:
- `sp`    - The `PetscDualSpace` object
- `f`     - The basis functional index
- `time`  - The time
- `cgeom` - A context with geometric information for this cell, we use v0 (the initial vertex) and J (the Jacobian)
- `Nc`    - The number of components for the function
- `func`  - The input function
- `ctx`   - A context for the function

Output Parameter:
- `value` - The output value

Calling sequence:
-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyDefault"))
"""
function PetscDualSpaceApplyDefault(petsclib::PetscLibType, sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFEGeom, Nc::PetscInt, func::external, ctx::Cvoid) end

@for_petsc function PetscDualSpaceApplyDefault(petsclib::$UnionPetscLib, sp::PetscDualSpace, f::$PetscInt, time::$PetscReal, cgeom::PetscFEGeom, Nc::$PetscInt, func::external, ctx::Cvoid )
	value_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyDefault, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, $PetscReal, Ptr{PetscFEGeom}, $PetscInt, external, Ptr{Cvoid}, Ptr{$PetscScalar}),
               sp, f, time, cgeom, Nc, func, ctx, value_,
              )

	value = value_[]

	return value
end 

"""
	spValue::PetscScalar = PetscDualSpaceApplyAllDefault(petsclib::PetscLibType,sp::PetscDualSpace, pointEval::PetscScalar) 
Apply all functionals from the dual space basis to the result of an evaluation at the points returned by `PetscDualSpaceGetAllData()`

Input Parameters:
- `sp`        - The `PetscDualSpace` object
- `pointEval` - Evaluation at the points returned by `PetscDualSpaceGetAllData()`

Output Parameter:
- `spValue` - The values of all dual space functionals

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyAllDefault"))
"""
function PetscDualSpaceApplyAllDefault(petsclib::PetscLibType, sp::PetscDualSpace, pointEval::PetscScalar) end

@for_petsc function PetscDualSpaceApplyAllDefault(petsclib::$UnionPetscLib, sp::PetscDualSpace, pointEval::$PetscScalar )
	spValue_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyAllDefault, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               sp, pointEval, spValue_,
              )

	spValue = spValue_[]

	return spValue
end 

"""
	spValue::PetscScalar = PetscDualSpaceApplyInteriorDefault(petsclib::PetscLibType,sp::PetscDualSpace, pointEval::PetscScalar) 
Apply interior functionals from the dual space basis to the result of an evaluation at the points returned by `PetscDualSpaceGetInteriorData()`

Input Parameters:
- `sp`        - The `PetscDualSpace` object
- `pointEval` - Evaluation at the points returned by `PetscDualSpaceGetInteriorData()`

Output Parameter:
- `spValue` - The values of interior dual space functionals

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyInteriorDefault"))
"""
function PetscDualSpaceApplyInteriorDefault(petsclib::PetscLibType, sp::PetscDualSpace, pointEval::PetscScalar) end

@for_petsc function PetscDualSpaceApplyInteriorDefault(petsclib::$UnionPetscLib, sp::PetscDualSpace, pointEval::$PetscScalar )
	spValue_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyInteriorDefault, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               sp, pointEval, spValue_,
              )

	spValue = spValue_[]

	return spValue
end 

"""
	PetscDualSpaceGetAllData(petsclib::PetscLibType,sp::PetscDualSpace, allNodes::PetscQuadrature, allMat::PetscMat) 
Get all quadrature nodes from this space, and the matrix that sends quadrature node values to degree

Input Parameter:
- `sp` - The dualspace

Output Parameters:
- `allNodes` - A `PetscQuadrature` object containing all evaluation nodes, pass `NULL` if not needed
- `allMat`   - A `Mat` for the node-to-dof transformation, pass `NULL` if not needed

Level: advanced

-seealso: `PetscQuadrature`, `PetscDualSpace`, `PetscDualSpaceCreate()`, `Mat`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetAllData"))
"""
function PetscDualSpaceGetAllData(petsclib::PetscLibType, sp::PetscDualSpace, allNodes::PetscQuadrature, allMat::PetscMat) end

@for_petsc function PetscDualSpaceGetAllData(petsclib::$UnionPetscLib, sp::PetscDualSpace, allNodes::PetscQuadrature, allMat::PetscMat )
	allMat_ = Ref(allMat.ptr)

    @chk ccall(
               (:PetscDualSpaceGetAllData, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{CMat}),
               sp, allNodes, allMat_,
              )

	allMat.ptr = C_NULL

	return nothing
end 

"""
	allNodes::PetscQuadrature,allMat::PetscMat = PetscDualSpaceCreateAllDataDefault(petsclib::PetscLibType,sp::PetscDualSpace) 
Create all evaluation nodes and the node

Input Parameter:
- `sp` - The dualspace

Output Parameters:
- `allNodes` - A `PetscQuadrature` object containing all evaluation nodes
- `allMat`   - A `Mat` for the node-to-dof transformation

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`, `Mat`, `PetscQuadrature`

# External Links
$(_doc_external("Dm/PetscDualSpaceCreateAllDataDefault"))
"""
function PetscDualSpaceCreateAllDataDefault(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceCreateAllDataDefault(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	allNodes_ = Ref{PetscQuadrature}()
	allMat_ = Ref{CMat}()

    @chk ccall(
               (:PetscDualSpaceCreateAllDataDefault, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{CMat}),
               sp, allNodes_, allMat_,
              )

	allNodes = allNodes_[]
	allMat = PetscMat(allMat_[], petsclib)

	return allNodes,allMat
end 

"""
	PetscDualSpaceGetInteriorData(petsclib::PetscLibType,sp::PetscDualSpace, intNodes::PetscQuadrature, intMat::PetscMat) 
Get all quadrature points necessary to compute the interior degrees of freedom from
this space, as well as the matrix that computes the degrees of freedom from the quadrature
values.

Input Parameter:
- `sp` - The dualspace

Output Parameters:
- `intNodes` - A `PetscQuadrature` object containing all evaluation points needed to evaluate interior degrees of freedom,
pass `NULL` if not needed
- `intMat`   - A matrix that computes dual space values from point values: size [spdim0 x (npoints * nc)], where spdim0 is
the size of the constrained layout (`PetscSectionGetConstrainStorageSize()`) of the dual space section,
npoints is the number of points in intNodes and nc is `PetscDualSpaceGetNumComponents()`.
Pass `NULL` if not needed

Level: advanced

-seealso: `PetscDualSpace`, `PetscQuadrature`, `Mat`, `PetscDualSpaceCreate()`, `PetscDualSpaceGetDimension()`, `PetscDualSpaceGetNumComponents()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetInteriorData"))
"""
function PetscDualSpaceGetInteriorData(petsclib::PetscLibType, sp::PetscDualSpace, intNodes::PetscQuadrature, intMat::PetscMat) end

@for_petsc function PetscDualSpaceGetInteriorData(petsclib::$UnionPetscLib, sp::PetscDualSpace, intNodes::PetscQuadrature, intMat::PetscMat )
	intMat_ = Ref(intMat.ptr)

    @chk ccall(
               (:PetscDualSpaceGetInteriorData, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{CMat}),
               sp, intNodes, intMat_,
              )

	intMat.ptr = C_NULL

	return nothing
end 

"""
	intNodes::PetscQuadrature,intMat::PetscMat = PetscDualSpaceCreateInteriorDataDefault(petsclib::PetscLibType,sp::PetscDualSpace) 
Create quadrature points by examining interior functionals and create the matrix mapping quadrature point values to interior dual space values

Input Parameter:
- `sp` - The dualspace

Output Parameters:
- `intNodes` - A `PetscQuadrature` object containing all evaluation points needed to evaluate interior degrees of freedom
- `intMat`   - A matrix that computes dual space values from point values: size [spdim0 x (npoints * nc)], where spdim0 is
the size of the constrained layout (`PetscSectionGetConstrainStorageSize()`) of the dual space section,
npoints is the number of points in allNodes and nc is `PetscDualSpaceGetNumComponents()`.

Level: advanced

-seealso: `PetscDualSpace`, `PetscQuadrature`, `Mat`, `PetscDualSpaceCreate()`, `PetscDualSpaceGetInteriorData()`

# External Links
$(_doc_external("Dm/PetscDualSpaceCreateInteriorDataDefault"))
"""
function PetscDualSpaceCreateInteriorDataDefault(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceCreateInteriorDataDefault(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	intNodes_ = Ref{PetscQuadrature}()
	intMat_ = Ref{CMat}()

    @chk ccall(
               (:PetscDualSpaceCreateInteriorDataDefault, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscQuadrature}, Ptr{CMat}),
               sp, intNodes_, intMat_,
              )

	intNodes = intNodes_[]
	intMat = PetscMat(intMat_[], petsclib)

	return intNodes,intMat
end 

"""
	equal::PetscBool = PetscDualSpaceEqual(petsclib::PetscLibType,A::PetscDualSpace, B::PetscDualSpace) 
Determine if two dual spaces are equivalent

Input Parameters:
- `A` - A `PetscDualSpace` object
- `B` - Another `PetscDualSpace` object

Output Parameter:
- `equal` - `PETSC_TRUE` if the dual spaces are equivalent

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceEqual"))
"""
function PetscDualSpaceEqual(petsclib::PetscLibType, A::PetscDualSpace, B::PetscDualSpace) end

@for_petsc function PetscDualSpaceEqual(petsclib::$UnionPetscLib, A::PetscDualSpace, B::PetscDualSpace )
	equal_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceEqual, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDualSpace, Ptr{PetscBool}),
               A, B, equal_,
              )

	equal = equal_[]

	return equal
end 

"""
	value::PetscScalar = PetscDualSpaceApplyFVM(petsclib::PetscLibType,sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFVCellGeom, Nc::PetscInt, func::external, ctx::Cvoid) 
Apply a functional from the dual space basis to an input function by assuming a point evaluation functional at the cell centroid.

Input Parameters:
- `sp`    - The `PetscDualSpace` object
- `f`     - The basis functional index
- `time`  - The time
- `cgeom` - A context with geometric information for this cell, we currently just use the centroid
- `Nc`    - The number of components for the function
- `func`  - The input function
- `ctx`   - A context for the function

Output Parameter:
- `value` - The output value (scalar)

Calling sequence:
-seealso: `PetscDualSpace`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceApplyFVM"))
"""
function PetscDualSpaceApplyFVM(petsclib::PetscLibType, sp::PetscDualSpace, f::PetscInt, time::PetscReal, cgeom::PetscFVCellGeom, Nc::PetscInt, func::external, ctx::Cvoid) end

@for_petsc function PetscDualSpaceApplyFVM(petsclib::$UnionPetscLib, sp::PetscDualSpace, f::$PetscInt, time::$PetscReal, cgeom::PetscFVCellGeom, Nc::$PetscInt, func::external, ctx::Cvoid )
	value_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceApplyFVM, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, $PetscReal, Ptr{PetscFVCellGeom}, $PetscInt, external, Ptr{Cvoid}, Ptr{$PetscScalar}),
               sp, f, time, cgeom, Nc, func, ctx, value_,
              )

	value = value_[]

	return value
end 

"""
	PetscDualSpaceGetHeightSubspace(petsclib::PetscLibType,sp::PetscDualSpace, height::PetscInt, subsp::PetscDualSpace) 
Get the subset of the dual space basis that is supported on a mesh point of a
given height.  This assumes that the reference cell is symmetric over points of this height.

Not Collective

Input Parameters:
- `sp`     - the `PetscDualSpace` object
- `height` - the height of the mesh point for which the subspace is desired

Output Parameter:
- `subsp` - the subspace.  Note that the functionals in the subspace are with respect to the intrinsic geometry of the
point, which will be of lesser dimension if height > 0.

Level: advanced

-seealso: `PetscDualSpace`, `PetscSpaceGetHeightSubspace()`, `PetscDualSpaceGetPointSubspace()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetHeightSubspace"))
"""
function PetscDualSpaceGetHeightSubspace(petsclib::PetscLibType, sp::PetscDualSpace, height::PetscInt, subsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetHeightSubspace(petsclib::$UnionPetscLib, sp::PetscDualSpace, height::$PetscInt, subsp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceGetHeightSubspace, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
               sp, height, subsp,
              )


	return nothing
end 

"""
	PetscDualSpaceGetPointSubspace(petsclib::PetscLibType,sp::PetscDualSpace, point::PetscInt, bdsp::PetscDualSpace) 
Get the subset of the dual space basis that is supported on a particular mesh point.

Not Collective

Input Parameters:
- `sp`    - the `PetscDualSpace` object
- `point` - the point (in the dual space's DM) for which the subspace is desired

Output Parameter:
- `bdsp` - the subspace.

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpaceGetHeightSubspace()`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetPointSubspace"))
"""
function PetscDualSpaceGetPointSubspace(petsclib::PetscLibType, sp::PetscDualSpace, point::PetscInt, bdsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetPointSubspace(petsclib::$UnionPetscLib, sp::PetscDualSpace, point::$PetscInt, bdsp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceGetPointSubspace, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
               sp, point, bdsp,
              )


	return nothing
end 

"""
	perms::PetscInt,flips::PetscScalar = PetscDualSpaceGetSymmetries(petsclib::PetscLibType,sp::PetscDualSpace) 
Returns a description of the symmetries of this basis

Not Collective

Input Parameter:
- `sp` - the `PetscDualSpace` object

Output Parameters:
- `perms` - Permutations of the interior degrees of freedom, parameterized by the point orientation
- `flips` - Sign reversal of the interior degrees of freedom, parameterized by the point orientation

Level: developer

-seealso: `PetscDualSpace`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetSymmetries"))
"""
function PetscDualSpaceGetSymmetries(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetSymmetries(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	perms_ = Ref{$PetscInt}()
	flips_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscDualSpaceGetSymmetries, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, $PetscScalar),
               sp, perms_, flips_,
              )

	perms = perms_[]
	flips = flips_[]

	return perms,flips
end 

"""
	k::PetscInt = PetscDualSpaceGetFormDegree(petsclib::PetscLibType,dsp::PetscDualSpace) 
Get the form degree k for the k
dual space's functionals.

Input Parameter:
- `dsp` - The `PetscDualSpace`

Output Parameter:
- `k` - The *signed* degree k of the k.  If k >= 0, this means that the degrees of freedom are k-forms, and are stored
in lexicographic order according to the basis of k-forms constructed from the wedge product of 1-forms.  So for example,
the 1-form basis in 3-D is (dx, dy, dz), and the 2-form basis in 3-D is (dx wedge dy, dx wedge dz, dy wedge dz).
If k < 0, this means that the degrees transform as k-forms, but are stored as (N-k) forms according to the
Hodge star map.  So for example if k = -2 and N = 3, this means that the degrees of freedom transform as 2-forms
but are stored as 1-forms.

Level: developer

-seealso: `PetscDualSpace`, `PetscDTAltV`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransform()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetFormDegree"))
"""
function PetscDualSpaceGetFormDegree(petsclib::PetscLibType, dsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetFormDegree(petsclib::$UnionPetscLib, dsp::PetscDualSpace )
	k_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetFormDegree, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               dsp, k_,
              )

	k = k_[]

	return k
end 

"""
	PetscDualSpaceSetFormDegree(petsclib::PetscLibType,dsp::PetscDualSpace, k::PetscInt) 
Set the form degree k for the k
dual space's functionals.

Input Parameters:
- `dsp` - The `PetscDualSpace`
- `k`   - The *signed* degree k of the k.  If k >= 0, this means that the degrees of freedom are k-forms, and are stored
in lexicographic order according to the basis of k-forms constructed from the wedge product of 1-forms.  So for example,
the 1-form basis in 3-D is (dx, dy, dz), and the 2-form basis in 3-D is (dx wedge dy, dx wedge dz, dy wedge dz).
If k < 0, this means that the degrees transform as k-forms, but are stored as (N-k) forms according to the
Hodge star map.  So for example if k = -2 and N = 3, this means that the degrees of freedom transform as 2-forms
but are stored as 1-forms.

Level: developer

-seealso: `PetscDualSpace`, `PetscDTAltV`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransform()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceSetFormDegree"))
"""
function PetscDualSpaceSetFormDegree(petsclib::PetscLibType, dsp::PetscDualSpace, k::PetscInt) end

@for_petsc function PetscDualSpaceSetFormDegree(petsclib::$UnionPetscLib, dsp::PetscDualSpace, k::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceSetFormDegree, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               dsp, k,
              )


	return nothing
end 

"""
	k::PetscInt = PetscDualSpaceGetDeRahm(petsclib::PetscLibType,dsp::PetscDualSpace) 
Get the k

Input Parameter:
- `dsp` - The `PetscDualSpace`

Output Parameter:
- `k` - The simplex dimension

Level: developer

-seealso: `PetscDualSpace`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransform()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceGetDeRahm"))
"""
function PetscDualSpaceGetDeRahm(petsclib::PetscLibType, dsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceGetDeRahm(petsclib::$UnionPetscLib, dsp::PetscDualSpace )
	k_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceGetDeRahm, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               dsp, k_,
              )

	k = k_[]

	return k
end 

"""
	PetscDualSpaceTransform(petsclib::PetscLibType,dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) 
Transform the function values

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `trans`     - The type of transform
- `isInverse` - Flag to invert the transform
- `fegeom`    - The cell geometry
- `Nv`        - The number of function samples
- `Nc`        - The number of function components
- `vals`      - The function values

Output Parameter:
- `vals` - The transformed function values

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceTransformGradient()`, `PetscDualSpaceTransformHessian()`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceTransform"))
"""
function PetscDualSpaceTransform(petsclib::PetscLibType, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscDualSpaceTransform(petsclib::$UnionPetscLib, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::$PetscInt, Nc::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpaceTransform, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDualSpaceTransformType, PetscBool, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, trans, isInverse, fegeom, Nv, Nc, vals,
              )


	return nothing
end 

"""
	PetscDualSpaceTransformGradient(petsclib::PetscLibType,dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) 
Transform the function gradient values

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `trans`     - The type of transform
- `isInverse` - Flag to invert the transform
- `fegeom`    - The cell geometry
- `Nv`        - The number of function gradient samples
- `Nc`        - The number of function components
- `vals`      - The function gradient values

Output Parameter:
- `vals` - The transformed function gradient values

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceTransform()`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceTransformGradient"))
"""
function PetscDualSpaceTransformGradient(petsclib::PetscLibType, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscDualSpaceTransformGradient(petsclib::$UnionPetscLib, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::$PetscInt, Nc::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpaceTransformGradient, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDualSpaceTransformType, PetscBool, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, trans, isInverse, fegeom, Nv, Nc, vals,
              )


	return nothing
end 

"""
	PetscDualSpaceTransformHessian(petsclib::PetscLibType,dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) 
Transform the function Hessian values

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `trans`     - The type of transform
- `isInverse` - Flag to invert the transform
- `fegeom`    - The cell geometry
- `Nv`        - The number of function Hessian samples
- `Nc`        - The number of function components
- `vals`      - The function gradient values

Output Parameter:
- `vals` - The transformed function Hessian values

Level: intermediate

-seealso: `PetscDualSpace`, `PetscDualSpaceTransform()`, `PetscDualSpacePullback()`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransformType`

# External Links
$(_doc_external("Dm/PetscDualSpaceTransformHessian"))
"""
function PetscDualSpaceTransformHessian(petsclib::PetscLibType, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::PetscInt, Nc::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscDualSpaceTransformHessian(petsclib::$UnionPetscLib, dsp::PetscDualSpace, trans::PetscDualSpaceTransformType, isInverse::PetscBool, fegeom::PetscFEGeom, Nv::$PetscInt, Nc::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpaceTransformHessian, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDualSpaceTransformType, PetscBool, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, trans, isInverse, fegeom, Nv, Nc, vals,
              )


	return nothing
end 

"""
	PetscDualSpacePullback(petsclib::PetscLibType,dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) 
Transform the given functional so that it operates on real space, rather than the reference element. Operationally, this means that we map the function evaluations depending on continuity requirements of our finite element method.

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `fegeom`    - The geometry for this cell
- `Nq`        - The number of function samples
- `Nc`        - The number of function components
- `pointEval` - The function values

Output Parameter:
- `pointEval` - The transformed function values

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpacePushforward()`, `PetscDualSpaceTransform()`, `PetscDualSpaceGetDeRahm()`

# External Links
$(_doc_external("Dm/PetscDualSpacePullback"))
"""
function PetscDualSpacePullback(petsclib::PetscLibType, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) end

@for_petsc function PetscDualSpacePullback(petsclib::$UnionPetscLib, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::$PetscInt, Nc::$PetscInt, pointEval::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpacePullback, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, fegeom, Nq, Nc, pointEval,
              )


	return nothing
end 

"""
	PetscDualSpacePushforward(petsclib::PetscLibType,dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) 
Transform the given function so that it operates on real space, rather than the reference element. Operationally, this means that we map the function evaluations depending on continuity requirements of our finite element method.

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `fegeom`    - The geometry for this cell
- `Nq`        - The number of function samples
- `Nc`        - The number of function components
- `pointEval` - The function values

Output Parameter:
- `pointEval` - The transformed function values

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpacePullback()`, `PetscDualSpaceTransform()`, `PetscDualSpaceGetDeRahm()`

# External Links
$(_doc_external("Dm/PetscDualSpacePushforward"))
"""
function PetscDualSpacePushforward(petsclib::PetscLibType, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) end

@for_petsc function PetscDualSpacePushforward(petsclib::$UnionPetscLib, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::$PetscInt, Nc::$PetscInt, pointEval::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpacePushforward, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, fegeom, Nq, Nc, pointEval,
              )


	return nothing
end 

"""
	PetscDualSpacePushforwardGradient(petsclib::PetscLibType,dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) 
Transform the given function gradient so that it operates on real space, rather than the reference element. Operationally, this means that we map the function evaluations depending on continuity requirements of our finite element method.

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `fegeom`    - The geometry for this cell
- `Nq`        - The number of function gradient samples
- `Nc`        - The number of function components
- `pointEval` - The function gradient values

Output Parameter:
- `pointEval` - The transformed function gradient values

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpacePushforward()`, `PetscDualSpacePullback()`, `PetscDualSpaceTransform()`, `PetscDualSpaceGetDeRahm()`

# External Links
$(_doc_external("Dm/PetscDualSpacePushforwardGradient"))
"""
function PetscDualSpacePushforwardGradient(petsclib::PetscLibType, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) end

@for_petsc function PetscDualSpacePushforwardGradient(petsclib::$UnionPetscLib, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::$PetscInt, Nc::$PetscInt, pointEval::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpacePushforwardGradient, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, fegeom, Nq, Nc, pointEval,
              )


	return nothing
end 

"""
	PetscDualSpacePushforwardHessian(petsclib::PetscLibType,dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) 
Transform the given function Hessian so that it operates on real space, rather than the reference element. Operationally, this means that we map the function evaluations depending on continuity requirements of our finite element method.

Input Parameters:
- `dsp`       - The `PetscDualSpace`
- `fegeom`    - The geometry for this cell
- `Nq`        - The number of function Hessian samples
- `Nc`        - The number of function components
- `pointEval` - The function gradient values

Output Parameter:
- `pointEval` - The transformed function Hessian values

Level: advanced

-seealso: `PetscDualSpace`, `PetscDualSpacePushforward()`, `PetscDualSpacePullback()`, `PetscDualSpaceTransform()`, `PetscDualSpaceGetDeRahm()`

# External Links
$(_doc_external("Dm/PetscDualSpacePushforwardHessian"))
"""
function PetscDualSpacePushforwardHessian(petsclib::PetscLibType, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::PetscInt, Nc::PetscInt, pointEval::Vector{PetscScalar}) end

@for_petsc function PetscDualSpacePushforwardHessian(petsclib::$UnionPetscLib, dsp::PetscDualSpace, fegeom::PetscFEGeom, Nq::$PetscInt, Nc::$PetscInt, pointEval::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDualSpacePushforwardHessian, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               dsp, fegeom, Nq, Nc, pointEval,
              )


	return nothing
end 

"""
	PetscDualSpaceSimpleSetDimension(petsclib::PetscLibType,sp::PetscDualSpace, dim::PetscInt) 
Set the number of functionals in the dual space basis

Logically Collective

Input Parameters:
- `sp`  - the `PetscDualSpace`
- `dim` - the basis dimension

Level: intermediate

-seealso: `PETSCDUALSPACESIMPLE`, `PetscDualSpace`, `PetscDualSpaceSimpleSetFunctional()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSimpleSetDimension"))
"""
function PetscDualSpaceSimpleSetDimension(petsclib::PetscLibType, sp::PetscDualSpace, dim::PetscInt) end

@for_petsc function PetscDualSpaceSimpleSetDimension(petsclib::$UnionPetscLib, sp::PetscDualSpace, dim::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceSimpleSetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               sp, dim,
              )


	return nothing
end 

"""
	PetscDualSpaceSimpleSetFunctional(petsclib::PetscLibType,sp::PetscDualSpace, func::PetscInt, q::PetscQuadrature) 
Set the given basis functional for this dual space

Not Collective

Input Parameters:
- `sp`   - the `PetscDualSpace`
- `func` - the basis index
- `q`    - the basis functional

Level: intermediate

-seealso: `PETSCDUALSPACESIMPLE`, `PetscDualSpace`, `PetscDualSpaceSimpleSetDimension()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSimpleSetFunctional"))
"""
function PetscDualSpaceSimpleSetFunctional(petsclib::PetscLibType, sp::PetscDualSpace, func::PetscInt, q::PetscQuadrature) end

@for_petsc function PetscDualSpaceSimpleSetFunctional(petsclib::$UnionPetscLib, sp::PetscDualSpace, func::$PetscInt, q::PetscQuadrature )

    @chk ccall(
               (:PetscDualSpaceSimpleSetFunctional, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, PetscQuadrature),
               sp, func, q,
              )


	return nothing
end 

"""
	continuous::PetscBool = PetscDualSpaceLagrangeGetContinuity(petsclib::PetscLibType,sp::PetscDualSpace) 
Retrieves the flag for element continuity

Not Collective

Input Parameter:
- `sp` - the `PetscDualSpace`

Output Parameter:
- `continuous` - flag for element continuity

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeSetContinuity()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetContinuity"))
"""
function PetscDualSpaceLagrangeGetContinuity(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetContinuity(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	continuous_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetContinuity, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, continuous_,
              )

	continuous = continuous_[]

	return continuous
end 

"""
	PetscDualSpaceLagrangeSetContinuity(petsclib::PetscLibType,sp::PetscDualSpace, continuous::PetscBool) 
Indicate whether the element is continuous

Logically Collective

Input Parameters:
- `sp`         - the `PetscDualSpace`
- `continuous` - flag for element continuity

Options Database Key:
- `-petscdualspace_lagrange_continuity <bool>` - use a continuous element

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeGetContinuity()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetContinuity"))
"""
function PetscDualSpaceLagrangeSetContinuity(petsclib::PetscLibType, sp::PetscDualSpace, continuous::PetscBool) end

@for_petsc function PetscDualSpaceLagrangeSetContinuity(petsclib::$UnionPetscLib, sp::PetscDualSpace, continuous::PetscBool )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetContinuity, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool),
               sp, continuous,
              )


	return nothing
end 

"""
	tensor::PetscBool = PetscDualSpaceLagrangeGetTensor(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the tensor nature of the dual space

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `tensor` - Whether the dual space has tensor layout (vs. simplicial)

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeSetTensor()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetTensor"))
"""
function PetscDualSpaceLagrangeGetTensor(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetTensor(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	tensor_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetTensor, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, tensor_,
              )

	tensor = tensor_[]

	return tensor
end 

"""
	PetscDualSpaceLagrangeSetTensor(petsclib::PetscLibType,sp::PetscDualSpace, tensor::PetscBool) 
Set the tensor nature of the dual space

Not Collective

Input Parameters:
- `sp`     - The `PetscDualSpace`
- `tensor` - Whether the dual space has tensor layout (vs. simplicial)

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeGetTensor()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetTensor"))
"""
function PetscDualSpaceLagrangeSetTensor(petsclib::PetscLibType, sp::PetscDualSpace, tensor::PetscBool) end

@for_petsc function PetscDualSpaceLagrangeSetTensor(petsclib::$UnionPetscLib, sp::PetscDualSpace, tensor::PetscBool )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetTensor, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool),
               sp, tensor,
              )


	return nothing
end 

"""
	trimmed::PetscBool = PetscDualSpaceLagrangeGetTrimmed(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the trimmed nature of the dual space

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `trimmed` - Whether the dual space represents to dual basis of a trimmed polynomial space (e.g. Raviart-Thomas and higher order / other form degree variants)

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeSetTrimmed()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetTrimmed"))
"""
function PetscDualSpaceLagrangeGetTrimmed(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetTrimmed(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	trimmed_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetTrimmed, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, trimmed_,
              )

	trimmed = trimmed_[]

	return trimmed
end 

"""
	PetscDualSpaceLagrangeSetTrimmed(petsclib::PetscLibType,sp::PetscDualSpace, trimmed::PetscBool) 
Set the trimmed nature of the dual space

Not Collective

Input Parameters:
- `sp`      - The `PetscDualSpace`
- `trimmed` - Whether the dual space represents to dual basis of a trimmed polynomial space (e.g. Raviart-Thomas and higher order / other form degree variants)

Level: intermediate

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeGetTrimmed()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetTrimmed"))
"""
function PetscDualSpaceLagrangeSetTrimmed(petsclib::PetscLibType, sp::PetscDualSpace, trimmed::PetscBool) end

@for_petsc function PetscDualSpaceLagrangeSetTrimmed(petsclib::$UnionPetscLib, sp::PetscDualSpace, trimmed::PetscBool )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetTrimmed, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool),
               sp, trimmed,
              )


	return nothing
end 

"""
	nodeType::PetscDTNodeType,boundary::PetscBool,exponent::PetscReal = PetscDualSpaceLagrangeGetNodeType(petsclib::PetscLibType,sp::PetscDualSpace) 
Get a description of how nodes are laid out for Lagrange polynomials in this
dual space

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameters:
- `nodeType` - The type of nodes
- `boundary` - Whether the node type is one that includes endpoints (if nodeType is `PETSCDTNODES_GAUSSJACOBI`, nodes that
include the boundary are Gauss-Lobatto-Jacobi nodes)
- `exponent` - If nodeType is `PETSCDTNODES_GAUSSJACOBI`, indicates the exponent used for both ends of the 1D Jacobi weight function
'0' is Gauss-Legendre, '-0.5' is Gauss-Chebyshev of the first type, '0.5' is Gauss-Chebyshev of the second type

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDTNodeType`, `PetscDualSpaceLagrangeSetNodeType()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetNodeType"))
"""
function PetscDualSpaceLagrangeGetNodeType(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetNodeType(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	nodeType_ = Ref{PetscDTNodeType}()
	boundary_ = Ref{PetscBool}()
	exponent_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetNodeType, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscDTNodeType}, Ptr{PetscBool}, Ptr{$PetscReal}),
               sp, nodeType_, boundary_, exponent_,
              )

	nodeType = unsafe_string(nodeType_[])
	boundary = boundary_[]
	exponent = exponent_[]

	return nodeType,boundary,exponent
end 

"""
	PetscDualSpaceLagrangeSetNodeType(petsclib::PetscLibType,sp::PetscDualSpace, nodeType::PetscDTNodeType, boundary::PetscBool, exponent::PetscReal) 
Set a description of how nodes are laid out for Lagrange polynomials in this
dual space

Logically Collective

Input Parameters:
- `sp`       - The `PetscDualSpace`
- `nodeType` - The type of nodes
- `boundary` - Whether the node type is one that includes endpoints (if nodeType is `PETSCDTNODES_GAUSSJACOBI`, nodes that
include the boundary are Gauss-Lobatto-Jacobi nodes)
- `exponent` - If nodeType is `PETSCDTNODES_GAUSSJACOBI`, indicates the exponent used for both ends of the 1D Jacobi weight function
'0' is Gauss-Legendre, '-0.5' is Gauss-Chebyshev of the first type, '0.5' is Gauss-Chebyshev of the second type

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDTNodeType`, `PetscDualSpaceLagrangeGetNodeType()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetNodeType"))
"""
function PetscDualSpaceLagrangeSetNodeType(petsclib::PetscLibType, sp::PetscDualSpace, nodeType::PetscDTNodeType, boundary::PetscBool, exponent::PetscReal) end

@for_petsc function PetscDualSpaceLagrangeSetNodeType(petsclib::$UnionPetscLib, sp::PetscDualSpace, nodeType::PetscDTNodeType, boundary::PetscBool, exponent::$PetscReal )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetNodeType, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscDTNodeType, PetscBool, $PetscReal),
               sp, nodeType, boundary, exponent,
              )


	return nothing
end 

"""
	useMoments::PetscBool = PetscDualSpaceLagrangeGetUseMoments(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the flag for using moment functionals

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `useMoments` - Moment flag

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeSetUseMoments()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetUseMoments"))
"""
function PetscDualSpaceLagrangeGetUseMoments(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetUseMoments(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	useMoments_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetUseMoments, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, useMoments_,
              )

	useMoments = useMoments_[]

	return useMoments
end 

"""
	PetscDualSpaceLagrangeSetUseMoments(petsclib::PetscLibType,sp::PetscDualSpace, useMoments::PetscBool) 
Set the flag for moment functionals

Logically Collective

Input Parameters:
- `sp`         - The `PetscDualSpace`
- `useMoments` - The flag for moment functionals

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeGetUseMoments()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetUseMoments"))
"""
function PetscDualSpaceLagrangeSetUseMoments(petsclib::PetscLibType, sp::PetscDualSpace, useMoments::PetscBool) end

@for_petsc function PetscDualSpaceLagrangeSetUseMoments(petsclib::$UnionPetscLib, sp::PetscDualSpace, useMoments::PetscBool )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetUseMoments, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool),
               sp, useMoments,
              )


	return nothing
end 

"""
	order::PetscInt = PetscDualSpaceLagrangeGetMomentOrder(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the order for moment integration

Not Collective

Input Parameter:
- `sp` - The `PetscDualSpace`

Output Parameter:
- `order` - Moment integration order

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeSetMomentOrder()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeGetMomentOrder"))
"""
function PetscDualSpaceLagrangeGetMomentOrder(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceLagrangeGetMomentOrder(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceLagrangeGetMomentOrder, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, order_,
              )

	order = order_[]

	return order
end 

"""
	PetscDualSpaceLagrangeSetMomentOrder(petsclib::PetscLibType,sp::PetscDualSpace, order::PetscInt) 
Set the order for moment integration

Logically Collective

Input Parameters:
- `sp`    - The `PetscDualSpace`
- `order` - The order for moment integration

Level: advanced

-seealso: `PETSCDUALSPACELAGRANGE`, `PetscDualSpace`, `PetscDualSpaceLagrangeGetMomentOrder()`

# External Links
$(_doc_external("Dm/PetscDualSpaceLagrangeSetMomentOrder"))
"""
function PetscDualSpaceLagrangeSetMomentOrder(petsclib::PetscLibType, sp::PetscDualSpace, order::PetscInt) end

@for_petsc function PetscDualSpaceLagrangeSetMomentOrder(petsclib::$UnionPetscLib, sp::PetscDualSpace, order::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceLagrangeSetMomentOrder, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               sp, order,
              )


	return nothing
end 

"""
	numSumSpaces::PetscInt = PetscDualSpaceSumGetNumSubspaces(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the number of spaces in the sum space

Input Parameter:
- `sp` - the dual space object

Output Parameter:
- `numSumSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumSetNumSubspaces()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumGetNumSubspaces"))
"""
function PetscDualSpaceSumGetNumSubspaces(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSumGetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	numSumSpaces_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDualSpaceSumGetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{$PetscInt}),
               sp, numSumSpaces_,
              )

	numSumSpaces = numSumSpaces_[]

	return numSumSpaces
end 

"""
	PetscDualSpaceSumSetNumSubspaces(petsclib::PetscLibType,sp::PetscDualSpace, numSumSpaces::PetscInt) 
Set the number of spaces in the sum space

Input Parameters:
- `sp`           - the dual space object
- `numSumSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumGetNumSubspaces()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumSetNumSubspaces"))
"""
function PetscDualSpaceSumSetNumSubspaces(petsclib::PetscLibType, sp::PetscDualSpace, numSumSpaces::PetscInt) end

@for_petsc function PetscDualSpaceSumSetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscDualSpace, numSumSpaces::$PetscInt )

    @chk ccall(
               (:PetscDualSpaceSumSetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt),
               sp, numSumSpaces,
              )


	return nothing
end 

"""
	concatenate::PetscBool = PetscDualSpaceSumGetConcatenate(petsclib::PetscLibType,sp::PetscDualSpace) 
Get the concatenate flag for this space.

Input Parameter:
- `sp` - the dual space object

Output Parameter:
- `concatenate` - flag indicating whether subspaces are concatenated.

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumSetConcatenate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumGetConcatenate"))
"""
function PetscDualSpaceSumGetConcatenate(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSumGetConcatenate(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	concatenate_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceSumGetConcatenate, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}),
               sp, concatenate_,
              )

	concatenate = concatenate_[]

	return concatenate
end 

"""
	PetscDualSpaceSumSetConcatenate(petsclib::PetscLibType,sp::PetscDualSpace, concatenate::PetscBool) 
Sets the concatenate flag for this space.

Input Parameters:
- `sp`          - the dual space object
- `concatenate` - are subspaces concatenated components (true) or direct summands (false)

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumGetConcatenate()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumSetConcatenate"))
"""
function PetscDualSpaceSumSetConcatenate(petsclib::PetscLibType, sp::PetscDualSpace, concatenate::PetscBool) end

@for_petsc function PetscDualSpaceSumSetConcatenate(petsclib::$UnionPetscLib, sp::PetscDualSpace, concatenate::PetscBool )

    @chk ccall(
               (:PetscDualSpaceSumSetConcatenate, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool),
               sp, concatenate,
              )


	return nothing
end 

"""
	PetscDualSpaceSumGetSubspace(petsclib::PetscLibType,sp::PetscDualSpace, s::PetscInt, subsp::PetscDualSpace) 
Get a space in the sum space

Input Parameters:
- `sp` - the dual space object
- `s`  - The space number

Output Parameter:
- `subsp` - the `PetscDualSpace`

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumSetSubspace()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumGetSubspace"))
"""
function PetscDualSpaceSumGetSubspace(petsclib::PetscLibType, sp::PetscDualSpace, s::PetscInt, subsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSumGetSubspace(petsclib::$UnionPetscLib, sp::PetscDualSpace, s::$PetscInt, subsp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceSumGetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, Ptr{PetscDualSpace}),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscDualSpaceSumSetSubspace(petsclib::PetscLibType,sp::PetscDualSpace, s::PetscInt, subsp::PetscDualSpace) 
Set a space in the sum space

Input Parameters:
- `sp`    - the dual space object
- `s`     - The space number
- `subsp` - the number of spaces

Level: intermediate

-seealso: `PETSCDUALSPACESUM`, `PetscDualSpace`, `PetscDualSpaceSumGetSubspace()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumSetSubspace"))
"""
function PetscDualSpaceSumSetSubspace(petsclib::PetscLibType, sp::PetscDualSpace, s::PetscInt, subsp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSumSetSubspace(petsclib::$UnionPetscLib, sp::PetscDualSpace, s::$PetscInt, subsp::PetscDualSpace )

    @chk ccall(
               (:PetscDualSpaceSumSetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, $PetscInt, PetscDualSpace),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscDualSpaceSumSetInterleave(petsclib::PetscLibType,sp::PetscDualSpace, interleave_basis::PetscBool, interleave_components::PetscBool) 
Set whether the basis functions and components of a uniform sum are interleaved

Logically collective

Input Parameters:
- `sp`                    - a `PetscDualSpace` of type `PETSCDUALSPACESUM`
- `interleave_basis`      - if `PETSC_TRUE`, the basis vectors of the subspaces are interleaved
- `interleave_components` - if `PETSC_TRUE` and the space concatenates components (`PetscDualSpaceSumGetConcatenate()`),
interleave the concatenated components

Level: developer

-seealso: `PetscDualSpace`, `PETSCDUALSPACESUM`, `PETSCFEVECTOR`, `PetscDualSpaceSumGetInterleave()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumSetInterleave"))
"""
function PetscDualSpaceSumSetInterleave(petsclib::PetscLibType, sp::PetscDualSpace, interleave_basis::PetscBool, interleave_components::PetscBool) end

@for_petsc function PetscDualSpaceSumSetInterleave(petsclib::$UnionPetscLib, sp::PetscDualSpace, interleave_basis::PetscBool, interleave_components::PetscBool )

    @chk ccall(
               (:PetscDualSpaceSumSetInterleave, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, PetscBool, PetscBool),
               sp, interleave_basis, interleave_components,
              )


	return nothing
end 

"""
	interleave_basis::PetscBool,interleave_components::PetscBool = PetscDualSpaceSumGetInterleave(petsclib::PetscLibType,sp::PetscDualSpace) 
Get whether the basis functions and components of a uniform sum are interleaved

Logically collective

Input Parameter:
- `sp` - a `PetscDualSpace` of type `PETSCDUALSPACESUM`

Output Parameters:
- `interleave_basis`      - if `PETSC_TRUE`, the basis vectors of the subspaces are interleaved
- `interleave_components` - if `PETSC_TRUE` and the space concatenates components (`PetscDualSpaceSumGetConcatenate()`),
interleave the concatenated components

Level: developer

-seealso: `PetscDualSpace`, `PETSCDUALSPACESUM`, `PETSCFEVECTOR`, `PetscDualSpaceSumSetInterleave()`

# External Links
$(_doc_external("Dm/PetscDualSpaceSumGetInterleave"))
"""
function PetscDualSpaceSumGetInterleave(petsclib::PetscLibType, sp::PetscDualSpace) end

@for_petsc function PetscDualSpaceSumGetInterleave(petsclib::$UnionPetscLib, sp::PetscDualSpace )
	interleave_basis_ = Ref{PetscBool}()
	interleave_components_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDualSpaceSumGetInterleave, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscBool}, Ptr{PetscBool}),
               sp, interleave_basis_, interleave_components_,
              )

	interleave_basis = interleave_basis_[]
	interleave_components = interleave_components_[]

	return interleave_basis,interleave_components
end 

"""
	sumSpace::PetscDualSpace = PetscDualSpaceCreateSum(petsclib::PetscLibType,numSubspaces::PetscInt, subspaces::Vector{PetscDualSpace}, concatenate::PetscBool) 
Create a finite element dual basis that is the sum of other dual bases

Collective

Input Parameters:
- `numSubspaces` - the number of spaces that will be added together
- `subspaces`    - an array of length `numSubspaces` of spaces
- `concatenate`  - if `PETSC_FALSE`, the sum-space has the same components as the individual dual spaces (`PetscDualSpaceGetNumComponents()`); if `PETSC_TRUE`, the individual components are concatenated to create a dual space with more components

Output Parameter:
- `sumSpace` - a `PetscDualSpace` of type `PETSCDUALSPACESUM`

Level: advanced

-seealso: `PetscDualSpace`, `PETSCDUALSPACESUM`, `PETSCSPACESUM`

# External Links
$(_doc_external("Dm/PetscDualSpaceCreateSum"))
"""
function PetscDualSpaceCreateSum(petsclib::PetscLibType, numSubspaces::PetscInt, subspaces::Vector{PetscDualSpace}, concatenate::PetscBool) end

@for_petsc function PetscDualSpaceCreateSum(petsclib::$UnionPetscLib, numSubspaces::$PetscInt, subspaces::Vector{PetscDualSpace}, concatenate::PetscBool )
	sumSpace_ = Ref{PetscDualSpace}()

    @chk ccall(
               (:PetscDualSpaceCreateSum, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscDualSpace}, PetscBool, Ptr{PetscDualSpace}),
               numSubspaces, subspaces, concatenate, sumSpace_,
              )

	sumSpace = sumSpace_[]

	return sumSpace
end 

"""
	PetscDualSpaceRefinedSetCellSpaces(petsclib::PetscLibType,sp::PetscDualSpace, cellSpaces::Vector{PetscDualSpace}) 
Set the dual spaces for the closures of each of the cells
in the multicell `DM` of a `PetscDualSpace`

Collective

Input Parameters:
- `sp`         - a `PetscDualSpace`
- `cellSpaces` - one `PetscDualSpace` for each of the cells.  The reference count of each cell space will be incremented,
so the user is still responsible for these spaces afterwards

Level: intermediate

-seealso: `PETSCDUALSPACEREFINED`, `PetscDualSpace`, `PetscFERefine()`

# External Links
$(_doc_external("Dm/PetscDualSpaceRefinedSetCellSpaces"))
"""
function PetscDualSpaceRefinedSetCellSpaces(petsclib::PetscLibType, sp::PetscDualSpace, cellSpaces::Vector{PetscDualSpace}) end

@for_petsc function PetscDualSpaceRefinedSetCellSpaces(petsclib::$UnionPetscLib, sp::PetscDualSpace, cellSpaces::Vector{PetscDualSpace} )

    @chk ccall(
               (:PetscDualSpaceRefinedSetCellSpaces, $petsc_library),
               PetscErrorCode,
               (PetscDualSpace, Ptr{PetscDualSpace}),
               sp, cellSpaces,
              )


	return nothing
end 

