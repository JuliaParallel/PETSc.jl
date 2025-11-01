# autodefined type arguments for class ------
mutable struct _n_PetscFV end
const PetscFV = Ptr{_n_PetscFV}

# -------------------------------------------------------
"""
	PetscFVFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `PetscFV` package. It is called
from `PetscFinalize()`.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("DM/PetscFVFinalizePackage"))
"""
function PetscFVFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscFVFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFVFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFVInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscFV` package. It is called
from `PetscDLLibraryRegister()` when using dynamic libraries, and on the first call to `PetscFVCreate()`
when using static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("DM/PetscFVInitializePackage"))
"""
function PetscFVInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscFVInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFVInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFVRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscFV` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine itself

-seealso: `PetscFV`, `PetscFVType`, `PetscFVRegisterAll()`, `PetscFVRegisterDestroy()`

# External Links
$(_doc_external("DM/PetscFVRegister"))
"""
function PetscFVRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscFVRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscFVRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscFVSetType(petsclib::PetscLibType,fvm::PetscFV, name::PetscFVType) 
Builds a particular `PetscFV`

Collective

Input Parameters:
- `fvm`  - The `PetscFV` object
- `name` - The type of FVM space

Options Database Key:
- `-petscfv_type <type>` - Sets the `PetscFVType`; use -help for a list of available types

Level: intermediate

-seealso: `PetscFV`, `PetscFVType`, `PetscFVGetType()`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVSetType"))
"""
function PetscFVSetType(petsclib::PetscLibType, fvm::PetscFV, name::PetscFVType) end

@for_petsc function PetscFVSetType(petsclib::$UnionPetscLib, fvm::PetscFV, name::PetscFVType )

    @chk ccall(
               (:PetscFVSetType, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscFVType),
               fvm, name,
              )


	return nothing
end 

"""
	name::PetscFVType = PetscFVGetType(petsclib::PetscLibType,fvm::PetscFV) 
Gets the `PetscFVType` (as a string) from a `PetscFV`.

Not Collective

Input Parameter:
- `fvm` - The `PetscFV`

Output Parameter:
- `name` - The `PetscFVType` name

Level: intermediate

-seealso: `PetscFV`, `PetscFVType`, `PetscFVSetType()`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVGetType"))
"""
function PetscFVGetType(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVGetType(petsclib::$UnionPetscLib, fvm::PetscFV )
	name_ = Ref{PetscFVType}()

    @chk ccall(
               (:PetscFVGetType, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscFVType}),
               fvm, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscFVViewFromOptions(petsclib::PetscLibType,A::PetscFV, obj::PetscObject, name::String) 
View a `PetscFV` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscFV` object
- `obj`  - Optional object that provides the options prefix
- `name` - command line option name

Level: intermediate

-seealso: `PetscFV`, `PetscFVView()`, `PetscObjectViewFromOptions()`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVViewFromOptions"))
"""
function PetscFVViewFromOptions(petsclib::PetscLibType, A::PetscFV, obj::PetscObject, name::String) end

@for_petsc function PetscFVViewFromOptions(petsclib::$UnionPetscLib, A::PetscFV, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscFVViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscFVView(petsclib::PetscLibType,fvm::PetscFV, v::PetscViewer) 
Views a `PetscFV`

Collective

Input Parameters:
- `fvm` - the `PetscFV` object to view
- `v`   - the viewer

Level: beginner

-seealso: `PetscFV`, `PetscViewer`, `PetscFVDestroy()`

# External Links
$(_doc_external("DM/PetscFVView"))
"""
function PetscFVView(petsclib::PetscLibType, fvm::PetscFV, v::PetscViewer) end

@for_petsc function PetscFVView(petsclib::$UnionPetscLib, fvm::PetscFV, v::PetscViewer )

    @chk ccall(
               (:PetscFVView, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscViewer),
               fvm, v,
              )


	return nothing
end 

"""
	PetscFVSetFromOptions(petsclib::PetscLibType,fvm::PetscFV) 
sets parameters in a `PetscFV` from the options database

Collective

Input Parameter:
- `fvm` - the `PetscFV` object to set options for

Options Database Key:
- `-petscfv_compute_gradients <bool>` - Determines whether cell gradients are calculated

Level: intermediate

-seealso: `PetscFV`, `PetscFVView()`

# External Links
$(_doc_external("DM/PetscFVSetFromOptions"))
"""
function PetscFVSetFromOptions(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVSetFromOptions(petsclib::$UnionPetscLib, fvm::PetscFV )

    @chk ccall(
               (:PetscFVSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscFV,),
               fvm,
              )


	return nothing
end 

"""
	PetscFVSetUp(petsclib::PetscLibType,fvm::PetscFV) 
Setup the data structures for the `PetscFV` based on the `PetscFVType` provided by `PetscFVSetType()`

Collective

Input Parameter:
- `fvm` - the `PetscFV` object to setup

Level: intermediate

-seealso: `PetscFV`, `PetscFVView()`, `PetscFVDestroy()`

# External Links
$(_doc_external("DM/PetscFVSetUp"))
"""
function PetscFVSetUp(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVSetUp(petsclib::$UnionPetscLib, fvm::PetscFV )

    @chk ccall(
               (:PetscFVSetUp, $petsc_library),
               PetscErrorCode,
               (PetscFV,),
               fvm,
              )


	return nothing
end 

"""
	PetscFVDestroy(petsclib::PetscLibType,fvm::PetscFV) 
Destroys a `PetscFV` object

Collective

Input Parameter:
- `fvm` - the `PetscFV` object to destroy

Level: beginner

-seealso: `PetscFV`, `PetscFVCreate()`, `PetscFVView()`

# External Links
$(_doc_external("DM/PetscFVDestroy"))
"""
function PetscFVDestroy(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVDestroy(petsclib::$UnionPetscLib, fvm::PetscFV )

    @chk ccall(
               (:PetscFVDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFV},),
               fvm,
              )


	return nothing
end 

"""
	fvm::PetscFV = PetscFVCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscFV` object. The type can then be set with `PetscFVSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscFV` object

Output Parameter:
- `fvm` - The `PetscFV` object

Level: beginner

-seealso: `PetscFVSetUp()`, `PetscFVSetType()`, `PETSCFVUPWIND`, `PetscFVDestroy()`

# External Links
$(_doc_external("DM/PetscFVCreate"))
"""
function PetscFVCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscFVCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	fvm_ = Ref{PetscFV}()

    @chk ccall(
               (:PetscFVCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscFV}),
               comm, fvm_,
              )

	fvm = fvm_[]

	return fvm
end 

"""
	PetscFVSetLimiter(petsclib::PetscLibType,fvm::PetscFV, lim::PetscLimiter) 
Set the `PetscLimiter` to the `PetscFV`

Logically Collective

Input Parameters:
- `fvm` - the `PetscFV` object
- `lim` - The `PetscLimiter`

Level: intermediate

-seealso: `PetscFV`, `PetscLimiter`, `PetscFVGetLimiter()`

# External Links
$(_doc_external("DM/PetscFVSetLimiter"))
"""
function PetscFVSetLimiter(petsclib::PetscLibType, fvm::PetscFV, lim::PetscLimiter) end

@for_petsc function PetscFVSetLimiter(petsclib::$UnionPetscLib, fvm::PetscFV, lim::PetscLimiter )

    @chk ccall(
               (:PetscFVSetLimiter, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscLimiter),
               fvm, lim,
              )


	return nothing
end 

"""
	PetscFVGetLimiter(petsclib::PetscLibType,fvm::PetscFV, lim::PetscLimiter) 
Get the `PetscLimiter` object from the `PetscFV`

Not Collective

Input Parameter:
- `fvm` - the `PetscFV` object

Output Parameter:
- `lim` - The `PetscLimiter`

Level: intermediate

-seealso: `PetscFV`, `PetscLimiter`, `PetscFVSetLimiter()`

# External Links
$(_doc_external("DM/PetscFVGetLimiter"))
"""
function PetscFVGetLimiter(petsclib::PetscLibType, fvm::PetscFV, lim::PetscLimiter) end

@for_petsc function PetscFVGetLimiter(petsclib::$UnionPetscLib, fvm::PetscFV, lim::PetscLimiter )

    @chk ccall(
               (:PetscFVGetLimiter, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscLimiter}),
               fvm, lim,
              )


	return nothing
end 

"""
	PetscFVSetNumComponents(petsclib::PetscLibType,fvm::PetscFV, comp::PetscInt) 
Set the number of field components in a `PetscFV`

Logically Collective

Input Parameters:
- `fvm`  - the `PetscFV` object
- `comp` - The number of components

Level: intermediate

-seealso: `PetscFV`, `PetscFVGetNumComponents()`

# External Links
$(_doc_external("DM/PetscFVSetNumComponents"))
"""
function PetscFVSetNumComponents(petsclib::PetscLibType, fvm::PetscFV, comp::PetscInt) end

@for_petsc function PetscFVSetNumComponents(petsclib::$UnionPetscLib, fvm::PetscFV, comp::$PetscInt )

    @chk ccall(
               (:PetscFVSetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt),
               fvm, comp,
              )


	return nothing
end 

"""
	comp::PetscInt = PetscFVGetNumComponents(petsclib::PetscLibType,fvm::PetscFV) 
Get the number of field components in a `PetscFV`

Not Collective

Input Parameter:
- `fvm` - the `PetscFV` object

Output Parameter:
- `comp` - The number of components

Level: intermediate

-seealso: `PetscFV`, `PetscFVSetNumComponents()`, `PetscFVSetComponentName()`

# External Links
$(_doc_external("DM/PetscFVGetNumComponents"))
"""
function PetscFVGetNumComponents(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVGetNumComponents(petsclib::$UnionPetscLib, fvm::PetscFV )
	comp_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFVGetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{$PetscInt}),
               fvm, comp_,
              )

	comp = comp_[]

	return comp
end 

"""
	PetscFVSetComponentName(petsclib::PetscLibType,fvm::PetscFV, comp::PetscInt, name::String) 
Set the name of a component (used in output and viewing) in a `PetscFV`

Logically Collective

Input Parameters:
- `fvm`  - the `PetscFV` object
- `comp` - the component number
- `name` - the component name

Level: intermediate

-seealso: `PetscFV`, `PetscFVGetComponentName()`

# External Links
$(_doc_external("DM/PetscFVSetComponentName"))
"""
function PetscFVSetComponentName(petsclib::PetscLibType, fvm::PetscFV, comp::PetscInt, name::String) end

@for_petsc function PetscFVSetComponentName(petsclib::$UnionPetscLib, fvm::PetscFV, comp::$PetscInt, name::String )

    @chk ccall(
               (:PetscFVSetComponentName, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt, Ptr{Cchar}),
               fvm, comp, name,
              )


	return nothing
end 

"""
	PetscFVGetComponentName(petsclib::PetscLibType,fvm::PetscFV, comp::PetscInt, name::String) 
Get the name of a component (used in output and viewing) in a `PetscFV`

Logically Collective

Input Parameters:
- `fvm`  - the `PetscFV` object
- `comp` - the component number

Output Parameter:
- `name` - the component name

Level: intermediate

-seealso: `PetscFV`, `PetscFVSetComponentName()`

# External Links
$(_doc_external("DM/PetscFVGetComponentName"))
"""
function PetscFVGetComponentName(petsclib::PetscLibType, fvm::PetscFV, comp::PetscInt, name::String) end

@for_petsc function PetscFVGetComponentName(petsclib::$UnionPetscLib, fvm::PetscFV, comp::$PetscInt, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscFVGetComponentName, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt, Ptr{Ptr{Cchar}}),
               fvm, comp, name_,
              )


	return nothing
end 

"""
	PetscFVSetSpatialDimension(petsclib::PetscLibType,fvm::PetscFV, dim::PetscInt) 
Set the spatial dimension of a `PetscFV`

Logically Collective

Input Parameters:
- `fvm` - the `PetscFV` object
- `dim` - The spatial dimension

Level: intermediate

-seealso: `PetscFV`, `PetscFVGetSpatialDimension()`

# External Links
$(_doc_external("DM/PetscFVSetSpatialDimension"))
"""
function PetscFVSetSpatialDimension(petsclib::PetscLibType, fvm::PetscFV, dim::PetscInt) end

@for_petsc function PetscFVSetSpatialDimension(petsclib::$UnionPetscLib, fvm::PetscFV, dim::$PetscInt )

    @chk ccall(
               (:PetscFVSetSpatialDimension, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt),
               fvm, dim,
              )


	return nothing
end 

"""
	dim::PetscInt = PetscFVGetSpatialDimension(petsclib::PetscLibType,fvm::PetscFV) 
Get the spatial dimension of a `PetscFV`

Not Collective

Input Parameter:
- `fvm` - the `PetscFV` object

Output Parameter:
- `dim` - The spatial dimension

Level: intermediate

-seealso: `PetscFV`, `PetscFVSetSpatialDimension()`

# External Links
$(_doc_external("DM/PetscFVGetSpatialDimension"))
"""
function PetscFVGetSpatialDimension(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVGetSpatialDimension(petsclib::$UnionPetscLib, fvm::PetscFV )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFVGetSpatialDimension, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{$PetscInt}),
               fvm, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	PetscFVSetComputeGradients(petsclib::PetscLibType,fvm::PetscFV, computeGradients::PetscBool) 
Toggle computation of cell gradients on a `PetscFV`

Logically Collective

Input Parameters:
- `fvm`              - the `PetscFV` object
- `computeGradients` - Flag to compute cell gradients

Level: intermediate

-seealso: `PetscFV`, `PetscFVGetComputeGradients()`

# External Links
$(_doc_external("DM/PetscFVSetComputeGradients"))
"""
function PetscFVSetComputeGradients(petsclib::PetscLibType, fvm::PetscFV, computeGradients::PetscBool) end

@for_petsc function PetscFVSetComputeGradients(petsclib::$UnionPetscLib, fvm::PetscFV, computeGradients::PetscBool )

    @chk ccall(
               (:PetscFVSetComputeGradients, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscBool),
               fvm, computeGradients,
              )


	return nothing
end 

"""
	computeGradients::PetscBool = PetscFVGetComputeGradients(petsclib::PetscLibType,fvm::PetscFV) 
Return flag for computation of cell gradients on a `PetscFV`

Not Collective

Input Parameter:
- `fvm` - the `PetscFV` object

Output Parameter:
- `computeGradients` - Flag to compute cell gradients

Level: intermediate

-seealso: `PetscFV`, `PetscFVSetComputeGradients()`

# External Links
$(_doc_external("DM/PetscFVGetComputeGradients"))
"""
function PetscFVGetComputeGradients(petsclib::PetscLibType, fvm::PetscFV) end

@for_petsc function PetscFVGetComputeGradients(petsclib::$UnionPetscLib, fvm::PetscFV )
	computeGradients_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscFVGetComputeGradients, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscBool}),
               fvm, computeGradients_,
              )

	computeGradients = computeGradients_[]

	return computeGradients
end 

"""
	PetscFVSetQuadrature(petsclib::PetscLibType,fvm::PetscFV, q::PetscQuadrature) 
Set the `PetscQuadrature` object for a `PetscFV`

Logically Collective

Input Parameters:
- `fvm` - the `PetscFV` object
- `q`   - The `PetscQuadrature`

Level: intermediate

-seealso: `PetscQuadrature`, `PetscFV`, `PetscFVGetQuadrature()`

# External Links
$(_doc_external("DM/PetscFVSetQuadrature"))
"""
function PetscFVSetQuadrature(petsclib::PetscLibType, fvm::PetscFV, q::PetscQuadrature) end

@for_petsc function PetscFVSetQuadrature(petsclib::$UnionPetscLib, fvm::PetscFV, q::PetscQuadrature )

    @chk ccall(
               (:PetscFVSetQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscQuadrature),
               fvm, q,
              )


	return nothing
end 

"""
	PetscFVGetQuadrature(petsclib::PetscLibType,fvm::PetscFV, q::PetscQuadrature) 
Get the `PetscQuadrature` from a `PetscFV`

Not Collective

Input Parameter:
- `fvm` - the `PetscFV` object

Output Parameter:
- `q` - The `PetscQuadrature`

Level: intermediate

-seealso: `PetscQuadrature`, `PetscFV`, `PetscFVSetQuadrature()`

# External Links
$(_doc_external("DM/PetscFVGetQuadrature"))
"""
function PetscFVGetQuadrature(petsclib::PetscLibType, fvm::PetscFV, q::PetscQuadrature) end

@for_petsc function PetscFVGetQuadrature(petsclib::$UnionPetscLib, fvm::PetscFV, q::PetscQuadrature )

    @chk ccall(
               (:PetscFVGetQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscQuadrature}),
               fvm, q,
              )


	return nothing
end 

"""
	PetscFVCreateDualSpace(petsclib::PetscLibType,fvm::PetscFV, ct::DMPolytopeType) 
Creates a `PetscDualSpace` appropriate for the `PetscFV`

Not Collective

Input Parameters:
- `fvm` - The `PetscFV` object
- `ct`  - The `DMPolytopeType` for the cell

Level: intermediate

-seealso: `PetscFVGetDualSpace()`, `PetscFVSetDualSpace()`, `PetscDualSpace`, `PetscFV`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVCreateDualSpace"))
"""
function PetscFVCreateDualSpace(petsclib::PetscLibType, fvm::PetscFV, ct::DMPolytopeType) end

@for_petsc function PetscFVCreateDualSpace(petsclib::$UnionPetscLib, fvm::PetscFV, ct::DMPolytopeType )

    @chk ccall(
               (:PetscFVCreateDualSpace, $petsc_library),
               PetscErrorCode,
               (PetscFV, DMPolytopeType),
               fvm, ct,
              )


	return nothing
end 

"""
	PetscFVGetDualSpace(petsclib::PetscLibType,fvm::PetscFV, sp::PetscDualSpace) 
Returns the `PetscDualSpace` used to define the inner product on a `PetscFV`

Not Collective

Input Parameter:
- `fvm` - The `PetscFV` object

Output Parameter:
- `sp` - The `PetscDualSpace` object

Level: intermediate

-seealso: `PetscFVSetDualSpace()`, `PetscFVCreateDualSpace()`, `PetscDualSpace`, `PetscFV`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVGetDualSpace"))
"""
function PetscFVGetDualSpace(petsclib::PetscLibType, fvm::PetscFV, sp::PetscDualSpace) end

@for_petsc function PetscFVGetDualSpace(petsclib::$UnionPetscLib, fvm::PetscFV, sp::PetscDualSpace )

    @chk ccall(
               (:PetscFVGetDualSpace, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscDualSpace}),
               fvm, sp,
              )


	return nothing
end 

"""
	PetscFVSetDualSpace(petsclib::PetscLibType,fvm::PetscFV, sp::PetscDualSpace) 
Sets the `PetscDualSpace` used to define the inner product

Not Collective

Input Parameters:
- `fvm` - The `PetscFV` object
- `sp`  - The `PetscDualSpace` object

Level: intermediate

-seealso: `PetscFVGetDualSpace()`, `PetscFVCreateDualSpace()`, `PetscDualSpace`, `PetscFV`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVSetDualSpace"))
"""
function PetscFVSetDualSpace(petsclib::PetscLibType, fvm::PetscFV, sp::PetscDualSpace) end

@for_petsc function PetscFVSetDualSpace(petsclib::$UnionPetscLib, fvm::PetscFV, sp::PetscDualSpace )

    @chk ccall(
               (:PetscFVSetDualSpace, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscDualSpace),
               fvm, sp,
              )


	return nothing
end 

"""
	PetscFVGetCellTabulation(petsclib::PetscLibType,fvm::PetscFV, T::PetscTabulation) 
Returns the tabulation of the basis functions at the quadrature points

Not Collective

Input Parameter:
- `fvm` - The `PetscFV` object

Output Parameter:
- `T` - The basis function values and derivatives at quadrature points

Level: intermediate

-seealso: `PetscFV`, `PetscTabulation`, `PetscFEGetCellTabulation()`, `PetscFVCreateTabulation()`, `PetscFVGetQuadrature()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("DM/PetscFVGetCellTabulation"))
"""
function PetscFVGetCellTabulation(petsclib::PetscLibType, fvm::PetscFV, T::PetscTabulation) end

@for_petsc function PetscFVGetCellTabulation(petsclib::$UnionPetscLib, fvm::PetscFV, T::PetscTabulation )

    @chk ccall(
               (:PetscFVGetCellTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscTabulation}),
               fvm, T,
              )


	return nothing
end 

"""
	T::PetscTabulation = PetscFVCreateTabulation(petsclib::PetscLibType,fvm::PetscFV, nrepl::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt) 
Tabulates the basis functions, and perhaps derivatives, at the points provided.

Not Collective

Input Parameters:
- `fvm`     - The `PetscFV` object
- `nrepl`   - The number of replicas
- `npoints` - The number of tabulation points in a replica
- `points`  - The tabulation point coordinates
- `K`       - The order of derivative to tabulate

Output Parameter:
- `T` - The basis function values and derivative at tabulation points

Level: intermediate

-seealso: `PetscFV`, `PetscTabulation`, `PetscFECreateTabulation()`, `PetscTabulationDestroy()`, `PetscFEGetCellTabulation()`

# External Links
$(_doc_external("DM/PetscFVCreateTabulation"))
"""
function PetscFVCreateTabulation(petsclib::PetscLibType, fvm::PetscFV, nrepl::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt) end

@for_petsc function PetscFVCreateTabulation(petsclib::$UnionPetscLib, fvm::PetscFV, nrepl::$PetscInt, npoints::$PetscInt, points::Vector{$PetscReal}, K::$PetscInt )
	T_ = Ref{PetscTabulation}()

    @chk ccall(
               (:PetscFVCreateTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{PetscTabulation}),
               fvm, nrepl, npoints, points, K, T_,
              )

	T = T_[]

	return T
end 

"""
	grad::Vector{PetscScalar} = PetscFVComputeGradient(petsclib::PetscLibType,fvm::PetscFV, numFaces::PetscInt, dx::Vector{PetscScalar}) 
Compute the gradient reconstruction matrix for a given cell

Input Parameters:
- `fvm`      - The `PetscFV` object
- `numFaces` - The number of cell faces which are not constrained
- `dx`       - The vector from the cell centroid to the neighboring cell centroid for each face

Output Parameter:
- `grad` - the gradient

Level: advanced

-seealso: `PetscFV`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVComputeGradient"))
"""
function PetscFVComputeGradient(petsclib::PetscLibType, fvm::PetscFV, numFaces::PetscInt, dx::Vector{PetscScalar}) end

@for_petsc function PetscFVComputeGradient(petsclib::$UnionPetscLib, fvm::PetscFV, numFaces::$PetscInt, dx::Vector{$PetscScalar} )
	grad = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFVComputeGradient, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               fvm, numFaces, dx, grad,
              )


	return grad
end 

"""
	fluxL::Vector{PetscScalar},fluxR::Vector{PetscScalar} = PetscFVIntegrateRHSFunction(petsclib::PetscLibType,fvm::PetscFV, prob::PetscDS, field::PetscInt, Nf::PetscInt, fgeom::PetscFVFaceGeom, neighborVol::PetscReal, uL::Vector{PetscScalar}, uR::Vector{PetscScalar}) 
Produce the cell residual vector for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `fvm`         - The `PetscFV` object for the field being integrated
- `prob`        - The `PetscDS` specifying the discretizations and continuum functions
- `field`       - The field being integrated
- `Nf`          - The number of faces in the chunk
- `fgeom`       - The face geometry for each face in the chunk
- `neighborVol` - The volume for each pair of cells in the chunk
- `uL`          - The state from the cell on the left
- `uR`          - The state from the cell on the right

Output Parameters:
- `fluxL` - the left fluxes for each face
- `fluxR` - the right fluxes for each face

Level: developer

-seealso: `PetscFV`, `PetscDS`, `PetscFVFaceGeom`, `PetscFVCreate()`

# External Links
$(_doc_external("DM/PetscFVIntegrateRHSFunction"))
"""
function PetscFVIntegrateRHSFunction(petsclib::PetscLibType, fvm::PetscFV, prob::PetscDS, field::PetscInt, Nf::PetscInt, fgeom::PetscFVFaceGeom, neighborVol::PetscReal, uL::Vector{PetscScalar}, uR::Vector{PetscScalar}) end

@for_petsc function PetscFVIntegrateRHSFunction(petsclib::$UnionPetscLib, fvm::PetscFV, prob::PetscDS, field::$PetscInt, Nf::$PetscInt, fgeom::PetscFVFaceGeom, neighborVol::$PetscReal, uL::Vector{$PetscScalar}, uR::Vector{$PetscScalar} )
	fluxL = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!
	fluxR = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFVIntegrateRHSFunction, $petsc_library),
               PetscErrorCode,
               (PetscFV, PetscDS, $PetscInt, $PetscInt, Ptr{PetscFVFaceGeom}, Ptr{$PetscReal}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               fvm, prob, field, Nf, fgeom, neighborVol, uL, uR, fluxL, fluxR,
              )


	return fluxL,fluxR
end 

"""
	PetscFVClone(petsclib::PetscLibType,fv::PetscFV, fvNew::PetscFV) 
Create a shallow copy of a `PetscFV` object that just references the internal objects.

Input Parameter:
- `fv` - The initial `PetscFV`

Output Parameter:
- `fvNew` - A clone of the `PetscFV`

Level: advanced

-seealso: `PetscFV`, `PetscFVType`, `PetscFVCreate()`, `PetscFVSetType()`

# External Links
$(_doc_external("DM/PetscFVClone"))
"""
function PetscFVClone(petsclib::PetscLibType, fv::PetscFV, fvNew::PetscFV) end

@for_petsc function PetscFVClone(petsclib::$UnionPetscLib, fv::PetscFV, fvNew::PetscFV )

    @chk ccall(
               (:PetscFVClone, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscFV}),
               fv, fvNew,
              )


	return nothing
end 

"""
	PetscFVRefine(petsclib::PetscLibType,fv::PetscFV, fvRef::PetscFV) 
Create a "refined" `PetscFV` object that refines the reference cell into
smaller copies.

Input Parameter:
- `fv` - The initial `PetscFV`

Output Parameter:
- `fvRef` - The refined `PetscFV`

Level: advanced

-seealso: `PetscFV`, `PetscFVType`, `PetscFVCreate()`, `PetscFVSetType()`

# External Links
$(_doc_external("DM/PetscFVRefine"))
"""
function PetscFVRefine(petsclib::PetscLibType, fv::PetscFV, fvRef::PetscFV) end

@for_petsc function PetscFVRefine(petsclib::$UnionPetscLib, fv::PetscFV, fvRef::PetscFV )

    @chk ccall(
               (:PetscFVRefine, $petsc_library),
               PetscErrorCode,
               (PetscFV, Ptr{PetscFV}),
               fv, fvRef,
              )


	return nothing
end 

"""
	PetscFVLeastSquaresSetMaxFaces(petsclib::PetscLibType,fvm::PetscFV, maxFaces::PetscInt) 
Set the maximum number of cell faces for gradient reconstruction

Not Collective

Input Parameters:
- `fvm`      - The `PetscFV` object
- `maxFaces` - The maximum number of cell faces

Level: intermediate

-seealso: `PetscFV`, `PetscFVCreate()`, `PETSCFVLEASTSQUARES`, `PetscFVComputeGradient()`

# External Links
$(_doc_external("DM/PetscFVLeastSquaresSetMaxFaces"))
"""
function PetscFVLeastSquaresSetMaxFaces(petsclib::PetscLibType, fvm::PetscFV, maxFaces::PetscInt) end

@for_petsc function PetscFVLeastSquaresSetMaxFaces(petsclib::$UnionPetscLib, fvm::PetscFV, maxFaces::$PetscInt )

    @chk ccall(
               (:PetscFVLeastSquaresSetMaxFaces, $petsc_library),
               PetscErrorCode,
               (PetscFV, $PetscInt),
               fvm, maxFaces,
              )


	return nothing
end 

