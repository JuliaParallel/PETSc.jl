mutable struct _n_PetscFE end
const PetscFE = Ptr{_n_PetscFE}

"""
	PetscFEFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `PetscFE` package. It is called
from `PetscFinalize()`.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Dm/PetscFEFinalizePackage"))
"""
function PetscFEFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscFEFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFEFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFEInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscFE` package. It is called
from `PetscDLLibraryRegister()` when using dynamic libraries, and on the first call to `PetscSpaceCreate()`
when using static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Dm/PetscFEInitializePackage"))
"""
function PetscFEInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscFEInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFEInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFERegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscFEType`

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: `PetscFE`, `PetscFEType`, `PetscFERegisterAll()`, `PetscFERegisterDestroy()`

# External Links
$(_doc_external("Dm/PetscFERegister"))
"""
function PetscFERegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscFERegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscFERegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscFESetType(petsclib::PetscLibType,fem::PetscFE, name::PetscFEType) 
Builds a particular `PetscFE`

Collective

Input Parameters:
- `fem`  - The `PetscFE` object
- `name` - The kind of FEM space

Options Database Key:
- `-petscfe_type <type>` - Sets the `PetscFE` type; use -help for a list of available types

Level: intermediate

-seealso: `PetscFEType`, `PetscFE`, `PetscFEGetType()`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFESetType"))
"""
function PetscFESetType(petsclib::PetscLibType, fem::PetscFE, name::PetscFEType) end

@for_petsc function PetscFESetType(petsclib::$UnionPetscLib, fem::PetscFE, name::PetscFEType )

    @chk ccall(
               (:PetscFESetType, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscFEType),
               fem, name,
              )


	return nothing
end 

"""
	name::PetscFEType = PetscFEGetType(petsclib::PetscLibType,fem::PetscFE) 
Gets the `PetscFEType` (as a string) from the `PetscFE` object.

Not Collective

Input Parameter:
- `fem` - The `PetscFE`

Output Parameter:
- `name` - The `PetscFEType` name

Level: intermediate

-seealso: `PetscFEType`, `PetscFE`, `PetscFESetType()`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetType"))
"""
function PetscFEGetType(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetType(petsclib::$UnionPetscLib, fem::PetscFE )
	name_ = Ref{PetscFEType}()

    @chk ccall(
               (:PetscFEGetType, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFEType}),
               fem, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscFEViewFromOptions(petsclib::PetscLibType,A::PetscFE, obj::PetscObject, name::String) 
View from a `PetscFE` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscFE` object
- `obj`  - Optional object that provides the options prefix, pass `NULL` to use the options prefix of `A`
- `name` - command line option name

Level: intermediate

-seealso: `PetscFE`, `PetscFEView()`, `PetscObjectViewFromOptions()`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEViewFromOptions"))
"""
function PetscFEViewFromOptions(petsclib::PetscLibType, A::PetscFE, obj::PetscObject, name::String) end

@for_petsc function PetscFEViewFromOptions(petsclib::$UnionPetscLib, A::PetscFE, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscFEViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscFEView(petsclib::PetscLibType,fem::PetscFE, viewer::PetscViewer) 
Views a `PetscFE`

Collective

Input Parameters:
- `fem`    - the `PetscFE` object to view
- `viewer` - the viewer

Level: beginner

-seealso: `PetscFE`, `PetscViewer`, `PetscFEDestroy()`, `PetscFEViewFromOptions()`

# External Links
$(_doc_external("Dm/PetscFEView"))
"""
function PetscFEView(petsclib::PetscLibType, fem::PetscFE, viewer::PetscViewer) end

@for_petsc function PetscFEView(petsclib::$UnionPetscLib, fem::PetscFE, viewer::PetscViewer )

    @chk ccall(
               (:PetscFEView, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscViewer),
               fem, viewer,
              )


	return nothing
end 

"""
	PetscFESetFromOptions(petsclib::PetscLibType,fem::PetscFE) 
sets parameters in a `PetscFE` from the options database

Collective

Input Parameter:
- `fem` - the `PetscFE` object to set options for

Options Database Keys:
- `-petscfe_num_blocks`  - the number of cell blocks to integrate concurrently
- `-petscfe_num_batches` - the number of cell batches to integrate serially

Level: intermediate

-seealso: `PetscFE`, `PetscFEView()`

# External Links
$(_doc_external("Dm/PetscFESetFromOptions"))
"""
function PetscFESetFromOptions(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFESetFromOptions(petsclib::$UnionPetscLib, fem::PetscFE )

    @chk ccall(
               (:PetscFESetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscFE,),
               fem,
              )


	return nothing
end 

"""
	PetscFESetUp(petsclib::PetscLibType,fem::PetscFE) 
Construct data structures for the `PetscFE` after the `PetscFEType` has been set

Collective

Input Parameter:
- `fem` - the `PetscFE` object to setup

Level: intermediate

-seealso: `PetscFE`, `PetscFEView()`, `PetscFEDestroy()`

# External Links
$(_doc_external("Dm/PetscFESetUp"))
"""
function PetscFESetUp(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFESetUp(petsclib::$UnionPetscLib, fem::PetscFE )

    @chk ccall(
               (:PetscFESetUp, $petsc_library),
               PetscErrorCode,
               (PetscFE,),
               fem,
              )


	return nothing
end 

"""
	PetscFEDestroy(petsclib::PetscLibType,fem::PetscFE) 
Destroys a `PetscFE` object

Collective

Input Parameter:
- `fem` - the `PetscFE` object to destroy

Level: beginner

-seealso: `PetscFE`, `PetscFEView()`

# External Links
$(_doc_external("Dm/PetscFEDestroy"))
"""
function PetscFEDestroy(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEDestroy(petsclib::$UnionPetscLib, fem::PetscFE )

    @chk ccall(
               (:PetscFEDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFE},),
               fem,
              )


	return nothing
end 

"""
	fem::PetscFE = PetscFECreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscFE` object. The type can then be set with `PetscFESetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscFE` object

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFE`, `PetscFEType`, `PetscFESetType()`, `PetscFECreateDefault()`, `PETSCFEGALERKIN`

# External Links
$(_doc_external("Dm/PetscFECreate"))
"""
function PetscFECreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscFECreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscFE}),
               comm, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	dim::PetscInt = PetscFEGetSpatialDimension(petsclib::PetscLibType,fem::PetscFE) 
Returns the spatial dimension of the element

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `dim` - The spatial dimension

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetSpatialDimension"))
"""
function PetscFEGetSpatialDimension(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetSpatialDimension(petsclib::$UnionPetscLib, fem::PetscFE )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFEGetSpatialDimension, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{$PetscInt}),
               fem, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	PetscFESetNumComponents(petsclib::PetscLibType,fem::PetscFE, comp::PetscInt) 
Sets the number of field components in the element

Not Collective

Input Parameters:
- `fem`  - The `PetscFE` object
- `comp` - The number of field components

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`, `PetscFEGetSpatialDimension()`, `PetscFEGetNumComponents()`

# External Links
$(_doc_external("Dm/PetscFESetNumComponents"))
"""
function PetscFESetNumComponents(petsclib::PetscLibType, fem::PetscFE, comp::PetscInt) end

@for_petsc function PetscFESetNumComponents(petsclib::$UnionPetscLib, fem::PetscFE, comp::$PetscInt )

    @chk ccall(
               (:PetscFESetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt),
               fem, comp,
              )


	return nothing
end 

"""
	comp::PetscInt = PetscFEGetNumComponents(petsclib::PetscLibType,fem::PetscFE) 
Returns the number of components in the element

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `comp` - The number of field components

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`, `PetscFEGetSpatialDimension()`

# External Links
$(_doc_external("Dm/PetscFEGetNumComponents"))
"""
function PetscFEGetNumComponents(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetNumComponents(petsclib::$UnionPetscLib, fem::PetscFE )
	comp_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFEGetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{$PetscInt}),
               fem, comp_,
              )

	comp = comp_[]

	return comp
end 

"""
	PetscFESetTileSizes(petsclib::PetscLibType,fem::PetscFE, blockSize::PetscInt, numBlocks::PetscInt, batchSize::PetscInt, numBatches::PetscInt) 
Sets the tile sizes for evaluation

Not Collective

Input Parameters:
- `fem`        - The `PetscFE` object
- `blockSize`  - The number of elements in a block
- `numBlocks`  - The number of blocks in a batch
- `batchSize`  - The number of elements in a batch
- `numBatches` - The number of batches in a chunk

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`, `PetscFEGetTileSizes()`

# External Links
$(_doc_external("Dm/PetscFESetTileSizes"))
"""
function PetscFESetTileSizes(petsclib::PetscLibType, fem::PetscFE, blockSize::PetscInt, numBlocks::PetscInt, batchSize::PetscInt, numBatches::PetscInt) end

@for_petsc function PetscFESetTileSizes(petsclib::$UnionPetscLib, fem::PetscFE, blockSize::$PetscInt, numBlocks::$PetscInt, batchSize::$PetscInt, numBatches::$PetscInt )

    @chk ccall(
               (:PetscFESetTileSizes, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
               fem, blockSize, numBlocks, batchSize, numBatches,
              )


	return nothing
end 

"""
	blockSize::PetscInt,numBlocks::PetscInt,batchSize::PetscInt,numBatches::PetscInt = PetscFEGetTileSizes(petsclib::PetscLibType,fem::PetscFE) 
Returns the tile sizes for evaluation

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameters:
- `blockSize`  - The number of elements in a block, pass `NULL` if not needed
- `numBlocks`  - The number of blocks in a batch, pass `NULL` if not needed
- `batchSize`  - The number of elements in a batch, pass `NULL` if not needed
- `numBatches` - The number of batches in a chunk, pass `NULL` if not needed

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`, `PetscFESetTileSizes()`

# External Links
$(_doc_external("Dm/PetscFEGetTileSizes"))
"""
function PetscFEGetTileSizes(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetTileSizes(petsclib::$UnionPetscLib, fem::PetscFE )
	blockSize_ = Ref{$PetscInt}()
	numBlocks_ = Ref{$PetscInt}()
	batchSize_ = Ref{$PetscInt}()
	numBatches_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFEGetTileSizes, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               fem, blockSize_, numBlocks_, batchSize_, numBatches_,
              )

	blockSize = blockSize_[]
	numBlocks = numBlocks_[]
	batchSize = batchSize_[]
	numBatches = numBatches_[]

	return blockSize,numBlocks,batchSize,numBatches
end 

"""
	PetscFEGetBasisSpace(petsclib::PetscLibType,fem::PetscFE, sp::PetscSpace) 
Returns the `PetscSpace` used for the approximation of the solution for the `PetscFE`

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `sp` - The `PetscSpace` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetBasisSpace"))
"""
function PetscFEGetBasisSpace(petsclib::PetscLibType, fem::PetscFE, sp::PetscSpace) end

@for_petsc function PetscFEGetBasisSpace(petsclib::$UnionPetscLib, fem::PetscFE, sp::PetscSpace )

    @chk ccall(
               (:PetscFEGetBasisSpace, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscSpace}),
               fem, sp,
              )


	return nothing
end 

"""
	PetscFESetBasisSpace(petsclib::PetscLibType,fem::PetscFE, sp::PetscSpace) 
Sets the `PetscSpace` used for the approximation of the solution

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `sp`  - The `PetscSpace` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscFECreate()`, `PetscFESetDualSpace()`

# External Links
$(_doc_external("Dm/PetscFESetBasisSpace"))
"""
function PetscFESetBasisSpace(petsclib::PetscLibType, fem::PetscFE, sp::PetscSpace) end

@for_petsc function PetscFESetBasisSpace(petsclib::$UnionPetscLib, fem::PetscFE, sp::PetscSpace )

    @chk ccall(
               (:PetscFESetBasisSpace, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscSpace),
               fem, sp,
              )


	return nothing
end 

"""
	PetscFEGetDualSpace(petsclib::PetscLibType,fem::PetscFE, sp::PetscDualSpace) 
Returns the `PetscDualSpace` used to define the inner product for a `PetscFE`

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `sp` - The `PetscDualSpace` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetDualSpace"))
"""
function PetscFEGetDualSpace(petsclib::PetscLibType, fem::PetscFE, sp::PetscDualSpace) end

@for_petsc function PetscFEGetDualSpace(petsclib::$UnionPetscLib, fem::PetscFE, sp::PetscDualSpace )

    @chk ccall(
               (:PetscFEGetDualSpace, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscDualSpace}),
               fem, sp,
              )


	return nothing
end 

"""
	PetscFESetDualSpace(petsclib::PetscLibType,fem::PetscFE, sp::PetscDualSpace) 
Sets the `PetscDualSpace` used to define the inner product

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `sp`  - The `PetscDualSpace` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscFECreate()`, `PetscFESetBasisSpace()`

# External Links
$(_doc_external("Dm/PetscFESetDualSpace"))
"""
function PetscFESetDualSpace(petsclib::PetscLibType, fem::PetscFE, sp::PetscDualSpace) end

@for_petsc function PetscFESetDualSpace(petsclib::$UnionPetscLib, fem::PetscFE, sp::PetscDualSpace )

    @chk ccall(
               (:PetscFESetDualSpace, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscDualSpace),
               fem, sp,
              )


	return nothing
end 

"""
	PetscFEGetQuadrature(petsclib::PetscLibType,fem::PetscFE, q::PetscQuadrature) 
Returns the `PetscQuadrature` used to calculate inner products

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `q` - The `PetscQuadrature` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetQuadrature"))
"""
function PetscFEGetQuadrature(petsclib::PetscLibType, fem::PetscFE, q::PetscQuadrature) end

@for_petsc function PetscFEGetQuadrature(petsclib::$UnionPetscLib, fem::PetscFE, q::PetscQuadrature )

    @chk ccall(
               (:PetscFEGetQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscQuadrature}),
               fem, q,
              )


	return nothing
end 

"""
	PetscFESetQuadrature(petsclib::PetscLibType,fem::PetscFE, q::PetscQuadrature) 
Sets the `PetscQuadrature` used to calculate inner products

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `q`   - The `PetscQuadrature` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`, `PetscFECreate()`, `PetscFEGetFaceQuadrature()`

# External Links
$(_doc_external("Dm/PetscFESetQuadrature"))
"""
function PetscFESetQuadrature(petsclib::PetscLibType, fem::PetscFE, q::PetscQuadrature) end

@for_petsc function PetscFESetQuadrature(petsclib::$UnionPetscLib, fem::PetscFE, q::PetscQuadrature )

    @chk ccall(
               (:PetscFESetQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscQuadrature),
               fem, q,
              )


	return nothing
end 

"""
	PetscFEGetFaceQuadrature(petsclib::PetscLibType,fem::PetscFE, q::PetscQuadrature) 
Returns the `PetscQuadrature` used to calculate inner products on faces

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `q` - The `PetscQuadrature` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`, `PetscFECreate()`, `PetscFESetQuadrature()`, `PetscFESetFaceQuadrature()`

# External Links
$(_doc_external("Dm/PetscFEGetFaceQuadrature"))
"""
function PetscFEGetFaceQuadrature(petsclib::PetscLibType, fem::PetscFE, q::PetscQuadrature) end

@for_petsc function PetscFEGetFaceQuadrature(petsclib::$UnionPetscLib, fem::PetscFE, q::PetscQuadrature )

    @chk ccall(
               (:PetscFEGetFaceQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscQuadrature}),
               fem, q,
              )


	return nothing
end 

"""
	PetscFESetFaceQuadrature(petsclib::PetscLibType,fem::PetscFE, q::PetscQuadrature) 
Sets the `PetscQuadrature` used to calculate inner products on faces

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `q`   - The `PetscQuadrature` object

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`, `PetscFECreate()`, `PetscFESetQuadrature()`

# External Links
$(_doc_external("Dm/PetscFESetFaceQuadrature"))
"""
function PetscFESetFaceQuadrature(petsclib::PetscLibType, fem::PetscFE, q::PetscQuadrature) end

@for_petsc function PetscFESetFaceQuadrature(petsclib::$UnionPetscLib, fem::PetscFE, q::PetscQuadrature )

    @chk ccall(
               (:PetscFESetFaceQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscQuadrature),
               fem, q,
              )


	return nothing
end 

"""
	PetscFECopyQuadrature(petsclib::PetscLibType,sfe::PetscFE, tfe::PetscFE) 
Copy both volumetric and surface quadrature to a new `PetscFE`

Not Collective

Input Parameters:
- `sfe` - The `PetscFE` source for the quadratures
- `tfe` - The `PetscFE` target for the quadratures

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`, `PetscFECreate()`, `PetscFESetQuadrature()`, `PetscFESetFaceQuadrature()`

# External Links
$(_doc_external("Dm/PetscFECopyQuadrature"))
"""
function PetscFECopyQuadrature(petsclib::PetscLibType, sfe::PetscFE, tfe::PetscFE) end

@for_petsc function PetscFECopyQuadrature(petsclib::$UnionPetscLib, sfe::PetscFE, tfe::PetscFE )

    @chk ccall(
               (:PetscFECopyQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscFE),
               sfe, tfe,
              )


	return nothing
end 

"""
	numDof::Vector{PetscInt} = PetscFEGetNumDof(petsclib::PetscLibType,fem::PetscFE) 
Returns the number of dofs (dual basis vectors) associated to mesh points on the reference cell of a given dimension

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `numDof` - Array of length `dim` with the number of dofs in each dimension

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFEGetNumDof"))
"""
function PetscFEGetNumDof(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetNumDof(petsclib::$UnionPetscLib, fem::PetscFE )
	numDof_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscFEGetNumDof, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{Ptr{$PetscInt}}),
               fem, numDof_,
              )

	numDof = unsafe_wrap(Array, numDof_[], VecGetLocalSize(petsclib, x); own = false)

	return numDof
end 

"""
	PetscFEGetCellTabulation(petsclib::PetscLibType,fem::PetscFE, k::PetscInt, T::PetscTabulation) 
Returns the tabulation of the basis functions at the quadrature points on the reference cell

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `k`   - The highest derivative we need to tabulate, very often 1

Output Parameter:
- `T` - The basis function values and derivatives at quadrature points

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscTabulation`, `PetscFECreateTabulation()`, `PetscTabulationDestroy()`

# External Links
$(_doc_external("Dm/PetscFEGetCellTabulation"))
"""
function PetscFEGetCellTabulation(petsclib::PetscLibType, fem::PetscFE, k::PetscInt, T::PetscTabulation) end

@for_petsc function PetscFEGetCellTabulation(petsclib::$UnionPetscLib, fem::PetscFE, k::$PetscInt, T::PetscTabulation )

    @chk ccall(
               (:PetscFEGetCellTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, Ptr{PetscTabulation}),
               fem, k, T,
              )


	return nothing
end 

"""
	PetscFEExpandFaceQuadrature(petsclib::PetscLibType,fe::PetscFE, fq::PetscQuadrature, efq::PetscQuadrature) 

# External Links
$(_doc_external("Dm/PetscFEExpandFaceQuadrature"))
"""
function PetscFEExpandFaceQuadrature(petsclib::PetscLibType, fe::PetscFE, fq::PetscQuadrature, efq::PetscQuadrature) end

@for_petsc function PetscFEExpandFaceQuadrature(petsclib::$UnionPetscLib, fe::PetscFE, fq::PetscQuadrature, efq::PetscQuadrature )

    @chk ccall(
               (:PetscFEExpandFaceQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscQuadrature, Ptr{PetscQuadrature}),
               fe, fq, efq,
              )


	return nothing
end 

"""
	PetscFEGetFaceTabulation(petsclib::PetscLibType,fem::PetscFE, k::PetscInt, Tf::PetscTabulation) 
Returns the tabulation of the basis functions at the face quadrature points for each face of the reference cell

Not Collective

Input Parameters:
- `fem` - The `PetscFE` object
- `k`   - The highest derivative we need to tabulate, very often 1

Output Parameter:
- `Tf` - The basis function values and derivatives at face quadrature points

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscTabulation`, `PetscFEGetCellTabulation()`, `PetscFECreateTabulation()`, `PetscTabulationDestroy()`

# External Links
$(_doc_external("Dm/PetscFEGetFaceTabulation"))
"""
function PetscFEGetFaceTabulation(petsclib::PetscLibType, fem::PetscFE, k::PetscInt, Tf::PetscTabulation) end

@for_petsc function PetscFEGetFaceTabulation(petsclib::$UnionPetscLib, fem::PetscFE, k::$PetscInt, Tf::PetscTabulation )

    @chk ccall(
               (:PetscFEGetFaceTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, Ptr{PetscTabulation}),
               fem, k, Tf,
              )


	return nothing
end 

"""
	PetscFEGetFaceCentroidTabulation(petsclib::PetscLibType,fem::PetscFE, Tc::PetscTabulation) 
Returns the tabulation of the basis functions at the face centroid points

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameter:
- `Tc` - The basis function values at face centroid points

Level: intermediate

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscTabulation`, `PetscFEGetFaceTabulation()`, `PetscFEGetCellTabulation()`, `PetscFECreateTabulation()`, `PetscTabulationDestroy()`

# External Links
$(_doc_external("Dm/PetscFEGetFaceCentroidTabulation"))
"""
function PetscFEGetFaceCentroidTabulation(petsclib::PetscLibType, fem::PetscFE, Tc::PetscTabulation) end

@for_petsc function PetscFEGetFaceCentroidTabulation(petsclib::$UnionPetscLib, fem::PetscFE, Tc::PetscTabulation )

    @chk ccall(
               (:PetscFEGetFaceCentroidTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscTabulation}),
               fem, Tc,
              )


	return nothing
end 

"""
	T::PetscTabulation = PetscFECreateTabulation(petsclib::PetscLibType,fem::PetscFE, nrepl::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt) 
Tabulates the basis functions, and perhaps derivatives, at the points provided.

Not Collective

Input Parameters:
- `fem`     - The `PetscFE` object
- `nrepl`   - The number of replicas
- `npoints` - The number of tabulation points in a replica
- `points`  - The tabulation point coordinates
- `K`       - The number of derivatives calculated

Output Parameter:
- `T` - The basis function values and derivatives at tabulation points

Level: intermediate

-seealso: `PetscTabulation`, `PetscFEGetCellTabulation()`, `PetscTabulationDestroy()`

# External Links
$(_doc_external("Dm/PetscFECreateTabulation"))
"""
function PetscFECreateTabulation(petsclib::PetscLibType, fem::PetscFE, nrepl::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt) end

@for_petsc function PetscFECreateTabulation(petsclib::$UnionPetscLib, fem::PetscFE, nrepl::$PetscInt, npoints::$PetscInt, points::Vector{$PetscReal}, K::$PetscInt )
	T_ = Ref{PetscTabulation}()

    @chk ccall(
               (:PetscFECreateTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{PetscTabulation}),
               fem, nrepl, npoints, points, K, T_,
              )

	T = T_[]

	return T
end 

"""
	PetscFEComputeTabulation(petsclib::PetscLibType,fem::PetscFE, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt, T::PetscTabulation) 
Tabulates the basis functions, and perhaps derivatives, at the points provided.

Not Collective

Input Parameters:
- `fem`     - The `PetscFE` object
- `npoints` - The number of tabulation points
- `points`  - The tabulation point coordinates
- `K`       - The number of derivatives calculated
- `T`       - An existing tabulation object with enough allocated space

Output Parameter:
- `T` - The basis function values and derivatives at tabulation points

Level: intermediate

-seealso: `PetscTabulation`, `PetscFEGetCellTabulation()`, `PetscTabulationDestroy()`

# External Links
$(_doc_external("Dm/PetscFEComputeTabulation"))
"""
function PetscFEComputeTabulation(petsclib::PetscLibType, fem::PetscFE, npoints::PetscInt, points::Vector{PetscReal}, K::PetscInt, T::PetscTabulation) end

@for_petsc function PetscFEComputeTabulation(petsclib::$UnionPetscLib, fem::PetscFE, npoints::$PetscInt, points::Vector{$PetscReal}, K::$PetscInt, T::PetscTabulation )

    @chk ccall(
               (:PetscFEComputeTabulation, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, Ptr{$PetscReal}, $PetscInt, PetscTabulation),
               fem, npoints, points, K, T,
              )


	return nothing
end 

"""
	trFE::PetscFE = PetscFECreateHeightTrace(petsclib::PetscLibType,fe::PetscFE, height::PetscInt) 

# External Links
$(_doc_external("Dm/PetscFECreateHeightTrace"))
"""
function PetscFECreateHeightTrace(petsclib::PetscLibType, fe::PetscFE, height::PetscInt) end

@for_petsc function PetscFECreateHeightTrace(petsclib::$UnionPetscLib, fe::PetscFE, height::$PetscInt )
	trFE_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateHeightTrace, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, Ptr{PetscFE}),
               fe, height, trFE_,
              )

	trFE = trFE_[]

	return trFE
end 

"""
	dim::PetscInt = PetscFEGetDimension(petsclib::PetscLibType,fem::PetscFE) 
Get the dimension of the finite element space on a cell

Not Collective

Input Parameter:
- `fem` - The `PetscFE`

Output Parameter:
- `dim` - The dimension

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`, `PetscSpaceGetDimension()`, `PetscDualSpaceGetDimension()`

# External Links
$(_doc_external("Dm/PetscFEGetDimension"))
"""
function PetscFEGetDimension(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEGetDimension(petsclib::$UnionPetscLib, fem::PetscFE )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFEGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{$PetscInt}),
               fem, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	PetscFEPushforward(petsclib::PetscLibType,fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) 
Map the reference element function to real space

Input Parameters:
- `fe`     - The `PetscFE`
- `fegeom` - The cell geometry
- `Nv`     - The number of function values
- `vals`   - The function values

Output Parameter:
- `vals` - The transformed function values

Level: advanced

-seealso: `PetscFE`, `PetscFEGeom`, `PetscDualSpace`, `PetscDualSpacePushforward()`

# External Links
$(_doc_external("Dm/PetscFEPushforward"))
"""
function PetscFEPushforward(petsclib::PetscLibType, fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscFEPushforward(petsclib::$UnionPetscLib, fe::PetscFE, fegeom::PetscFEGeom, Nv::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscFEPushforward, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFEGeom}, $PetscInt, Ptr{$PetscScalar}),
               fe, fegeom, Nv, vals,
              )


	return nothing
end 

"""
	PetscFEPushforwardGradient(petsclib::PetscLibType,fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) 
Map the reference element function gradient to real space

Input Parameters:
- `fe`     - The `PetscFE`
- `fegeom` - The cell geometry
- `Nv`     - The number of function gradient values
- `vals`   - The function gradient values

Output Parameter:
- `vals` - The transformed function gradient values

Level: advanced

-seealso: `PetscFE`, `PetscFEGeom`, `PetscDualSpace`, `PetscFEPushforward()`, `PetscDualSpacePushforwardGradient()`, `PetscDualSpacePushforward()`

# External Links
$(_doc_external("Dm/PetscFEPushforwardGradient"))
"""
function PetscFEPushforwardGradient(petsclib::PetscLibType, fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscFEPushforwardGradient(petsclib::$UnionPetscLib, fe::PetscFE, fegeom::PetscFEGeom, Nv::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscFEPushforwardGradient, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFEGeom}, $PetscInt, Ptr{$PetscScalar}),
               fe, fegeom, Nv, vals,
              )


	return nothing
end 

"""
	PetscFEPushforwardHessian(petsclib::PetscLibType,fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) 
Map the reference element function Hessian to real space

Input Parameters:
- `fe`     - The `PetscFE`
- `fegeom` - The cell geometry
- `Nv`     - The number of function Hessian values
- `vals`   - The function Hessian values

Output Parameter:
- `vals` - The transformed function Hessian values

Level: advanced

-seealso: `PetscFE`, `PetscFEGeom`, `PetscDualSpace`, `PetscFEPushforward()`, `PetscDualSpacePushforwardHessian()`, `PetscDualSpacePushforward()`

# External Links
$(_doc_external("Dm/PetscFEPushforwardHessian"))
"""
function PetscFEPushforwardHessian(petsclib::PetscLibType, fe::PetscFE, fegeom::PetscFEGeom, Nv::PetscInt, vals::Vector{PetscScalar}) end

@for_petsc function PetscFEPushforwardHessian(petsclib::$UnionPetscLib, fe::PetscFE, fegeom::PetscFEGeom, Nv::$PetscInt, vals::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscFEPushforwardHessian, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFEGeom}, $PetscInt, Ptr{$PetscScalar}),
               fe, fegeom, Nv, vals,
              )


	return nothing
end 

"""
	integral::Vector{PetscScalar} = PetscFEIntegrate(petsclib::PetscLibType,prob::PetscDS, field::PetscInt, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}) 
Produce the integral for the given field for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `prob`            - The `PetscDS` specifying the discretizations and continuum functions
- `field`           - The field being integrated
- `Ne`              - The number of elements in the chunk
- `cgeom`           - The cell geometry for each cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements

Output Parameter:
- `integral` - the integral for this field

Level: intermediate

-seealso: `PetscFE`, `PetscDS`, `PetscFEIntegrateResidual()`, `PetscFEIntegrateBd()`

# External Links
$(_doc_external("Dm/PetscFEIntegrate"))
"""
function PetscFEIntegrate(petsclib::PetscLibType, prob::PetscDS, field::PetscInt, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}) end

@for_petsc function PetscFEIntegrate(petsclib::$UnionPetscLib, prob::PetscDS, field::$PetscInt, Ne::$PetscInt, cgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar} )
	integral = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrate, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               prob, field, Ne, cgeom, coefficients, probAux, coefficientsAux, integral,
              )


	return integral
end 

#=
"""
	PetscFEIntegrateBd(petsclib::PetscLibType,prob::PetscDS, field::PetscInt, ::Cvoid(obj_func) 
Produce the integral for the given field for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `prob`            - The `PetscDS` specifying the discretizations and continuum functions
- `field`           - The field being integrated
- `obj_func`        - The function to be integrated
- `Ne`              - The number of elements in the chunk
- `geom`            - The face geometry for each face in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements

Output Parameter:
- `integral` - the integral for this field

Level: intermediate

-seealso: `PetscFE`, `PetscDS`, `PetscFEIntegrateResidual()`, `PetscFEIntegrate()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateBd"))
"""
function PetscFEIntegrateBd(petsclib::PetscLibType, prob::PetscDS, field::PetscInt, ::Cvoid(obj_func) end

@for_petsc function PetscFEIntegrateBd(petsclib::$UnionPetscLib, prob::PetscDS, field::$PetscInt, ::Cvoid(obj_func )

    @chk ccall(
               (:PetscFEIntegrateBd, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid(obj_func}),
               prob, field, g0,
              )


	return nothing
end 
=#

"""
	elemVec::Vector{PetscScalar} = PetscFEIntegrateResidual(petsclib::PetscLibType,ds::PetscDS, key::PetscFormKey, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) 
Produce the element residual vector for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `ds`              - The `PetscDS` specifying the discretizations and continuum functions
- `key`             - The (label+value, field) being integrated
- `Ne`              - The number of elements in the chunk
- `cgeom`           - The cell geometry for each cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time

Output Parameter:
- `elemVec` - the element residual vectors from each element

Level: intermediate

-seealso: `PetscFEIntegrateBdResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateResidual"))
"""
function PetscFEIntegrateResidual(petsclib::PetscLibType, ds::PetscDS, key::PetscFormKey, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) end

@for_petsc function PetscFEIntegrateResidual(petsclib::$UnionPetscLib, ds::PetscDS, key::PetscFormKey, Ne::$PetscInt, cgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal )
	elemVec = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscFormKey, $PetscInt, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, Ptr{$PetscScalar}),
               ds, key, Ne, cgeom, coefficients, coefficients_t, probAux, coefficientsAux, t, elemVec,
              )


	return elemVec
end 

"""
	elemVec::Vector{PetscScalar} = PetscFEIntegrateBdResidual(petsclib::PetscLibType,ds::PetscDS, wf::PetscWeakForm, key::PetscFormKey, Ne::PetscInt, fgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) 
Produce the element residual vector for a chunk of elements by quadrature integration over a boundary

Not Collective

Input Parameters:
- `ds`              - The `PetscDS` specifying the discretizations and continuum functions
- `wf`              - The PetscWeakForm object holding the pointwise functions
- `key`             - The (label+value, field) being integrated
- `Ne`              - The number of elements in the chunk
- `fgeom`           - The face geometry for each cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time

Output Parameter:
- `elemVec` - the element residual vectors from each element

Level: intermediate

-seealso: `PetscFEIntegrateResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateBdResidual"))
"""
function PetscFEIntegrateBdResidual(petsclib::PetscLibType, ds::PetscDS, wf::PetscWeakForm, key::PetscFormKey, Ne::PetscInt, fgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) end

@for_petsc function PetscFEIntegrateBdResidual(petsclib::$UnionPetscLib, ds::PetscDS, wf::PetscWeakForm, key::PetscFormKey, Ne::$PetscInt, fgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal )
	elemVec = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscWeakForm, PetscFormKey, $PetscInt, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, Ptr{$PetscScalar}),
               ds, wf, key, Ne, fgeom, coefficients, coefficients_t, probAux, coefficientsAux, t, elemVec,
              )


	return elemVec
end 

"""
	elemVec::Vector{PetscScalar} = PetscFEIntegrateHybridResidual(petsclib::PetscLibType,ds::PetscDS, dsIn::PetscDS, key::PetscFormKey, s::PetscInt, Ne::PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) 
Produce the element residual vector for a chunk of hybrid element faces by quadrature integration

Not Collective

Input Parameters:
- `ds`              - The `PetscDS` specifying the discretizations and continuum functions
- `dsIn`            - The `PetscDS` specifying the discretizations and continuum functions for input
- `key`             - The (label+value, field) being integrated
- `s`               - The side of the cell being integrated, 0 for negative and 1 for positive
- `Ne`              - The number of elements in the chunk
- `fgeom`           - The face geometry for each cell in the chunk
- `cgeom`           - The cell geometry for each neighbor cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time

Output Parameter:
- `elemVec` - the element residual vectors from each element

Level: developer

-seealso: `PetscFEIntegrateResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateHybridResidual"))
"""
function PetscFEIntegrateHybridResidual(petsclib::PetscLibType, ds::PetscDS, dsIn::PetscDS, key::PetscFormKey, s::PetscInt, Ne::PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal) end

@for_petsc function PetscFEIntegrateHybridResidual(petsclib::$UnionPetscLib, ds::PetscDS, dsIn::PetscDS, key::PetscFormKey, s::$PetscInt, Ne::$PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal )
	elemVec = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateHybridResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS, PetscFormKey, $PetscInt, $PetscInt, Ptr{PetscFEGeom}, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, Ptr{$PetscScalar}),
               ds, dsIn, key, s, Ne, fgeom, cgeom, coefficients, coefficients_t, probAux, coefficientsAux, t, elemVec,
              )


	return elemVec
end 

"""
	elemMat::Vector{PetscScalar} = PetscFEIntegrateJacobian(petsclib::PetscLibType,rds::PetscDS, cds::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, dsAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) 
Produce the element Jacobian for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `rds`             - The `PetscDS` specifying the row discretizations and continuum functions
- `cds`             - The `PetscDS` specifying the column discretizations
- `jtype`           - The type of matrix pointwise functions that should be used
- `key`             - The (label+value, fieldI*Nf + fieldJ) being integrated
- `Ne`              - The number of elements in the chunk
- `cgeom`           - The cell geometry for each cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements for the Jacobian evaluation point
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `dsAux`           - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time
- `u_tshift`        - A multiplier for the dF/du_t term (as opposed to the dF/du term)

Output Parameter:
- `elemMat` - the element matrices for the Jacobian from each element

Level: intermediate

-seealso: `PetscFEIntegrateResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateJacobian"))
"""
function PetscFEIntegrateJacobian(petsclib::PetscLibType, rds::PetscDS, cds::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::PetscInt, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, dsAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) end

@for_petsc function PetscFEIntegrateJacobian(petsclib::$UnionPetscLib, rds::PetscDS, cds::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::$PetscInt, cgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, dsAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal, u_tshift::$PetscReal )
	elemMat = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS, PetscFEJacobianType, PetscFormKey, $PetscInt, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, $PetscReal, Ptr{$PetscScalar}),
               rds, cds, jtype, key, Ne, cgeom, coefficients, coefficients_t, dsAux, coefficientsAux, t, u_tshift, elemMat,
              )


	return elemMat
end 

"""
	elemMat::Vector{PetscScalar} = PetscFEIntegrateBdJacobian(petsclib::PetscLibType,ds::PetscDS, wf::PetscWeakForm, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::PetscInt, fgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) 
Produce the boundary element Jacobian for a chunk of elements by quadrature integration

Not Collective

Input Parameters:
- `ds`              - The `PetscDS` specifying the discretizations and continuum functions
- `wf`              - The PetscWeakForm holding the pointwise functions
- `jtype`           - The type of matrix pointwise functions that should be used
- `key`             - The (label+value, fieldI*Nf + fieldJ) being integrated
- `Ne`              - The number of elements in the chunk
- `fgeom`           - The face geometry for each cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements for the Jacobian evaluation point
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time
- `u_tshift`        - A multiplier for the dF/du_t term (as opposed to the dF/du term)

Output Parameter:
- `elemMat` - the element matrices for the Jacobian from each element

Level: intermediate

-seealso: `PetscFEIntegrateJacobian()`, `PetscFEIntegrateResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateBdJacobian"))
"""
function PetscFEIntegrateBdJacobian(petsclib::PetscLibType, ds::PetscDS, wf::PetscWeakForm, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::PetscInt, fgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) end

@for_petsc function PetscFEIntegrateBdJacobian(petsclib::$UnionPetscLib, ds::PetscDS, wf::PetscWeakForm, jtype::PetscFEJacobianType, key::PetscFormKey, Ne::$PetscInt, fgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal, u_tshift::$PetscReal )
	elemMat = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscWeakForm, PetscFEJacobianType, PetscFormKey, $PetscInt, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, $PetscReal, Ptr{$PetscScalar}),
               ds, wf, jtype, key, Ne, fgeom, coefficients, coefficients_t, probAux, coefficientsAux, t, u_tshift, elemMat,
              )


	return elemMat
end 

"""
	elemMat::Vector{PetscScalar} = PetscFEIntegrateHybridJacobian(petsclib::PetscLibType,ds::PetscDS, dsIn::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, s::PetscInt, Ne::PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) 
Produce the boundary element Jacobian for a chunk of hybrid elements by quadrature integration

Not Collective

Input Parameters:
- `ds`              - The `PetscDS` specifying the discretizations and continuum functions for the output
- `dsIn`            - The `PetscDS` specifying the discretizations and continuum functions for the input
- `jtype`           - The type of matrix pointwise functions that should be used
- `key`             - The (label+value, fieldI*Nf + fieldJ) being integrated
- `s`               - The side of the cell being integrated, 0 for negative and 1 for positive
- `Ne`              - The number of elements in the chunk
- `fgeom`           - The face geometry for each cell in the chunk
- `cgeom`           - The cell geometry for each neighbor cell in the chunk
- `coefficients`    - The array of FEM basis coefficients for the elements for the Jacobian evaluation point
- `coefficients_t`  - The array of FEM basis time derivative coefficients for the elements
- `probAux`         - The `PetscDS` specifying the auxiliary discretizations
- `coefficientsAux` - The array of FEM auxiliary basis coefficients for the elements
- `t`               - The time
- `u_tshift`        - A multiplier for the dF/du_t term (as opposed to the dF/du term)

Output Parameter:
- `elemMat` - the element matrices for the Jacobian from each element

Level: developer

-seealso: `PetscFEIntegrateJacobian()`, `PetscFEIntegrateResidual()`

# External Links
$(_doc_external("Dm/PetscFEIntegrateHybridJacobian"))
"""
function PetscFEIntegrateHybridJacobian(petsclib::PetscLibType, ds::PetscDS, dsIn::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, s::PetscInt, Ne::PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{PetscScalar}, coefficients_t::Vector{PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{PetscScalar}, t::PetscReal, u_tshift::PetscReal) end

@for_petsc function PetscFEIntegrateHybridJacobian(petsclib::$UnionPetscLib, ds::PetscDS, dsIn::PetscDS, jtype::PetscFEJacobianType, key::PetscFormKey, s::$PetscInt, Ne::$PetscInt, fgeom::PetscFEGeom, cgeom::PetscFEGeom, coefficients::Vector{$PetscScalar}, coefficients_t::Vector{$PetscScalar}, probAux::PetscDS, coefficientsAux::Vector{$PetscScalar}, t::$PetscReal, u_tshift::$PetscReal )
	elemMat = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscFEIntegrateHybridJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS, PetscFEJacobianType, PetscFormKey, $PetscInt, $PetscInt, Ptr{PetscFEGeom}, Ptr{PetscFEGeom}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, PetscDS, Ptr{$PetscScalar}, $PetscReal, $PetscReal, Ptr{$PetscScalar}),
               ds, dsIn, jtype, key, s, Ne, fgeom, cgeom, coefficients, coefficients_t, probAux, coefficientsAux, t, u_tshift, elemMat,
              )


	return elemMat
end 

"""
	PetscFEGetHeightSubspace(petsclib::PetscLibType,fe::PetscFE, height::PetscInt, subfe::PetscFE) 
Get the subspace of this space for a mesh point of a given height

Input Parameters:
- `fe`     - The finite element space
- `height` - The height of the `DMPLEX` point

Output Parameter:
- `subfe` - The subspace of this `PetscFE` space

Level: advanced

-seealso: `PetscFECreateDefault()`

# External Links
$(_doc_external("Dm/PetscFEGetHeightSubspace"))
"""
function PetscFEGetHeightSubspace(petsclib::PetscLibType, fe::PetscFE, height::PetscInt, subfe::PetscFE) end

@for_petsc function PetscFEGetHeightSubspace(petsclib::$UnionPetscLib, fe::PetscFE, height::$PetscInt, subfe::PetscFE )

    @chk ccall(
               (:PetscFEGetHeightSubspace, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, Ptr{PetscFE}),
               fe, height, subfe,
              )


	return nothing
end 

"""
	PetscFERefine(petsclib::PetscLibType,fe::PetscFE, feRef::PetscFE) 
Create a "refined" `PetscFE` object that refines the reference cell into
smaller copies.

Collective

Input Parameter:
- `fe` - The initial `PetscFE`

Output Parameter:
- `feRef` - The refined `PetscFE`

Level: advanced

-seealso: `PetscFEType`, `PetscFECreate()`, `PetscFESetType()`

# External Links
$(_doc_external("Dm/PetscFERefine"))
"""
function PetscFERefine(petsclib::PetscLibType, fe::PetscFE, feRef::PetscFE) end

@for_petsc function PetscFERefine(petsclib::$UnionPetscLib, fe::PetscFE, feRef::PetscFE )

    @chk ccall(
               (:PetscFERefine, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFE}),
               fe, feRef,
              )


	return nothing
end 

"""
	fem::PetscFE = PetscFECreateFromSpaces(petsclib::PetscLibType,P::PetscSpace, Q::PetscDualSpace, q::PetscQuadrature, fq::PetscQuadrature) 
Create a `PetscFE` from the basis and dual spaces

Collective

Input Parameters:
- `P`  - The basis space
- `Q`  - The dual space
- `q`  - The cell quadrature
- `fq` - The face quadrature

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFE`, `PetscSpace`, `PetscDualSpace`, `PetscQuadrature`,
`PetscFECreateLagrangeByCell()`, `PetscFECreateDefault()`, `PetscFECreateByCell()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFECreateFromSpaces"))
"""
function PetscFECreateFromSpaces(petsclib::PetscLibType, P::PetscSpace, Q::PetscDualSpace, q::PetscQuadrature, fq::PetscQuadrature) end

@for_petsc function PetscFECreateFromSpaces(petsclib::$UnionPetscLib, P::PetscSpace, Q::PetscDualSpace, q::PetscQuadrature, fq::PetscQuadrature )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateFromSpaces, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscDualSpace, PetscQuadrature, PetscQuadrature, Ptr{PetscFE}),
               P, Q, q, fq, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	fem::PetscFE = PetscFECreateDefault(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, isSimplex::PetscBool, prefix::String, qorder::PetscInt) 
Create a `PetscFE` for basic FEM computation

Collective

Input Parameters:
- `comm`      - The MPI comm
- `dim`       - The spatial dimension
- `Nc`        - The number of components
- `isSimplex` - Flag for simplex reference cell, otherwise its a tensor product
- `prefix`    - The options prefix, or `NULL`
- `qorder`    - The quadrature order or `PETSC_DETERMINE` to use `PetscSpace` polynomial degree

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFECreateLagrange()`, `PetscFECreateByCell()`, `PetscSpaceSetFromOptions()`, `PetscDualSpaceSetFromOptions()`, `PetscFESetFromOptions()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFECreateDefault"))
"""
function PetscFECreateDefault(petsclib::PetscLibType, comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, isSimplex::PetscBool, prefix::String, qorder::PetscInt) end

@for_petsc function PetscFECreateDefault(petsclib::$UnionPetscLib, comm::MPI_Comm, dim::$PetscInt, Nc::$PetscInt, isSimplex::PetscBool, prefix::String, qorder::$PetscInt )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateDefault, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, PetscBool, Ptr{Cchar}, $PetscInt, Ptr{PetscFE}),
               comm, dim, Nc, isSimplex, prefix, qorder, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	fem::PetscFE = PetscFECreateByCell(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, ct::DMPolytopeType, prefix::String, qorder::PetscInt) 
Create a `PetscFE` for basic FEM computation

Collective

Input Parameters:
- `comm`   - The MPI comm
- `dim`    - The spatial dimension
- `Nc`     - The number of components
- `ct`     - The celltype of the reference cell
- `prefix` - The options prefix, or `NULL`
- `qorder` - The quadrature order or `PETSC_DETERMINE` to use `PetscSpace` polynomial degree

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFECreateDefault()`, `PetscFECreateLagrange()`, `PetscSpaceSetFromOptions()`, `PetscDualSpaceSetFromOptions()`, `PetscFESetFromOptions()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFECreateByCell"))
"""
function PetscFECreateByCell(petsclib::PetscLibType, comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, ct::DMPolytopeType, prefix::String, qorder::PetscInt) end

@for_petsc function PetscFECreateByCell(petsclib::$UnionPetscLib, comm::MPI_Comm, dim::$PetscInt, Nc::$PetscInt, ct::DMPolytopeType, prefix::String, qorder::$PetscInt )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateByCell, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, DMPolytopeType, Ptr{Cchar}, $PetscInt, Ptr{PetscFE}),
               comm, dim, Nc, ct, prefix, qorder, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	fem::PetscFE = PetscFECreateLagrange(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, isSimplex::PetscBool, k::PetscInt, qorder::PetscInt) 
Create a `PetscFE` for the basic Lagrange space of degree k

Collective

Input Parameters:
- `comm`      - The MPI comm
- `dim`       - The spatial dimension
- `Nc`        - The number of components
- `isSimplex` - Flag for simplex reference cell, otherwise its a tensor product
- `k`         - The degree k of the space
- `qorder`    - The quadrature order or `PETSC_DETERMINE` to use `PetscSpace` polynomial degree

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFECreateLagrangeByCell()`, `PetscFECreateDefault()`, `PetscFECreateByCell()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFECreateLagrange"))
"""
function PetscFECreateLagrange(petsclib::PetscLibType, comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, isSimplex::PetscBool, k::PetscInt, qorder::PetscInt) end

@for_petsc function PetscFECreateLagrange(petsclib::$UnionPetscLib, comm::MPI_Comm, dim::$PetscInt, Nc::$PetscInt, isSimplex::PetscBool, k::$PetscInt, qorder::$PetscInt )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateLagrange, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, PetscBool, $PetscInt, $PetscInt, Ptr{PetscFE}),
               comm, dim, Nc, isSimplex, k, qorder, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	fem::PetscFE = PetscFECreateLagrangeByCell(petsclib::PetscLibType,comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, ct::DMPolytopeType, k::PetscInt, qorder::PetscInt) 
Create a `PetscFE` for the basic Lagrange space of degree k

Collective

Input Parameters:
- `comm`   - The MPI comm
- `dim`    - The spatial dimension
- `Nc`     - The number of components
- `ct`     - The celltype of the reference cell
- `k`      - The degree k of the space
- `qorder` - The quadrature order or `PETSC_DETERMINE` to use `PetscSpace` polynomial degree

Output Parameter:
- `fem` - The `PetscFE` object

Level: beginner

-seealso: `PetscFECreateLagrange()`, `PetscFECreateDefault()`, `PetscFECreateByCell()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFECreateLagrangeByCell"))
"""
function PetscFECreateLagrangeByCell(petsclib::PetscLibType, comm::MPI_Comm, dim::PetscInt, Nc::PetscInt, ct::DMPolytopeType, k::PetscInt, qorder::PetscInt) end

@for_petsc function PetscFECreateLagrangeByCell(petsclib::$UnionPetscLib, comm::MPI_Comm, dim::$PetscInt, Nc::$PetscInt, ct::DMPolytopeType, k::$PetscInt, qorder::$PetscInt )
	fem_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateLagrangeByCell, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, DMPolytopeType, $PetscInt, $PetscInt, Ptr{PetscFE}),
               comm, dim, Nc, ct, k, qorder, fem_,
              )

	fem = fem_[]

	return fem
end 

"""
	PetscFELimitDegree(petsclib::PetscLibType,fe::PetscFE, minDegree::PetscInt, maxDegree::PetscInt, newfe::PetscFE) 
Copy a `PetscFE` but limit the degree to be in the given range

Collective

Input Parameters:
- `fe`        - The `PetscFE`
- `minDegree` - The minimum degree, or `PETSC_DETERMINE` for no limit
- `maxDegree` - The maximum degree, or `PETSC_DETERMINE` for no limit

Output Parameter:
- `newfe` - The `PetscFE` object

Level: advanced

-seealso: `PetscFECreateLagrange()`, `PetscFECreateDefault()`, `PetscFECreateByCell()`, `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFELimitDegree"))
"""
function PetscFELimitDegree(petsclib::PetscLibType, fe::PetscFE, minDegree::PetscInt, maxDegree::PetscInt, newfe::PetscFE) end

@for_petsc function PetscFELimitDegree(petsclib::$UnionPetscLib, fe::PetscFE, minDegree::$PetscInt, maxDegree::$PetscInt, newfe::PetscFE )

    @chk ccall(
               (:PetscFELimitDegree, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, $PetscInt, Ptr{PetscFE}),
               fe, minDegree, maxDegree, newfe,
              )


	return nothing
end 

"""
	dgfe::PetscFE = PetscFECreateBrokenElement(petsclib::PetscLibType,cgfe::PetscFE) 
Create a discontinuous version of the input `PetscFE`

Collective

Input Parameters:
- `cgfe` - The continuous `PetscFE` object

Output Parameter:
- `dgfe` - The discontinuous `PetscFE` object

Level: advanced

-seealso: `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`, `PetscFECreateLagrange()`, `PetscFECreateLagrangeByCell()`, `PetscDualSpaceLagrangeSetContinuity()`

# External Links
$(_doc_external("Dm/PetscFECreateBrokenElement"))
"""
function PetscFECreateBrokenElement(petsclib::PetscLibType, cgfe::PetscFE) end

@for_petsc function PetscFECreateBrokenElement(petsclib::$UnionPetscLib, cgfe::PetscFE )
	dgfe_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateBrokenElement, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFE}),
               cgfe, dgfe_,
              )

	dgfe = dgfe_[]

	return dgfe
end 

"""
	PetscFESetName(petsclib::PetscLibType,fe::PetscFE, name::String) 
Names the `PetscFE` and its subobjects

Not Collective

Input Parameters:
- `fe`   - The `PetscFE`
- `name` - The name

Level: intermediate

-seealso: `PetscFECreate()`, `PetscSpaceCreate()`, `PetscDualSpaceCreate()`

# External Links
$(_doc_external("Dm/PetscFESetName"))
"""
function PetscFESetName(petsclib::PetscLibType, fe::PetscFE, name::String) end

@for_petsc function PetscFESetName(petsclib::$UnionPetscLib, fe::PetscFE, name::String )

    @chk ccall(
               (:PetscFESetName, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{Cchar}),
               fe, name,
              )


	return nothing
end 

"""
	cgeom::PetscFEGeom = PetscFECreateCellGeometry(petsclib::PetscLibType,fe::PetscFE, quad::PetscQuadrature) 

# External Links
$(_doc_external("Dm/PetscFECreateCellGeometry"))
"""
function PetscFECreateCellGeometry(petsclib::PetscLibType, fe::PetscFE, quad::PetscQuadrature) end

@for_petsc function PetscFECreateCellGeometry(petsclib::$UnionPetscLib, fe::PetscFE, quad::PetscQuadrature )
	cgeom_ = Ref{PetscFEGeom}()

    @chk ccall(
               (:PetscFECreateCellGeometry, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscQuadrature, Ptr{PetscFEGeom}),
               fe, quad, cgeom_,
              )

	cgeom = cgeom_[]

	return cgeom
end 

"""
	PetscFEDestroyCellGeometry(petsclib::PetscLibType,fe::PetscFE, cgeom::PetscFEGeom) 

# External Links
$(_doc_external("Dm/PetscFEDestroyCellGeometry"))
"""
function PetscFEDestroyCellGeometry(petsclib::PetscLibType, fe::PetscFE, cgeom::PetscFEGeom) end

@for_petsc function PetscFEDestroyCellGeometry(petsclib::$UnionPetscLib, fe::PetscFE, cgeom::PetscFEGeom )

    @chk ccall(
               (:PetscFEDestroyCellGeometry, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscFEGeom}),
               fe, cgeom,
              )


	return nothing
end 

"""
	geom::PetscFEGeom = PetscFEGeomCreate(petsclib::PetscLibType,quad::PetscQuadrature, numCells::PetscInt, dimEmbed::PetscInt, mode::PetscFEGeomMode) 
Create a `PetscFEGeom` object to manage geometry for a group of cells

Input Parameters:
- `quad`     - A `PetscQuadrature` determining the tabulation
- `numCells` - The number of cells in the group
- `dimEmbed` - The coordinate dimension
- `mode`     - Type of geometry data to store

Output Parameter:
- `geom` - The `PetscFEGeom` object, which is a struct not a `PetscObject`

Level: beginner

-seealso: `PetscFEGeom`, `PetscQuadrature`, `PetscFEGeomDestroy()`, `PetscFEGeomComplete()`

# External Links
$(_doc_external("Dm/PetscFEGeomCreate"))
"""
function PetscFEGeomCreate(petsclib::PetscLibType, quad::PetscQuadrature, numCells::PetscInt, dimEmbed::PetscInt, mode::PetscFEGeomMode) end

@for_petsc function PetscFEGeomCreate(petsclib::$UnionPetscLib, quad::PetscQuadrature, numCells::$PetscInt, dimEmbed::$PetscInt, mode::PetscFEGeomMode )
	geom_ = Ref{PetscFEGeom}()

    @chk ccall(
               (:PetscFEGeomCreate, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt, $PetscInt, PetscFEGeomMode, PetscFEGeom),
               quad, numCells, dimEmbed, mode, geom_,
              )

	geom = geom_[]

	return geom
end 

"""
	PetscFEGeomDestroy(petsclib::PetscLibType,geom::PetscFEGeom) 
Destroy a `PetscFEGeom` object

Input Parameter:
- `geom` - `PetscFEGeom` object

Level: beginner

-seealso: `PetscFEGeom`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomDestroy"))
"""
function PetscFEGeomDestroy(petsclib::PetscLibType, geom::PetscFEGeom) end

@for_petsc function PetscFEGeomDestroy(petsclib::$UnionPetscLib, geom::PetscFEGeom )

    @chk ccall(
               (:PetscFEGeomDestroy, $petsc_library),
               PetscErrorCode,
               (PetscFEGeom,),
               geom,
              )


	return nothing
end 

"""
	PetscFEGeomGetChunk(petsclib::PetscLibType,geom::PetscFEGeom, cStart::PetscInt, cEnd::PetscInt, chunkGeom::Vector{PetscFEGeom}) 
Get a chunk of cells in the group as a `PetscFEGeom`

Input Parameters:
- `geom`   - `PetscFEGeom` object
- `cStart` - The first cell in the chunk
- `cEnd`   - The first cell not in the chunk

Output Parameter:
- `chunkGeom` - an array of cells of length `cEnd` - `cStart`

Level: intermediate

-seealso: `PetscFEGeom`, `PetscFEGeomRestoreChunk()`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomGetChunk"))
"""
function PetscFEGeomGetChunk(petsclib::PetscLibType, geom::PetscFEGeom, cStart::PetscInt, cEnd::PetscInt, chunkGeom::Vector{PetscFEGeom}) end

@for_petsc function PetscFEGeomGetChunk(petsclib::$UnionPetscLib, geom::PetscFEGeom, cStart::$PetscInt, cEnd::$PetscInt, chunkGeom::Vector{PetscFEGeom} )
	chunkGeom_ = Ref(pointer(chunkGeom))

    @chk ccall(
               (:PetscFEGeomGetChunk, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{Ptr{PetscFEGeom}}),
               geom, cStart, cEnd, chunkGeom_,
              )


	return nothing
end 

"""
	PetscFEGeomRestoreChunk(petsclib::PetscLibType,geom::PetscFEGeom, cStart::PetscInt, cEnd::PetscInt, chunkGeom::PetscFEGeom) 
Restore the chunk obtained with `PetscFEGeomCreateChunk()`

Input Parameters:
- `geom`      - `PetscFEGeom` object
- `cStart`    - The first cell in the chunk
- `cEnd`      - The first cell not in the chunk
- `chunkGeom` - The chunk of cells

Level: intermediate

-seealso: `PetscFEGeom`, `PetscFEGeomGetChunk()`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomRestoreChunk"))
"""
function PetscFEGeomRestoreChunk(petsclib::PetscLibType, geom::PetscFEGeom, cStart::PetscInt, cEnd::PetscInt, chunkGeom::PetscFEGeom) end

@for_petsc function PetscFEGeomRestoreChunk(petsclib::$UnionPetscLib, geom::PetscFEGeom, cStart::$PetscInt, cEnd::$PetscInt, chunkGeom::PetscFEGeom )

    @chk ccall(
               (:PetscFEGeomRestoreChunk, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFEGeom}, $PetscInt, $PetscInt, PetscFEGeom),
               geom, cStart, cEnd, chunkGeom,
              )


	return nothing
end 

"""
	PetscFEGeomGetPoint(petsclib::PetscLibType,geom::PetscFEGeom, c::PetscInt, p::PetscInt, pcoords::Vector{PetscReal}, pgeom::PetscFEGeom) 
Get the geometry for cell `c` at point `p` as a `PetscFEGeom`

Input Parameters:
- `geom`    - `PetscFEGeom` object
- `c`       - The cell
- `p`       - The point
- `pcoords` - The reference coordinates of point `p`, or `NULL`

Output Parameter:
- `pgeom` - The geometry of cell `c` at point `p`

Level: intermediate

-seealso: `PetscFEGeom`, `PetscFEGeomRestoreChunk()`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomGetPoint"))
"""
function PetscFEGeomGetPoint(petsclib::PetscLibType, geom::PetscFEGeom, c::PetscInt, p::PetscInt, pcoords::Vector{PetscReal}, pgeom::PetscFEGeom) end

@for_petsc function PetscFEGeomGetPoint(petsclib::$UnionPetscLib, geom::PetscFEGeom, c::$PetscInt, p::$PetscInt, pcoords::Vector{$PetscReal}, pgeom::PetscFEGeom )

    @chk ccall(
               (:PetscFEGeomGetPoint, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{PetscFEGeom}),
               geom, c, p, pcoords, pgeom,
              )


	return nothing
end 

"""
	PetscFEGeomGetCellPoint(petsclib::PetscLibType,geom::PetscFEGeom, c::PetscInt, p::PetscInt, pgeom::PetscFEGeom) 
Get the cell geometry for cell `c` at point `p` as a `PetscFEGeom`

Input Parameters:
- `geom` - `PetscFEGeom` object
- `c`    - The cell
- `p`    - The point

Output Parameter:
- `pgeom` - The cell geometry of cell `c` at point `p`

Level: intermediate

-seealso: `PetscFEGeom`, `PetscFEGeomMode`, `PetscFEGeomRestoreChunk()`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomGetCellPoint"))
"""
function PetscFEGeomGetCellPoint(petsclib::PetscLibType, geom::PetscFEGeom, c::PetscInt, p::PetscInt, pgeom::PetscFEGeom) end

@for_petsc function PetscFEGeomGetCellPoint(petsclib::$UnionPetscLib, geom::PetscFEGeom, c::$PetscInt, p::$PetscInt, pgeom::PetscFEGeom )

    @chk ccall(
               (:PetscFEGeomGetCellPoint, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFEGeom}, $PetscInt, $PetscInt, Ptr{PetscFEGeom}),
               geom, c, p, pgeom,
              )


	return nothing
end 

"""
	PetscFEGeomComplete(petsclib::PetscLibType,geom::PetscFEGeom) 
Calculate derived quantities from a base geometry specification

Input Parameter:
- `geom` - `PetscFEGeom` object

Level: intermediate

-seealso: `PetscFEGeom`, `PetscFEGeomCreate()`

# External Links
$(_doc_external("Dm/PetscFEGeomComplete"))
"""
function PetscFEGeomComplete(petsclib::PetscLibType, geom::PetscFEGeom) end

@for_petsc function PetscFEGeomComplete(petsclib::$UnionPetscLib, geom::PetscFEGeom )

    @chk ccall(
               (:PetscFEGeomComplete, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFEGeom},),
               geom,
              )


	return nothing
end 

"""
	numSubelements::PetscInt,v0::Vector{PetscReal},jac::Vector{PetscReal},invjac::Vector{PetscReal} = PetscFECompositeGetMapping(petsclib::PetscLibType,fem::PetscFE) 
Returns the mappings from the reference element to each subelement

Not Collective

Input Parameter:
- `fem` - The `PetscFE` object

Output Parameters:
- `numSubelements` - The number of sub elements
- `v0`             - The affine transformation for each element, an array of length dim * Nc. Pass `NULL` to ignore.
- `jac`            - The Jacobian for each element, an array of length dim^2 * Nc. Pass `NULL` to ignore.
- `invjac`         - The inverse of the Jacobian, an array of length dim^2 * Nc. Pass `NULL` to ignore.

Level: intermediate

-seealso: `PetscFE`, `PetscFECreate()`

# External Links
$(_doc_external("Dm/PetscFECompositeGetMapping"))
"""
function PetscFECompositeGetMapping(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFECompositeGetMapping(petsclib::$UnionPetscLib, fem::PetscFE )
	numSubelements_ = Ref{$PetscInt}()
	v0_ = Ref{Ptr{$PetscReal}}()
	jac_ = Ref{Ptr{$PetscReal}}()
	invjac_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:PetscFECompositeGetMapping, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               fem, numSubelements_, v0_, jac_, invjac_,
              )

	numSubelements = numSubelements_[]
	v0 = unsafe_wrap(Array, v0_[], VecGetLocalSize(petsclib, x); own = false)
	jac = unsafe_wrap(Array, jac_[], VecGetLocalSize(petsclib, x); own = false)
	invjac = unsafe_wrap(Array, invjac_[], VecGetLocalSize(petsclib, x); own = false)

	return numSubelements,v0,jac,invjac
end 

"""
	PetscFEOpenCLSetRealType(petsclib::PetscLibType,fem::PetscFE, realType::PetscDataType) 
Set the scalar type for running on the OpenCL accelerator

Input Parameters:
- `fem`      - The `PetscFE`
- `realType` - The scalar type

Level: developer

-seealso: `PetscFE`, `PetscFEOpenCLGetRealType()`

# External Links
$(_doc_external("Dm/PetscFEOpenCLSetRealType"))
"""
function PetscFEOpenCLSetRealType(petsclib::PetscLibType, fem::PetscFE, realType::PetscDataType) end

@for_petsc function PetscFEOpenCLSetRealType(petsclib::$UnionPetscLib, fem::PetscFE, realType::PetscDataType )

    @chk ccall(
               (:PetscFEOpenCLSetRealType, $petsc_library),
               PetscErrorCode,
               (PetscFE, PetscDataType),
               fem, realType,
              )


	return nothing
end 

"""
	realType::PetscDataType = PetscFEOpenCLGetRealType(petsclib::PetscLibType,fem::PetscFE) 
Get the scalar type for running on the OpenCL accelerator

Input Parameter:
- `fem` - The `PetscFE`

Output Parameter:
- `realType` - The scalar type

Level: developer

-seealso: `PetscFE`, `PetscFEOpenCLSetRealType()`

# External Links
$(_doc_external("Dm/PetscFEOpenCLGetRealType"))
"""
function PetscFEOpenCLGetRealType(petsclib::PetscLibType, fem::PetscFE) end

@for_petsc function PetscFEOpenCLGetRealType(petsclib::$UnionPetscLib, fem::PetscFE )
	realType_ = Ref{PetscDataType}()

    @chk ccall(
               (:PetscFEOpenCLGetRealType, $petsc_library),
               PetscErrorCode,
               (PetscFE, Ptr{PetscDataType}),
               fem, realType_,
              )

	realType = unsafe_string(realType_[])

	return realType
end 

"""
	vector_fe::PetscFE = PetscFECreateVector(petsclib::PetscLibType,scalar_fe::PetscFE, num_copies::PetscInt, interleave_basis::PetscBool, interleave_components::PetscBool) 
Create a vector
`PetscFE`.

Collective

Input Parameters:
- `scalar_fe`             - a `PetscFE` finite element
- `num_copies`            - a positive integer
- `interleave_basis`      - if `PETSC_TRUE`, the first `num_copies` basis vectors
of the output finite element will be copies of the
first basis vector of `scalar_fe`, and so on for the
other basis vectors; otherwise all of the first-copy
basis vectors will come first, followed by all of the
second-copy, and so on.
- `interleave_components` - if `PETSC_TRUE`, the first `num_copies` components
of the output finite element will be copies of the
first component of `scalar_fe`, and so on for the
other components; otherwise all of the first-copy
components will come first, followed by all of the
second-copy, and so on.

Output Parameter:
- `vector_fe` - a `PetscFE` of type `PETSCFEVECTOR` that represent a discretization space with `num_copies` copies of `scalar_fe`

Level: intermediate

-seealso: `PetscFE`, `PetscFEType`, `PetscFECreate()`, `PetscFESetType()`, `PETSCFEBASIC`, `PETSCFEVECTOR`

# External Links
$(_doc_external("Dm/PetscFECreateVector"))
"""
function PetscFECreateVector(petsclib::PetscLibType, scalar_fe::PetscFE, num_copies::PetscInt, interleave_basis::PetscBool, interleave_components::PetscBool) end

@for_petsc function PetscFECreateVector(petsclib::$UnionPetscLib, scalar_fe::PetscFE, num_copies::$PetscInt, interleave_basis::PetscBool, interleave_components::PetscBool )
	vector_fe_ = Ref{PetscFE}()

    @chk ccall(
               (:PetscFECreateVector, $petsc_library),
               PetscErrorCode,
               (PetscFE, $PetscInt, PetscBool, PetscBool, Ptr{PetscFE}),
               scalar_fe, num_copies, interleave_basis, interleave_components, vector_fe_,
              )

	vector_fe = vector_fe_[]

	return vector_fe
end 

"""
	q::PetscQuadrature = PetscQuadratureCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a `PetscQuadrature` object

Collective

Input Parameter:
- `comm` - The communicator for the `PetscQuadrature` object

Output Parameter:
- `q` - The `PetscQuadrature` object

Level: beginner

-seealso: `PetscQuadrature`, `Petscquadraturedestroy()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureCreate"))
"""
function PetscQuadratureCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscQuadratureCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	q_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscQuadratureCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscQuadrature}),
               comm, q_,
              )

	q = q_[]

	return q
end 

"""
	r::PetscQuadrature = PetscQuadratureDuplicate(petsclib::PetscLibType,q::PetscQuadrature) 
Create a deep copy of the `PetscQuadrature` object

Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Output Parameter:
- `r` - The new `PetscQuadrature` object

Level: beginner

-seealso: `PetscQuadrature`, `PetscQuadratureCreate()`, `PetscQuadratureDestroy()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureDuplicate"))
"""
function PetscQuadratureDuplicate(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureDuplicate(petsclib::$UnionPetscLib, q::PetscQuadrature )
	r_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscQuadratureDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{PetscQuadrature}),
               q, r_,
              )

	r = r_[]

	return r
end 

"""
	PetscQuadratureDestroy(petsclib::PetscLibType,q::PetscQuadrature) 
Destroys a `PetscQuadrature` object

Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Level: beginner

-seealso: `PetscQuadrature`, `PetscQuadratureCreate()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureDestroy"))
"""
function PetscQuadratureDestroy(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureDestroy(petsclib::$UnionPetscLib, q::PetscQuadrature )

    @chk ccall(
               (:PetscQuadratureDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscQuadrature},),
               q,
              )


	return nothing
end 

"""
	ct::DMPolytopeType = PetscQuadratureGetCellType(petsclib::PetscLibType,q::PetscQuadrature) 
Return the cell type of the integration domain

Not Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Output Parameter:
- `ct` - The cell type of the integration domain

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureSetCellType()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureGetCellType"))
"""
function PetscQuadratureGetCellType(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureGetCellType(petsclib::$UnionPetscLib, q::PetscQuadrature )
	ct_ = Ref{DMPolytopeType}()

    @chk ccall(
               (:PetscQuadratureGetCellType, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{DMPolytopeType}),
               q, ct_,
              )

	ct = unsafe_string(ct_[])

	return ct
end 

"""
	PetscQuadratureSetCellType(petsclib::PetscLibType,q::PetscQuadrature, ct::DMPolytopeType) 
Set the cell type of the integration domain

Not Collective

Input Parameters:
- `q`  - The `PetscQuadrature` object
- `ct` - The cell type of the integration domain

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureGetCellType()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureSetCellType"))
"""
function PetscQuadratureSetCellType(petsclib::PetscLibType, q::PetscQuadrature, ct::DMPolytopeType) end

@for_petsc function PetscQuadratureSetCellType(petsclib::$UnionPetscLib, q::PetscQuadrature, ct::DMPolytopeType )

    @chk ccall(
               (:PetscQuadratureSetCellType, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, DMPolytopeType),
               q, ct,
              )


	return nothing
end 

"""
	order::PetscInt = PetscQuadratureGetOrder(petsclib::PetscLibType,q::PetscQuadrature) 
Return the order of the method in the `PetscQuadrature`

Not Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Output Parameter:
- `order` - The order of the quadrature, i.e. the highest degree polynomial that is exactly integrated

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureSetOrder()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureGetOrder"))
"""
function PetscQuadratureGetOrder(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureGetOrder(petsclib::$UnionPetscLib, q::PetscQuadrature )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscQuadratureGetOrder, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{$PetscInt}),
               q, order_,
              )

	order = order_[]

	return order
end 

"""
	PetscQuadratureSetOrder(petsclib::PetscLibType,q::PetscQuadrature, order::PetscInt) 
Set the order of the method in the `PetscQuadrature`

Not Collective

Input Parameters:
- `q`     - The `PetscQuadrature` object
- `order` - The order of the quadrature, i.e. the highest degree polynomial that is exactly integrated

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureGetOrder()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureSetOrder"))
"""
function PetscQuadratureSetOrder(petsclib::PetscLibType, q::PetscQuadrature, order::PetscInt) end

@for_petsc function PetscQuadratureSetOrder(petsclib::$UnionPetscLib, q::PetscQuadrature, order::$PetscInt )

    @chk ccall(
               (:PetscQuadratureSetOrder, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt),
               q, order,
              )


	return nothing
end 

"""
	Nc::PetscInt = PetscQuadratureGetNumComponents(petsclib::PetscLibType,q::PetscQuadrature) 
Return the number of components for functions to be integrated

Not Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Output Parameter:
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureSetNumComponents()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureGetNumComponents"))
"""
function PetscQuadratureGetNumComponents(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureGetNumComponents(petsclib::$UnionPetscLib, q::PetscQuadrature )
	Nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscQuadratureGetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{$PetscInt}),
               q, Nc_,
              )

	Nc = Nc_[]

	return Nc
end 

"""
	PetscQuadratureSetNumComponents(petsclib::PetscLibType,q::PetscQuadrature, Nc::PetscInt) 
Sets the number of components for functions to be integrated

Not Collective

Input Parameters:
- `q`  - The `PetscQuadrature` object
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureGetNumComponents()`, `PetscQuadratureGetData()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureSetNumComponents"))
"""
function PetscQuadratureSetNumComponents(petsclib::PetscLibType, q::PetscQuadrature, Nc::PetscInt) end

@for_petsc function PetscQuadratureSetNumComponents(petsclib::$UnionPetscLib, q::PetscQuadrature, Nc::$PetscInt )

    @chk ccall(
               (:PetscQuadratureSetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt),
               q, Nc,
              )


	return nothing
end 

"""
	dim::PetscInt,Nc::PetscInt,npoints::PetscInt,points::Vector{PetscReal},weights::Vector{PetscReal} = PetscQuadratureGetData(petsclib::PetscLibType,q::PetscQuadrature) 
Returns the data defining the `PetscQuadrature`

Not Collective

Input Parameter:
- `q` - The `PetscQuadrature` object

Output Parameters:
- `dim`     - The spatial dimension
- `Nc`      - The number of components
- `npoints` - The number of quadrature points
- `points`  - The coordinates of each quadrature point
- `weights` - The weight of each quadrature point

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureCreate()`, `PetscQuadratureSetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureGetData"))
"""
function PetscQuadratureGetData(petsclib::PetscLibType, q::PetscQuadrature) end

@for_petsc function PetscQuadratureGetData(petsclib::$UnionPetscLib, q::PetscQuadrature )
	dim_ = Ref{$PetscInt}()
	Nc_ = Ref{$PetscInt}()
	npoints_ = Ref{$PetscInt}()
	points_ = Ref{Ptr{$PetscReal}}()
	weights_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:PetscQuadratureGetData, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               q, dim_, Nc_, npoints_, points_, weights_,
              )

	dim = dim_[]
	Nc = Nc_[]
	npoints = npoints_[]
	points = unsafe_wrap(Array, points_[], VecGetLocalSize(petsclib, x); own = false)
	weights = unsafe_wrap(Array, weights_[], VecGetLocalSize(petsclib, x); own = false)

	return dim,Nc,npoints,points,weights
end 

"""
	equal::PetscBool = PetscQuadratureEqual(petsclib::PetscLibType,A::PetscQuadrature, B::PetscQuadrature) 
determine whether two quadratures are equivalent

Input Parameters:
- `A` - A `PetscQuadrature` object
- `B` - Another `PetscQuadrature` object

Output Parameter:
- `equal` - `PETSC_TRUE` if the quadratures are the same

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureCreate()`

# External Links
$(_doc_external("Dm/PetscQuadratureEqual"))
"""
function PetscQuadratureEqual(petsclib::PetscLibType, A::PetscQuadrature, B::PetscQuadrature) end

@for_petsc function PetscQuadratureEqual(petsclib::$UnionPetscLib, A::PetscQuadrature, B::PetscQuadrature )
	equal_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscQuadratureEqual, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, PetscQuadrature, Ptr{PetscBool}),
               A, B, equal_,
              )

	equal = equal_[]

	return equal
end 

"""
	PetscQuadraturePushForward(petsclib::PetscLibType,q::PetscQuadrature, imageDim::PetscInt, origin::Vector{PetscReal}, originImage::Vector{PetscReal}, J::Vector{PetscReal}, formDegree::PetscInt, Jinvstarq::PetscQuadrature) 
Push forward a quadrature functional under an affine transformation.

Collective

Input Parameters:
- `q`           - the quadrature functional
- `imageDim`    - the dimension of the image of the transformation
- `origin`      - a point in the original space
- `originImage` - the image of the origin under the transformation
- `J`           - the Jacobian of the image: an [imageDim x dim] matrix in row major order
- `formDegree`  - transform the quadrature weights as k-forms of this form degree (if the number of components is a multiple of (dim choose `formDegree`),
it is assumed that they represent multiple k-forms) [see `PetscDTAltVPullback()` for interpretation of `formDegree`]

Output Parameter:
- `Jinvstarq` - a quadrature rule where each point is the image of a point in the original quadrature rule, and where the k-form weights have
been pulled-back by the pseudoinverse of `J` to the k-form weights in the image space.

Level: intermediate

-seealso: `PetscQuadrature`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("Dm/PetscQuadraturePushForward"))
"""
function PetscQuadraturePushForward(petsclib::PetscLibType, q::PetscQuadrature, imageDim::PetscInt, origin::Vector{PetscReal}, originImage::Vector{PetscReal}, J::Vector{PetscReal}, formDegree::PetscInt, Jinvstarq::PetscQuadrature) end

@for_petsc function PetscQuadraturePushForward(petsclib::$UnionPetscLib, q::PetscQuadrature, imageDim::$PetscInt, origin::Vector{$PetscReal}, originImage::Vector{$PetscReal}, J::Vector{$PetscReal}, formDegree::$PetscInt, Jinvstarq::PetscQuadrature )

    @chk ccall(
               (:PetscQuadraturePushForward, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{PetscQuadrature}),
               q, imageDim, origin, originImage, J, formDegree, Jinvstarq,
              )


	return nothing
end 

"""
	PetscQuadratureSetData(petsclib::PetscLibType,q::PetscQuadrature, dim::PetscInt, Nc::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, weights::Vector{PetscReal}) 
Sets the data defining the quadrature

Not Collective

Input Parameters:
- `q`       - The `PetscQuadrature` object
- `dim`     - The spatial dimension
- `Nc`      - The number of components
- `npoints` - The number of quadrature points
- `points`  - The coordinates of each quadrature point
- `weights` - The weight of each quadrature point

Level: intermediate

-seealso: `PetscQuadrature`, `PetscQuadratureCreate()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureSetData"))
"""
function PetscQuadratureSetData(petsclib::PetscLibType, q::PetscQuadrature, dim::PetscInt, Nc::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, weights::Vector{PetscReal}) end

@for_petsc function PetscQuadratureSetData(petsclib::$UnionPetscLib, q::PetscQuadrature, dim::$PetscInt, Nc::$PetscInt, npoints::$PetscInt, points::Vector{$PetscReal}, weights::Vector{$PetscReal} )

    @chk ccall(
               (:PetscQuadratureSetData, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               q, dim, Nc, npoints, points, weights,
              )


	return nothing
end 

"""
	PetscQuadratureView(petsclib::PetscLibType,quad::PetscQuadrature, viewer::PetscViewer) 
View a `PetscQuadrature` object

Collective

Input Parameters:
- `quad`   - The `PetscQuadrature` object
- `viewer` - The `PetscViewer` object

Level: beginner

-seealso: `PetscQuadrature`, `PetscViewer`, `PetscQuadratureCreate()`, `PetscQuadratureGetData()`

# External Links
$(_doc_external("Dm/PetscQuadratureView"))
"""
function PetscQuadratureView(petsclib::PetscLibType, quad::PetscQuadrature, viewer::PetscViewer) end

@for_petsc function PetscQuadratureView(petsclib::$UnionPetscLib, quad::PetscQuadrature, viewer::PetscViewer )

    @chk ccall(
               (:PetscQuadratureView, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, PetscViewer),
               quad, viewer,
              )


	return nothing
end 

"""
	PetscQuadratureExpandComposite(petsclib::PetscLibType,q::PetscQuadrature, numSubelements::PetscInt, v0::Vector{PetscReal}, jac::Vector{PetscReal}, qref::PetscQuadrature) 
Return a quadrature over the composite element, which has the original quadrature in each subelement

Not Collective; No Fortran Support

Input Parameters:
- `q`              - The original `PetscQuadrature`
- `numSubelements` - The number of subelements the original element is divided into
- `v0`             - An array of the initial points for each subelement
- `jac`            - An array of the Jacobian mappings from the reference to each subelement

Output Parameter:
- `qref` - The dimension

Level: intermediate

-seealso: `PetscQuadrature`, `PetscFECreate()`, `PetscSpaceGetDimension()`, `PetscDualSpaceGetDimension()`

# External Links
$(_doc_external("Dm/PetscQuadratureExpandComposite"))
"""
function PetscQuadratureExpandComposite(petsclib::PetscLibType, q::PetscQuadrature, numSubelements::PetscInt, v0::Vector{PetscReal}, jac::Vector{PetscReal}, qref::PetscQuadrature) end

@for_petsc function PetscQuadratureExpandComposite(petsclib::$UnionPetscLib, q::PetscQuadrature, numSubelements::$PetscInt, v0::Vector{$PetscReal}, jac::Vector{$PetscReal}, qref::PetscQuadrature )

    @chk ccall(
               (:PetscQuadratureExpandComposite, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{PetscQuadrature}),
               q, numSubelements, v0, jac, qref,
              )


	return nothing
end 

"""
	Np::PetscInt = PetscQuadratureComputePermutations(petsclib::PetscLibType,quad::PetscQuadrature, perm::Vector{IS}) 
Compute permutations of quadrature points corresponding to domain orientations

Input Parameter:
- `quad` - The `PetscQuadrature`

Output Parameters:
- `Np`   - The number of domain orientations
- `perm` - An array of `IS` permutations, one for ech orientation,

Level: developer

-seealso: `PetscQuadratureSetCellType()`, `PetscQuadrature`

# External Links
$(_doc_external("Dm/PetscQuadratureComputePermutations"))
"""
function PetscQuadratureComputePermutations(petsclib::PetscLibType, quad::PetscQuadrature, perm::Vector{IS}) end

@for_petsc function PetscQuadratureComputePermutations(petsclib::$UnionPetscLib, quad::PetscQuadrature, perm::Vector{IS} )
	Np_ = Ref{$PetscInt}()
	perm_ = Ref(pointer(perm))

    @chk ccall(
               (:PetscQuadratureComputePermutations, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
               quad, Np_, perm_,
              )

	Np = Np_[]

	return Np
end 

"""
	PetscTabulationDestroy(petsclib::PetscLibType,T::PetscTabulation) 
Frees memory from the associated tabulation.

Not Collective

Input Parameter:
- `T` - The tabulation

Level: intermediate

-seealso: `PetscTabulation`, `PetscFECreateTabulation()`, `PetscFEGetCellTabulation()`

# External Links
$(_doc_external("Dm/PetscTabulationDestroy"))
"""
function PetscTabulationDestroy(petsclib::PetscLibType, T::PetscTabulation) end

@for_petsc function PetscTabulationDestroy(petsclib::$UnionPetscLib, T::PetscTabulation )

    @chk ccall(
               (:PetscTabulationDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscTabulation},),
               T,
              )


	return nothing
end 

"""
	PetscWeakFormCopy(petsclib::PetscLibType,wf::PetscWeakForm, wfNew::PetscWeakForm) 
Copy the pointwise functions to another `PetscWeakForm`

Not Collective

Input Parameter:
- `wf` - The original `PetscWeakForm`

Output Parameter:
- `wfNew` - The copy of the `PetscWeakForm`

Level: intermediate

-seealso: `PetscWeakForm`, `PetscWeakFormCreate()`, `PetscWeakFormDestroy()`

# External Links
$(_doc_external("Dm/PetscWeakFormCopy"))
"""
function PetscWeakFormCopy(petsclib::PetscLibType, wf::PetscWeakForm, wfNew::PetscWeakForm) end

@for_petsc function PetscWeakFormCopy(petsclib::$UnionPetscLib, wf::PetscWeakForm, wfNew::PetscWeakForm )

    @chk ccall(
               (:PetscWeakFormCopy, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, PetscWeakForm),
               wf, wfNew,
              )


	return nothing
end 

"""
	PetscWeakFormClear(petsclib::PetscLibType,wf::PetscWeakForm) 
Clear all functions from the `PetscWeakForm`

Not Collective

Input Parameter:
- `wf` - The original `PetscWeakForm`

Level: intermediate

-seealso: `PetscWeakForm`, `PetscWeakFormCopy()`, `PetscWeakFormCreate()`, `PetscWeakFormDestroy()`

# External Links
$(_doc_external("Dm/PetscWeakFormClear"))
"""
function PetscWeakFormClear(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormClear(petsclib::$UnionPetscLib, wf::PetscWeakForm )

    @chk ccall(
               (:PetscWeakFormClear, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm,),
               wf,
              )


	return nothing
end 

"""
	PetscWeakFormRewriteKeys(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, Nv::PetscInt, values::Vector{PetscInt}) 
Change any key on the given label to use the new set of label values

Not Collective

Input Parameters:
- `wf`     - The original `PetscWeakForm`
- `label`  - The label to change keys for
- `Nv`     - The number of new label values
- `values` - The set of new values to relabel keys with

Level: intermediate

-seealso: `PetscWeakForm`, `DMLabel`, `PetscWeakFormReplaceLabel()`, `PetscWeakFormCreate()`, `PetscWeakFormDestroy()`

# External Links
$(_doc_external("Dm/PetscWeakFormRewriteKeys"))
"""
function PetscWeakFormRewriteKeys(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, Nv::PetscInt, values::Vector{PetscInt}) end

@for_petsc function PetscWeakFormRewriteKeys(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, Nv::$PetscInt, values::Vector{$PetscInt} )

    @chk ccall(
               (:PetscWeakFormRewriteKeys, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, Ptr{$PetscInt}),
               wf, label, Nv, values,
              )


	return nothing
end 

"""
	PetscWeakFormReplaceLabel(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel) 
Change any key on a label of the same name to use the new label

Not Collective

Input Parameters:
- `wf`    - The original `PetscWeakForm`
- `label` - The label to change keys for

Level: intermediate

-seealso: `PetscWeakForm`, `DMLabel`, `PetscWeakFormRewriteKeys()`, `PetscWeakFormCreate()`, `PetscWeakFormDestroy()`

# External Links
$(_doc_external("Dm/PetscWeakFormReplaceLabel"))
"""
function PetscWeakFormReplaceLabel(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel) end

@for_petsc function PetscWeakFormReplaceLabel(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel )

    @chk ccall(
               (:PetscWeakFormReplaceLabel, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel),
               wf, label,
              )


	return nothing
end 

"""
	PetscWeakFormClearIndex(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, kind::PetscWeakFormKind, ind::PetscInt) 

# External Links
$(_doc_external("Dm/PetscWeakFormClearIndex"))
"""
function PetscWeakFormClearIndex(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, kind::PetscWeakFormKind, ind::PetscInt) end

@for_petsc function PetscWeakFormClearIndex(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, kind::PetscWeakFormKind, ind::$PetscInt )

    @chk ccall(
               (:PetscWeakFormClearIndex, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, PetscWeakFormKind, $PetscInt),
               wf, label, val, f, part, kind, ind,
              )


	return nothing
end 

#=
"""
	n::PetscInt = PetscWeakFormGetObjective(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(obj) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetObjective"))
"""
function PetscWeakFormGetObjective(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(obj) end

@for_petsc function PetscWeakFormGetObjective(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, ::Cvoid(obj )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetObjective, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid(obj),
               wf, label, val, f, part, n_, g0,
              )

	n = n_[]

	return n
end 


"""
	PetscWeakFormSetObjective(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n::PetscInt, ::Cvoid(obj) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetObjective"))
"""
function PetscWeakFormSetObjective(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n::PetscInt, ::Cvoid(obj) end

@for_petsc function PetscWeakFormSetObjective(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, n::$PetscInt, ::Cvoid(obj )

    @chk ccall(
               (:PetscWeakFormSetObjective, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid(obj),
               wf, label, val, f, part, n, g0,
              )


	return nothing
end 
=#

"""
	PetscWeakFormAddObjective(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, obj::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddObjective"))
"""
function PetscWeakFormAddObjective(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, obj::external) end

@for_petsc function PetscWeakFormAddObjective(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, obj::external )

    @chk ccall(
               (:PetscWeakFormAddObjective, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, external),
               wf, label, val, f, part, obj,
              )


	return nothing
end 

#=
"""
	PetscWeakFormGetIndexObjective(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ind::PetscInt, ::Cvoid(obj) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetIndexObjective"))
"""
function PetscWeakFormGetIndexObjective(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ind::PetscInt, ::Cvoid(obj) end

@for_petsc function PetscWeakFormGetIndexObjective(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, ind::$PetscInt, ::Cvoid(obj )

    @chk ccall(
               (:PetscWeakFormGetIndexObjective, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid(obj),
               wf, label, val, f, part, ind, g0,
              )


	return nothing
end 
=#
"""
	PetscWeakFormSetIndexObjective(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ind::PetscInt, obj::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexObjective"))
"""
function PetscWeakFormSetIndexObjective(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ind::PetscInt, obj::external) end

@for_petsc function PetscWeakFormSetIndexObjective(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, ind::$PetscInt, obj::external )

    @chk ccall(
               (:PetscWeakFormSetIndexObjective, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external),
               wf, label, val, f, part, ind, obj,
              )


	return nothing
end 

#=
"""
	n0::PetscInt = PetscWeakFormGetResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(f0) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetResidual"))
"""
function PetscWeakFormGetResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(f0) end

@for_petsc function PetscWeakFormGetResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, ::Cvoid(f0 )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid(f0),
               wf, label, val, f, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 
=#
"""
	PetscWeakFormAddResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, f0::external, f1::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddResidual"))
"""
function PetscWeakFormAddResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, f0::external, f1::external) end

@for_petsc function PetscWeakFormAddResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, f0::external, f1::external )

    @chk ccall(
               (:PetscWeakFormAddResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, external, external),
               wf, label, val, f, part, f0, f1,
              )


	return nothing
end 

#=
"""
	PetscWeakFormSetResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n0::PetscInt, ::Cvoid(f0) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetResidual"))
"""
function PetscWeakFormSetResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n0::PetscInt, ::Cvoid(f0) end

@for_petsc function PetscWeakFormSetResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, n0::$PetscInt, ::Cvoid(f0 )

    @chk ccall(
               (:PetscWeakFormSetResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid(f0),
               wf, label, val, f, part, n0, g0,
              )


	return nothing
end 
=#
"""
	PetscWeakFormSetIndexResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i0::PetscInt, f0::external, i1::PetscInt, f1::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexResidual"))
"""
function PetscWeakFormSetIndexResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i0::PetscInt, f0::external, i1::PetscInt, f1::external) end

@for_petsc function PetscWeakFormSetIndexResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, i0::$PetscInt, f0::external, i1::$PetscInt, f1::external )

    @chk ccall(
               (:PetscWeakFormSetIndexResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, part, i0, f0, i1, f1,
              )


	return nothing
end 

#=
"""
	n0::PetscInt = PetscWeakFormGetBdResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(f0) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetBdResidual"))
"""
function PetscWeakFormGetBdResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, ::Cvoid(f0) end

@for_petsc function PetscWeakFormGetBdResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, ::Cvoid(f0 )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid(f0),
               wf, label, val, f, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 
=#
"""
	PetscWeakFormAddBdResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, f0::external, f1::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddBdResidual"))
"""
function PetscWeakFormAddBdResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, f0::external, f1::external) end

@for_petsc function PetscWeakFormAddBdResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, f0::external, f1::external )

    @chk ccall(
               (:PetscWeakFormAddBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, external, external),
               wf, label, val, f, part, f0, f1,
              )


	return nothing
end 

#=
"""
	PetscWeakFormSetBdResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n0::PetscInt, ::Cvoid(f0) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetBdResidual"))
"""
function PetscWeakFormSetBdResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n0::PetscInt, ::Cvoid(f0) end

@for_petsc function PetscWeakFormSetBdResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, n0::$PetscInt, ::Cvoid(f0 )

    @chk ccall(
               (:PetscWeakFormSetBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid(f0),
               wf, label, val, f, part, n0, g0,
              )


	return nothing
end 
=#

"""
	PetscWeakFormSetIndexBdResidual(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i0::PetscInt, f0::external, i1::PetscInt, f1::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexBdResidual"))
"""
function PetscWeakFormSetIndexBdResidual(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i0::PetscInt, f0::external, i1::PetscInt, f1::external) end

@for_petsc function PetscWeakFormSetIndexBdResidual(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, i0::$PetscInt, f0::external, i1::$PetscInt, f1::external )

    @chk ccall(
               (:PetscWeakFormSetIndexBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, part, i0, f0, i1, f1,
              )


	return nothing
end 

"""
	hasJac::PetscBool = PetscWeakFormHasJacobian(petsclib::PetscLibType,wf::PetscWeakForm) 

# External Links
$(_doc_external("Dm/PetscWeakFormHasJacobian"))
"""
function PetscWeakFormHasJacobian(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormHasJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	hasJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscWeakFormHasJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{PetscBool}),
               wf, hasJac_,
              )

	hasJac = hasJac_[]

	return hasJac
end 

#=
"""
	n0::PetscInt = PetscWeakFormGetJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetJacobian"))
"""
function PetscWeakFormGetJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormGetJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::Cvoid )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, g, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 
=#

"""
	PetscWeakFormAddJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddJacobian"))
"""
function PetscWeakFormAddJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) end

@for_petsc function PetscWeakFormAddJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::external, g1::external, g2::external, g3::external )

    @chk ccall(
               (:PetscWeakFormAddJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, external, external, external),
               wf, label, val, f, g, part, g0, g1, g2, g3,
              )


	return nothing
end 

#=
"""
	PetscWeakFormSetJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetJacobian"))
"""
function PetscWeakFormSetJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormSetJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, n0::$PetscInt, g0::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, g, part, n0, g0,
              )


	return nothing
end 
=#

"""
	PetscWeakFormSetIndexJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexJacobian"))
"""
function PetscWeakFormSetIndexJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) end

@for_petsc function PetscWeakFormSetIndexJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, i0::$PetscInt, g0::external, i1::$PetscInt, g1::external, i2::$PetscInt, g2::external, i3::$PetscInt, g3::external )

    @chk ccall(
               (:PetscWeakFormSetIndexJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, g, part, i0, g0, i1, g1, i2, g2, i3, g3,
              )


	return nothing
end 

"""
	hasJacPre::PetscBool = PetscWeakFormHasJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm) 

# External Links
$(_doc_external("Dm/PetscWeakFormHasJacobianPreconditioner"))
"""
function PetscWeakFormHasJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormHasJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	hasJacPre_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscWeakFormHasJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{PetscBool}),
               wf, hasJacPre_,
              )

	hasJacPre = hasJacPre_[]

	return hasJacPre
end 

#=
"""
	n0::PetscInt = PetscWeakFormGetJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetJacobianPreconditioner"))
"""
function PetscWeakFormGetJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormGetJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::Cvoid )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, g, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 
=#

"""
	PetscWeakFormAddJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddJacobianPreconditioner"))
"""
function PetscWeakFormAddJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) end

@for_petsc function PetscWeakFormAddJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::external, g1::external, g2::external, g3::external )

    @chk ccall(
               (:PetscWeakFormAddJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, external, external, external),
               wf, label, val, f, g, part, g0, g1, g2, g3,
              )


	return nothing
end 

#=
"""
	PetscWeakFormSetJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetJacobianPreconditioner"))
"""
function PetscWeakFormSetJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormSetJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, n0::$PetscInt, g0::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, g, part, n0, g0,
              )


	return nothing
end 
=#

"""
	PetscWeakFormSetIndexJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexJacobianPreconditioner"))
"""
function PetscWeakFormSetIndexJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) end

@for_petsc function PetscWeakFormSetIndexJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, i0::$PetscInt, g0::external, i1::$PetscInt, g1::external, i2::$PetscInt, g2::external, i3::$PetscInt, g3::external )

    @chk ccall(
               (:PetscWeakFormSetIndexJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, g, part, i0, g0, i1, g1, i2, g2, i3, g3,
              )


	return nothing
end 

"""
	hasJac::PetscBool = PetscWeakFormHasBdJacobian(petsclib::PetscLibType,wf::PetscWeakForm) 

# External Links
$(_doc_external("Dm/PetscWeakFormHasBdJacobian"))
"""
function PetscWeakFormHasBdJacobian(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormHasBdJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	hasJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscWeakFormHasBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{PetscBool}),
               wf, hasJac_,
              )

	hasJac = hasJac_[]

	return hasJac
end 

"""
	n0::PetscInt = PetscWeakFormGetBdJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetBdJacobian"))
"""
function PetscWeakFormGetBdJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormGetBdJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::Cvoid )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, g, part, n0_, g0
              )

	n0 = n0_[]

	return n0
end 

"""
	PetscWeakFormAddBdJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddBdJacobian"))
"""
function PetscWeakFormAddBdJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) end

@for_petsc function PetscWeakFormAddBdJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::external, g1::external, g2::external, g3::external )

    @chk ccall(
               (:PetscWeakFormAddBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, external, external, external),
               wf, label, val, f, g, part, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscWeakFormSetBdJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetBdJacobian"))
"""
function PetscWeakFormSetBdJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormSetBdJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, n0::$PetscInt, g0::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, g, part, n0, g0,
              )


	return nothing
end 

"""
	PetscWeakFormSetIndexBdJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexBdJacobian"))
"""
function PetscWeakFormSetIndexBdJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) end

@for_petsc function PetscWeakFormSetIndexBdJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, i0::$PetscInt, g0::external, i1::$PetscInt, g1::external, i2::$PetscInt, g2::external, i3::$PetscInt, g3::external )

    @chk ccall(
               (:PetscWeakFormSetIndexBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, g, part, i0, g0, i1, g1, i2, g2, i3, g3,
              )


	return nothing
end 

"""
	hasJacPre::PetscBool = PetscWeakFormHasBdJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm) 

# External Links
$(_doc_external("Dm/PetscWeakFormHasBdJacobianPreconditioner"))
"""
function PetscWeakFormHasBdJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormHasBdJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	hasJacPre_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscWeakFormHasBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{PetscBool}),
               wf, hasJacPre_,
              )

	hasJacPre = hasJacPre_[]

	return hasJacPre
end 

"""
	n0::PetscInt = PetscWeakFormGetBdJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetBdJacobianPreconditioner"))
"""
function PetscWeakFormGetBdJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormGetBdJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::Cvoid )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, g, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 

"""
	PetscWeakFormAddBdJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddBdJacobianPreconditioner"))
"""
function PetscWeakFormAddBdJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) end

@for_petsc function PetscWeakFormAddBdJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::external, g1::external, g2::external, g3::external )

    @chk ccall(
               (:PetscWeakFormAddBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, external, external, external),
               wf, label, val, f, g, part, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscWeakFormSetBdJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetBdJacobianPreconditioner"))
"""
function PetscWeakFormSetBdJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormSetBdJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, n0::$PetscInt, g0::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, g, part, n0, g0,
              )


	return nothing
end 

"""
	PetscWeakFormSetIndexBdJacobianPreconditioner(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexBdJacobianPreconditioner"))
"""
function PetscWeakFormSetIndexBdJacobianPreconditioner(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) end

@for_petsc function PetscWeakFormSetIndexBdJacobianPreconditioner(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, i0::$PetscInt, g0::external, i1::$PetscInt, g1::external, i2::$PetscInt, g2::external, i3::$PetscInt, g3::external )

    @chk ccall(
               (:PetscWeakFormSetIndexBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, g, part, i0, g0, i1, g1, i2, g2, i3, g3,
              )


	return nothing
end 

"""
	hasDynJac::PetscBool = PetscWeakFormHasDynamicJacobian(petsclib::PetscLibType,wf::PetscWeakForm) 

# External Links
$(_doc_external("Dm/PetscWeakFormHasDynamicJacobian"))
"""
function PetscWeakFormHasDynamicJacobian(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormHasDynamicJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	hasDynJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscWeakFormHasDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{PetscBool}),
               wf, hasDynJac_,
              )

	hasDynJac = hasDynJac_[]

	return hasDynJac
end 

"""
	n0::PetscInt = PetscWeakFormGetDynamicJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetDynamicJacobian"))
"""
function PetscWeakFormGetDynamicJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormGetDynamicJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::Cvoid )
	n0_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, g, part, n0_, g0,
              )

	n0 = n0_[]

	return n0
end 

"""
	PetscWeakFormAddDynamicJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormAddDynamicJacobian"))
"""
function PetscWeakFormAddDynamicJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, g0::external, g1::external, g2::external, g3::external) end

@for_petsc function PetscWeakFormAddDynamicJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, g0::external, g1::external, g2::external, g3::external )

    @chk ccall(
               (:PetscWeakFormAddDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, external, external, external),
               wf, label, val, f, g, part, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscWeakFormSetDynamicJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetDynamicJacobian"))
"""
function PetscWeakFormSetDynamicJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, n0::PetscInt, g0::Cvoid) end

@for_petsc function PetscWeakFormSetDynamicJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, n0::$PetscInt, g0::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, g, part, n0, g0
              )


	return nothing
end 

"""
	PetscWeakFormSetIndexDynamicJacobian(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexDynamicJacobian"))
"""
function PetscWeakFormSetIndexDynamicJacobian(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, g::PetscInt, part::PetscInt, i0::PetscInt, g0::external, i1::PetscInt, g1::external, i2::PetscInt, g2::external, i3::PetscInt, g3::external) end

@for_petsc function PetscWeakFormSetIndexDynamicJacobian(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, g::$PetscInt, part::$PetscInt, i0::$PetscInt, g0::external, i1::$PetscInt, g1::external, i2::$PetscInt, g2::external, i3::$PetscInt, g3::external )

    @chk ccall(
               (:PetscWeakFormSetIndexDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external, $PetscInt, external, $PetscInt, external, $PetscInt, external),
               wf, label, val, f, g, part, i0, g0, i1, g1, i2, g2, i3, g3,
              )


	return nothing
end 

"""
	n::PetscInt = PetscWeakFormGetRiemannSolver(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, r::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormGetRiemannSolver"))
"""
function PetscWeakFormGetRiemannSolver(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, r::Cvoid) end

@for_petsc function PetscWeakFormGetRiemannSolver(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, r::Cvoid)
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetRiemannSolver, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Cvoid),
               wf, label, val, f, part, n_, r,
              )

	n = n_[]

	return n
end 

"""
	PetscWeakFormSetRiemannSolver(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n::PetscInt, r::Cvoid) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetRiemannSolver"))
"""
function PetscWeakFormSetRiemannSolver(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, n::PetscInt, r::Cvoid) end

@for_petsc function PetscWeakFormSetRiemannSolver(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, n::$PetscInt, r::Cvoid )

    @chk ccall(
               (:PetscWeakFormSetRiemannSolver, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Cvoid),
               wf, label, val, f, part, n, r,
              )


	return nothing
end 

"""
	PetscWeakFormSetIndexRiemannSolver(petsclib::PetscLibType,wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i::PetscInt, r::external) 

# External Links
$(_doc_external("Dm/PetscWeakFormSetIndexRiemannSolver"))
"""
function PetscWeakFormSetIndexRiemannSolver(petsclib::PetscLibType, wf::PetscWeakForm, label::DMLabel, val::PetscInt, f::PetscInt, part::PetscInt, i::PetscInt, r::external) end

@for_petsc function PetscWeakFormSetIndexRiemannSolver(petsclib::$UnionPetscLib, wf::PetscWeakForm, label::DMLabel, val::$PetscInt, f::$PetscInt, part::$PetscInt, i::$PetscInt, r::external )

    @chk ccall(
               (:PetscWeakFormSetIndexRiemannSolver, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, DMLabel, $PetscInt, $PetscInt, $PetscInt, $PetscInt, external),
               wf, label, val, f, part, i, r,
              )


	return nothing
end 

"""
	Nf::PetscInt = PetscWeakFormGetNumFields(petsclib::PetscLibType,wf::PetscWeakForm) 
Returns the number of fields in a `PetscWeakForm`

Not Collective

Input Parameter:
- `wf` - The `PetscWeakForm` object

Output Parameter:
- `Nf` - The number of fields

Level: beginner

-seealso: `PetscWeakForm`, `PetscWeakFormSetNumFields()`, `PetscWeakFormCreate()`

# External Links
$(_doc_external("Dm/PetscWeakFormGetNumFields"))
"""
function PetscWeakFormGetNumFields(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormGetNumFields(petsclib::$UnionPetscLib, wf::PetscWeakForm )
	Nf_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscWeakFormGetNumFields, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, Ptr{$PetscInt}),
               wf, Nf_,
              )

	Nf = Nf_[]

	return Nf
end 

"""
	PetscWeakFormSetNumFields(petsclib::PetscLibType,wf::PetscWeakForm, Nf::PetscInt) 
Sets the number of fields

Not Collective

Input Parameters:
- `wf` - The `PetscWeakForm` object
- `Nf` - The number of fields

Level: beginner

-seealso: `PetscWeakForm`, `PetscWeakFormGetNumFields()`, `PetscWeakFormCreate()`

# External Links
$(_doc_external("Dm/PetscWeakFormSetNumFields"))
"""
function PetscWeakFormSetNumFields(petsclib::PetscLibType, wf::PetscWeakForm, Nf::PetscInt) end

@for_petsc function PetscWeakFormSetNumFields(petsclib::$UnionPetscLib, wf::PetscWeakForm, Nf::$PetscInt )

    @chk ccall(
               (:PetscWeakFormSetNumFields, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, $PetscInt),
               wf, Nf,
              )


	return nothing
end 

"""
	PetscWeakFormDestroy(petsclib::PetscLibType,wf::PetscWeakForm) 
Destroys a `PetscWeakForm` object

Collective

Input Parameter:
- `wf` - the `PetscWeakForm` object to destroy

Level: developer

-seealso: `PetscWeakForm`, `PetscWeakFormCreate()`, `PetscWeakFormView()`

# External Links
$(_doc_external("Dm/PetscWeakFormDestroy"))
"""
function PetscWeakFormDestroy(petsclib::PetscLibType, wf::PetscWeakForm) end

@for_petsc function PetscWeakFormDestroy(petsclib::$UnionPetscLib, wf::PetscWeakForm )

    @chk ccall(
               (:PetscWeakFormDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscWeakForm},),
               wf,
              )


	return nothing
end 

"""
	PetscWeakFormView(petsclib::PetscLibType,wf::PetscWeakForm, v::PetscViewer) 
Views a `PetscWeakForm`

Collective

Input Parameters:
- `wf` - the `PetscWeakForm` object to view
- `v`  - the viewer

Level: developer

-seealso: `PetscViewer`, `PetscWeakForm`, `PetscWeakFormDestroy()`, `PetscWeakFormCreate()`

# External Links
$(_doc_external("Dm/PetscWeakFormView"))
"""
function PetscWeakFormView(petsclib::PetscLibType, wf::PetscWeakForm, v::PetscViewer) end

@for_petsc function PetscWeakFormView(petsclib::$UnionPetscLib, wf::PetscWeakForm, v::PetscViewer )

    @chk ccall(
               (:PetscWeakFormView, $petsc_library),
               PetscErrorCode,
               (PetscWeakForm, PetscViewer),
               wf, v,
              )


	return nothing
end 

"""
	wf::PetscWeakForm = PetscWeakFormCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscWeakForm` object.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscWeakForm` object

Output Parameter:
- `wf` - The `PetscWeakForm` object

Level: beginner

-seealso: `PetscWeakForm`, `PetscDS`, `PetscWeakFormDestroy()`

# External Links
$(_doc_external("Dm/PetscWeakFormCreate"))
"""
function PetscWeakFormCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscWeakFormCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	wf_ = Ref{PetscWeakForm}()

    @chk ccall(
               (:PetscWeakFormCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscWeakForm}),
               comm, wf_,
              )

	wf = wf_[]

	return wf
end 

