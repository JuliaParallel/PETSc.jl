# autodefined type arguments for class ------
mutable struct _n_PetscSpace end
const PetscSpace = Ptr{_n_PetscSpace}

mutable struct _n_PetscDualSpace end
const PetscDualSpace = Ptr{_n_PetscDualSpace}

# -------------------------------------------------------
"""
	PetscSpaceRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscSpace` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine for the implementation type

-seealso: `PetscSpace`, `PetscSpaceRegisterAll()`, `PetscSpaceRegisterDestroy()`

# External Links
$(_doc_external("DM/PetscSpaceRegister"))
"""
function PetscSpaceRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscSpaceRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscSpaceRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscSpaceSetType(petsclib::PetscLibType,sp::PetscSpace, name::PetscSpaceType) 
Builds a particular `PetscSpace`

Collective

Input Parameters:
- `sp`   - The `PetscSpace` object
- `name` - The kind of space

Options Database Key:
- `-petscspace_type <type>` - Sets the `PetscSpace` type; use -help for a list of available types

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceType`, `PetscSpaceGetType()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceSetType"))
"""
function PetscSpaceSetType(petsclib::PetscLibType, sp::PetscSpace, name::PetscSpaceType) end

@for_petsc function PetscSpaceSetType(petsclib::$UnionPetscLib, sp::PetscSpace, name::PetscSpaceType )

    @chk ccall(
               (:PetscSpaceSetType, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscSpaceType),
               sp, name,
              )


	return nothing
end 

"""
	name::PetscSpaceType = PetscSpaceGetType(petsclib::PetscLibType,sp::PetscSpace) 
Gets the `PetscSpaceType` (as a string) from the object.

Not Collective

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameter:
- `name` - The `PetscSpace` type name

Level: intermediate

-seealso: `PetscSpaceType`, `PetscSpace`, `PetscSpaceSetType()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceGetType"))
"""
function PetscSpaceGetType(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceGetType(petsclib::$UnionPetscLib, sp::PetscSpace )
	name_ = Ref{PetscSpaceType}()

    @chk ccall(
               (:PetscSpaceGetType, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{PetscSpaceType}),
               sp, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscSpaceViewFromOptions(petsclib::PetscLibType,A::PetscSpace, obj::PetscObject, name::String) 
View a `PetscSpace` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscSpace` object
- `obj`  - Optional object that provides the options name prefix
- `name` - command line option name

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceView()`, `PetscObjectViewFromOptions()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceViewFromOptions"))
"""
function PetscSpaceViewFromOptions(petsclib::PetscLibType, A::PetscSpace, obj::PetscObject, name::String) end

@for_petsc function PetscSpaceViewFromOptions(petsclib::$UnionPetscLib, A::PetscSpace, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscSpaceViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscSpaceView(petsclib::PetscLibType,sp::PetscSpace, v::PetscViewer) 
Views a `PetscSpace`

Collective

Input Parameters:
- `sp` - the `PetscSpace` object to view
- `v`  - the viewer

Level: beginner

-seealso: `PetscSpace`, `PetscViewer`, `PetscSpaceViewFromOptions()`, `PetscSpaceDestroy()`

# External Links
$(_doc_external("DM/PetscSpaceView"))
"""
function PetscSpaceView(petsclib::PetscLibType, sp::PetscSpace, v::PetscViewer) end

@for_petsc function PetscSpaceView(petsclib::$UnionPetscLib, sp::PetscSpace, v::PetscViewer )

    @chk ccall(
               (:PetscSpaceView, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscViewer),
               sp, v,
              )


	return nothing
end 

"""
	PetscSpaceSetFromOptions(petsclib::PetscLibType,sp::PetscSpace) 
sets parameters in a `PetscSpace` from the options database

Collective

Input Parameter:
- `sp` - the `PetscSpace` object to set options for

Options Database Keys:
- `-petscspace_degree <deg>`   - the approximation order of the space
- `-petscspace_variables <n>`  - the number of different variables, e.g. x and y
- `-petscspace_components <c>` - the number of components, say d for a vector field

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceView()`

# External Links
$(_doc_external("DM/PetscSpaceSetFromOptions"))
"""
function PetscSpaceSetFromOptions(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceSetFromOptions(petsclib::$UnionPetscLib, sp::PetscSpace )

    @chk ccall(
               (:PetscSpaceSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSpace,),
               sp,
              )


	return nothing
end 

"""
	PetscSpaceSetUp(petsclib::PetscLibType,sp::PetscSpace) 
Construct data structures for the `PetscSpace`

Collective

Input Parameter:
- `sp` - the `PetscSpace` object to setup

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceView()`, `PetscSpaceDestroy()`

# External Links
$(_doc_external("DM/PetscSpaceSetUp"))
"""
function PetscSpaceSetUp(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceSetUp(petsclib::$UnionPetscLib, sp::PetscSpace )

    @chk ccall(
               (:PetscSpaceSetUp, $petsc_library),
               PetscErrorCode,
               (PetscSpace,),
               sp,
              )


	return nothing
end 

"""
	PetscSpaceDestroy(petsclib::PetscLibType,sp::PetscSpace) 
Destroys a `PetscSpace` object

Collective

Input Parameter:
- `sp` - the `PetscSpace` object to destroy

Level: beginner

-seealso: `PetscSpace`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceDestroy"))
"""
function PetscSpaceDestroy(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceDestroy(petsclib::$UnionPetscLib, sp::PetscSpace )

    @chk ccall(
               (:PetscSpaceDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSpace},),
               sp,
              )


	return nothing
end 

"""
	sp::PetscSpace = PetscSpaceCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscSpace` object. The type can then be set with `PetscSpaceSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscSpace` object

Output Parameter:
- `sp` - The `PetscSpace` object

Level: beginner

-seealso: `PetscSpace`, `PetscSpaceSetType()`, `PETSCSPACEPOLYNOMIAL`

# External Links
$(_doc_external("DM/PetscSpaceCreate"))
"""
function PetscSpaceCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSpaceCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	sp_ = Ref{PetscSpace}()

    @chk ccall(
               (:PetscSpaceCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSpace}),
               comm, sp_,
              )

	sp = sp_[]

	return sp
end 

"""
	dim::PetscInt = PetscSpaceGetDimension(petsclib::PetscLibType,sp::PetscSpace) 
Return the dimension of this space, i.e. the number of basis vectors

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameter:
- `dim` - The dimension

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceGetDegree()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceGetDimension"))
"""
function PetscSpaceGetDimension(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceGetDimension(petsclib::$UnionPetscLib, sp::PetscSpace )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	minDegree::PetscInt,maxDegree::PetscInt = PetscSpaceGetDegree(petsclib::PetscLibType,sp::PetscSpace) 
Return the polynomial degrees that characterize this space

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameters:
- `minDegree` - The degree of the largest polynomial space contained in the space, pass `NULL` if not needed
- `maxDegree` - The degree of the smallest polynomial space containing the space, pass `NULL` if not needed

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceSetDegree()`, `PetscSpaceGetDimension()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceGetDegree"))
"""
function PetscSpaceGetDegree(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceGetDegree(petsclib::$UnionPetscLib, sp::PetscSpace )
	minDegree_ = Ref{$PetscInt}()
	maxDegree_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceGetDegree, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}, Ptr{$PetscInt}),
               sp, minDegree_, maxDegree_,
              )

	minDegree = minDegree_[]
	maxDegree = maxDegree_[]

	return minDegree,maxDegree
end 

"""
	PetscSpaceSetDegree(petsclib::PetscLibType,sp::PetscSpace, degree::PetscInt, maxDegree::PetscInt) 
Set the degree of approximation for this space.

Input Parameters:
- `sp`        - The `PetscSpace`
- `degree`    - The degree of the largest polynomial space contained in the space
- `maxDegree` - The degree of the largest polynomial space containing the space.  One of degree and maxDegree can be `PETSC_DETERMINE`.

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceGetDegree()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceSetDegree"))
"""
function PetscSpaceSetDegree(petsclib::PetscLibType, sp::PetscSpace, degree::PetscInt, maxDegree::PetscInt) end

@for_petsc function PetscSpaceSetDegree(petsclib::$UnionPetscLib, sp::PetscSpace, degree::$PetscInt, maxDegree::$PetscInt )

    @chk ccall(
               (:PetscSpaceSetDegree, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, $PetscInt),
               sp, degree, maxDegree,
              )


	return nothing
end 

"""
	Nc::PetscInt = PetscSpaceGetNumComponents(petsclib::PetscLibType,sp::PetscSpace) 
Return the number of components for this space

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameter:
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceSetNumComponents()`, `PetscSpaceGetNumVariables()`, `PetscSpaceGetDimension()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceGetNumComponents"))
"""
function PetscSpaceGetNumComponents(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceGetNumComponents(petsclib::$UnionPetscLib, sp::PetscSpace )
	Nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceGetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, Nc_,
              )

	Nc = Nc_[]

	return Nc
end 

"""
	PetscSpaceSetNumComponents(petsclib::PetscLibType,sp::PetscSpace, Nc::PetscInt) 
Set the number of components for this space

Input Parameters:
- `sp` - The `PetscSpace`
- `Nc` - The number of components

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceGetNumComponents()`, `PetscSpaceSetNumVariables()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceSetNumComponents"))
"""
function PetscSpaceSetNumComponents(petsclib::PetscLibType, sp::PetscSpace, Nc::PetscInt) end

@for_petsc function PetscSpaceSetNumComponents(petsclib::$UnionPetscLib, sp::PetscSpace, Nc::$PetscInt )

    @chk ccall(
               (:PetscSpaceSetNumComponents, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt),
               sp, Nc,
              )


	return nothing
end 

"""
	PetscSpaceSetNumVariables(petsclib::PetscLibType,sp::PetscSpace, n::PetscInt) 
Set the number of variables for this space

Input Parameters:
- `sp` - The `PetscSpace`
- `n`  - The number of variables, e.g. x, y, z...

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceGetNumVariables()`, `PetscSpaceSetNumComponents()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceSetNumVariables"))
"""
function PetscSpaceSetNumVariables(petsclib::PetscLibType, sp::PetscSpace, n::PetscInt) end

@for_petsc function PetscSpaceSetNumVariables(petsclib::$UnionPetscLib, sp::PetscSpace, n::$PetscInt )

    @chk ccall(
               (:PetscSpaceSetNumVariables, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt),
               sp, n,
              )


	return nothing
end 

"""
	n::PetscInt = PetscSpaceGetNumVariables(petsclib::PetscLibType,sp::PetscSpace) 
Return the number of variables for this space

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameter:
- `n` - The number of variables, e.g. x, y, z...

Level: intermediate

-seealso: `PetscSpace`, `PetscSpaceSetNumVariables()`, `PetscSpaceGetNumComponents()`, `PetscSpaceGetDimension()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceGetNumVariables"))
"""
function PetscSpaceGetNumVariables(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceGetNumVariables(petsclib::$UnionPetscLib, sp::PetscSpace )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceGetNumVariables, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, n_,
              )

	n = n_[]

	return n
end 

"""
	B::Vector{PetscReal},D::Vector{PetscReal},H::Vector{PetscReal} = PetscSpaceEvaluate(petsclib::PetscLibType,sp::PetscSpace, npoints::PetscInt, points::Vector{PetscReal}) 
Evaluate the basis functions and their derivatives (jet) at each point

Input Parameters:
- `sp`      - The `PetscSpace`
- `npoints` - The number of evaluation points, in reference coordinates
- `points`  - The point coordinates

Output Parameters:
- `B` - The function evaluations in a `npoints` x `nfuncs` array
- `D` - The derivative evaluations in a `npoints` x `nfuncs` x `dim` array
- `H` - The second derivative evaluations in a `npoints` x `nfuncs` x `dim` x `dim` array

Level: beginner

-seealso: `PetscSpace`, `PetscFECreateTabulation()`, `PetscFEGetCellTabulation()`, `PetscSpaceCreate()`

# External Links
$(_doc_external("DM/PetscSpaceEvaluate"))
"""
function PetscSpaceEvaluate(petsclib::PetscLibType, sp::PetscSpace, npoints::PetscInt, points::Vector{PetscReal}) end

@for_petsc function PetscSpaceEvaluate(petsclib::$UnionPetscLib, sp::PetscSpace, npoints::$PetscInt, points::Vector{$PetscReal} )
	B = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	D = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	H = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscSpaceEvaluate, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sp, npoints, points, B, D, H,
              )


	return B,D,H
end 

"""
	PetscSpaceGetHeightSubspace(petsclib::PetscLibType,sp::PetscSpace, height::PetscInt, subsp::PetscSpace) 
Get the subset of the primal space basis that is supported on a mesh point of a given height.

Not Collective

Input Parameters:
- `sp`     - the `PetscSpace` object
- `height` - the height of the mesh point for which the subspace is desired

Output Parameter:
- `subsp` - the subspace

Level: advanced

-seealso: `PetscDualSpaceGetHeightSubspace()`, `PetscSpace`

# External Links
$(_doc_external("DM/PetscSpaceGetHeightSubspace"))
"""
function PetscSpaceGetHeightSubspace(petsclib::PetscLibType, sp::PetscSpace, height::PetscInt, subsp::PetscSpace) end

@for_petsc function PetscSpaceGetHeightSubspace(petsclib::$UnionPetscLib, sp::PetscSpace, height::$PetscInt, subsp::PetscSpace )

    @chk ccall(
               (:PetscSpaceGetHeightSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, Ptr{PetscSpace}),
               sp, height, subsp,
              )


	return nothing
end 

"""
	subspace::PetscSpace = PetscSpaceCreateSubspace(petsclib::PetscLibType,origSpace::PetscSpace, dualSubspace::PetscDualSpace, x::PetscReal, Jx::PetscReal, u::PetscReal, Ju::PetscReal, copymode::PetscCopyMode) 
creates a subspace from a an `origSpace` and its dual `dualSubspace`

Input Parameters:
- `origSpace`    - the original `PetscSpace`
- `dualSubspace` - no idea
- `x`            - no idea
- `Jx`           - no idea
- `u`            - no idea
- `Ju`           - no idea
- `copymode`     - whether to copy, borrow, or own some of the input arrays I guess

Output Parameter:
- `subspace` - the subspace

Level: advanced

-seealso: `PetscSpace`, `PetscDualSpace`, `PetscCopyMode`, `PetscSpaceType`

# External Links
$(_doc_external("DM/PetscSpaceCreateSubspace"))
"""
function PetscSpaceCreateSubspace(petsclib::PetscLibType, origSpace::PetscSpace, dualSubspace::PetscDualSpace, x::PetscReal, Jx::PetscReal, u::PetscReal, Ju::PetscReal, copymode::PetscCopyMode) end

@for_petsc function PetscSpaceCreateSubspace(petsclib::$UnionPetscLib, origSpace::PetscSpace, dualSubspace::PetscDualSpace, x::$PetscReal, Jx::$PetscReal, u::$PetscReal, Ju::$PetscReal, copymode::PetscCopyMode )
	subspace_ = Ref{PetscSpace}()

    @chk ccall(
               (:PetscSpaceCreateSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscDualSpace, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, PetscCopyMode, Ptr{PetscSpace}),
               origSpace, dualSubspace, x, Jx, u, Ju, copymode, subspace_,
              )

	subspace = subspace_[]

	return subspace
end 

"""
	PetscSpacePointSetPoints(petsclib::PetscLibType,sp::PetscSpace, q::PetscQuadrature) 
Sets the evaluation points for the space to coincide with the points of a quadrature rule

Logically Collective

Input Parameters:
- `sp` - The `PetscSpace`
- `q`  - The `PetscQuadrature` defining the points

Level: intermediate

-seealso: `PetscSpace`, `PetscQuadrature`, `PetscSpaceCreate()`, `PetscSpaceSetType()`

# External Links
$(_doc_external("DM/PetscSpacePointSetPoints"))
"""
function PetscSpacePointSetPoints(petsclib::PetscLibType, sp::PetscSpace, q::PetscQuadrature) end

@for_petsc function PetscSpacePointSetPoints(petsclib::$UnionPetscLib, sp::PetscSpace, q::PetscQuadrature )

    @chk ccall(
               (:PetscSpacePointSetPoints, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscQuadrature),
               sp, q,
              )


	return nothing
end 

"""
	PetscSpacePointGetPoints(petsclib::PetscLibType,sp::PetscSpace, q::PetscQuadrature) 
Gets the evaluation points for the space as the points of a quadrature rule

Logically Collective

Input Parameter:
- `sp` - The `PetscSpace`

Output Parameter:
- `q` - The `PetscQuadrature` defining the points

Level: intermediate

-seealso: `PetscSpace`, `PetscQuadrature`, `PetscSpaceCreate()`, `PetscSpaceSetType()`

# External Links
$(_doc_external("DM/PetscSpacePointGetPoints"))
"""
function PetscSpacePointGetPoints(petsclib::PetscLibType, sp::PetscSpace, q::PetscQuadrature) end

@for_petsc function PetscSpacePointGetPoints(petsclib::$UnionPetscLib, sp::PetscSpace, q::PetscQuadrature )

    @chk ccall(
               (:PetscSpacePointGetPoints, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{PetscQuadrature}),
               sp, q,
              )


	return nothing
end 

"""
	numSumSpaces::PetscInt = PetscSpaceSumGetNumSubspaces(petsclib::PetscLibType,sp::PetscSpace) 
Get the number of spaces in the sum space

Input Parameter:
- `sp` - the function space object

Output Parameter:
- `numSumSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumSetNumSubspaces()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceSumGetNumSubspaces"))
"""
function PetscSpaceSumGetNumSubspaces(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceSumGetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscSpace )
	numSumSpaces_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceSumGetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, numSumSpaces_,
              )

	numSumSpaces = numSumSpaces_[]

	return numSumSpaces
end 

"""
	PetscSpaceSumSetNumSubspaces(petsclib::PetscLibType,sp::PetscSpace, numSumSpaces::PetscInt) 
Set the number of spaces in the sum space

Input Parameters:
- `sp`           - the function space object
- `numSumSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumGetNumSubspaces()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceSumSetNumSubspaces"))
"""
function PetscSpaceSumSetNumSubspaces(petsclib::PetscLibType, sp::PetscSpace, numSumSpaces::PetscInt) end

@for_petsc function PetscSpaceSumSetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscSpace, numSumSpaces::$PetscInt )

    @chk ccall(
               (:PetscSpaceSumSetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt),
               sp, numSumSpaces,
              )


	return nothing
end 

"""
	concatenate::PetscBool = PetscSpaceSumGetConcatenate(petsclib::PetscLibType,sp::PetscSpace) 
Get the concatenate flag for this space.

Input Parameter:
- `sp` - the function space object

Output Parameter:
- `concatenate` - flag indicating whether subspaces are concatenated.

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumSetConcatenate()`

# External Links
$(_doc_external("DM/PetscSpaceSumGetConcatenate"))
"""
function PetscSpaceSumGetConcatenate(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceSumGetConcatenate(petsclib::$UnionPetscLib, sp::PetscSpace )
	concatenate_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSpaceSumGetConcatenate, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{PetscBool}),
               sp, concatenate_,
              )

	concatenate = concatenate_[]

	return concatenate
end 

"""
	PetscSpaceSumSetConcatenate(petsclib::PetscLibType,sp::PetscSpace, concatenate::PetscBool) 
Sets the concatenate flag for this space.

Input Parameters:
- `sp`          - the function space object
- `concatenate` - are subspaces concatenated components (true) or direct summands (false)

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumGetConcatenate()`

# External Links
$(_doc_external("DM/PetscSpaceSumSetConcatenate"))
"""
function PetscSpaceSumSetConcatenate(petsclib::PetscLibType, sp::PetscSpace, concatenate::PetscBool) end

@for_petsc function PetscSpaceSumSetConcatenate(petsclib::$UnionPetscLib, sp::PetscSpace, concatenate::PetscBool )

    @chk ccall(
               (:PetscSpaceSumSetConcatenate, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscBool),
               sp, concatenate,
              )


	return nothing
end 

"""
	PetscSpaceSumGetSubspace(petsclib::PetscLibType,sp::PetscSpace, s::PetscInt, subsp::PetscSpace) 
Get a space in the sum space

Input Parameters:
- `sp` - the function space object
- `s`  - The space number

Output Parameter:
- `subsp` - the `PetscSpace`

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumSetSubspace()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceSumGetSubspace"))
"""
function PetscSpaceSumGetSubspace(petsclib::PetscLibType, sp::PetscSpace, s::PetscInt, subsp::PetscSpace) end

@for_petsc function PetscSpaceSumGetSubspace(petsclib::$UnionPetscLib, sp::PetscSpace, s::$PetscInt, subsp::PetscSpace )

    @chk ccall(
               (:PetscSpaceSumGetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, Ptr{PetscSpace}),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscSpaceSumSetSubspace(petsclib::PetscLibType,sp::PetscSpace, s::PetscInt, subsp::PetscSpace) 
Set a space in the sum space

Input Parameters:
- `sp`    - the function space object
- `s`     - The space number
- `subsp` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACESUM`, `PetscSpace`, `PetscSpaceSumGetSubspace()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceSumSetSubspace"))
"""
function PetscSpaceSumSetSubspace(petsclib::PetscLibType, sp::PetscSpace, s::PetscInt, subsp::PetscSpace) end

@for_petsc function PetscSpaceSumSetSubspace(petsclib::$UnionPetscLib, sp::PetscSpace, s::$PetscInt, subsp::PetscSpace )

    @chk ccall(
               (:PetscSpaceSumSetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, PetscSpace),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscSpaceSumSetInterleave(petsclib::PetscLibType,sp::PetscSpace, interleave_basis::PetscBool, interleave_components::PetscBool) 
Set whether the basis functions and components of a uniform sum are interleaved

Logically collective

Input Parameters:
- `sp`                    - a `PetscSpace` of type `PETSCSPACESUM`
- `interleave_basis`      - if `PETSC_TRUE`, the basis vectors of the subspaces are interleaved
- `interleave_components` - if `PETSC_TRUE` and the space concatenates components (`PetscSpaceSumGetConcatenate()`),
interleave the concatenated components

Level: developer

-seealso: `PetscSpace`, `PETSCSPACESUM`, `PETSCFEVECTOR`, `PetscSpaceSumGetInterleave()`

# External Links
$(_doc_external("DM/PetscSpaceSumSetInterleave"))
"""
function PetscSpaceSumSetInterleave(petsclib::PetscLibType, sp::PetscSpace, interleave_basis::PetscBool, interleave_components::PetscBool) end

@for_petsc function PetscSpaceSumSetInterleave(petsclib::$UnionPetscLib, sp::PetscSpace, interleave_basis::PetscBool, interleave_components::PetscBool )

    @chk ccall(
               (:PetscSpaceSumSetInterleave, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscBool, PetscBool),
               sp, interleave_basis, interleave_components,
              )


	return nothing
end 

"""
	interleave_basis::PetscBool,interleave_components::PetscBool = PetscSpaceSumGetInterleave(petsclib::PetscLibType,sp::PetscSpace) 
Get whether the basis functions and components of a uniform sum are interleaved

Logically collective

Input Parameter:
- `sp` - a `PetscSpace` of type `PETSCSPACESUM`

Output Parameters:
- `interleave_basis`      - if `PETSC_TRUE`, the basis vectors of the subspaces are interleaved
- `interleave_components` - if `PETSC_TRUE` and the space concatenates components (`PetscSpaceSumGetConcatenate()`),
interleave the concatenated components

Level: developer

-seealso: `PetscSpace`, `PETSCSPACESUM`, `PETSCFEVECTOR`, `PetscSpaceSumSetInterleave()`

# External Links
$(_doc_external("DM/PetscSpaceSumGetInterleave"))
"""
function PetscSpaceSumGetInterleave(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceSumGetInterleave(petsclib::$UnionPetscLib, sp::PetscSpace )
	interleave_basis_ = Ref{PetscBool}()
	interleave_components_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSpaceSumGetInterleave, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{PetscBool}, Ptr{PetscBool}),
               sp, interleave_basis_, interleave_components_,
              )

	interleave_basis = interleave_basis_[]
	interleave_components = interleave_components_[]

	return interleave_basis,interleave_components
end 

"""
	PetscSpaceTensorSetNumSubspaces(petsclib::PetscLibType,sp::PetscSpace, numTensSpaces::PetscInt) 
Set the number of spaces in the tensor product space

Input Parameters:
- `sp`            - the function space object
- `numTensSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACETENSOR`, `PetscSpace`, `PetscSpaceTensorGetNumSubspaces()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceTensorSetNumSubspaces"))
"""
function PetscSpaceTensorSetNumSubspaces(petsclib::PetscLibType, sp::PetscSpace, numTensSpaces::PetscInt) end

@for_petsc function PetscSpaceTensorSetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscSpace, numTensSpaces::$PetscInt )

    @chk ccall(
               (:PetscSpaceTensorSetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt),
               sp, numTensSpaces,
              )


	return nothing
end 

"""
	numTensSpaces::PetscInt = PetscSpaceTensorGetNumSubspaces(petsclib::PetscLibType,sp::PetscSpace) 
Get the number of spaces in the tensor product space

Input Parameter:
- `sp` - the function space object

Output Parameter:
- `numTensSpaces` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACETENSOR`, `PetscSpace`, `PetscSpaceTensorSetNumSubspaces()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceTensorGetNumSubspaces"))
"""
function PetscSpaceTensorGetNumSubspaces(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpaceTensorGetNumSubspaces(petsclib::$UnionPetscLib, sp::PetscSpace )
	numTensSpaces_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpaceTensorGetNumSubspaces, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, numTensSpaces_,
              )

	numTensSpaces = numTensSpaces_[]

	return numTensSpaces
end 

"""
	PetscSpaceTensorSetSubspace(petsclib::PetscLibType,sp::PetscSpace, s::PetscInt, subsp::PetscSpace) 
Set a space in the tensor product space

Input Parameters:
- `sp`    - the function space object
- `s`     - The space number
- `subsp` - the number of spaces

Level: intermediate

-seealso: `PETSCSPACETENSOR`, `PetscSpace`, `PetscSpaceTensorGetSubspace()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceTensorSetSubspace"))
"""
function PetscSpaceTensorSetSubspace(petsclib::PetscLibType, sp::PetscSpace, s::PetscInt, subsp::PetscSpace) end

@for_petsc function PetscSpaceTensorSetSubspace(petsclib::$UnionPetscLib, sp::PetscSpace, s::$PetscInt, subsp::PetscSpace )

    @chk ccall(
               (:PetscSpaceTensorSetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, PetscSpace),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscSpaceTensorGetSubspace(petsclib::PetscLibType,sp::PetscSpace, s::PetscInt, subsp::PetscSpace) 
Get a space in the tensor product space

Input Parameters:
- `sp` - the function space object
- `s`  - The space number

Output Parameter:
- `subsp` - the `PetscSpace`

Level: intermediate

-seealso: `PETSCSPACETENSOR`, `PetscSpace`, `PetscSpaceTensorSetSubspace()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpaceTensorGetSubspace"))
"""
function PetscSpaceTensorGetSubspace(petsclib::PetscLibType, sp::PetscSpace, s::PetscInt, subsp::PetscSpace) end

@for_petsc function PetscSpaceTensorGetSubspace(petsclib::$UnionPetscLib, sp::PetscSpace, s::$PetscInt, subsp::PetscSpace )

    @chk ccall(
               (:PetscSpaceTensorGetSubspace, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt, Ptr{PetscSpace}),
               sp, s, subsp,
              )


	return nothing
end 

"""
	PetscSpacePTrimmedSetFormDegree(petsclib::PetscLibType,sp::PetscSpace, formDegree::PetscInt) 
Set the form degree of the trimmed polynomials.

Input Parameters:
- `sp`         - the function space object
- `formDegree` - the form degree

Options Database Key:
- `-petscspace_ptrimmed_form_degree <int>` - The trimmed polynomial form degree

Level: intermediate

-seealso: `PetscSpace`, `PetscDTAltV`, `PetscDTPTrimmedEvalJet()`, `PetscSpacePTrimmedGetFormDegree()`

# External Links
$(_doc_external("DM/PetscSpacePTrimmedSetFormDegree"))
"""
function PetscSpacePTrimmedSetFormDegree(petsclib::PetscLibType, sp::PetscSpace, formDegree::PetscInt) end

@for_petsc function PetscSpacePTrimmedSetFormDegree(petsclib::$UnionPetscLib, sp::PetscSpace, formDegree::$PetscInt )

    @chk ccall(
               (:PetscSpacePTrimmedSetFormDegree, $petsc_library),
               PetscErrorCode,
               (PetscSpace, $PetscInt),
               sp, formDegree,
              )


	return nothing
end 

"""
	formDegree::PetscInt = PetscSpacePTrimmedGetFormDegree(petsclib::PetscLibType,sp::PetscSpace) 
Get the form degree of the trimmed polynomials.

Input Parameter:
- `sp` - the function space object

Output Parameter:
- `formDegree` - the form degree

Level: intermediate

-seealso: `PetscSpace`, `PetscDTAltV`, `PetscDTPTrimmedEvalJet()`, `PetscSpacePTrimmedSetFormDegree()`

# External Links
$(_doc_external("DM/PetscSpacePTrimmedGetFormDegree"))
"""
function PetscSpacePTrimmedGetFormDegree(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpacePTrimmedGetFormDegree(petsclib::$UnionPetscLib, sp::PetscSpace )
	formDegree_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSpacePTrimmedGetFormDegree, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{$PetscInt}),
               sp, formDegree_,
              )

	formDegree = formDegree_[]

	return formDegree
end 

"""
	PetscSpacePolynomialSetTensor(petsclib::PetscLibType,sp::PetscSpace, tensor::PetscBool) 
Set whether a function space is a space of tensor polynomials.

Input Parameters:
- `sp`     - the function space object
- `tensor` - `PETSC_TRUE` for a tensor polynomial space, `PETSC_FALSE` for a polynomial space

Options Database Key:
- `-petscspace_poly_tensor <bool>` - Whether to use tensor product polynomials in higher dimension

Level: intermediate

-seealso: `PetscSpace`, `PetscSpacePolynomialGetTensor()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpacePolynomialSetTensor"))
"""
function PetscSpacePolynomialSetTensor(petsclib::PetscLibType, sp::PetscSpace, tensor::PetscBool) end

@for_petsc function PetscSpacePolynomialSetTensor(petsclib::$UnionPetscLib, sp::PetscSpace, tensor::PetscBool )

    @chk ccall(
               (:PetscSpacePolynomialSetTensor, $petsc_library),
               PetscErrorCode,
               (PetscSpace, PetscBool),
               sp, tensor,
              )


	return nothing
end 

"""
	tensor::PetscBool = PetscSpacePolynomialGetTensor(petsclib::PetscLibType,sp::PetscSpace) 
Get whether a function space is a space of tensor
polynomials.

Input Parameter:
- `sp` - the function space object

Output Parameter:
- `tensor` - `PETSC_TRUE` for a tensor polynomial space, `PETSC_FALSE` for a polynomial space

Level: intermediate

-seealso: `PetscSpace`, `PetscSpacePolynomialSetTensor()`, `PetscSpaceSetDegree()`, `PetscSpaceSetNumVariables()`

# External Links
$(_doc_external("DM/PetscSpacePolynomialGetTensor"))
"""
function PetscSpacePolynomialGetTensor(petsclib::PetscLibType, sp::PetscSpace) end

@for_petsc function PetscSpacePolynomialGetTensor(petsclib::$UnionPetscLib, sp::PetscSpace )
	tensor_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSpacePolynomialGetTensor, $petsc_library),
               PetscErrorCode,
               (PetscSpace, Ptr{PetscBool}),
               sp, tensor_,
              )

	tensor = tensor_[]

	return tensor
end 

