# autodefined type arguments for class ------
mutable struct _n_MatColoring end
const MatColoring = Ptr{_n_MatColoring}

mutable struct _n_MatFDColoring end
const MatFDColoring = Ptr{_n_MatFDColoring}

mutable struct MatFDColoringFn end

mutable struct _n_MatPartitioning end
const MatPartitioning = Ptr{_n_MatPartitioning}

mutable struct MatNullSpaceRemoveFn end

mutable struct MatMFFDFn end

mutable struct MatMFFDiFn end

mutable struct MatMFFDiBaseFn end

mutable struct MatMFFDCheckhFn end

mutable struct _n_MatCoarsen end
const MatCoarsen = Ptr{_n_MatCoarsen}

mutable struct _n_PetscCoarsenData end
const PetscCoarsenData = Ptr{_n_PetscCoarsenData}

# -------------------------------------------------------
"""
	MatColoringPatch(petsclib::PetscLibType,mat::PetscMat, ncolors::PetscInt, n::PetscInt, colorarray::Vector{ISColoringValue}, iscoloring::ISColoring) 
Used inside matrix coloring routines that use `MatGetRowIJ()` and/or
`MatGetColumnIJ()`.

Collective

Input Parameters:
- `mat`        - the matrix
- `ncolors`    - maximum color value
- `n`          - number of entries in colorarray
- `colorarray` - array indicating color for each column

Output Parameter:
- `iscoloring` - coloring generated using colorarray information

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetRowIJ()`, `MatGetColumnIJ()`

# External Links
$(_doc_external("Mat/MatColoringPatch"))
"""
function MatColoringPatch(petsclib::PetscLibType, mat::PetscMat, ncolors::PetscInt, n::PetscInt, colorarray::Vector{ISColoringValue}, iscoloring::ISColoring) end

@for_petsc function MatColoringPatch(petsclib::$UnionPetscLib, mat::PetscMat, ncolors::$PetscInt, n::$PetscInt, colorarray::Vector{ISColoringValue}, iscoloring::ISColoring )

    @chk ccall(
               (:MatColoringPatch, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{ISColoringValue}, Ptr{ISColoring}),
               mat, ncolors, n, colorarray, iscoloring,
              )


	return nothing
end 

"""
	MatColoringRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new sparse matrix coloring to the  matrix package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of Coloring (for example `MATCOLORINGSL`)
- `function` - function pointer that creates the coloring

Level: developer

-seealso: `MatColoringType`, `MatColoringRegisterDestroy()`, `MatColoringRegisterAll()`

# External Links
$(_doc_external("Mat/MatColoringRegister"))
"""
function MatColoringRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatColoringRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatColoringRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	mcptr::MatColoring = MatColoringCreate(petsclib::PetscLibType,m::PetscMat) 
Creates a matrix coloring context.

Collective

Input Parameter:
- `m` - a `Mat` from which a coloring is derived

Output Parameter:
- `mcptr` - the new `MatColoring` context

Options Database Keys:
- `-mat_coloring_type`      - the type of coloring algorithm used. See `MatColoringType`.
- `-mat_coloring_maxcolors` - the maximum number of relevant colors, all nodes not in a color are in maxcolors+1
- `-mat_coloring_distance`  - compute a distance 1,2,... coloring.
- `-mat_coloring_view`      - print information about the coloring and the produced index sets
- `-mat_coloring_test`      - debugging option that prints all coloring incompatibilities
- `-mat_is_coloring_test`   - debugging option that throws an error if MatColoringApply() generates an incorrect iscoloring

Level: beginner

-seealso: `MatColoringSetFromOptions()`, `MatColoring`, `MatColoringApply()`, `MatFDColoringCreate()`, `DMCreateColoring()`, `MatColoringType`

# External Links
$(_doc_external("Mat/MatColoringCreate"))
"""
function MatColoringCreate(petsclib::PetscLibType, m::PetscMat) end

@for_petsc function MatColoringCreate(petsclib::$UnionPetscLib, m::PetscMat )
	mcptr_ = Ref{MatColoring}()

    @chk ccall(
               (:MatColoringCreate, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatColoring}),
               m, mcptr_,
              )

	mcptr = mcptr_[]

	return mcptr
end 

"""
	MatColoringDestroy(petsclib::PetscLibType,mc::MatColoring) 
Destroys the matrix coloring context

Collective

Input Parameter:
- `mc` - the `MatColoring` context

Level: beginner

-seealso: `MatColoring`, `MatColoringCreate()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringDestroy"))
"""
function MatColoringDestroy(petsclib::PetscLibType, mc::MatColoring) end

@for_petsc function MatColoringDestroy(petsclib::$UnionPetscLib, mc::MatColoring )

    @chk ccall(
               (:MatColoringDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatColoring},),
               mc,
              )


	return nothing
end 

"""
	MatColoringSetType(petsclib::PetscLibType,mc::MatColoring, type::MatColoringType) 
Sets the type of coloring algorithm used

Collective

Input Parameters:
- `mc`   - the `MatColoring` context
- `type` - the type of coloring

Options Database Key:
- `-mat_coloring_type type` - the name of the type

Level: beginner

-seealso: `MatColoring`, `MatColoringSetFromOptions()`, `MatColoringType`, `MatColoringCreate()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringSetType"))
"""
function MatColoringSetType(petsclib::PetscLibType, mc::MatColoring, type::MatColoringType) end

@for_petsc function MatColoringSetType(petsclib::$UnionPetscLib, mc::MatColoring, type::MatColoringType )

    @chk ccall(
               (:MatColoringSetType, $petsc_library),
               PetscErrorCode,
               (MatColoring, MatColoringType),
               mc, type,
              )


	return nothing
end 

"""
	MatColoringSetFromOptions(petsclib::PetscLibType,mc::MatColoring) 
Sets `MatColoring` options from options database

Collective

Input Parameter:
- `mc` - `MatColoring` context

Options Database Keys:
- `-mat_coloring_type`      - the type of coloring algorithm used. See `MatColoringType`.
- `-mat_coloring_maxcolors` - the maximum number of relevant colors, all nodes not in a color are in maxcolors+1
- `-mat_coloring_distance`  - compute a distance 1,2,... coloring.
- `-mat_coloring_view`      - print information about the coloring and the produced index sets
- `-snes_fd_color`          - instruct SNES to using coloring and then `MatFDColoring` to compute the Jacobians
- `-snes_fd_color_use_mat`  - instruct `SNES` to color the matrix directly instead of the `DM` from which the matrix comes (the default)

Level: beginner

-seealso: `MatColoring`, `MatColoringApply()`, `MatColoringSetDistance()`, `MatColoringSetType()`, `SNESComputeJacobianDefaultColor()`, `MatColoringType`

# External Links
$(_doc_external("Mat/MatColoringSetFromOptions"))
"""
function MatColoringSetFromOptions(petsclib::PetscLibType, mc::MatColoring) end

@for_petsc function MatColoringSetFromOptions(petsclib::$UnionPetscLib, mc::MatColoring )

    @chk ccall(
               (:MatColoringSetFromOptions, $petsc_library),
               PetscErrorCode,
               (MatColoring,),
               mc,
              )


	return nothing
end 

"""
	MatColoringSetDistance(petsclib::PetscLibType,mc::MatColoring, dist::PetscInt) 
Sets the distance of the coloring

Logically Collective

Input Parameters:
- `mc`   - the `MatColoring` context
- `dist` - the distance the coloring should compute

Options Database Key:
- `-mat_coloring_type` - the type of coloring algorithm used. See `MatColoringType`.

Level: beginner

-seealso: `MatColoring`, `MatColoringSetFromOptions()`, `MatColoringGetDistance()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringSetDistance"))
"""
function MatColoringSetDistance(petsclib::PetscLibType, mc::MatColoring, dist::PetscInt) end

@for_petsc function MatColoringSetDistance(petsclib::$UnionPetscLib, mc::MatColoring, dist::$PetscInt )

    @chk ccall(
               (:MatColoringSetDistance, $petsc_library),
               PetscErrorCode,
               (MatColoring, $PetscInt),
               mc, dist,
              )


	return nothing
end 

"""
	dist::PetscInt = MatColoringGetDistance(petsclib::PetscLibType,mc::MatColoring) 
Gets the distance of the coloring

Logically Collective

Input Parameter:
- `mc` - the `MatColoring` context

Output Parameter:
- `dist` - the current distance being used for the coloring.

Level: beginner

-seealso: `MatColoring`, `MatColoringSetDistance()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringGetDistance"))
"""
function MatColoringGetDistance(petsclib::PetscLibType, mc::MatColoring) end

@for_petsc function MatColoringGetDistance(petsclib::$UnionPetscLib, mc::MatColoring )
	dist_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatColoringGetDistance, $petsc_library),
               PetscErrorCode,
               (MatColoring, Ptr{$PetscInt}),
               mc, dist_,
              )

	dist = dist_[]

	return dist
end 

"""
	MatColoringSetMaxColors(petsclib::PetscLibType,mc::MatColoring, maxcolors::PetscInt) 
Sets the maximum number of colors to produce

Logically Collective

Input Parameters:
- `mc`        - the `MatColoring` context
- `maxcolors` - the maximum number of colors to produce

Level: beginner

-seealso: `MatColoring`, `MatColoringGetMaxColors()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringSetMaxColors"))
"""
function MatColoringSetMaxColors(petsclib::PetscLibType, mc::MatColoring, maxcolors::PetscInt) end

@for_petsc function MatColoringSetMaxColors(petsclib::$UnionPetscLib, mc::MatColoring, maxcolors::$PetscInt )

    @chk ccall(
               (:MatColoringSetMaxColors, $petsc_library),
               PetscErrorCode,
               (MatColoring, $PetscInt),
               mc, maxcolors,
              )


	return nothing
end 

"""
	maxcolors::PetscInt = MatColoringGetMaxColors(petsclib::PetscLibType,mc::MatColoring) 
Gets the maximum number of colors

Logically Collective

Input Parameter:
- `mc` - the `MatColoring` context

Output Parameter:
- `maxcolors` - the current maximum number of colors to produce

Level: beginner

-seealso: `MatColoring`, `MatColoringSetMaxColors()`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringGetMaxColors"))
"""
function MatColoringGetMaxColors(petsclib::PetscLibType, mc::MatColoring) end

@for_petsc function MatColoringGetMaxColors(petsclib::$UnionPetscLib, mc::MatColoring )
	maxcolors_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatColoringGetMaxColors, $petsc_library),
               PetscErrorCode,
               (MatColoring, Ptr{$PetscInt}),
               mc, maxcolors_,
              )

	maxcolors = maxcolors_[]

	return maxcolors
end 

"""
	MatColoringApply(petsclib::PetscLibType,mc::MatColoring, coloring::ISColoring) 
Apply the coloring to the matrix, producing index
sets corresponding to a number of independent sets in the induced
graph.

Collective

Input Parameter:
- `mc` - the `MatColoring` context

Output Parameter:
- `coloring` - the `ISColoring` instance containing the coloring

Level: beginner

-seealso: `ISColoring`, `MatColoring`, `MatColoringCreate()`

# External Links
$(_doc_external("Mat/MatColoringApply"))
"""
function MatColoringApply(petsclib::PetscLibType, mc::MatColoring, coloring::ISColoring) end

@for_petsc function MatColoringApply(petsclib::$UnionPetscLib, mc::MatColoring, coloring::ISColoring )

    @chk ccall(
               (:MatColoringApply, $petsc_library),
               PetscErrorCode,
               (MatColoring, Ptr{ISColoring}),
               mc, coloring,
              )


	return nothing
end 

"""
	MatColoringView(petsclib::PetscLibType,mc::MatColoring, viewer::PetscViewer) 
Output details about the `MatColoring`.

Collective

Input Parameters:
- `mc`     - the `MatColoring` context
- `viewer` - the Viewer context

Level: beginner

-seealso: `PetscViewer`, `MatColoring`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringView"))
"""
function MatColoringView(petsclib::PetscLibType, mc::MatColoring, viewer::PetscViewer) end

@for_petsc function MatColoringView(petsclib::$UnionPetscLib, mc::MatColoring, viewer::PetscViewer )

    @chk ccall(
               (:MatColoringView, $petsc_library),
               PetscErrorCode,
               (MatColoring, PetscViewer),
               mc, viewer,
              )


	return nothing
end 

"""
	MatColoringSetWeightType(petsclib::PetscLibType,mc::MatColoring, wt::MatColoringWeightType) 
Set the type of weight computation used while computing the coloring

Logically Collective

Input Parameters:
- `mc` - the `MatColoring` context
- `wt` - the weight type

Level: beginner

-seealso: `MatColoring`, `MatColoringWeightType`, `MatColoringApply()`

# External Links
$(_doc_external("Mat/MatColoringSetWeightType"))
"""
function MatColoringSetWeightType(petsclib::PetscLibType, mc::MatColoring, wt::MatColoringWeightType) end

@for_petsc function MatColoringSetWeightType(petsclib::$UnionPetscLib, mc::MatColoring, wt::MatColoringWeightType )

    @chk ccall(
               (:MatColoringSetWeightType, $petsc_library),
               PetscErrorCode,
               (MatColoring, MatColoringWeightType),
               mc, wt,
              )


	return nothing
end 

"""
	degrees::PetscInt = MatColoringGetDegrees(petsclib::PetscLibType,G::PetscMat, distance::PetscInt) 

# External Links
$(_doc_external("Mat/MatColoringGetDegrees"))
"""
function MatColoringGetDegrees(petsclib::PetscLibType, G::PetscMat, distance::PetscInt) end

@for_petsc function MatColoringGetDegrees(petsclib::$UnionPetscLib, G::PetscMat, distance::$PetscInt )
	degrees_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatColoringGetDegrees, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               G, distance, degrees_,
              )

	degrees = degrees_[]

	return degrees
end 

"""
	MatColoringCreateWeights(petsclib::PetscLibType,mc::MatColoring, weights::PetscReal, lperm::PetscInt) 

# External Links
$(_doc_external("Mat/MatColoringCreateWeights"))
"""
function MatColoringCreateWeights(petsclib::PetscLibType, mc::MatColoring, weights::PetscReal, lperm::PetscInt) end

@for_petsc function MatColoringCreateWeights(petsclib::$UnionPetscLib, mc::MatColoring, weights::$PetscReal, lperm::$PetscInt )

    @chk ccall(
               (:MatColoringCreateWeights, $petsc_library),
               PetscErrorCode,
               (MatColoring, $PetscReal, $PetscInt),
               mc, weights, lperm,
              )


	return nothing
end 

"""
	weights::PetscReal,lperm::PetscInt = MatColoringSetWeights(petsclib::PetscLibType,mc::MatColoring) 

# External Links
$(_doc_external("Mat/MatColoringSetWeights"))
"""
function MatColoringSetWeights(petsclib::PetscLibType, mc::MatColoring) end

@for_petsc function MatColoringSetWeights(petsclib::$UnionPetscLib, mc::MatColoring )
	weights_ = Ref{$PetscReal}()
	lperm_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatColoringSetWeights, $petsc_library),
               PetscErrorCode,
               (MatColoring, Ptr{$PetscReal}, Ptr{$PetscInt}),
               mc, weights_, lperm_,
              )

	weights = weights_[]
	lperm = lperm_[]

	return weights,lperm
end 

"""
	MatFDColoringSetF(petsclib::PetscLibType,fd::MatFDColoring, F::PetscVec) 

# External Links
$(_doc_external("Mat/MatFDColoringSetF"))
"""
function MatFDColoringSetF(petsclib::PetscLibType, fd::MatFDColoring, F::PetscVec) end

@for_petsc function MatFDColoringSetF(petsclib::$UnionPetscLib, fd::MatFDColoring, F::PetscVec )

    @chk ccall(
               (:MatFDColoringSetF, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, CVec),
               fd, F,
              )


	return nothing
end 

"""
	MatFDColoringView(petsclib::PetscLibType,c::MatFDColoring, viewer::PetscViewer) 
Views a finite difference coloring context.

Collective

Input Parameters:
- `c`      - the coloring context
- `viewer` - visualization context

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`

# External Links
$(_doc_external("Mat/MatFDColoringView"))
"""
function MatFDColoringView(petsclib::PetscLibType, c::MatFDColoring, viewer::PetscViewer) end

@for_petsc function MatFDColoringView(petsclib::$UnionPetscLib, c::MatFDColoring, viewer::PetscViewer )

    @chk ccall(
               (:MatFDColoringView, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, PetscViewer),
               c, viewer,
              )


	return nothing
end 

"""
	MatFDColoringSetParameters(petsclib::PetscLibType,matfd::MatFDColoring, error::PetscReal, umin::PetscReal) 
Sets the parameters for the approximation of
a sparse Jacobian matrix using finite differences and matrix coloring

Logically Collective

Input Parameters:
- `matfd` - the coloring context
- `error` - relative error
- `umin`  - minimum allowable u-value magnitude

Level: advanced

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringSetFromOptions()`

# External Links
$(_doc_external("Mat/MatFDColoringSetParameters"))
"""
function MatFDColoringSetParameters(petsclib::PetscLibType, matfd::MatFDColoring, error::PetscReal, umin::PetscReal) end

@for_petsc function MatFDColoringSetParameters(petsclib::$UnionPetscLib, matfd::MatFDColoring, error::$PetscReal, umin::$PetscReal )

    @chk ccall(
               (:MatFDColoringSetParameters, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, $PetscReal, $PetscReal),
               matfd, error, umin,
              )


	return nothing
end 

"""
	MatFDColoringSetBlockSize(petsclib::PetscLibType,matfd::MatFDColoring, brows::PetscInt, bcols::PetscInt) 
Sets block size for efficient inserting entries of Jacobian matrix.

Logically Collective

Input Parameters:
- `matfd` - the coloring context
- `brows` - number of rows in the block
- `bcols` - number of columns in the block

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringSetFromOptions()`

# External Links
$(_doc_external("Mat/MatFDColoringSetBlockSize"))
"""
function MatFDColoringSetBlockSize(petsclib::PetscLibType, matfd::MatFDColoring, brows::PetscInt, bcols::PetscInt) end

@for_petsc function MatFDColoringSetBlockSize(petsclib::$UnionPetscLib, matfd::MatFDColoring, brows::$PetscInt, bcols::$PetscInt )

    @chk ccall(
               (:MatFDColoringSetBlockSize, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, $PetscInt, $PetscInt),
               matfd, brows, bcols,
              )


	return nothing
end 

"""
	MatFDColoringSetUp(petsclib::PetscLibType,mat::PetscMat, iscoloring::ISColoring, color::MatFDColoring) 
Sets up the internal data structures of matrix coloring context for the later use.

Collective

Input Parameters:
- `mat`        - the matrix containing the nonzero structure of the Jacobian
- `iscoloring` - the coloring of the matrix; usually obtained with `MatGetColoring()` or `DMCreateColoring()`
- `color`      - the matrix coloring context

Level: beginner

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringDestroy()`

# External Links
$(_doc_external("Mat/MatFDColoringSetUp"))
"""
function MatFDColoringSetUp(petsclib::PetscLibType, mat::PetscMat, iscoloring::ISColoring, color::MatFDColoring) end

@for_petsc function MatFDColoringSetUp(petsclib::$UnionPetscLib, mat::PetscMat, iscoloring::ISColoring, color::MatFDColoring )

    @chk ccall(
               (:MatFDColoringSetUp, $petsc_library),
               PetscErrorCode,
               (CMat, ISColoring, MatFDColoring),
               mat, iscoloring, color,
              )


	return nothing
end 

"""
	MatFDColoringGetFunction(petsclib::PetscLibType,matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid) 
Gets the function to use for computing the Jacobian.

Not Collective

Input Parameter:
- `matfd` - the coloring context

Output Parameters:
- `f`    - the function, see `MatFDColoringFn` for the calling sequence
- `fctx` - the optional user-defined function context

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringSetFunction()`, `MatFDColoringSetFromOptions()`, `MatFDColoringFn`

# External Links
$(_doc_external("Mat/MatFDColoringGetFunction"))
"""
function MatFDColoringGetFunction(petsclib::PetscLibType, matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid) end

@for_petsc function MatFDColoringGetFunction(petsclib::$UnionPetscLib, matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid )

    @chk ccall(
               (:MatFDColoringGetFunction, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, MatFDColoringFn, Cvoid),
               matfd, f, fctx,
              )


	return nothing
end 

"""
	MatFDColoringSetFunction(petsclib::PetscLibType,matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid) 
Sets the function to use for computing the Jacobian.

Logically Collective

Input Parameters:
- `matfd` - the coloring context
- `f`     - the function, see `MatFDColoringFn` for the calling sequence
- `fctx`  - the optional user-defined function context

Level: advanced

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringGetFunction()`, `MatFDColoringSetFromOptions()`, `MatFDColoringFn`

# External Links
$(_doc_external("Mat/MatFDColoringSetFunction"))
"""
function MatFDColoringSetFunction(petsclib::PetscLibType, matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid) end

@for_petsc function MatFDColoringSetFunction(petsclib::$UnionPetscLib, matfd::MatFDColoring, f::MatFDColoringFn, fctx::Cvoid )

    @chk ccall(
               (:MatFDColoringSetFunction, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, Ptr{MatFDColoringFn}, Ptr{Cvoid}),
               matfd, f, fctx,
              )


	return nothing
end 

"""
	MatFDColoringSetFromOptions(petsclib::PetscLibType,matfd::MatFDColoring) 
Sets coloring finite difference parameters from
the options database.

Collective

The Jacobian, F'(u), is estimated with the differencing approximation
-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringView()`, `MatFDColoringSetParameters()`

# External Links
$(_doc_external("Mat/MatFDColoringSetFromOptions"))
"""
function MatFDColoringSetFromOptions(petsclib::PetscLibType, matfd::MatFDColoring) end

@for_petsc function MatFDColoringSetFromOptions(petsclib::$UnionPetscLib, matfd::MatFDColoring )

    @chk ccall(
               (:MatFDColoringSetFromOptions, $petsc_library),
               PetscErrorCode,
               (MatFDColoring,),
               matfd,
              )


	return nothing
end 

"""
	MatFDColoringSetType(petsclib::PetscLibType,matfd::MatFDColoring, type::MatMFFDType) 
Sets the approach for computing the finite difference parameter

Collective

Input Parameters:
- `matfd` - the coloring context
- `type`  - either `MATMFFD_WP` or `MATMFFD_DS`

Options Database Key:
- `-mat_fd_type` - "wp" or "ds"

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringView()`, `MatFDColoringSetParameters()`

# External Links
$(_doc_external("Mat/MatFDColoringSetType"))
"""
function MatFDColoringSetType(petsclib::PetscLibType, matfd::MatFDColoring, type::MatMFFDType) end

@for_petsc function MatFDColoringSetType(petsclib::$UnionPetscLib, matfd::MatFDColoring, type::MatMFFDType )

    @chk ccall(
               (:MatFDColoringSetType, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, MatMFFDType),
               matfd, type,
              )


	return nothing
end 

"""
	color::MatFDColoring = MatFDColoringCreate(petsclib::PetscLibType,mat::PetscMat, iscoloring::ISColoring) 
Creates a matrix coloring context for finite difference
computation of Jacobians.

Collective

Input Parameters:
- `mat`        - the matrix containing the nonzero structure of the Jacobian
- `iscoloring` - the coloring of the matrix; usually obtained with `MatColoringCreate()` or `DMCreateColoring()`

Output Parameter:
- `color` - the new coloring context

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringDestroy()`, `SNESComputeJacobianDefaultColor()`, `ISColoringCreate()`,
`MatFDColoringSetFunction()`, `MatFDColoringSetFromOptions()`, `MatFDColoringApply()`,
`MatFDColoringView()`, `MatFDColoringSetParameters()`, `MatColoringCreate()`, `DMCreateColoring()`, `MatFDColoringSetValues()`

# External Links
$(_doc_external("Mat/MatFDColoringCreate"))
"""
function MatFDColoringCreate(petsclib::PetscLibType, mat::PetscMat, iscoloring::ISColoring) end

@for_petsc function MatFDColoringCreate(petsclib::$UnionPetscLib, mat::PetscMat, iscoloring::ISColoring )
	color_ = Ref{MatFDColoring}()

    @chk ccall(
               (:MatFDColoringCreate, $petsc_library),
               PetscErrorCode,
               (CMat, ISColoring, Ptr{MatFDColoring}),
               mat, iscoloring, color_,
              )

	color = color_[]

	return color
end 

"""
	MatFDColoringDestroy(petsclib::PetscLibType,c::MatFDColoring) 
Destroys a matrix coloring context that was created
via `MatFDColoringCreate()`.

Collective

Input Parameter:
- `c` - coloring context

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`

# External Links
$(_doc_external("Mat/MatFDColoringDestroy"))
"""
function MatFDColoringDestroy(petsclib::PetscLibType, c::MatFDColoring) end

@for_petsc function MatFDColoringDestroy(petsclib::$UnionPetscLib, c::MatFDColoring )

    @chk ccall(
               (:MatFDColoringDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatFDColoring},),
               c,
              )


	return nothing
end 

"""
	n::PetscInt,cols::Vector{PetscInt} = MatFDColoringGetPerturbedColumns(petsclib::PetscLibType,coloring::MatFDColoring) 
Returns the indices of the columns that
that are currently being perturbed.

Not Collective

Input Parameter:
- `coloring` - coloring context created with `MatFDColoringCreate()`

Output Parameters:
- `n`    - the number of local columns being perturbed
- `cols` - the column indices, in global numbering

Level: advanced

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringDestroy()`, `MatFDColoringView()`, `MatFDColoringApply()`

# External Links
$(_doc_external("Mat/MatFDColoringGetPerturbedColumns"))
"""
function MatFDColoringGetPerturbedColumns(petsclib::PetscLibType, coloring::MatFDColoring) end

@for_petsc function MatFDColoringGetPerturbedColumns(petsclib::$UnionPetscLib, coloring::MatFDColoring )
	n_ = Ref{$PetscInt}()
	cols_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatFDColoringGetPerturbedColumns, $petsc_library),
               PetscErrorCode,
               (MatFDColoring, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               coloring, n_, cols_,
              )

	n = n_[]
	cols = unsafe_wrap(Array, cols_[], VecGetLocalSize(petsclib, x); own = false)

	return n,cols
end 

"""
	MatFDColoringApply(petsclib::PetscLibType,J::PetscMat, coloring::MatFDColoring, x1::PetscVec, sctx::Cvoid) 
Given a matrix for which a `MatFDColoring` context
has been created, computes the Jacobian for a function via finite differences.

Collective

Input Parameters:
- `J`        - matrix to store Jacobian entries into
- `coloring` - coloring context created with `MatFDColoringCreate()`
- `x1`       - location at which Jacobian is to be computed
- `sctx`     - context required by function, if this is being used with the `SNES` solver then it is `SNES` object, otherwise it is `NULL`

Options Database Keys:
- `-mat_fd_type`                       - "wp" or "ds"  (see `MATMFFD_WP` or `MATMFFD_DS`)
- `-mat_fd_coloring_view`              - Activates basic viewing or coloring
- `-mat_fd_coloring_view draw`         - Activates drawing of coloring
- `-mat_fd_coloring_view ::ascii_info` - Activates viewing of coloring info

Level: intermediate

-seealso: `Mat`, `MatFDColoring`, `MatFDColoringCreate()`, `MatFDColoringDestroy()`, `MatFDColoringView()`, `MatFDColoringSetFunction()`, `MatFDColoringSetValues()`

# External Links
$(_doc_external("Mat/MatFDColoringApply"))
"""
function MatFDColoringApply(petsclib::PetscLibType, J::PetscMat, coloring::MatFDColoring, x1::PetscVec, sctx::Cvoid) end

@for_petsc function MatFDColoringApply(petsclib::$UnionPetscLib, J::PetscMat, coloring::MatFDColoring, x1::PetscVec, sctx::Cvoid )

    @chk ccall(
               (:MatFDColoringApply, $petsc_library),
               PetscErrorCode,
               (CMat, MatFDColoring, CVec, Ptr{Cvoid}),
               J, coloring, x1, sctx,
              )


	return nothing
end 

"""
	MatFDColoringSetValues(petsclib::PetscLibType,J::PetscMat, coloring::MatFDColoring, y::Vector{PetscScalar}) 
takes a matrix in compressed color format and enters the matrix into a PETSc `Mat`

Collective

Input Parameters:
- `J`        - the sparse matrix
- `coloring` - created with `MatFDColoringCreate()` and a local coloring
- `y`        - column major storage of matrix values with one color of values per column, the number of rows of `y` should match
the number of local rows of `J` and the number of columns is the number of colors.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatFDColoringCreate()`, `ISColoring`, `ISColoringCreate()`, `ISColoringSetType()`, `IS_COLORING_LOCAL`, `MatFDColoringSetBlockSize()`

# External Links
$(_doc_external("Mat/MatFDColoringSetValues"))
"""
function MatFDColoringSetValues(petsclib::PetscLibType, J::PetscMat, coloring::MatFDColoring, y::Vector{PetscScalar}) end

@for_petsc function MatFDColoringSetValues(petsclib::$UnionPetscLib, J::PetscMat, coloring::MatFDColoring, y::Vector{$PetscScalar} )

    @chk ccall(
               (:MatFDColoringSetValues, $petsc_library),
               PetscErrorCode,
               (CMat, MatFDColoring, Ptr{$PetscScalar}),
               J, coloring, y,
              )


	return nothing
end 

"""
	MatFDColoringUseDM(petsclib::PetscLibType,coloring::PetscMat, fdcoloring::MatFDColoring) 
allows a `MatFDColoring` object to use the `DM` associated with the matrix to compute a `IS_COLORING_LOCAL` coloring

Input Parameters:
- `coloring`   - The matrix to get the `DM` from
- `fdcoloring` - the `MatFDColoring` object

Level: advanced

Developer Note:
This routine exists because the PETSc `Mat` library does not know about the `DM` objects

See also: 
=== 
`DM`, `MatFDColoring`, `MatFDColoringCreate()`, `ISColoringType`

# External Links
$(_doc_external("Dm/MatFDColoringUseDM"))
"""
function MatFDColoringUseDM(petsclib::PetscLibType, coloring::PetscMat, fdcoloring::MatFDColoring) end

@for_petsc function MatFDColoringUseDM(petsclib::$UnionPetscLib, coloring::PetscMat, fdcoloring::MatFDColoring )

    @chk ccall(
               (:MatFDColoringUseDM, $petsc_library),
               PetscErrorCode,
               (CMat, MatFDColoring),
               coloring, fdcoloring,
              )


	return nothing
end 

"""
	MatTransposeColoringDestroy(petsclib::PetscLibType,c::MatTransposeColoring) 
Destroys a coloring context for matrix product C = A*B^T that was created
via `MatTransposeColoringCreate()`.

Collective

Input Parameter:
- `c` - coloring context

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTransposeColoringCreate()`

# External Links
$(_doc_external("Mat/MatTransposeColoringDestroy"))
"""
function MatTransposeColoringDestroy(petsclib::PetscLibType, c::MatTransposeColoring) end

@for_petsc function MatTransposeColoringDestroy(petsclib::$UnionPetscLib, c::MatTransposeColoring )

    @chk ccall(
               (:MatTransposeColoringDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatTransposeColoring},),
               c,
              )


	return nothing
end 

"""
	color::MatTransposeColoring = MatTransposeColoringCreate(petsclib::PetscLibType,mat::PetscMat, iscoloring::ISColoring) 
Creates a matrix coloring context for the matrix product C = A*B^T.

Collective

Input Parameters:
- `mat`        - the matrix product C
- `iscoloring` - the coloring of the matrix; usually obtained with `MatColoringCreate()` or `DMCreateColoring()`

Output Parameter:
- `color` - the new coloring context

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTransposeColoringDestroy()`, `MatTransColoringApplySpToDen()`,
`MatTransColoringApplyDenToSp()`

# External Links
$(_doc_external("Mat/MatTransposeColoringCreate"))
"""
function MatTransposeColoringCreate(petsclib::PetscLibType, mat::PetscMat, iscoloring::ISColoring) end

@for_petsc function MatTransposeColoringCreate(petsclib::$UnionPetscLib, mat::PetscMat, iscoloring::ISColoring )
	color_ = Ref{MatTransposeColoring}()

    @chk ccall(
               (:MatTransposeColoringCreate, $petsc_library),
               PetscErrorCode,
               (CMat, ISColoring, Ptr{MatTransposeColoring}),
               mat, iscoloring, color_,
              )

	color = color_[]

	return color
end 

"""
	MatPartitioningRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new sparse matrix partitioning to the  matrix package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of partitioning (for example `MATPARTITIONINGCURRENT`) or `MATPARTITIONINGPARMETIS`
- `function` - function pointer that creates the partitioning type

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningCreate()`, `MatPartitioningRegisterDestroy()`, `MatPartitioningRegisterAll()`

# External Links
$(_doc_external("Mat/MatPartitioningRegister"))
"""
function MatPartitioningRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatPartitioningRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatPartitioningRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	type::MatPartitioningType = MatPartitioningGetType(petsclib::PetscLibType,partitioning::MatPartitioning) 
Gets the Partitioning method type and name (as a string)
from the partitioning context.

Not Collective

Input Parameter:
- `partitioning` - the partitioning context

Output Parameter:
- `type` - partitioner type

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningCreate()`, `MatPartitioningRegisterDestroy()`, `MatPartitioningRegisterAll()`

# External Links
$(_doc_external("Mat/MatPartitioningGetType"))
"""
function MatPartitioningGetType(petsclib::PetscLibType, partitioning::MatPartitioning) end

@for_petsc function MatPartitioningGetType(petsclib::$UnionPetscLib, partitioning::MatPartitioning )
	type_ = Ref{MatPartitioningType}()

    @chk ccall(
               (:MatPartitioningGetType, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{MatPartitioningType}),
               partitioning, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	MatPartitioningSetNParts(petsclib::PetscLibType,part::MatPartitioning, n::PetscInt) 
Set how many partitions need to be created;
by default this is one per processor. Certain partitioning schemes may
in fact only support that option.

Collective

Input Parameters:
- `part` - the partitioning context
- `n`    - the number of partitions

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningCreate()`, `MatPartitioningApply()`

# External Links
$(_doc_external("Mat/MatPartitioningSetNParts"))
"""
function MatPartitioningSetNParts(petsclib::PetscLibType, part::MatPartitioning, n::PetscInt) end

@for_petsc function MatPartitioningSetNParts(petsclib::$UnionPetscLib, part::MatPartitioning, n::$PetscInt )

    @chk ccall(
               (:MatPartitioningSetNParts, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscInt),
               part, n,
              )


	return nothing
end 

"""
	MatPartitioningApplyND(petsclib::PetscLibType,matp::MatPartitioning, partitioning::IS) 
Gets a nested dissection partitioning for a matrix.

Collective

Input Parameter:
- `matp` - the matrix partitioning object

Output Parameter:
- `partitioning` - the partitioning. For each local node, a positive value indicates the processor
number the node has been assigned to. Negative x values indicate the separator level -(x+1).

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioningRegister()`, `MatPartitioningCreate()`,
`MatPartitioningDestroy()`, `MatPartitioningSetAdjacency()`, `ISPartitioningToNumbering()`,
`ISPartitioningCount()`

# External Links
$(_doc_external("Mat/MatPartitioningApplyND"))
"""
function MatPartitioningApplyND(petsclib::PetscLibType, matp::MatPartitioning, partitioning::IS) end

@for_petsc function MatPartitioningApplyND(petsclib::$UnionPetscLib, matp::MatPartitioning, partitioning::IS )
	partitioning_ = Ref(partitioning.ptr)

    @chk ccall(
               (:MatPartitioningApplyND, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{CIS}),
               matp, partitioning_,
              )

	partitioning.ptr = C_NULL

	return nothing
end 

"""
	MatPartitioningApply(petsclib::PetscLibType,matp::MatPartitioning, partitioning::IS) 
Gets a partitioning for the graph represented by a sparse matrix.

Collective

Input Parameter:
- `matp` - the matrix partitioning object

Output Parameter:
- `partitioning` - the partitioning. For each local node this tells the MPI rank that that node is assigned to.

Options Database Keys:
- `-mat_partitioning_type <type>` - set the partitioning package or algorithm to use
- `-mat_partitioning_view`        - display information about the partitioning object

Level: beginner

The user can define additional partitionings; see `MatPartitioningRegister()`.

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningRegister()`, `MatPartitioningCreate()`,
`MatPartitioningDestroy()`, `MatPartitioningSetAdjacency()`, `ISPartitioningToNumbering()`,
`ISPartitioningCount()`

# External Links
$(_doc_external("Mat/MatPartitioningApply"))
"""
function MatPartitioningApply(petsclib::PetscLibType, matp::MatPartitioning, partitioning::IS) end

@for_petsc function MatPartitioningApply(petsclib::$UnionPetscLib, matp::MatPartitioning, partitioning::IS )
	partitioning_ = Ref(partitioning.ptr)

    @chk ccall(
               (:MatPartitioningApply, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{CIS}),
               matp, partitioning_,
              )

	partitioning.ptr = C_NULL

	return nothing
end 

"""
	MatPartitioningImprove(petsclib::PetscLibType,matp::MatPartitioning, partitioning::IS) 
Improves the quality of a given partition.

Collective

Input Parameters:
- `matp`         - the matrix partitioning object
- `partitioning` - the original partitioning. For each local node this tells the processor
number that that node is assigned to.

Options Database Key:
- `-mat_partitioning_improve` - improve the quality of the given partition

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningApply()`, `MatPartitioningCreate()`,
`MatPartitioningDestroy()`, `MatPartitioningSetAdjacency()`, `ISPartitioningToNumbering()`,
`ISPartitioningCount()`

# External Links
$(_doc_external("Mat/MatPartitioningImprove"))
"""
function MatPartitioningImprove(petsclib::PetscLibType, matp::MatPartitioning, partitioning::IS) end

@for_petsc function MatPartitioningImprove(petsclib::$UnionPetscLib, matp::MatPartitioning, partitioning::IS )
	partitioning_ = Ref(partitioning.ptr)

    @chk ccall(
               (:MatPartitioningImprove, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{CIS}),
               matp, partitioning_,
              )

	partitioning.ptr = C_NULL

	return nothing
end 

"""
	MatPartitioningViewImbalance(petsclib::PetscLibType,matp::MatPartitioning, partitioning::IS) 
Display partitioning imbalance information.

Collective

Input Parameters:
- `matp`         - the matrix partitioning object
- `partitioning` - the partitioning. For each local node this tells the MPI rank that that node is assigned to.

Options Database Key:
- `-mat_partitioning_view_balance` - view the balance information from the last partitioning

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningApply()`, `MatPartitioningView()`

# External Links
$(_doc_external("Mat/MatPartitioningViewImbalance"))
"""
function MatPartitioningViewImbalance(petsclib::PetscLibType, matp::MatPartitioning, partitioning::IS) end

@for_petsc function MatPartitioningViewImbalance(petsclib::$UnionPetscLib, matp::MatPartitioning, partitioning::IS )

    @chk ccall(
               (:MatPartitioningViewImbalance, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, CIS),
               matp, partitioning,
              )


	return nothing
end 

"""
	MatPartitioningSetAdjacency(petsclib::PetscLibType,part::MatPartitioning, adj::PetscMat) 
Sets the adjacency graph (matrix) of the thing to be
partitioned.

Collective

Input Parameters:
- `part` - the partitioning context
- `adj`  - the adjacency matrix, this can be any `MatType` but the natural representation is `MATMPIADJ`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningCreate()`

# External Links
$(_doc_external("Mat/MatPartitioningSetAdjacency"))
"""
function MatPartitioningSetAdjacency(petsclib::PetscLibType, part::MatPartitioning, adj::PetscMat) end

@for_petsc function MatPartitioningSetAdjacency(petsclib::$UnionPetscLib, part::MatPartitioning, adj::PetscMat )

    @chk ccall(
               (:MatPartitioningSetAdjacency, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, CMat),
               part, adj,
              )


	return nothing
end 

"""
	MatPartitioningDestroy(petsclib::PetscLibType,part::MatPartitioning) 
Destroys the partitioning context.

Collective

Input Parameter:
- `part` - the partitioning context

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningType`, `MatPartitioningCreate()`

# External Links
$(_doc_external("Mat/MatPartitioningDestroy"))
"""
function MatPartitioningDestroy(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningDestroy(petsclib::$UnionPetscLib, part::MatPartitioning )

    @chk ccall(
               (:MatPartitioningDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatPartitioning},),
               part,
              )


	return nothing
end 

"""
	MatPartitioningSetVertexWeights(petsclib::PetscLibType,part::MatPartitioning, weights::Vector{PetscInt}) 
Sets the weights for vertices for a partitioning.

Logically Collective

Input Parameters:
- `part`    - the partitioning context
- `weights` - the weights, on each process this array must have the same size as the number of local rows times the value passed with `MatPartitioningSetNumberVertexWeights()` or
1 if that is not provided

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningCreate()`, `MatPartitioningSetType()`, `MatPartitioningSetPartitionWeights()`, `MatPartitioningSetNumberVertexWeights()`

# External Links
$(_doc_external("Mat/MatPartitioningSetVertexWeights"))
"""
function MatPartitioningSetVertexWeights(petsclib::PetscLibType, part::MatPartitioning, weights::Vector{PetscInt}) end

@for_petsc function MatPartitioningSetVertexWeights(petsclib::$UnionPetscLib, part::MatPartitioning, weights::Vector{$PetscInt} )

    @chk ccall(
               (:MatPartitioningSetVertexWeights, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscInt}),
               part, weights,
              )


	return nothing
end 

"""
	MatPartitioningSetPartitionWeights(petsclib::PetscLibType,part::MatPartitioning, weights::Vector{PetscReal}) 
Sets the weights for each partition.

Logically Collective

Input Parameters:
- `part`    - the partitioning context
- `weights` - An array of size nparts that is used to specify the fraction of
vertex weight that should be distributed to each sub-domain for
the balance constraint. If all of the sub-domains are to be of
the same size, then each of the nparts elements should be set
to a value of 1/nparts. Note that the sum of all of the weights
should be one.

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningSetVertexWeights()`, `MatPartitioningCreate()`, `MatPartitioningSetType()`

# External Links
$(_doc_external("Mat/MatPartitioningSetPartitionWeights"))
"""
function MatPartitioningSetPartitionWeights(petsclib::PetscLibType, part::MatPartitioning, weights::Vector{PetscReal}) end

@for_petsc function MatPartitioningSetPartitionWeights(petsclib::$UnionPetscLib, part::MatPartitioning, weights::Vector{$PetscReal} )

    @chk ccall(
               (:MatPartitioningSetPartitionWeights, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscReal}),
               part, weights,
              )


	return nothing
end 

"""
	MatPartitioningSetUseEdgeWeights(petsclib::PetscLibType,part::MatPartitioning, use_edge_weights::PetscBool) 
Set a flag to indicate whether or not to use edge weights.

Logically Collective

Input Parameters:
- `part`             - the partitioning context
- `use_edge_weights` - the flag indicateing whether or not to use edge weights. By default no edge weights will be used,
that is, use_edge_weights is set to FALSE. If set use_edge_weights to TRUE, users need to make sure legal
edge weights are stored in an ADJ matrix.

Options Database Key:
- `-mat_partitioning_use_edge_weights` - (true or false)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningCreate()`, `MatPartitioningSetType()`, `MatPartitioningSetVertexWeights()`, `MatPartitioningSetPartitionWeights()`

# External Links
$(_doc_external("Mat/MatPartitioningSetUseEdgeWeights"))
"""
function MatPartitioningSetUseEdgeWeights(petsclib::PetscLibType, part::MatPartitioning, use_edge_weights::PetscBool) end

@for_petsc function MatPartitioningSetUseEdgeWeights(petsclib::$UnionPetscLib, part::MatPartitioning, use_edge_weights::PetscBool )

    @chk ccall(
               (:MatPartitioningSetUseEdgeWeights, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, PetscBool),
               part, use_edge_weights,
              )


	return nothing
end 

"""
	use_edge_weights::PetscBool = MatPartitioningGetUseEdgeWeights(petsclib::PetscLibType,part::MatPartitioning) 
Get a flag that indicates whether or not to edge weights are used.

Logically Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `use_edge_weights` - the flag indicateing whether or not to edge weights are used.

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningCreate()`, `MatPartitioningSetType()`, `MatPartitioningSetVertexWeights()`, `MatPartitioningSetPartitionWeights()`,
`MatPartitioningSetUseEdgeWeights`

# External Links
$(_doc_external("Mat/MatPartitioningGetUseEdgeWeights"))
"""
function MatPartitioningGetUseEdgeWeights(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningGetUseEdgeWeights(petsclib::$UnionPetscLib, part::MatPartitioning )
	use_edge_weights_ = Ref{PetscBool}()

    @chk ccall(
               (:MatPartitioningGetUseEdgeWeights, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{PetscBool}),
               part, use_edge_weights_,
              )

	use_edge_weights = use_edge_weights_[]

	return use_edge_weights
end 

"""
	newp::MatPartitioning = MatPartitioningCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a partitioning context.

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newp` - location to put the context

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningSetType()`, `MatPartitioningApply()`, `MatPartitioningDestroy()`,
`MatPartitioningSetAdjacency()`

# External Links
$(_doc_external("Mat/MatPartitioningCreate"))
"""
function MatPartitioningCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function MatPartitioningCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newp_ = Ref{MatPartitioning}()

    @chk ccall(
               (:MatPartitioningCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{MatPartitioning}),
               comm, newp_,
              )

	newp = newp_[]

	return newp
end 

"""
	MatPartitioningViewFromOptions(petsclib::PetscLibType,A::MatPartitioning, obj::PetscObject, name::String) 
View a partitioning context from the options database

Collective

Input Parameters:
- `A`    - the partitioning context
- `obj`  - Optional object that provides the prefix used in the options database check
- `name` - command line option

Options Database Key:
- `-mat_partitioning_view [viewertype]:...` - the viewer and its options

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningView()`, `PetscObjectViewFromOptions()`, `MatPartitioningCreate()`

# External Links
$(_doc_external("Mat/MatPartitioningViewFromOptions"))
"""
function MatPartitioningViewFromOptions(petsclib::PetscLibType, A::MatPartitioning, obj::PetscObject, name::String) end

@for_petsc function MatPartitioningViewFromOptions(petsclib::$UnionPetscLib, A::MatPartitioning, obj::PetscObject, name::String )

    @chk ccall(
               (:MatPartitioningViewFromOptions, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	MatPartitioningView(petsclib::PetscLibType,part::MatPartitioning, viewer::PetscViewer) 
Prints the partitioning data structure.

Collective

Input Parameters:
- `part`   - the partitioning context
- `viewer` - optional visualization context

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `PetscViewer`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Mat/MatPartitioningView"))
"""
function MatPartitioningView(petsclib::PetscLibType, part::MatPartitioning, viewer::PetscViewer) end

@for_petsc function MatPartitioningView(petsclib::$UnionPetscLib, part::MatPartitioning, viewer::PetscViewer )

    @chk ccall(
               (:MatPartitioningView, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, PetscViewer),
               part, viewer,
              )


	return nothing
end 

"""
	MatPartitioningSetType(petsclib::PetscLibType,part::MatPartitioning, type::MatPartitioningType) 
Sets the type of partitioner to use

Collective

Input Parameters:
- `part` - the partitioning context.
- `type` - a known method

Options Database Key:
- `-mat_partitioning_type  <type>` - (for instance, parmetis), use -help for a list of available methods or see  `MatPartitioningType`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningCreate()`, `MatPartitioningApply()`, `MatPartitioningType`

# External Links
$(_doc_external("Mat/MatPartitioningSetType"))
"""
function MatPartitioningSetType(petsclib::PetscLibType, part::MatPartitioning, type::MatPartitioningType) end

@for_petsc function MatPartitioningSetType(petsclib::$UnionPetscLib, part::MatPartitioning, type::MatPartitioningType )

    @chk ccall(
               (:MatPartitioningSetType, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, MatPartitioningType),
               part, type,
              )


	return nothing
end 

"""
	MatPartitioningSetFromOptions(petsclib::PetscLibType,part::MatPartitioning) 
Sets various partitioning options from the
options database for the partitioning object

Collective

Input Parameter:
- `part` - the partitioning context.

Options Database Keys:
- `-mat_partitioning_type  <type>` - (for instance, parmetis), use -help for a list of available methods
- `-mat_partitioning_nparts`       - number of subgraphs

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`

# External Links
$(_doc_external("Mat/MatPartitioningSetFromOptions"))
"""
function MatPartitioningSetFromOptions(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningSetFromOptions(petsclib::$UnionPetscLib, part::MatPartitioning )

    @chk ccall(
               (:MatPartitioningSetFromOptions, $petsc_library),
               PetscErrorCode,
               (MatPartitioning,),
               part,
              )


	return nothing
end 

"""
	MatPartitioningSetNumberVertexWeights(petsclib::PetscLibType,partitioning::MatPartitioning, ncon::PetscInt) 
Sets the number of weights per vertex

Not Collective

Input Parameters:
- `partitioning` - the partitioning context
- `ncon`         - the number of weights

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatPartitioning`, `MatPartitioningSetVertexWeights()`

# External Links
$(_doc_external("Mat/MatPartitioningSetNumberVertexWeights"))
"""
function MatPartitioningSetNumberVertexWeights(petsclib::PetscLibType, partitioning::MatPartitioning, ncon::PetscInt) end

@for_petsc function MatPartitioningSetNumberVertexWeights(petsclib::$UnionPetscLib, partitioning::MatPartitioning, ncon::$PetscInt )

    @chk ccall(
               (:MatPartitioningSetNumberVertexWeights, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscInt),
               partitioning, ncon,
              )


	return nothing
end 

"""
	MatPartitioningHierarchicalGetFineparts(petsclib::PetscLibType,part::MatPartitioning, fineparts::IS) 

# External Links
$(_doc_external("Mat/MatPartitioningHierarchicalGetFineparts"))
"""
function MatPartitioningHierarchicalGetFineparts(petsclib::PetscLibType, part::MatPartitioning, fineparts::IS) end

@for_petsc function MatPartitioningHierarchicalGetFineparts(petsclib::$UnionPetscLib, part::MatPartitioning, fineparts::IS )
	fineparts_ = Ref(fineparts.ptr)

    @chk ccall(
               (:MatPartitioningHierarchicalGetFineparts, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{CIS}),
               part, fineparts_,
              )

	fineparts.ptr = C_NULL

	return nothing
end 

"""
	MatPartitioningHierarchicalGetCoarseparts(petsclib::PetscLibType,part::MatPartitioning, coarseparts::IS) 

# External Links
$(_doc_external("Mat/MatPartitioningHierarchicalGetCoarseparts"))
"""
function MatPartitioningHierarchicalGetCoarseparts(petsclib::PetscLibType, part::MatPartitioning, coarseparts::IS) end

@for_petsc function MatPartitioningHierarchicalGetCoarseparts(petsclib::$UnionPetscLib, part::MatPartitioning, coarseparts::IS )
	coarseparts_ = Ref(coarseparts.ptr)

    @chk ccall(
               (:MatPartitioningHierarchicalGetCoarseparts, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{CIS}),
               part, coarseparts_,
              )

	coarseparts.ptr = C_NULL

	return nothing
end 

"""
	MatPartitioningHierarchicalSetNcoarseparts(petsclib::PetscLibType,part::MatPartitioning, ncoarseparts::PetscInt) 

# External Links
$(_doc_external("Mat/MatPartitioningHierarchicalSetNcoarseparts"))
"""
function MatPartitioningHierarchicalSetNcoarseparts(petsclib::PetscLibType, part::MatPartitioning, ncoarseparts::PetscInt) end

@for_petsc function MatPartitioningHierarchicalSetNcoarseparts(petsclib::$UnionPetscLib, part::MatPartitioning, ncoarseparts::$PetscInt )

    @chk ccall(
               (:MatPartitioningHierarchicalSetNcoarseparts, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscInt),
               part, ncoarseparts,
              )


	return nothing
end 

"""
	MatPartitioningHierarchicalSetNfineparts(petsclib::PetscLibType,part::MatPartitioning, nfineparts::PetscInt) 

# External Links
$(_doc_external("Mat/MatPartitioningHierarchicalSetNfineparts"))
"""
function MatPartitioningHierarchicalSetNfineparts(petsclib::PetscLibType, part::MatPartitioning, nfineparts::PetscInt) end

@for_petsc function MatPartitioningHierarchicalSetNfineparts(petsclib::$UnionPetscLib, part::MatPartitioning, nfineparts::$PetscInt )

    @chk ccall(
               (:MatPartitioningHierarchicalSetNfineparts, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscInt),
               part, nfineparts,
              )


	return nothing
end 

"""
	MatPartitioningPartySetGlobal(petsclib::PetscLibType,part::MatPartitioning, glob::String) 
Set global method for Party partitioner.

Collective

Input Parameters:
- `part`   - the partitioning context
- `global` - a string representing the method

Options Database Key:
- `-mat_partitioning_party_global <method>` - the global method

Level: advanced

-seealso: `MATPARTITIONINGPARTY`, `MatPartitioningPartySetLocal()`

# External Links
$(_doc_external("Mat/MatPartitioningPartySetGlobal"))
"""
function MatPartitioningPartySetGlobal(petsclib::PetscLibType, part::MatPartitioning, glob::String) end

@for_petsc function MatPartitioningPartySetGlobal(petsclib::$UnionPetscLib, part::MatPartitioning, glob::String )

    @chk ccall(
               (:MatPartitioningPartySetGlobal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{Cchar}),
               part, glob,
              )


	return nothing
end 

"""
	MatPartitioningPartySetLocal(petsclib::PetscLibType,part::MatPartitioning, loc::String) 
Set local method used by the Party partitioner.

Collective

Input Parameters:
- `part`  - the partitioning context
- `local` - a string representing the method

Options Database Key:
- `-mat_partitioning_party_local <method>` - the local method

Level: advanced

-seealso: `MATPARTITIONINGPARTY`, `MatPartitioningPartySetGlobal()`

# External Links
$(_doc_external("Mat/MatPartitioningPartySetLocal"))
"""
function MatPartitioningPartySetLocal(petsclib::PetscLibType, part::MatPartitioning, loc::String) end

@for_petsc function MatPartitioningPartySetLocal(petsclib::$UnionPetscLib, part::MatPartitioning, loc::String )

    @chk ccall(
               (:MatPartitioningPartySetLocal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{Cchar}),
               part, loc,
              )


	return nothing
end 

"""
	MatPartitioningPartySetCoarseLevel(petsclib::PetscLibType,part::MatPartitioning, level::PetscReal) 
Set the coarse level parameter for the
Party partitioner.

Collective

Input Parameters:
- `part`  - the partitioning context
- `level` - the coarse level in range [0.0,1.0]

Options Database Key:
- `-mat_partitioning_party_coarse <l>` - Coarse level

Level: advanced

-seealso: `MATPARTITIONINGPARTY`

# External Links
$(_doc_external("Mat/MatPartitioningPartySetCoarseLevel"))
"""
function MatPartitioningPartySetCoarseLevel(petsclib::PetscLibType, part::MatPartitioning, level::PetscReal) end

@for_petsc function MatPartitioningPartySetCoarseLevel(petsclib::$UnionPetscLib, part::MatPartitioning, level::$PetscReal )

    @chk ccall(
               (:MatPartitioningPartySetCoarseLevel, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscReal),
               part, level,
              )


	return nothing
end 

"""
	MatPartitioningPartySetMatchOptimization(petsclib::PetscLibType,part::MatPartitioning, opt::PetscBool) 
Activate matching optimization for
graph reduction.

Collective

Input Parameters:
- `part` - the partitioning context
- `opt`  - boolean flag

Options Database Key:
- `-mat_partitioning_party_match_optimization` - Matching optimization on/off

Level: advanced

-seealso: `MATPARTITIONINGPARTY`

# External Links
$(_doc_external("Mat/MatPartitioningPartySetMatchOptimization"))
"""
function MatPartitioningPartySetMatchOptimization(petsclib::PetscLibType, part::MatPartitioning, opt::PetscBool) end

@for_petsc function MatPartitioningPartySetMatchOptimization(petsclib::$UnionPetscLib, part::MatPartitioning, opt::PetscBool )

    @chk ccall(
               (:MatPartitioningPartySetMatchOptimization, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, PetscBool),
               part, opt,
              )


	return nothing
end 

"""
	MatPartitioningPartySetBipart(petsclib::PetscLibType,part::MatPartitioning, bp::PetscBool) 
Activate or deactivate recursive bisection in the Party partitioner

Collective

Input Parameters:
- `part` - the partitioning context
- `bp`   - boolean flag

Options Database Key:
- `-mat_partitioning_party_bipart` - Bipartitioning option on/off

Level: advanced

-seealso: `MATPARTITIONINGPARTY`

# External Links
$(_doc_external("Mat/MatPartitioningPartySetBipart"))
"""
function MatPartitioningPartySetBipart(petsclib::PetscLibType, part::MatPartitioning, bp::PetscBool) end

@for_petsc function MatPartitioningPartySetBipart(petsclib::$UnionPetscLib, part::MatPartitioning, bp::PetscBool )

    @chk ccall(
               (:MatPartitioningPartySetBipart, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, PetscBool),
               part, bp,
              )


	return nothing
end 

"""
	MatPartitioningParmetisSetCoarseSequential(petsclib::PetscLibType,part::MatPartitioning) 
Use the sequential code to
do the partitioning of the coarse grid.

Logically Collective

Input Parameter:
- `part` - the partitioning context

Level: advanced

-seealso: `MATPARTITIONINGPARMETIS`

# External Links
$(_doc_external("Mat/MatPartitioningParmetisSetCoarseSequential"))
"""
function MatPartitioningParmetisSetCoarseSequential(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningParmetisSetCoarseSequential(petsclib::$UnionPetscLib, part::MatPartitioning )

    @chk ccall(
               (:MatPartitioningParmetisSetCoarseSequential, $petsc_library),
               PetscErrorCode,
               (MatPartitioning,),
               part,
              )


	return nothing
end 

"""
	MatPartitioningParmetisSetRepartition(petsclib::PetscLibType,part::MatPartitioning) 
Repartition
current mesh to rebalance computation.

Logically Collective

Input Parameter:
- `part` - the partitioning context

Level: advanced

-seealso: `MATPARTITIONINGPARMETIS`

# External Links
$(_doc_external("Mat/MatPartitioningParmetisSetRepartition"))
"""
function MatPartitioningParmetisSetRepartition(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningParmetisSetRepartition(petsclib::$UnionPetscLib, part::MatPartitioning )

    @chk ccall(
               (:MatPartitioningParmetisSetRepartition, $petsc_library),
               PetscErrorCode,
               (MatPartitioning,),
               part,
              )


	return nothing
end 

"""
	cut::PetscInt = MatPartitioningParmetisGetEdgeCut(petsclib::PetscLibType,part::MatPartitioning) 
Returns the number of edge cuts in the vertex partition.

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `cut` - the edge cut

Level: advanced

-seealso: `MATPARTITIONINGPARMETIS`

# External Links
$(_doc_external("Mat/MatPartitioningParmetisGetEdgeCut"))
"""
function MatPartitioningParmetisGetEdgeCut(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningParmetisGetEdgeCut(petsclib::$UnionPetscLib, part::MatPartitioning )
	cut_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatPartitioningParmetisGetEdgeCut, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscInt}),
               part, cut_,
              )

	cut = cut_[]

	return cut
end 

"""
	MatPartitioningPTScotchSetImbalance(petsclib::PetscLibType,part::MatPartitioning, imb::PetscReal) 
Sets the value of the load imbalance
ratio to be used during strategy selection.

Collective

Input Parameters:
- `part` - the partitioning context
- `imb`  - the load imbalance ratio

Options Database Key:
- `-mat_partitioning_ptscotch_imbalance <imb>` - set load imbalance ratio

-seealso: `MATPARTITIONINGSCOTCH`, `MatPartitioningPTScotchSetStrategy()`, `MatPartitioningPTScotchGetImbalance()`

# External Links
$(_doc_external("Mat/MatPartitioningPTScotchSetImbalance"))
"""
function MatPartitioningPTScotchSetImbalance(petsclib::PetscLibType, part::MatPartitioning, imb::PetscReal) end

@for_petsc function MatPartitioningPTScotchSetImbalance(petsclib::$UnionPetscLib, part::MatPartitioning, imb::$PetscReal )

    @chk ccall(
               (:MatPartitioningPTScotchSetImbalance, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscReal),
               part, imb,
              )


	return nothing
end 

"""
	imb::PetscReal = MatPartitioningPTScotchGetImbalance(petsclib::PetscLibType,part::MatPartitioning) 
Gets the value of the load imbalance
ratio used during strategy selection.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `imb` - the load imbalance ratio

Level: advanced

-seealso: `MATPARTITIONINGSCOTCH`, `MatPartitioningPTScotchSetImbalance()`

# External Links
$(_doc_external("Mat/MatPartitioningPTScotchGetImbalance"))
"""
function MatPartitioningPTScotchGetImbalance(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningPTScotchGetImbalance(petsclib::$UnionPetscLib, part::MatPartitioning )
	imb_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatPartitioningPTScotchGetImbalance, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscReal}),
               part, imb_,
              )

	imb = imb_[]

	return imb
end 

"""
	MatPartitioningPTScotchSetStrategy(petsclib::PetscLibType,part::MatPartitioning, strategy::MPPTScotchStrategyType) 
Sets the strategy to be used in PTScotch.

Collective

Input Parameters:
- `part`     - the partitioning context
- `strategy` - the strategy, one of
-seealso: `MATPARTITIONINGSCOTCH`, `MatPartitioningPTScotchSetImbalance()`, `MatPartitioningPTScotchGetStrategy()`

# External Links
$(_doc_external("Mat/MatPartitioningPTScotchSetStrategy"))
"""
function MatPartitioningPTScotchSetStrategy(petsclib::PetscLibType, part::MatPartitioning, strategy::MPPTScotchStrategyType) end

@for_petsc function MatPartitioningPTScotchSetStrategy(petsclib::$UnionPetscLib, part::MatPartitioning, strategy::MPPTScotchStrategyType )

    @chk ccall(
               (:MatPartitioningPTScotchSetStrategy, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, MPPTScotchStrategyType),
               part, strategy,
              )


	return nothing
end 

"""
	MatPartitioningPTScotchGetStrategy(petsclib::PetscLibType,part::MatPartitioning, strategy::MPPTScotchStrategyType) 
Gets the strategy used in PTScotch.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `strategy` - the strategy

Level: advanced

-seealso: `MATPARTITIONINGSCOTCH`, `MatPartitioningPTScotchSetStrategy()`

# External Links
$(_doc_external("Mat/MatPartitioningPTScotchGetStrategy"))
"""
function MatPartitioningPTScotchGetStrategy(petsclib::PetscLibType, part::MatPartitioning, strategy::MPPTScotchStrategyType) end

@for_petsc function MatPartitioningPTScotchGetStrategy(petsclib::$UnionPetscLib, part::MatPartitioning, strategy::MPPTScotchStrategyType )

    @chk ccall(
               (:MatPartitioningPTScotchGetStrategy, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{MPPTScotchStrategyType}),
               part, strategy,
              )


	return nothing
end 

"""
	MatPartitioningChacoSetGlobal(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoGlobalType) 
Set the global method for Chaco partitioner.

Collective

Input Parameters:
- `part`   - the partitioning context
- `method` - one of `MP_CHACO_MULTILEVEL`, `MP_CHACO_SPECTRAL`, `MP_CHACO_LINEAR`,
`MP_CHACO_RANDOM` or `MP_CHACO_SCATTERED`

Options Database Key:
- `-mat_partitioning_chaco_global <method>` - the global method

Level: advanced

-seealso: `MatPartitioning`, `MatPartioningSetType()`, `MatPartitioningType`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetLocal()`, `MatPartitioningChacoGetGlobal()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetGlobal"))
"""
function MatPartitioningChacoSetGlobal(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoGlobalType) end

@for_petsc function MatPartitioningChacoSetGlobal(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoGlobalType )

    @chk ccall(
               (:MatPartitioningChacoSetGlobal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, MPChacoGlobalType),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoGetGlobal(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoGlobalType) 
Get the global method used by the Chaco partitioner.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `method` - the method

Level: advanced

-seealso: `MatPartitioningType`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetGlobal()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoGetGlobal"))
"""
function MatPartitioningChacoGetGlobal(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoGlobalType) end

@for_petsc function MatPartitioningChacoGetGlobal(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoGlobalType )

    @chk ccall(
               (:MatPartitioningChacoGetGlobal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{MPChacoGlobalType}),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoSetLocal(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoLocalType) 
Set the local method for the Chaco partitioner.

Collective

Input Parameters:
- `part`   - the partitioning context
- `method` - one of `MP_CHACO_KERNIGHAN` or `MP_CHACO_NONE`

Options Database Key:
- `-mat_partitioning_chaco_local <method>` - the local method

Level: advanced

-seealso: `MatPartitioningType`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetGlobal()`, `MatPartitioningChacoGetLocal()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetLocal"))
"""
function MatPartitioningChacoSetLocal(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoLocalType) end

@for_petsc function MatPartitioningChacoSetLocal(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoLocalType )

    @chk ccall(
               (:MatPartitioningChacoSetLocal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, MPChacoLocalType),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoGetLocal(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoLocalType) 
Get local method used by the Chaco partitioner.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `method` - the method

Level: advanced

-seealso: `MatPartitioningType`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetLocal()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoGetLocal"))
"""
function MatPartitioningChacoGetLocal(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoLocalType) end

@for_petsc function MatPartitioningChacoGetLocal(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoLocalType )

    @chk ccall(
               (:MatPartitioningChacoGetLocal, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{MPChacoLocalType}),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoSetCoarseLevel(petsclib::PetscLibType,part::MatPartitioning, level::PetscReal) 
Set the coarse level parameter for the
Chaco partitioner.

Collective

Input Parameters:
- `part`  - the partitioning context
- `level` - the coarse level in range [0.0,1.0]

Options Database Key:
- `-mat_partitioning_chaco_coarse <l>` - Coarse level

Level: advanced

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetCoarseLevel"))
"""
function MatPartitioningChacoSetCoarseLevel(petsclib::PetscLibType, part::MatPartitioning, level::PetscReal) end

@for_petsc function MatPartitioningChacoSetCoarseLevel(petsclib::$UnionPetscLib, part::MatPartitioning, level::$PetscReal )

    @chk ccall(
               (:MatPartitioningChacoSetCoarseLevel, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscReal),
               part, level,
              )


	return nothing
end 

"""
	MatPartitioningChacoSetEigenSolver(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoEigenType) 
Set the eigensolver method for Chaco partitioner.

Collective

Input Parameters:
- `part`   - the partitioning context
- `method` - one of `MP_CHACO_LANCZOS` or `MP_CHACO_RQI`

Options Database Key:
- `-mat_partitioning_chaco_eigen_solver <method>` - the eigensolver

Level: advanced

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenTol()`, `MatPartitioningChacoSetEigenNumber()`,
`MatPartitioningChacoGetEigenSolver()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetEigenSolver"))
"""
function MatPartitioningChacoSetEigenSolver(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoEigenType) end

@for_petsc function MatPartitioningChacoSetEigenSolver(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoEigenType )

    @chk ccall(
               (:MatPartitioningChacoSetEigenSolver, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, MPChacoEigenType),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoGetEigenSolver(petsclib::PetscLibType,part::MatPartitioning, method::MPChacoEigenType) 
Get the eigensolver used by the Chaco partitioner.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `method` - the method

Level: advanced

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenSolver()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoGetEigenSolver"))
"""
function MatPartitioningChacoGetEigenSolver(petsclib::PetscLibType, part::MatPartitioning, method::MPChacoEigenType) end

@for_petsc function MatPartitioningChacoGetEigenSolver(petsclib::$UnionPetscLib, part::MatPartitioning, method::MPChacoEigenType )

    @chk ccall(
               (:MatPartitioningChacoGetEigenSolver, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{MPChacoEigenType}),
               part, method,
              )


	return nothing
end 

"""
	MatPartitioningChacoSetEigenTol(petsclib::PetscLibType,part::MatPartitioning, tol::PetscReal) 
Sets the tolerance for the eigensolver used by Chaco

Collective

Input Parameters:
- `part` - the partitioning context
- `tol`  - the tolerance

Options Database Key:
- `-mat_partitioning_chaco_eigen_tol <tol>` - Tolerance for eigensolver

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenSolver()`, `MatPartitioningChacoGetEigenTol()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetEigenTol"))
"""
function MatPartitioningChacoSetEigenTol(petsclib::PetscLibType, part::MatPartitioning, tol::PetscReal) end

@for_petsc function MatPartitioningChacoSetEigenTol(petsclib::$UnionPetscLib, part::MatPartitioning, tol::$PetscReal )

    @chk ccall(
               (:MatPartitioningChacoSetEigenTol, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscReal),
               part, tol,
              )


	return nothing
end 

"""
	tol::PetscReal = MatPartitioningChacoGetEigenTol(petsclib::PetscLibType,part::MatPartitioning) 
Gets the eigensolver tolerance used by Chaco

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `tol` - the tolerance

Level: advanced

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenTol()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoGetEigenTol"))
"""
function MatPartitioningChacoGetEigenTol(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningChacoGetEigenTol(petsclib::$UnionPetscLib, part::MatPartitioning )
	tol_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatPartitioningChacoGetEigenTol, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscReal}),
               part, tol_,
              )

	tol = tol_[]

	return tol
end 

"""
	MatPartitioningChacoSetEigenNumber(petsclib::PetscLibType,part::MatPartitioning, num::PetscInt) 
Sets the number of eigenvectors to compute by Chaco during partitioning
during partitioning.

Collective

Input Parameters:
- `part` - the partitioning context
- `num`  - the number of eigenvectors

Options Database Key:
- `-mat_partitioning_chaco_eigen_number <n>` - Number of eigenvectors

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenSolver()`, `MatPartitioningChacoGetEigenTol()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoSetEigenNumber"))
"""
function MatPartitioningChacoSetEigenNumber(petsclib::PetscLibType, part::MatPartitioning, num::PetscInt) end

@for_petsc function MatPartitioningChacoSetEigenNumber(petsclib::$UnionPetscLib, part::MatPartitioning, num::$PetscInt )

    @chk ccall(
               (:MatPartitioningChacoSetEigenNumber, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, $PetscInt),
               part, num,
              )


	return nothing
end 

"""
	num::PetscInt = MatPartitioningChacoGetEigenNumber(petsclib::PetscLibType,part::MatPartitioning) 
Gets the number of eigenvectors used by Chaco.

Not Collective

Input Parameter:
- `part` - the partitioning context

Output Parameter:
- `num` - number of eigenvectors

Level: advanced

-seealso: `MatPartitioningType`, `MatPartitioning`, `MATPARTITIONINGCHACO`, `MatPartitioningChacoSetEigenNumber()`

# External Links
$(_doc_external("Mat/MatPartitioningChacoGetEigenNumber"))
"""
function MatPartitioningChacoGetEigenNumber(petsclib::PetscLibType, part::MatPartitioning) end

@for_petsc function MatPartitioningChacoGetEigenNumber(petsclib::$UnionPetscLib, part::MatPartitioning )
	num_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatPartitioningChacoGetEigenNumber, $petsc_library),
               PetscErrorCode,
               (MatPartitioning, Ptr{$PetscInt}),
               part, num_,
              )

	num = num_[]

	return num
end 

"""
	MatNullSpaceSetFunction(petsclib::PetscLibType,sp::MatNullSpace, rem::MatNullSpaceRemoveFn, ctx::Cvoid) 
set a function that removes a null space from a vector
out of null spaces.

Logically Collective

Input Parameters:
- `sp`  - the `MatNullSpace` null space object
- `rem` - the function that removes the null space
- `ctx` - context for the remove function

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceDestroy()`, `MatNullSpaceRemove()`, `MatSetNullSpace()`, `MatNullSpaceCreate()`, `MatNullSpaceRemoveFn`

# External Links
$(_doc_external("Mat/MatNullSpaceSetFunction"))
"""
function MatNullSpaceSetFunction(petsclib::PetscLibType, sp::MatNullSpace, rem::MatNullSpaceRemoveFn, ctx::Cvoid) end

@for_petsc function MatNullSpaceSetFunction(petsclib::$UnionPetscLib, sp::MatNullSpace, rem::MatNullSpaceRemoveFn, ctx::Cvoid )

    @chk ccall(
               (:MatNullSpaceSetFunction, $petsc_library),
               PetscErrorCode,
               (MatNullSpace, Ptr{MatNullSpaceRemoveFn}, Ptr{Cvoid}),
               sp, rem, ctx,
              )


	return nothing
end 

"""
	has_const::PetscBool,n::PetscInt = MatNullSpaceGetVecs(petsclib::PetscLibType,sp::MatNullSpace, vecs::Vector{PetscVec}) 
get the vectors defining the null space

Not Collective

Input Parameter:
- `sp` - null space object

Output Parameters:
- `has_const` - `PETSC_TRUE` if the null space contains the constant vector, otherwise `PETSC_FALSE`
- `n`         - number of vectors (excluding constant vector) in the null space
- `vecs`      - returns array of length `n` containing the orthonormal vectors that span the null space (excluding the constant vector), `NULL` if `n` is 0

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceCreate()`, `MatGetNullSpace()`, `MatGetNearNullSpace()`

# External Links
$(_doc_external("Mat/MatNullSpaceGetVecs"))
"""
function MatNullSpaceGetVecs(petsclib::PetscLibType, sp::MatNullSpace, vecs::Vector{PetscVec}) end

@for_petsc function MatNullSpaceGetVecs(petsclib::$UnionPetscLib, sp::MatNullSpace, vecs::Vector{PetscVec} )
	has_const_ = Ref{PetscBool}()
	n_ = Ref{$PetscInt}()
	vecs_ = Ref(pointer(vecs))

    @chk ccall(
               (:MatNullSpaceGetVecs, $petsc_library),
               PetscErrorCode,
               (MatNullSpace, Ptr{PetscBool}, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
               sp, has_const_, n_, vecs_,
              )

	has_const = has_const_[]
	n = n_[]

	return has_const,n
end 

"""
	sp::MatNullSpace = MatNullSpaceCreateRigidBody(petsclib::PetscLibType,coords::PetscVec) 
create rigid body modes from coordinates

Collective

Input Parameter:
- `coords` - block of coordinates of each node, must have block size set

Output Parameter:
- `sp` - the null space

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatSetNullSpace()`, `PCGAMG`

# External Links
$(_doc_external("Mat/MatNullSpaceCreateRigidBody"))
"""
function MatNullSpaceCreateRigidBody(petsclib::PetscLibType, coords::PetscVec) end

@for_petsc function MatNullSpaceCreateRigidBody(petsclib::$UnionPetscLib, coords::PetscVec )
	sp_ = Ref{MatNullSpace}()

    @chk ccall(
               (:MatNullSpaceCreateRigidBody, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{MatNullSpace}),
               coords, sp_,
              )

	sp = sp_[]

	return sp
end 

"""
	MatNullSpaceView(petsclib::PetscLibType,sp::MatNullSpace, viewer::PetscViewer) 
Visualizes a null space object.

Collective

Input Parameters:
- `sp`     - the null space
- `viewer` - visualization context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `PetscViewer`, `MatNullSpaceCreate()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Mat/MatNullSpaceView"))
"""
function MatNullSpaceView(petsclib::PetscLibType, sp::MatNullSpace, viewer::PetscViewer) end

@for_petsc function MatNullSpaceView(petsclib::$UnionPetscLib, sp::MatNullSpace, viewer::PetscViewer )

    @chk ccall(
               (:MatNullSpaceView, $petsc_library),
               PetscErrorCode,
               (MatNullSpace, PetscViewer),
               sp, viewer,
              )


	return nothing
end 

"""
	SP::MatNullSpace = MatNullSpaceCreate(petsclib::PetscLibType,comm::MPI_Comm, has_cnst::PetscBool, n::PetscInt, vecs::Vector{PetscVec}) 
Creates a `MatNullSpace` data structure used to project vectors out of null spaces.

Collective

Input Parameters:
- `comm`     - the MPI communicator associated with the object
- `has_cnst` - `PETSC_TRUE` if the null space contains the constant vector; otherwise `PETSC_FALSE`
- `n`        - number of vectors (excluding constant vector) in null space
- `vecs`     - the vectors that span the null space (excluding the constant vector);
these vectors must be orthonormal. These vectors are NOT copied, so do not change them
after this call. You should free the array that you pass in and destroy the vectors (this will reduce the reference count
for them by one).

Output Parameter:
- `SP` - the null space context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceDestroy()`, `MatNullSpaceRemove()`, `MatSetNullSpace()`, `MatNullSpaceSetFunction()`

# External Links
$(_doc_external("Mat/MatNullSpaceCreate"))
"""
function MatNullSpaceCreate(petsclib::PetscLibType, comm::MPI_Comm, has_cnst::PetscBool, n::PetscInt, vecs::Vector{PetscVec}) end

@for_petsc function MatNullSpaceCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, has_cnst::PetscBool, n::$PetscInt, vecs::Vector{PetscVec} )
	SP_ = Ref{MatNullSpace}()

    @chk ccall(
               (:MatNullSpaceCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscBool, $PetscInt, Ptr{CVec}, Ptr{MatNullSpace}),
               comm, has_cnst, n, vecs, SP_,
              )

	SP = SP_[]

	return SP
end 

"""
	MatNullSpaceDestroy(petsclib::PetscLibType,sp::MatNullSpace) 
Destroys a data structure used to project vectors out of null spaces.

Collective

Input Parameter:
- `sp` - the null space context to be destroyed

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceCreate()`, `MatNullSpaceRemove()`, `MatNullSpaceSetFunction()`

# External Links
$(_doc_external("Mat/MatNullSpaceDestroy"))
"""
function MatNullSpaceDestroy(petsclib::PetscLibType, sp::MatNullSpace) end

@for_petsc function MatNullSpaceDestroy(petsclib::$UnionPetscLib, sp::MatNullSpace )

    @chk ccall(
               (:MatNullSpaceDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatNullSpace},),
               Ref(sp),
              )


	return nothing
end 

"""
	MatNullSpaceRemove(petsclib::PetscLibType,sp::MatNullSpace, vec::PetscVec) 
Removes all the components of a null space from a vector.

Collective

Input Parameters:
- `sp`  - the null space context (if this is `NULL` then no null space is removed)
- `vec` - the vector from which the null space is to be removed

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceCreate()`, `MatNullSpaceDestroy()`, `MatNullSpaceSetFunction()`

# External Links
$(_doc_external("Mat/MatNullSpaceRemove"))
"""
function MatNullSpaceRemove(petsclib::PetscLibType, sp::MatNullSpace, vec::PetscVec) end

@for_petsc function MatNullSpaceRemove(petsclib::$UnionPetscLib, sp::MatNullSpace, vec::PetscVec )

    @chk ccall(
               (:MatNullSpaceRemove, $petsc_library),
               PetscErrorCode,
               (MatNullSpace, CVec),
               sp, vec,
              )


	return nothing
end 

"""
	isNull::PetscBool = MatNullSpaceTest(petsclib::PetscLibType,sp::MatNullSpace, mat::PetscMat) 
Tests if the claimed null space is really a null space of a matrix

Collective

Input Parameters:
- `sp`  - the null space context
- `mat` - the matrix

Output Parameter:
- `isNull` - `PETSC_TRUE` if the nullspace is valid for this matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatNullSpaceCreate()`, `MatNullSpaceDestroy()`, `MatNullSpaceSetFunction()`

# External Links
$(_doc_external("Mat/MatNullSpaceTest"))
"""
function MatNullSpaceTest(petsclib::PetscLibType, sp::MatNullSpace, mat::PetscMat) end

@for_petsc function MatNullSpaceTest(petsclib::$UnionPetscLib, sp::MatNullSpace, mat::PetscMat )
	isNull_ = Ref{PetscBool}()

    @chk ccall(
               (:MatNullSpaceTest, $petsc_library),
               PetscErrorCode,
               (MatNullSpace, CMat, Ptr{PetscBool}),
               sp, mat, isNull_,
              )

	isNull = isNull_[]

	return isNull
end 

"""
	MatMFFDWPSetComputeNormU(petsclib::PetscLibType,A::PetscMat, flag::PetscBool) 
Sets whether it computes the ||U|| used by the Walker
PETSc routine for computing h. With any Krylov solver this need only
be computed during the first iteration and kept for later.

Input Parameters:
- `A`    - the `MATMFFD` matrix
- `flag` - `PETSC_TRUE` causes it to compute ||U||, `PETSC_FALSE` uses the previous value

Options Database Key:
- `-mat_mffd_compute_normu <true,false>` - true by default, false can save calculations but you
must be sure that ||U|| has not changed in the mean time.

Level: advanced

-seealso: `MATMFFD_WP`, `MATMFFD`, `MatMFFDSetFunctionError()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Mat/MatMFFDWPSetComputeNormU"))
"""
function MatMFFDWPSetComputeNormU(petsclib::PetscLibType, A::PetscMat, flag::PetscBool) end

@for_petsc function MatMFFDWPSetComputeNormU(petsclib::$UnionPetscLib, A::PetscMat, flag::PetscBool )

    @chk ccall(
               (:MatMFFDWPSetComputeNormU, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flag,
              )


	return nothing
end 

"""
	MatMFFDFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the MATMFFD` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `PetscFinalize()`, `MatCreateMFFD()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Mat/MatMFFDFinalizePackage"))
"""
function MatMFFDFinalizePackage(petsclib::PetscLibType) end

@for_petsc function MatMFFDFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:MatMFFDFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	MatMFFDInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the MATMFFD` package. It is called
from `MatInitializePackage()`.

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `PetscInitialize()`

# External Links
$(_doc_external("Mat/MatMFFDInitializePackage"))
"""
function MatMFFDInitializePackage(petsclib::PetscLibType) end

@for_petsc function MatMFFDInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:MatMFFDInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	MatMFFDSetType(petsclib::PetscLibType,mat::PetscMat, ftype::MatMFFDType) 
Sets the method that is used to compute the
differencing parameter for finite difference matrix-free formulations.

Input Parameters:
- `mat`   - the "matrix-free" matrix created via `MatCreateSNESMF()`, or `MatCreateMFFD()`
or `MatSetType`(mat,`MATMFFD`);
- `ftype` - the type requested, either `MATMFFD_WP` or `MATMFFD_DS`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MATMFFD_WP`, `MATMFFD_DS`, `MatCreateSNESMF()`, `MatMFFDRegister()`, `MatMFFDSetFunction()`, `MatCreateMFFD()`

# External Links
$(_doc_external("Mat/MatMFFDSetType"))
"""
function MatMFFDSetType(petsclib::PetscLibType, mat::PetscMat, ftype::MatMFFDType) end

@for_petsc function MatMFFDSetType(petsclib::$UnionPetscLib, mat::PetscMat, ftype::MatMFFDType )

    @chk ccall(
               (:MatMFFDSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatMFFDType),
               mat, ftype,
              )


	return nothing
end 

"""
	MatMFFDRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method to the `MATMFFD` registry.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined compute-h module
- `function` - routine to create method context

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDRegisterAll()`, `MatMFFDRegisterDestroy()`

# External Links
$(_doc_external("Mat/MatMFFDRegister"))
"""
function MatMFFDRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatMFFDRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatMFFDRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatMFFDSetOptionsPrefix(petsclib::PetscLibType,mat::PetscMat, prefix::String) 
Sets the prefix used for searching for all
MATMFFD` options in the database.

Collective

Input Parameters:
- `mat`    - the `MATMFFD` context
- `prefix` - the prefix to prepend to all option names

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatSetFromOptions()`, `MatCreateSNESMF()`, `MatCreateMFFD()`

# External Links
$(_doc_external("Mat/MatMFFDSetOptionsPrefix"))
"""
function MatMFFDSetOptionsPrefix(petsclib::PetscLibType, mat::PetscMat, prefix::String) end

@for_petsc function MatMFFDSetOptionsPrefix(petsclib::$UnionPetscLib, mat::PetscMat, prefix::String )

    @chk ccall(
               (:MatMFFDSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               mat, prefix,
              )


	return nothing
end 

"""
	h::PetscScalar = MatMFFDGetH(petsclib::PetscLibType,mat::PetscMat) 
Gets the last value that was used as the differencing for a `MATMFFD` matrix
parameter.

Not Collective

Input Parameters:
- `mat` - the `MATMFFD` matrix

Output Parameter:
- `h` - the differencing step size

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatCreateSNESMF()`, `MatMFFDSetHHistory()`, `MatCreateMFFD()`, `MatMFFDResetHHistory()`

# External Links
$(_doc_external("Mat/MatMFFDGetH"))
"""
function MatMFFDGetH(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatMFFDGetH(petsclib::$UnionPetscLib, mat::PetscMat )
	h_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatMFFDGetH, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, h_,
              )

	h = h_[]

	return h
end 

"""
	MatMFFDSetFunction(petsclib::PetscLibType,mat::PetscMat, func::MatMFFDFn, funcctx::Cvoid) 
Sets the function used in applying the matrix

Logically Collective

Input Parameters:
- `mat`     - the matrix-free matrix `MATMFFD` created via `MatCreateSNESMF()` or `MatCreateMFFD()`
- `func`    - the function to use
- `funcctx` - optional function context passed to function

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDFn`, `MatCreateSNESMF()`, `MatMFFDGetH()`, `MatCreateMFFD()`,
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`, `SNESSetFunction()`

# External Links
$(_doc_external("Mat/MatMFFDSetFunction"))
"""
function MatMFFDSetFunction(petsclib::PetscLibType, mat::PetscMat, func::MatMFFDFn, funcctx::Cvoid) end

@for_petsc function MatMFFDSetFunction(petsclib::$UnionPetscLib, mat::PetscMat, func::MatMFFDFn, funcctx::Cvoid )

    @chk ccall(
               (:MatMFFDSetFunction, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatMFFDFn}, Ptr{Cvoid}),
               mat, func, funcctx,
              )


	return nothing
end 

"""
	MatMFFDSetFunctioni(petsclib::PetscLibType,mat::PetscMat, funci::MatMFFDiFn) 
Sets the function for computing a single component for a `MATMFFD` matrix

Logically Collective

Input Parameters:
- `mat`   - the matrix-free matrix `MATMFFD`
- `funci` - the function to use

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDiFn`, `MatCreateSNESMF()`, `MatMFFDGetH()`, `MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`,
`SNESSetFunction()`, `MatGetDiagonal()`

# External Links
$(_doc_external("Mat/MatMFFDSetFunctioni"))
"""
function MatMFFDSetFunctioni(petsclib::PetscLibType, mat::PetscMat, funci::MatMFFDiFn) end

@for_petsc function MatMFFDSetFunctioni(petsclib::$UnionPetscLib, mat::PetscMat, funci::MatMFFDiFn )

    @chk ccall(
               (:MatMFFDSetFunctioni, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatMFFDiFn}),
               mat, funci,
              )


	return nothing
end 

"""
	MatMFFDSetFunctioniBase(petsclib::PetscLibType,mat::PetscMat, func::MatMFFDiBaseFn) 
Sets the function to compute the base vector for a single component function evaluation for a `MATMFFD` matrix

Logically Collective

Input Parameters:
- `mat`  - the `MATMFFD` matrix-free matrix
- `func` - the function to use

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatCreateSNESMF()`, `MatMFFDGetH()`, `MatCreateMFFD()`,
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`, `SNESSetFunction()`, `MatGetDiagonal()`

# External Links
$(_doc_external("Mat/MatMFFDSetFunctioniBase"))
"""
function MatMFFDSetFunctioniBase(petsclib::PetscLibType, mat::PetscMat, func::MatMFFDiBaseFn) end

@for_petsc function MatMFFDSetFunctioniBase(petsclib::$UnionPetscLib, mat::PetscMat, func::MatMFFDiBaseFn )

    @chk ccall(
               (:MatMFFDSetFunctioniBase, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatMFFDiBaseFn}),
               mat, func,
              )


	return nothing
end 

"""
	MatMFFDSetPeriod(petsclib::PetscLibType,mat::PetscMat, period::PetscInt) 
Sets how often the step

Logically Collective

Input Parameters:
- `mat`    - the `MATMFFD` matrix-free matrix
- `period` - 1 for every time, 2 for every second etc

Options Database Key:
- `-mat_mffd_period <period>` - Sets how often `h` is recomputed

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatCreateSNESMF()`, `MatMFFDGetH()`,
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`

# External Links
$(_doc_external("Mat/MatMFFDSetPeriod"))
"""
function MatMFFDSetPeriod(petsclib::PetscLibType, mat::PetscMat, period::PetscInt) end

@for_petsc function MatMFFDSetPeriod(petsclib::$UnionPetscLib, mat::PetscMat, period::$PetscInt )

    @chk ccall(
               (:MatMFFDSetPeriod, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               mat, period,
              )


	return nothing
end 

"""
	MatMFFDSetFunctionError(petsclib::PetscLibType,mat::PetscMat, error::PetscReal) 
Sets the error_rel for the approximation of matrix

Logically Collective

Input Parameters:
- `mat`   - the `MATMFFD` matrix-free matrix
- `error` - relative error (should be set to the square root of the relative error in the function evaluations)

Options Database Key:
- `-mat_mffd_err <error_rel>` - Sets error_rel

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatCreateSNESMF()`, `MatMFFDGetH()`, `MatCreateMFFD()`,
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`

# External Links
$(_doc_external("Mat/MatMFFDSetFunctionError"))
"""
function MatMFFDSetFunctionError(petsclib::PetscLibType, mat::PetscMat, error::PetscReal) end

@for_petsc function MatMFFDSetFunctionError(petsclib::$UnionPetscLib, mat::PetscMat, error::$PetscReal )

    @chk ccall(
               (:MatMFFDSetFunctionError, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               mat, error,
              )


	return nothing
end 

"""
	MatMFFDSetHHistory(petsclib::PetscLibType,J::PetscMat, history::Vector{PetscScalar}, nhistory::PetscInt) 
Sets an array to collect a history of the
differencing values (h) computed for the matrix-free product `MATMFFD` matrix

Logically Collective

Input Parameters:
- `J`        - the `MATMFFD` matrix-free matrix
- `history`  - space to hold the history
- `nhistory` - number of entries in history, if more entries are generated than
nhistory, then the later ones are discarded

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDGetH()`, `MatCreateSNESMF()`,
`MatMFFDResetHHistory()`, `MatMFFDSetFunctionError()`

# External Links
$(_doc_external("Mat/MatMFFDSetHHistory"))
"""
function MatMFFDSetHHistory(petsclib::PetscLibType, J::PetscMat, history::Vector{PetscScalar}, nhistory::PetscInt) end

@for_petsc function MatMFFDSetHHistory(petsclib::$UnionPetscLib, J::PetscMat, history::Vector{$PetscScalar}, nhistory::$PetscInt )

    @chk ccall(
               (:MatMFFDSetHHistory, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}, $PetscInt),
               J, history, nhistory,
              )


	return nothing
end 

"""
	MatMFFDResetHHistory(petsclib::PetscLibType,J::PetscMat) 
Resets the counter to zero to begin
collecting a new set of differencing histories for the `MATMFFD` matrix

Logically Collective

Input Parameter:
- `J` - the matrix-free matrix context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDGetH()`, `MatCreateSNESMF()`,
`MatMFFDSetHHistory()`, `MatMFFDSetFunctionError()`

# External Links
$(_doc_external("Mat/MatMFFDResetHHistory"))
"""
function MatMFFDResetHHistory(petsclib::PetscLibType, J::PetscMat) end

@for_petsc function MatMFFDResetHHistory(petsclib::$UnionPetscLib, J::PetscMat )

    @chk ccall(
               (:MatMFFDResetHHistory, $petsc_library),
               PetscErrorCode,
               (CMat,),
               J,
              )


	return nothing
end 

"""
	MatMFFDSetBase(petsclib::PetscLibType,J::PetscMat, U::PetscVec, F::PetscVec) 
Sets the vector `U` at which matrix vector products of the
Jacobian are computed for the `MATMFFD` matrix

Logically Collective

Input Parameters:
- `J` - the `MATMFFD` matrix
- `U` - the vector
- `F` - (optional) vector that contains F(u) if it has been already computed

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMult()`

# External Links
$(_doc_external("Mat/MatMFFDSetBase"))
"""
function MatMFFDSetBase(petsclib::PetscLibType, J::PetscMat, U::PetscVec, F::PetscVec) end

@for_petsc function MatMFFDSetBase(petsclib::$UnionPetscLib, J::PetscMat, U::PetscVec, F::PetscVec )

    @chk ccall(
               (:MatMFFDSetBase, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               J, U, F,
              )


	return nothing
end 

"""
	MatMFFDSetCheckh(petsclib::PetscLibType,J::PetscMat, fun::MatMFFDCheckhFn, ctx::Cvoid) 
Sets a function that checks the computed `h` and adjusts
it to satisfy some criteria for the `MATMFFD` matrix

Logically Collective

Input Parameters:
- `J`   - the `MATMFFD` matrix
- `fun` - the function that checks `h`, see `MatMFFDCheckhFn`
- `ctx` - any context needed by the function

Options Database Keys:
- `-mat_mffd_check_positivity <bool>` - Ensure that U + h*a  is non-negative

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDCheckhFn`, `MatMFFDCheckPositivity()`

# External Links
$(_doc_external("Mat/MatMFFDSetCheckh"))
"""
function MatMFFDSetCheckh(petsclib::PetscLibType, J::PetscMat, fun::MatMFFDCheckhFn, ctx::Cvoid) end

@for_petsc function MatMFFDSetCheckh(petsclib::$UnionPetscLib, J::PetscMat, fun::MatMFFDCheckhFn, ctx::Cvoid )

    @chk ccall(
               (:MatMFFDSetCheckh, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatMFFDCheckhFn}, Ptr{Cvoid}),
               J, fun, ctx,
              )


	return nothing
end 

"""
	MatMFFDCheckPositivity(petsclib::PetscLibType,dummy::Cvoid, U::PetscVec, a::PetscVec, h::PetscScalar) 
Checks that all entries in U + h*a  are positive or
zero, decreases `h` until this is satisfied for a `MATMFFD` matrix

Logically Collective

Input Parameters:
- `dummy` - context variable (unused)
- `U`     - base vector that is added to
- `a`     - vector that is added
- `h`     - scaling factor on `a`, may be changed on output

Options Database Keys:
- `-mat_mffd_check_positivity <bool>` - Ensure that U + h*a is nonnegative

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatMFFDSetCheckh()`

# External Links
$(_doc_external("Mat/MatMFFDCheckPositivity"))
"""
function MatMFFDCheckPositivity(petsclib::PetscLibType, dummy::Cvoid, U::PetscVec, a::PetscVec, h::PetscScalar) end

@for_petsc function MatMFFDCheckPositivity(petsclib::$UnionPetscLib, dummy::Cvoid, U::PetscVec, a::PetscVec, h::$PetscScalar )

    @chk ccall(
               (:MatMFFDCheckPositivity, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, CVec, CVec, Ptr{$PetscScalar}),
               dummy, U, a, h,
              )


	return nothing
end 

"""
	MatMFFDDSSetUmin(petsclib::PetscLibType,A::PetscMat, umin::PetscReal) 
Sets the "umin" parameter used by the
PETSc routine for computing the differencing parameter, h, which is used
for matrix-free Jacobian-vector products for a `MATMFFD` matrix.

Input Parameters:
- `A`    - the `MATMFFD` matrix
- `umin` - the parameter

Level: advanced

-seealso: `MATMFFD`, `MatMFFDSetFunctionError()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Mat/MatMFFDDSSetUmin"))
"""
function MatMFFDDSSetUmin(petsclib::PetscLibType, A::PetscMat, umin::PetscReal) end

@for_petsc function MatMFFDDSSetUmin(petsclib::$UnionPetscLib, A::PetscMat, umin::$PetscReal )

    @chk ccall(
               (:MatMFFDDSSetUmin, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               A, umin,
              )


	return nothing
end 

"""
	MatMFFDComputeJacobian(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, jac::PetscMat, B::PetscMat, dummy::Cvoid) 
Tells the matrix
Jacobian matrix-vector products will be computed at, i.e. J(x) * a. The x is obtained
from the `SNES` object (using `SNESGetSolution()`).

Collective

Input Parameters:
- `snes`  - the nonlinear solver context
- `x`     - the point at which the Jacobian-vector products will be performed
- `jac`   - the matrix-free Jacobian object of `MatType` `MATMFFD`, likely obtained with `MatCreateSNESMF()`
- `B`     - either the same as `jac` or another matrix type (ignored)
- `dummy` - the user context (ignored)

Options Database Key:
- `-snes_mf` - use the matrix created with `MatSNESMFCreate()` to setup the Jacobian for each new solution in the Newton process

Level: developer

-seealso: [](ch_snes), `MatMFFDGetH()`, `MatCreateSNESMF()`, `MatMFFDSetBase()`, `MatCreateMFFD()`, `MATMFFD`,
`MatMFFDSetHHistory()`, `MatMFFDSetFunctionError()`, `SNESSetJacobian()`

# External Links
$(_doc_external("Snes/MatMFFDComputeJacobian"))
"""
function MatMFFDComputeJacobian(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, jac::PetscMat, B::PetscMat, dummy::Cvoid) end

@for_petsc function MatMFFDComputeJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, jac::PetscMat, B::PetscMat, dummy::Cvoid )

    @chk ccall(
               (:MatMFFDComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x, jac, B, dummy,
              )


	return nothing
end 

"""
	MatCoarsenRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new sparse matrix coarsening algorithm to the matrix package.

Logically Collective, No Fortran Support

Input Parameters:
- `sname`    - name of coarsen (for example `MATCOARSENMIS`)
- `function` - function pointer that creates the coarsen type

Level: developer

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenSetType()`, `MatCoarsenCreate()`, `MatCoarsenRegisterDestroy()`, `MatCoarsenRegisterAll()`

# External Links
$(_doc_external("Mat/MatCoarsenRegister"))
"""
function MatCoarsenRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatCoarsenRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatCoarsenRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	type::MatCoarsenType = MatCoarsenGetType(petsclib::PetscLibType,coarsen::MatCoarsen) 
Gets the Coarsen method type and name (as a string)
from the coarsen context.

Not Collective

Input Parameter:
- `coarsen` - the coarsen context

Output Parameter:
- `type` - coarsener type

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenCreate()`, `MatCoarsenType`, `MatCoarsenSetType()`, `MatCoarsenRegister()`

# External Links
$(_doc_external("Mat/MatCoarsenGetType"))
"""
function MatCoarsenGetType(petsclib::PetscLibType, coarsen::MatCoarsen) end

@for_petsc function MatCoarsenGetType(petsclib::$UnionPetscLib, coarsen::MatCoarsen )
	type_ = Ref{MatCoarsenType}()

    @chk ccall(
               (:MatCoarsenGetType, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, Ptr{MatCoarsenType}),
               coarsen, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	MatCoarsenApply(petsclib::PetscLibType,coarser::MatCoarsen) 
Gets a coarsen for a matrix.

Collective

Input Parameter:
- `coarser` - the coarsen

Options Database Keys:
- `-mat_coarsen_type mis|hem|misk` - mis: maximal independent set based; misk: distance k MIS; hem: heavy edge matching
- `-mat_coarsen_view`              - view the coarsening object

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenSetFromOptions()`, `MatCoarsenSetType()`, `MatCoarsenRegister()`, `MatCoarsenCreate()`,
`MatCoarsenDestroy()`, `MatCoarsenSetAdjacency()`
`MatCoarsenGetData()`

# External Links
$(_doc_external("Mat/MatCoarsenApply"))
"""
function MatCoarsenApply(petsclib::PetscLibType, coarser::MatCoarsen) end

@for_petsc function MatCoarsenApply(petsclib::$UnionPetscLib, coarser::MatCoarsen )

    @chk ccall(
               (:MatCoarsenApply, $petsc_library),
               PetscErrorCode,
               (MatCoarsen,),
               coarser,
              )


	return nothing
end 

"""
	MatCoarsenSetAdjacency(petsclib::PetscLibType,agg::MatCoarsen, adj::PetscMat) 
Sets the adjacency graph (matrix) of the thing to be coarsened.

Collective

Input Parameters:
- `agg` - the coarsen context
- `adj` - the adjacency matrix

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenSetFromOptions()`, `Mat`, `MatCoarsenCreate()`, `MatCoarsenApply()`

# External Links
$(_doc_external("Mat/MatCoarsenSetAdjacency"))
"""
function MatCoarsenSetAdjacency(petsclib::PetscLibType, agg::MatCoarsen, adj::PetscMat) end

@for_petsc function MatCoarsenSetAdjacency(petsclib::$UnionPetscLib, agg::MatCoarsen, adj::PetscMat )

    @chk ccall(
               (:MatCoarsenSetAdjacency, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, CMat),
               agg, adj,
              )


	return nothing
end 

"""
	MatCoarsenSetStrictAggs(petsclib::PetscLibType,agg::MatCoarsen, str::PetscBool) 
Set whether to keep strict (non overlapping) aggregates in the linked list of aggregates for a coarsen context

Logically Collective

Input Parameters:
- `agg` - the coarsen context
- `str` - `PETSC_TRUE` keep strict aggregates, `PETSC_FALSE` allow overlap

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenCreate()`, `MatCoarsenSetFromOptions()`

# External Links
$(_doc_external("Mat/MatCoarsenSetStrictAggs"))
"""
function MatCoarsenSetStrictAggs(petsclib::PetscLibType, agg::MatCoarsen, str::PetscBool) end

@for_petsc function MatCoarsenSetStrictAggs(petsclib::$UnionPetscLib, agg::MatCoarsen, str::PetscBool )

    @chk ccall(
               (:MatCoarsenSetStrictAggs, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, PetscBool),
               agg, str,
              )


	return nothing
end 

"""
	MatCoarsenDestroy(petsclib::PetscLibType,agg::MatCoarsen) 
Destroys the coarsen context.

Collective

Input Parameter:
- `agg` - the coarsen context

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenCreate()`

# External Links
$(_doc_external("Mat/MatCoarsenDestroy"))
"""
function MatCoarsenDestroy(petsclib::PetscLibType, agg::MatCoarsen) end

@for_petsc function MatCoarsenDestroy(petsclib::$UnionPetscLib, agg::MatCoarsen )

    @chk ccall(
               (:MatCoarsenDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MatCoarsen},),
               agg,
              )


	return nothing
end 

"""
	MatCoarsenViewFromOptions(petsclib::PetscLibType,A::MatCoarsen, obj::PetscObject, name::String) 
View the coarsener from the options database

Collective

Input Parameters:
- `A`    - the coarsen context
- `obj`  - Optional object that provides the prefix for the option name
- `name` - command line option (usually `-mat_coarsen_view`)

Options Database Key:
- `-mat_coarsen_view [viewertype]:...` - the viewer and its options

-seealso: `MatCoarsen`, `MatCoarsenView`, `PetscObjectViewFromOptions()`, `MatCoarsenCreate()`

# External Links
$(_doc_external("Mat/MatCoarsenViewFromOptions"))
"""
function MatCoarsenViewFromOptions(petsclib::PetscLibType, A::MatCoarsen, obj::PetscObject, name::String) end

@for_petsc function MatCoarsenViewFromOptions(petsclib::$UnionPetscLib, A::MatCoarsen, obj::PetscObject, name::String )

    @chk ccall(
               (:MatCoarsenViewFromOptions, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	MatCoarsenView(petsclib::PetscLibType,agg::MatCoarsen, viewer::PetscViewer) 
Prints the coarsen data structure.

Collective

Input Parameters:
- `agg`    - the coarsen context
- `viewer` - optional visualization context

For viewing the options database see `MatCoarsenViewFromOptions()`

Level: advanced

-seealso: `MatCoarsen`, `PetscViewer`, `PetscViewerASCIIOpen()`, `MatCoarsenViewFromOptions`

# External Links
$(_doc_external("Mat/MatCoarsenView"))
"""
function MatCoarsenView(petsclib::PetscLibType, agg::MatCoarsen, viewer::PetscViewer) end

@for_petsc function MatCoarsenView(petsclib::$UnionPetscLib, agg::MatCoarsen, viewer::PetscViewer )

    @chk ccall(
               (:MatCoarsenView, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, PetscViewer),
               agg, viewer,
              )


	return nothing
end 

"""
	MatCoarsenSetType(petsclib::PetscLibType,coarser::MatCoarsen, type::MatCoarsenType) 
Sets the type of aggregator to use

Collective

Input Parameters:
- `coarser` - the coarsen context.
- `type`    - a known coarsening method

Options Database Key:
- `-mat_coarsen_type  <type>` - maximal independent set based; distance k MIS; heavy edge matching

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenCreate()`, `MatCoarsenApply()`, `MatCoarsenType`, `MatCoarsenGetType()`

# External Links
$(_doc_external("Mat/MatCoarsenSetType"))
"""
function MatCoarsenSetType(petsclib::PetscLibType, coarser::MatCoarsen, type::MatCoarsenType) end

@for_petsc function MatCoarsenSetType(petsclib::$UnionPetscLib, coarser::MatCoarsen, type::MatCoarsenType )

    @chk ccall(
               (:MatCoarsenSetType, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, MatCoarsenType),
               coarser, type,
              )


	return nothing
end 

"""
	MatCoarsenSetGreedyOrdering(petsclib::PetscLibType,coarser::MatCoarsen, perm::IS) 
Sets the ordering of the vertices to use with a greedy coarsening method

Logically Collective

Input Parameters:
- `coarser` - the coarsen context
- `perm`    - vertex ordering of (greedy) algorithm

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenCreate()`, `MatCoarsenSetType()`

# External Links
$(_doc_external("Mat/MatCoarsenSetGreedyOrdering"))
"""
function MatCoarsenSetGreedyOrdering(petsclib::PetscLibType, coarser::MatCoarsen, perm::IS) end

@for_petsc function MatCoarsenSetGreedyOrdering(petsclib::$UnionPetscLib, coarser::MatCoarsen, perm::IS )

    @chk ccall(
               (:MatCoarsenSetGreedyOrdering, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, CIS),
               coarser, perm,
              )


	return nothing
end 

"""
	MatCoarsenGetData(petsclib::PetscLibType,coarser::MatCoarsen, llist::PetscCoarsenData) 
Gets the weights for vertices for a coarsener.

Logically Collective, No Fortran Support

Input Parameter:
- `coarser` - the coarsen context

Output Parameter:
- `llist` - linked list of aggregates

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenApply()`, `MatCoarsenCreate()`, `MatCoarsenSetType()`, `PetscCoarsenData`

# External Links
$(_doc_external("Mat/MatCoarsenGetData"))
"""
function MatCoarsenGetData(petsclib::PetscLibType, coarser::MatCoarsen, llist::PetscCoarsenData) end

@for_petsc function MatCoarsenGetData(petsclib::$UnionPetscLib, coarser::MatCoarsen, llist::PetscCoarsenData )

    @chk ccall(
               (:MatCoarsenGetData, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, PetscCoarsenData),
               coarser, llist,
              )


	return nothing
end 

"""
	MatCoarsenSetFromOptions(petsclib::PetscLibType,coarser::MatCoarsen) 
Sets various coarsen options from the options database.

Collective

Input Parameter:
- `coarser` - the coarsen context.

Options Database Key:
- `-mat_coarsen_type  <type>`                                                       - mis: maximal independent set based; misk: distance k MIS; hem: heavy edge matching
- `-mat_coarsen_max_it <its> number of iterations to use in the coarsening process` - see `MatCoarsenSetMaximumIterations()`

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenApply()`, `MatCoarsenCreate()`, `MatCoarsenSetType()`,
`MatCoarsenSetMaximumIterations()`

# External Links
$(_doc_external("Mat/MatCoarsenSetFromOptions"))
"""
function MatCoarsenSetFromOptions(petsclib::PetscLibType, coarser::MatCoarsen) end

@for_petsc function MatCoarsenSetFromOptions(petsclib::$UnionPetscLib, coarser::MatCoarsen )

    @chk ccall(
               (:MatCoarsenSetFromOptions, $petsc_library),
               PetscErrorCode,
               (MatCoarsen,),
               coarser,
              )


	return nothing
end 

"""
	MatCoarsenSetMaximumIterations(petsclib::PetscLibType,coarse::MatCoarsen, n::PetscInt) 
Maximum `MATCOARSENHEM` iterations to use

Logically Collective

Input Parameters:
- `coarse` - the coarsen context
- `n`      - number of HEM iterations

Options Database Key:
- `-mat_coarsen_max_it <default=4>` - Maximum `MATCOARSENHEM` iterations to use

Level: intermediate

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenApply()`, `MatCoarsenCreate()`, `MatCoarsenSetType()`

# External Links
$(_doc_external("Mat/MatCoarsenSetMaximumIterations"))
"""
function MatCoarsenSetMaximumIterations(petsclib::PetscLibType, coarse::MatCoarsen, n::PetscInt) end

@for_petsc function MatCoarsenSetMaximumIterations(petsclib::$UnionPetscLib, coarse::MatCoarsen, n::$PetscInt )

    @chk ccall(
               (:MatCoarsenSetMaximumIterations, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, $PetscInt),
               coarse, n,
              )


	return nothing
end 

"""
	MatCoarsenSetStrengthIndex(petsclib::PetscLibType,coarse::MatCoarsen, n::PetscInt, idx::Vector{PetscInt}) 
Index array to use for index to use for strength of connection

Logically Collective

Input Parameters:
- `coarse` - the coarsen context
- `n`      - number of indices
- `idx`    - array of indices

Options Database Key:
- `-mat_coarsen_strength_index` - array of subset of variables per vertex to use for strength norm, -1 for using all (default)

Level: intermediate

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenApply()`, `MatCoarsenCreate()`, `MatCoarsenSetType()`

# External Links
$(_doc_external("Mat/MatCoarsenSetStrengthIndex"))
"""
function MatCoarsenSetStrengthIndex(petsclib::PetscLibType, coarse::MatCoarsen, n::PetscInt, idx::Vector{PetscInt}) end

@for_petsc function MatCoarsenSetStrengthIndex(petsclib::$UnionPetscLib, coarse::MatCoarsen, n::$PetscInt, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatCoarsenSetStrengthIndex, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, $PetscInt, Ptr{$PetscInt}),
               coarse, n, idx,
              )


	return nothing
end 

"""
	MatCoarsenSetThreshold(petsclib::PetscLibType,coarse::MatCoarsen, b::PetscReal) 
Set the threshold for HEM

Logically Collective

Input Parameters:
- `coarse` - the coarsen context
- `b`      - threshold value

Options Database Key:
- `-mat_coarsen_threshold <-1>` - threshold

Level: intermediate

-seealso: `MatCoarsen`, `MatCoarsenType`, `MatCoarsenApply()`, `MatCoarsenCreate()`, `MatCoarsenSetType()`

# External Links
$(_doc_external("Mat/MatCoarsenSetThreshold"))
"""
function MatCoarsenSetThreshold(petsclib::PetscLibType, coarse::MatCoarsen, b::PetscReal) end

@for_petsc function MatCoarsenSetThreshold(petsclib::$UnionPetscLib, coarse::MatCoarsen, b::$PetscReal )

    @chk ccall(
               (:MatCoarsenSetThreshold, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, $PetscReal),
               coarse, b,
              )


	return nothing
end 

"""
	newcrs::MatCoarsen = MatCoarsenCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a coarsen context.

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newcrs` - location to put the context

Level: advanced

-seealso: `MatCoarsen`, `MatCoarsenSetType()`, `MatCoarsenApply()`, `MatCoarsenDestroy()`,
`MatCoarsenSetAdjacency()`, `MatCoarsenGetData()`

# External Links
$(_doc_external("Mat/MatCoarsenCreate"))
"""
function MatCoarsenCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function MatCoarsenCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newcrs_ = Ref{MatCoarsen}()

    @chk ccall(
               (:MatCoarsenCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{MatCoarsen}),
               comm, newcrs_,
              )

	newcrs = newcrs_[]

	return newcrs
end 

"""
	MatCoarsenMISKSetDistance(petsclib::PetscLibType,crs::MatCoarsen, k::PetscInt) 
the distance to be used by MISK

Collective

Input Parameters:
- `crs` - the coarsen
- `k`   - the distance

Options Database Key:
- `-mat_coarsen_misk_distance <k>` - distance for MIS

Level: advanced

-seealso: `MATCOARSENMISK`, `MatCoarsen`, `MatCoarsenSetFromOptions()`, `MatCoarsenSetType()`, `MatCoarsenRegister()`, `MatCoarsenCreate()`,
`MatCoarsenDestroy()`, `MatCoarsenSetAdjacency()`, `MatCoarsenMISKGetDistance()`
`MatCoarsenGetData()`

# External Links
$(_doc_external("Mat/MatCoarsenMISKSetDistance"))
"""
function MatCoarsenMISKSetDistance(petsclib::PetscLibType, crs::MatCoarsen, k::PetscInt) end

@for_petsc function MatCoarsenMISKSetDistance(petsclib::$UnionPetscLib, crs::MatCoarsen, k::$PetscInt )

    @chk ccall(
               (:MatCoarsenMISKSetDistance, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, $PetscInt),
               crs, k,
              )


	return nothing
end 

"""
	k::PetscInt = MatCoarsenMISKGetDistance(petsclib::PetscLibType,crs::MatCoarsen) 
gets the distance to be used by MISK

Collective

Input Parameter:
- `crs` - the coarsen

Output Parameter:
- `k` - the distance

Level: advanced

-seealso: `MATCOARSENMISK`, `MatCoarsen`, `MatCoarsenSetFromOptions()`, `MatCoarsenSetType()`,
`MatCoarsenRegister()`, `MatCoarsenCreate()`, `MatCoarsenDestroy()`,
`MatCoarsenSetAdjacency()`, `MatCoarsenGetData()`

# External Links
$(_doc_external("Mat/MatCoarsenMISKGetDistance"))
"""
function MatCoarsenMISKGetDistance(petsclib::PetscLibType, crs::MatCoarsen) end

@for_petsc function MatCoarsenMISKGetDistance(petsclib::$UnionPetscLib, crs::MatCoarsen )
	k_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatCoarsenMISKGetDistance, $petsc_library),
               PetscErrorCode,
               (MatCoarsen, Ptr{$PetscInt}),
               crs, k_,
              )

	k = k_[]

	return k
end 

