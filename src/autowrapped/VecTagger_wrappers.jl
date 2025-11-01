# autodefined type arguments for class ------
mutable struct _n_VecTagger end
const VecTagger = Ptr{_n_VecTagger}

# -------------------------------------------------------
"""
	tagger::VecTagger = VecTaggerCreate(petsclib::PetscLibType,comm::MPI_Comm) 
create a `VecTagger` context.

Collective

Input Parameter:
- `comm` - communicator on which the `VecTagger` will operate

Output Parameter:
- `tagger` - new Vec tagger context

Level: advanced

-seealso: `VecTagger`, `VecTaggerSetBlockSize()`, `VecTaggerSetFromOptions()`, `VecTaggerSetUp()`, `VecTaggerComputeIS()`, `VecTaggerComputeBoxes()`, `VecTaggerDestroy()`

# External Links
$(_doc_external("Vec/VecTaggerCreate"))
"""
function VecTaggerCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function VecTaggerCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	tagger_ = Ref{VecTagger}()

    @chk ccall(
               (:VecTaggerCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{VecTagger}),
               comm, tagger_,
              )

	tagger = tagger_[]

	return tagger
end 

"""
	VecTaggerSetType(petsclib::PetscLibType,tagger::VecTagger, type::VecTaggerType) 
set the Vec tagger implementation

Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `type`   - a known method

Options Database Key:
- `-vec_tagger_type <type>` - Sets the method; use -help for a list
of available methods (for instance, absolute, relative, cdf, or, and)

Level: advanced

-seealso: `VecTaggerType`, `VecTaggerCreate()`, `VecTagger`

# External Links
$(_doc_external("Vec/VecTaggerSetType"))
"""
function VecTaggerSetType(petsclib::PetscLibType, tagger::VecTagger, type::VecTaggerType) end

@for_petsc function VecTaggerSetType(petsclib::$UnionPetscLib, tagger::VecTagger, type::VecTaggerType )

    @chk ccall(
               (:VecTaggerSetType, $petsc_library),
               PetscErrorCode,
               (VecTagger, VecTaggerType),
               tagger, type,
              )


	return nothing
end 

"""
	type::VecTaggerType = VecTaggerGetType(petsclib::PetscLibType,tagger::VecTagger) 
Gets the `VecTaggerType` name (as a string) from the `VecTagger`.

Not Collective

Input Parameter:
- `tagger` - The `VecTagger` context

Output Parameter:
- `type` - The `VecTagger` type name

Level: advanced

-seealso: `VecTaggerSetType()`, `VecTaggerCreate()`, `VecTaggerSetFromOptions()`, `VecTagger`, `VecTaggerType`

# External Links
$(_doc_external("Vec/VecTaggerGetType"))
"""
function VecTaggerGetType(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerGetType(petsclib::$UnionPetscLib, tagger::VecTagger )
	type_ = Ref{VecTaggerType}()

    @chk ccall(
               (:VecTaggerGetType, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{VecTaggerType}),
               tagger, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	VecTaggerDestroy(petsclib::PetscLibType,tagger::VecTagger) 
destroy a `VecTagger` context

Collective

Input Parameter:
- `tagger` - address of tagger

Level: advanced

-seealso: `VecTaggerCreate()`, `VecTaggerSetType()`, `VecTagger`

# External Links
$(_doc_external("Vec/VecTaggerDestroy"))
"""
function VecTaggerDestroy(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerDestroy(petsclib::$UnionPetscLib, tagger::VecTagger )

    @chk ccall(
               (:VecTaggerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{VecTagger},),
               tagger,
              )


	return nothing
end 

"""
	VecTaggerSetUp(petsclib::PetscLibType,tagger::VecTagger) 
set up a `VecTagger` context

Collective

Input Parameter:
- `tagger` - Vec tagger object

Level: advanced

-seealso: `VecTaggerSetFromOptions()`, `VecTaggerSetType()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerSetUp"))
"""
function VecTaggerSetUp(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerSetUp(petsclib::$UnionPetscLib, tagger::VecTagger )

    @chk ccall(
               (:VecTaggerSetUp, $petsc_library),
               PetscErrorCode,
               (VecTagger,),
               tagger,
              )


	return nothing
end 

"""
	VecTaggerSetFromOptions(petsclib::PetscLibType,tagger::VecTagger) 
set `VecTagger` options using the options database

Logically Collective

Input Parameter:
- `tagger` - vec tagger

Options Database Keys:
- `-vec_tagger_type`       - implementation type, see `VecTaggerSetType()`
- `-vec_tagger_block_size` - set the block size, see `VecTaggerSetBlockSize()`
- `-vec_tagger_invert`     - invert the index set returned by `VecTaggerComputeIS()`

Level: advanced

-seealso: `VecTagger`, `VecTaggerCreate()`, `VecTaggerSetUp()`


# External Links
$(_doc_external("Vec/VecTaggerSetFromOptions"))
"""
function VecTaggerSetFromOptions(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerSetFromOptions(petsclib::$UnionPetscLib, tagger::VecTagger )

    @chk ccall(
               (:VecTaggerSetFromOptions, $petsc_library),
               PetscErrorCode,
               (VecTagger,),
               tagger,
              )


	return nothing
end 

"""
	VecTaggerSetBlockSize(petsclib::PetscLibType,tagger::VecTagger, blocksize::PetscInt) 
set the block size of the set of indices returned by `VecTaggerComputeIS()`.

Logically Collective

Input Parameters:
- `tagger`    - vec tagger
- `blocksize` - block size of the criteria used to tagger vectors

Level: advanced

-seealso: `VecTaggerComputeIS()`, `VecTaggerGetBlockSize()`, `VecSetBlockSize()`, `VecGetBlockSize()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerSetBlockSize"))
"""
function VecTaggerSetBlockSize(petsclib::PetscLibType, tagger::VecTagger, blocksize::PetscInt) end

@for_petsc function VecTaggerSetBlockSize(petsclib::$UnionPetscLib, tagger::VecTagger, blocksize::$PetscInt )

    @chk ccall(
               (:VecTaggerSetBlockSize, $petsc_library),
               PetscErrorCode,
               (VecTagger, $PetscInt),
               tagger, blocksize,
              )


	return nothing
end 

"""
	blocksize::PetscInt = VecTaggerGetBlockSize(petsclib::PetscLibType,tagger::VecTagger) 
get the block size of the indices created by `VecTaggerComputeIS()`.

Logically Collective

Input Parameter:
- `tagger` - vec tagger

Output Parameter:
- `blocksize` - block size of the vectors the tagger operates on

Level: advanced

-seealso: `VecTaggerComputeIS()`, `VecTaggerSetBlockSize()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerGetBlockSize"))
"""
function VecTaggerGetBlockSize(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerGetBlockSize(petsclib::$UnionPetscLib, tagger::VecTagger )
	blocksize_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecTaggerGetBlockSize, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{$PetscInt}),
               tagger, blocksize_,
              )

	blocksize = blocksize_[]

	return blocksize
end 

"""
	VecTaggerSetInvert(petsclib::PetscLibType,tagger::VecTagger, invert::PetscBool) 
If the tagged index sets are based on boxes that can be returned by `VecTaggerComputeBoxes()`,
then this option inverts values used to compute the IS, i.e., from being in the union of the boxes to being in the
intersection of their exteriors.

Logically Collective

Input Parameters:
- `tagger` - vec tagger
- `invert` - `PETSC_TRUE` to invert, `PETSC_FALSE` to use the indices as is

Level: advanced

-seealso: `VecTaggerComputeIS()`, `VecTaggerGetInvert()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerSetInvert"))
"""
function VecTaggerSetInvert(petsclib::PetscLibType, tagger::VecTagger, invert::PetscBool) end

@for_petsc function VecTaggerSetInvert(petsclib::$UnionPetscLib, tagger::VecTagger, invert::PetscBool )

    @chk ccall(
               (:VecTaggerSetInvert, $petsc_library),
               PetscErrorCode,
               (VecTagger, PetscBool),
               tagger, invert,
              )


	return nothing
end 

"""
	invert::PetscBool = VecTaggerGetInvert(petsclib::PetscLibType,tagger::VecTagger) 
get whether the set of indices returned by `VecTaggerComputeIS()` are inverted

Logically Collective

Input Parameter:
- `tagger` - vec tagger

Output Parameter:
- `invert` - `PETSC_TRUE` to invert, `PETSC_FALSE` to use the indices as is

Level: advanced

-seealso: `VecTaggerComputeIS()`, `VecTaggerSetInvert()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerGetInvert"))
"""
function VecTaggerGetInvert(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerGetInvert(petsclib::$UnionPetscLib, tagger::VecTagger )
	invert_ = Ref{PetscBool}()

    @chk ccall(
               (:VecTaggerGetInvert, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{PetscBool}),
               tagger, invert_,
              )

	invert = invert_[]

	return invert
end 

"""
	VecTaggerView(petsclib::PetscLibType,tagger::VecTagger, viewer::PetscViewer) 
view a `VecTagger` context

Collective

Input Parameters:
- `tagger` - vec tagger
- `viewer` - viewer to display tagger, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: advanced

-seealso: `VecTaggerCreate()`, `VecTagger`

# External Links
$(_doc_external("Vec/VecTaggerView"))
"""
function VecTaggerView(petsclib::PetscLibType, tagger::VecTagger, viewer::PetscViewer) end

@for_petsc function VecTaggerView(petsclib::$UnionPetscLib, tagger::VecTagger, viewer::PetscViewer )

    @chk ccall(
               (:VecTaggerView, $petsc_library),
               PetscErrorCode,
               (VecTagger, PetscViewer),
               tagger, viewer,
              )


	return nothing
end 

"""
	numBoxes::PetscInt,listed::PetscBool = VecTaggerComputeBoxes(petsclib::PetscLibType,tagger::VecTagger, vec::PetscVec, boxes::Vector{VecTaggerBox}) 
If the tagged index set can be summarized as a list of boxes of values, returns that list, otherwise returns
in listed `PETSC_FALSE`

Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `vec`    - the vec to tag

Output Parameters:
- `numBoxes` - the number of boxes in the tag definition
- `boxes`    - a newly allocated list of boxes.  This is a flat array of (BlockSize * `numBoxe`s) pairs that the user can free with `PetscFree()`.
- `listed`   - `PETSC_TRUE` if a list was created, pass in `NULL` if not needed

Level: advanced

-seealso: `VecTaggerComputeIS()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerComputeBoxes"))
"""
function VecTaggerComputeBoxes(petsclib::PetscLibType, tagger::VecTagger, vec::PetscVec, boxes::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerComputeBoxes(petsclib::$UnionPetscLib, tagger::VecTagger, vec::PetscVec, boxes::Vector{VecTaggerBox} )
	numBoxes_ = Ref{$PetscInt}()
	boxes_ = Ref(pointer(boxes))
	listed_ = Ref{PetscBool}()

    @chk ccall(
               (:VecTaggerComputeBoxes, $petsc_library),
               PetscErrorCode,
               (VecTagger, CVec, Ptr{$PetscInt}, Ptr{Ptr{VecTaggerBox}}, Ptr{PetscBool}),
               tagger, vec, numBoxes_, boxes_, listed_,
              )

	numBoxes = numBoxes_[]
	listed = listed_[]

	return numBoxes,listed
end 

"""
	listed::PetscBool = VecTaggerComputeIS(petsclib::PetscLibType,tagger::VecTagger, vec::PetscVec, is::Vector{IS}) 
Use a `VecTagger` context to tag a set of indices based on a vector's values

Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `vec`    - the vec to tag

Output Parameters:
- `is`     - a list of the indices tagged by the tagger, i.e., if the number of local indices will be n / bs, where n is `VecGetLocalSize()` and bs is `VecTaggerGetBlockSize()`.
- `listed` - routine was able to compute the `IS`, pass in `NULL` if not needed

Level: advanced

-seealso: `VecTaggerComputeBoxes()`, `VecTagger`, `VecTaggerCreate()`

# External Links
$(_doc_external("Vec/VecTaggerComputeIS"))
"""
function VecTaggerComputeIS(petsclib::PetscLibType, tagger::VecTagger, vec::PetscVec, is::Vector{IS}) end

@for_petsc function VecTaggerComputeIS(petsclib::$UnionPetscLib, tagger::VecTagger, vec::PetscVec, is::Vector{IS} )
	listed_ = Ref{PetscBool}()

    @chk ccall(
               (:VecTaggerComputeIS, $petsc_library),
               PetscErrorCode,
               (VecTagger, CVec, Ptr{CIS}, Ptr{PetscBool}),
               tagger, vec, is, listed_,
              )

	listed = listed_[]

	return listed
end 

"""
	VecTaggerInitializePackage(petsclib::PetscLibType) 
Initialize VecTagger package

Logically Collective

Level: developer

-seealso: `VecTaggerFinalizePackage()`

# External Links
$(_doc_external("Vec/VecTaggerInitializePackage"))
"""
function VecTaggerInitializePackage(petsclib::PetscLibType) end

@for_petsc function VecTaggerInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecTaggerInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecTaggerFinalizePackage(petsclib::PetscLibType) 
Finalize VecTagger package, it is called from PetscFinalize()

Logically Collective

Level: developer

-seealso: `VecTaggerInitializePackage()`

# External Links
$(_doc_external("Vec/VecTaggerFinalizePackage"))
"""
function VecTaggerFinalizePackage(petsclib::PetscLibType) end

@for_petsc function VecTaggerFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecTaggerFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecTaggerRegisterAll(petsclib::PetscLibType) 
Registers all the `VecTagger` communication implementations

Not Collective

Level: advanced

-seealso: `VecTaggerRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecTaggerRegisterAll"))
"""
function VecTaggerRegisterAll(petsclib::PetscLibType) end

@for_petsc function VecTaggerRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecTaggerRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecTaggerRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds an implementation of the `VecTagger` communication protocol.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined implementation
- `function` - routine to create method context

Level: advanced

-seealso: `VecTaggerType`, `VecTaggerCreate()`, `VecTagger`, `VecTaggerRegisterAll()`, `VecTaggerRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecTaggerRegister"))
"""
function VecTaggerRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function VecTaggerRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:VecTaggerRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	nsubs::PetscInt = VecTaggerOrGetSubs(petsclib::PetscLibType,tagger::VecTagger, subs::Vector{VecTagger}) 
Get the sub `VecTagger`s whose union defines the outer `VecTagger`

Not Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameters:
- `nsubs` - the number of sub `VecTagger`s
- `subs`  - the sub `VecTagger`s

Level: advanced

-seealso: `VecTaggerOrSetSubs()`

# External Links
$(_doc_external("Vec/VecTaggerOrGetSubs"))
"""
function VecTaggerOrGetSubs(petsclib::PetscLibType, tagger::VecTagger, subs::Vector{VecTagger}) end

@for_petsc function VecTaggerOrGetSubs(petsclib::$UnionPetscLib, tagger::VecTagger, subs::Vector{VecTagger} )
	nsubs_ = Ref{$PetscInt}()
	subs_ = Ref(pointer(subs))

    @chk ccall(
               (:VecTaggerOrGetSubs, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{$PetscInt}, Ptr{Ptr{VecTagger}}),
               tagger, nsubs_, subs_,
              )

	nsubs = nsubs_[]

	return nsubs
end 

"""
	VecTaggerOrSetSubs(petsclib::PetscLibType,tagger::VecTagger, nsubs::PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode) 
Set the sub `VecTagger`s whose union defines the outer `VecTagger`

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `nsubs`  - the number of sub `VecTagger`s
- `subs`   - the sub `VecTagger`s
- `mode`   - the copy mode to use for `subs`

Level: advanced

-seealso: `VecTaggetOrGetStubs()`

# External Links
$(_doc_external("Vec/VecTaggerOrSetSubs"))
"""
function VecTaggerOrSetSubs(petsclib::PetscLibType, tagger::VecTagger, nsubs::PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode) end

@for_petsc function VecTaggerOrSetSubs(petsclib::$UnionPetscLib, tagger::VecTagger, nsubs::$PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode )

    @chk ccall(
               (:VecTaggerOrSetSubs, $petsc_library),
               PetscErrorCode,
               (VecTagger, $PetscInt, Ptr{VecTagger}, PetscCopyMode),
               tagger, nsubs, subs, mode,
              )


	return nothing
end 

"""
	VecTaggerCDFSetMethod(petsclib::PetscLibType,tagger::VecTagger, method::VecTaggerCDFMethod) 
Set the method used to compute absolute boxes from CDF boxes

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `method` - the method

Level: advanced

-seealso: `Vec`, `VecTagger`, `VecTaggerCDFMethod`

# External Links
$(_doc_external("Vec/VecTaggerCDFSetMethod"))
"""
function VecTaggerCDFSetMethod(petsclib::PetscLibType, tagger::VecTagger, method::VecTaggerCDFMethod) end

@for_petsc function VecTaggerCDFSetMethod(petsclib::$UnionPetscLib, tagger::VecTagger, method::VecTaggerCDFMethod )

    @chk ccall(
               (:VecTaggerCDFSetMethod, $petsc_library),
               PetscErrorCode,
               (VecTagger, VecTaggerCDFMethod),
               tagger, method,
              )


	return nothing
end 

"""
	VecTaggerCDFGetMethod(petsclib::PetscLibType,tagger::VecTagger, method::VecTaggerCDFMethod) 
Get the method used to compute absolute boxes from CDF boxes

Logically Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameter:
- `method` - the method

Level: advanced

-seealso: `Vec`, `VecTagger`, `VecTaggerCDFMethod`

# External Links
$(_doc_external("Vec/VecTaggerCDFGetMethod"))
"""
function VecTaggerCDFGetMethod(petsclib::PetscLibType, tagger::VecTagger, method::VecTaggerCDFMethod) end

@for_petsc function VecTaggerCDFGetMethod(petsclib::$UnionPetscLib, tagger::VecTagger, method::VecTaggerCDFMethod )

    @chk ccall(
               (:VecTaggerCDFGetMethod, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{VecTaggerCDFMethod}),
               tagger, method,
              )


	return nothing
end 

"""
	VecTaggerCDFIterativeSetTolerances(petsclib::PetscLibType,tagger::VecTagger, maxit::PetscInt, rtol::PetscReal, atol::PetscReal) 
Set the tolerances for iterative computation of absolute boxes from CDF boxes.

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `maxit`  - the maximum number of iterations: 0 indicates the absolute values will be estimated from an initial guess based only on the minimum, maximum, mean,
and standard deviation of the box endpoints.
- `rtol`   - the acceptable relative tolerance in the absolute values from the initial guess
- `atol`   - the acceptable absolute tolerance in the absolute values from the initial guess

Level: advanced

-seealso: `VecTagger`, `VecTaggerCDFSetMethod()`

# External Links
$(_doc_external("Vec/VecTaggerCDFIterativeSetTolerances"))
"""
function VecTaggerCDFIterativeSetTolerances(petsclib::PetscLibType, tagger::VecTagger, maxit::PetscInt, rtol::PetscReal, atol::PetscReal) end

@for_petsc function VecTaggerCDFIterativeSetTolerances(petsclib::$UnionPetscLib, tagger::VecTagger, maxit::$PetscInt, rtol::$PetscReal, atol::$PetscReal )

    @chk ccall(
               (:VecTaggerCDFIterativeSetTolerances, $petsc_library),
               PetscErrorCode,
               (VecTagger, $PetscInt, $PetscReal, $PetscReal),
               tagger, maxit, rtol, atol,
              )


	return nothing
end 

"""
	maxit::PetscInt,rtol::PetscReal,atol::PetscReal = VecTaggerCDFIterativeGetTolerances(petsclib::PetscLibType,tagger::VecTagger) 
Get the tolerances for iterative computation of absolute boxes from CDF boxes.

Logically Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameters:
- `maxit` - the maximum number of iterations: 0 indicates the absolute values will be estimated from an initial guess based only on the minimum, maximum,
mean, and standard deviation of the box endpoints.
- `rtol`  - the acceptable relative tolerance in the absolute values from the initial guess
- `atol`  - the acceptable absolute tolerance in the absolute values from the initial guess

Level: advanced

-seealso: `VecTagger`, `VecTaggerCDFSetMethod()`

# External Links
$(_doc_external("Vec/VecTaggerCDFIterativeGetTolerances"))
"""
function VecTaggerCDFIterativeGetTolerances(petsclib::PetscLibType, tagger::VecTagger) end

@for_petsc function VecTaggerCDFIterativeGetTolerances(petsclib::$UnionPetscLib, tagger::VecTagger )
	maxit_ = Ref{$PetscInt}()
	rtol_ = Ref{$PetscReal}()
	atol_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecTaggerCDFIterativeGetTolerances, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               tagger, maxit_, rtol_, atol_,
              )

	maxit = maxit_[]
	rtol = rtol_[]
	atol = atol_[]

	return maxit,rtol,atol
end 

"""
	VecTaggerCDFSetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Set the cumulative box defining the values to be tagged by the tagger, where cumulative boxes are subsets of [0,1], where 0 indicates the smallest value present in the vector and 1 indicates the largest.

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `box`    - a blocksize array of `VecTaggerBox` boxes

Level: advanced

-seealso: `VecTagger`, `VecTaggerCDFGetBox()`, `VecTaggerBox`

# External Links
$(_doc_external("Vec/VecTaggerCDFSetBox"))
"""
function VecTaggerCDFSetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerCDFSetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )

    @chk ccall(
               (:VecTaggerCDFSetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{VecTaggerBox}),
               tagger, box,
              )


	return nothing
end 

"""
	VecTaggerCDFGetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Get the cumulative box (multi
are subsets of [0,1], where 0 indicates the smallest value present in the vector and 1 indicates the largest.

Logically Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameter:
- `box` - a blocksize array of `VecTaggerBox` boxes

Level: advanced

-seealso: `VecTagger`, `VecTaggerCDFSetBox()`, `VecTaggerBox`

# External Links
$(_doc_external("Vec/VecTaggerCDFGetBox"))
"""
function VecTaggerCDFGetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerCDFGetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )
	box_ = Ref(pointer(box))

    @chk ccall(
               (:VecTaggerCDFGetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{Ptr{VecTaggerBox}}),
               tagger, box_,
              )


	return nothing
end 

"""
	VecTaggerRelativeSetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Set the relative box defining the values to be tagged by the tagger, where relative boxes are subsets of [0,1] (or [0,1]+[0,1]i for complex scalars), where 0 indicates the smallest value present in the vector and 1 indicates the largest.

Logically Collective

Input Parameters:
- `tagger` - the VecTagger context
- `box`    - a blocksize list of VecTaggerBox boxes

Level: advanced

-seealso: `VecTaggerRelativeGetBox()`

# External Links
$(_doc_external("Vec/VecTaggerRelativeSetBox"))
"""
function VecTaggerRelativeSetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerRelativeSetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )

    @chk ccall(
               (:VecTaggerRelativeSetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{VecTaggerBox}),
               tagger, box,
              )


	return nothing
end 

"""
	VecTaggerRelativeGetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Get the relative box defining the values to be tagged by the tagger, where relative boxess are subsets of [0,1] (or [0,1]+[0,1]i for complex scalars), where 0 indicates the smallest value present in the vector and 1 indicates the largest.

Logically Collective

Input Parameter:
- `tagger` - the VecTagger context

Output Parameter:
- `box` - a blocksize list of VecTaggerBox boxes

Level: advanced

-seealso: `VecTaggerRelativeSetBox()`

# External Links
$(_doc_external("Vec/VecTaggerRelativeGetBox"))
"""
function VecTaggerRelativeGetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerRelativeGetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )
	box_ = Ref(pointer(box))

    @chk ccall(
               (:VecTaggerRelativeGetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{Ptr{VecTaggerBox}}),
               tagger, box_,
              )


	return nothing
end 

"""
	nsubs::PetscInt = VecTaggerAndGetSubs(petsclib::PetscLibType,tagger::VecTagger, subs::Vector{VecTagger}) 
Get the sub `VecTagger`s whose intersection defines the outer `VecTagger`

Not Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameters:
- `nsubs` - the number of sub `VecTagger`s
- `subs`  - the sub `VecTagger`s

Level: advanced

-seealso: `VecTagger`, `VecTaggerAndSetSubs()`

# External Links
$(_doc_external("Vec/VecTaggerAndGetSubs"))
"""
function VecTaggerAndGetSubs(petsclib::PetscLibType, tagger::VecTagger, subs::Vector{VecTagger}) end

@for_petsc function VecTaggerAndGetSubs(petsclib::$UnionPetscLib, tagger::VecTagger, subs::Vector{VecTagger} )
	nsubs_ = Ref{$PetscInt}()
	subs_ = Ref(pointer(subs))

    @chk ccall(
               (:VecTaggerAndGetSubs, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{$PetscInt}, Ptr{Ptr{VecTagger}}),
               tagger, nsubs_, subs_,
              )

	nsubs = nsubs_[]

	return nsubs
end 

"""
	VecTaggerAndSetSubs(petsclib::PetscLibType,tagger::VecTagger, nsubs::PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode) 
Set the sub `VecTagger`s whose intersection defines the outer `VecTagger`

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `nsubs`  - the number of sub `VecTagger`s
- `subs`   - the sub `VecTagger`s
- `mode`   - the copy mode to use for `subs`

Level: advanced

-seealso: `VecTagger`

# External Links
$(_doc_external("Vec/VecTaggerAndSetSubs"))
"""
function VecTaggerAndSetSubs(petsclib::PetscLibType, tagger::VecTagger, nsubs::PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode) end

@for_petsc function VecTaggerAndSetSubs(petsclib::$UnionPetscLib, tagger::VecTagger, nsubs::$PetscInt, subs::Vector{VecTagger}, mode::PetscCopyMode )

    @chk ccall(
               (:VecTaggerAndSetSubs, $petsc_library),
               PetscErrorCode,
               (VecTagger, $PetscInt, Ptr{VecTagger}, PetscCopyMode),
               tagger, nsubs, subs, mode,
              )


	return nothing
end 

"""
	VecTaggerAbsoluteSetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Set the box defining the values to be tagged by the tagger.

Logically Collective

Input Parameters:
- `tagger` - the `VecTagger` context
- `box`    - the box: a blocksize array of `VecTaggerBox` boxes

Level: advanced

-seealso: `VecTagger`, `VecTaggerBox`, `VecTaggerAbsoluteGetBox()`

# External Links
$(_doc_external("Vec/VecTaggerAbsoluteSetBox"))
"""
function VecTaggerAbsoluteSetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerAbsoluteSetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )

    @chk ccall(
               (:VecTaggerAbsoluteSetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{VecTaggerBox}),
               tagger, box,
              )


	return nothing
end 

"""
	VecTaggerAbsoluteGetBox(petsclib::PetscLibType,tagger::VecTagger, box::Vector{VecTaggerBox}) 
Get the box defining the values to be tagged by the tagger.

Logically Collective

Input Parameter:
- `tagger` - the `VecTagger` context

Output Parameter:
- `box` - the box: a blocksize array of `VecTaggerBox` boxes

Level: advanced

-seealso: `VecTagger`, `VecTaggerBox`, `VecTaggerAbsoluteSetBox()`

# External Links
$(_doc_external("Vec/VecTaggerAbsoluteGetBox"))
"""
function VecTaggerAbsoluteGetBox(petsclib::PetscLibType, tagger::VecTagger, box::Vector{VecTaggerBox}) end

@for_petsc function VecTaggerAbsoluteGetBox(petsclib::$UnionPetscLib, tagger::VecTagger, box::Vector{VecTaggerBox} )
	box_ = Ref(pointer(box))

    @chk ccall(
               (:VecTaggerAbsoluteGetBox, $petsc_library),
               PetscErrorCode,
               (VecTagger, Ptr{Ptr{VecTaggerBox}}),
               tagger, box_,
              )


	return nothing
end 

