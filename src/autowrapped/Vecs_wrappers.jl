# autodefined type arguments for class ------
mutable struct _n_Vecs end
const Vecs = Ptr{_n_Vecs}

mutable struct n_PetscRandom end
const PetscRandom = Ptr{n_PetscRandom}
# -------------------------------------------------------
"""
	VecScatterPetscToFFTW(petsclib::PetscLibType,A::PetscMat, x::PetscVec, y::PetscVec) 
Copies a PETSc vector to the vector that goes into `MATFFTW` calls.

Collective

Input Parameters:
- `A` - FFTW matrix
- `x` - the PETSc vector

Output Parameter:
- `y` - the FFTW vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATFFTW`, `VecScatterFFTWToPetsc()`, `MatCreateVecsFFTW()`

# External Links
$(_doc_external("Mat/VecScatterPetscToFFTW"))
"""
function VecScatterPetscToFFTW(petsclib::PetscLibType, A::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function VecScatterPetscToFFTW(petsclib::$UnionPetscLib, A::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecScatterPetscToFFTW, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	VecScatterFFTWToPetsc(petsclib::PetscLibType,A::PetscMat, x::PetscVec, y::PetscVec) 
Converts `MATFFTW` output vector to a PETSc vector.

Collective

Input Parameters:
- `A` - `MATFFTW` matrix
- `x` - FFTW vector

Output Parameter:
- `y` - PETSc vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `VecScatterPetscToFFTW()`, `MATFFTW`, `MatCreateVecsFFTW()`

# External Links
$(_doc_external("Mat/VecScatterFFTWToPetsc"))
"""
function VecScatterFFTWToPetsc(petsclib::PetscLibType, A::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function VecScatterFFTWToPetsc(petsclib::$UnionPetscLib, A::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecScatterFFTWToPetsc, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	VecScale(petsclib::PetscLibType,x::PetscVec, alpha::PetscScalar) 
Scales a vector.

Logically Collective

Input Parameters:
- `x`     - the vector
- `alpha` - the scalar

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecSet()`

# External Links
$(_doc_external("Vec/VecScale"))
"""
function VecScale(petsclib::PetscLibType, x::PetscVec, alpha::PetscScalar) end

@for_petsc function VecScale(petsclib::$UnionPetscLib, x::PetscVec, alpha::$PetscScalar )

    @chk ccall(
               (:VecScale, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar),
               x, alpha,
              )


	return nothing
end 

"""
	VecSet(petsclib::PetscLibType,x::PetscVec, alpha::PetscScalar) 
Sets all components of a vector to a single scalar value.

Logically Collective

Input Parameters:
- `x`     - the vector
- `alpha` - the scalar

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecSetValues()`, `VecSetValuesBlocked()`, `VecSetRandom()`

# External Links
$(_doc_external("Vec/VecSet"))
"""
function VecSet(petsclib::PetscLibType, x::PetscVec, alpha::PetscScalar) end

@for_petsc function VecSet(petsclib::$UnionPetscLib, x::PetscVec, alpha::$PetscScalar )

    @chk ccall(
               (:VecSet, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar),
               x, alpha,
              )


	return nothing
end 

"""
	VecSetValues(petsclib::PetscLibType,x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) 
Inserts or adds values into certain locations of a vector.

Not Collective

Input Parameters:
- `x`    - vector to insert in
- `ni`   - number of elements to add
- `ix`   - indices where to add
- `y`    - array of values. Pass `NULL` to set all zeroes.
- `iora` - either `INSERT_VALUES` to replace the current values or `ADD_VALUES` to add values to any existing entries

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValuesLocal()`,
`VecSetValue()`, `VecSetValuesBlocked()`, `InsertMode`, `INSERT_VALUES`, `ADD_VALUES`, `VecGetValues()`

# External Links
$(_doc_external("Vec/VecSetValues"))
"""
VecSetValues(petsclib::PetscLibType, x::PetscVec, ni::PetscInt, ix::Vector, y::Vector, iora::InsertMode)

@for_petsc function VecSetValues(petsclib::$UnionPetscLib, x::AbstractPetscVec, ni::$PetscInt, ix::Vector{$PetscInt}, y::Vector{$PetscScalar}, iora::InsertMode )

    @chk ccall(
               (:VecSetValues, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               x, ni, ix, y, iora,
              )


	return nothing
end 

"""
	VecSetValuesBlocked(petsclib::PetscLibType,x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) 
Inserts or adds blocks of values into certain locations of a vector.

Not Collective

Input Parameters:
- `x`    - vector to insert in
- `ni`   - number of blocks to add
- `ix`   - indices where to add in block count, rather than element count
- `y`    - array of values. Pass `NULL` to set all zeroes.
- `iora` - either `INSERT_VALUES` replaces existing entries with new values, `ADD_VALUES`, adds values to any existing entries

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValuesBlockedLocal()`,
`VecSetValues()`

# External Links
$(_doc_external("Vec/VecSetValuesBlocked"))
"""
function VecSetValuesBlocked(petsclib::PetscLibType, x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) end

@for_petsc function VecSetValuesBlocked(petsclib::$UnionPetscLib, x::PetscVec, ni::$PetscInt, ix::Vector{$PetscInt}, y::Vector{$PetscScalar}, iora::InsertMode )

    @chk ccall(
               (:VecSetValuesBlocked, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               x, ni, ix, y, iora,
              )


	return nothing
end 

"""
	VecSetValuesLocal(petsclib::PetscLibType,x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) 
Inserts or adds values into certain locations of a vector,
using a local ordering of the nodes.

Not Collective

Input Parameters:
- `x`    - vector to insert in
- `ni`   - number of elements to add
- `ix`   - indices where to add
- `y`    - array of values. Pass `NULL` to set all zeroes.
- `iora` - either `INSERT_VALUES` replaces existing entries with new values, `ADD_VALUES` adds values to any existing entries

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValues()`, `VecSetLocalToGlobalMapping()`,
`VecSetValuesBlockedLocal()`

# External Links
$(_doc_external("Vec/VecSetValuesLocal"))
"""
function VecSetValuesLocal(petsclib::PetscLibType, x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) end

@for_petsc function VecSetValuesLocal(petsclib::$UnionPetscLib, x::PetscVec, ni::$PetscInt, ix::Vector{$PetscInt}, y::Vector{$PetscScalar}, iora::InsertMode )

    @chk ccall(
               (:VecSetValuesLocal, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               x, ni, ix, y, iora,
              )


	return nothing
end 

"""
	VecSetValuesBlockedLocal(petsclib::PetscLibType,x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) 
Inserts or adds values into certain locations of a vector,
using a local ordering of the nodes.

Not Collective

Input Parameters:
- `x`    - vector to insert in
- `ni`   - number of blocks to add
- `ix`   - indices where to add in block count, not element count
- `y`    - array of values. Pass `NULL` to set all zeroes.
- `iora` - either `INSERT_VALUES` replaces existing entries with new values, `ADD_VALUES` adds values to any existing entries

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValues()`, `VecSetValuesBlocked()`,
`VecSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Vec/VecSetValuesBlockedLocal"))
"""
function VecSetValuesBlockedLocal(petsclib::PetscLibType, x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}, y::Vector{PetscScalar}, iora::InsertMode) end

@for_petsc function VecSetValuesBlockedLocal(petsclib::$UnionPetscLib, x::PetscVec, ni::$PetscInt, ix::Vector{$PetscInt}, y::Vector{$PetscScalar}, iora::InsertMode )

    @chk ccall(
               (:VecSetValuesBlockedLocal, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               x, ni, ix, y, iora,
              )


	return nothing
end 

"""
	nstash::PetscInt,reallocs::PetscInt,bnstash::PetscInt,breallocs::PetscInt = VecStashGetInfo(petsclib::PetscLibType,vec::PetscVec) 
Gets how many values are currently in the vector stash, i.e. need
to be communicated to other processors during the `VecAssemblyBegin()`/`VecAssemblyEnd()` process

Not Collective

Input Parameter:
- `vec` - the vector

Output Parameters:
- `nstash`    - the size of the stash
- `reallocs`  - the number of additional mallocs incurred in building the stash
- `bnstash`   - the size of the block stash
- `breallocs` - the number of additional mallocs incurred in building the block stash (from `VecSetValuesBlocked()`)

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecStashSetInitialSize()`, `VecStashView()`

# External Links
$(_doc_external("Vec/VecStashGetInfo"))
"""
function VecStashGetInfo(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecStashGetInfo(petsclib::$UnionPetscLib, vec::PetscVec )
	nstash_ = Ref{$PetscInt}()
	reallocs_ = Ref{$PetscInt}()
	bnstash_ = Ref{$PetscInt}()
	breallocs_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecStashGetInfo, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               vec, nstash_, reallocs_, bnstash_, breallocs_,
              )

	nstash = nstash_[]
	reallocs = reallocs_[]
	bnstash = bnstash_[]
	breallocs = breallocs_[]

	return nstash,reallocs,bnstash,breallocs
end 

"""
	VecSetLocalToGlobalMapping(petsclib::PetscLibType,x::PetscVec, mapping::ISLocalToGlobalMapping) 
Sets a local numbering to global numbering used
by the routine `VecSetValuesLocal()` to allow users to insert vector entries
using a local (per-processor) numbering.

Logically Collective

Input Parameters:
- `x`       - vector
- `mapping` - mapping created with `ISLocalToGlobalMappingCreate()` or `ISLocalToGlobalMappingCreateIS()`

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValues()`, `VecSetValuesLocal()`,
`VecGetLocalToGlobalMapping()`, `VecSetValuesBlockedLocal()`

# External Links
$(_doc_external("Vec/VecSetLocalToGlobalMapping"))
"""
function VecSetLocalToGlobalMapping(petsclib::PetscLibType, x::PetscVec, mapping::ISLocalToGlobalMapping) end

@for_petsc function VecSetLocalToGlobalMapping(petsclib::$UnionPetscLib, x::PetscVec, mapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:VecSetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CVec, ISLocalToGlobalMapping),
               x, mapping,
              )


	return nothing
end 

"""
	VecSetPreallocationCOO(petsclib::PetscLibType,x::PetscVec, ncoo::PetscCount, coo_i::Vector{PetscInt}) 
set preallocation for a vector using a coordinate format of the entries with global indices

Collective

Input Parameters:
- `x`     - vector being preallocated
- `ncoo`  - number of entries
- `coo_i` - entry indices

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecSetValuesCOO()`, `VecSetPreallocationCOOLocal()`

# External Links
$(_doc_external("Vec/VecSetPreallocationCOO"))
"""
function VecSetPreallocationCOO(petsclib::PetscLibType, x::PetscVec, ncoo::PetscCount, coo_i::Vector{PetscInt}) end

@for_petsc function VecSetPreallocationCOO(petsclib::$UnionPetscLib, x::PetscVec, ncoo::PetscCount, coo_i::Vector{$PetscInt} )

    @chk ccall(
               (:VecSetPreallocationCOO, $petsc_library),
               PetscErrorCode,
               (CVec, PetscCount, Ptr{$PetscInt}),
               x, ncoo, coo_i,
              )


	return nothing
end 

"""
	VecSetPreallocationCOOLocal(petsclib::PetscLibType,x::PetscVec, ncoo::PetscCount, coo_i::Vector{PetscInt}) 
set preallocation for vectors using a coordinate format of the entries with local indices

Collective

Input Parameters:
- `x`     - vector being preallocated
- `ncoo`  - number of entries
- `coo_i` - row indices (local numbering; may be modified)

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecSetPreallocationCOO()`, `VecSetValuesCOO()`

# External Links
$(_doc_external("Vec/VecSetPreallocationCOOLocal"))
"""
function VecSetPreallocationCOOLocal(petsclib::PetscLibType, x::PetscVec, ncoo::PetscCount, coo_i::Vector{PetscInt}) end

@for_petsc function VecSetPreallocationCOOLocal(petsclib::$UnionPetscLib, x::PetscVec, ncoo::PetscCount, coo_i::Vector{$PetscInt} )

    @chk ccall(
               (:VecSetPreallocationCOOLocal, $petsc_library),
               PetscErrorCode,
               (CVec, PetscCount, Ptr{$PetscInt}),
               x, ncoo, coo_i,
              )


	return nothing
end 

"""
	VecSetValuesCOO(petsclib::PetscLibType,x::PetscVec, coo_v::Vector{PetscScalar}, imode::InsertMode) 
set values at once in a vector preallocated using `VecSetPreallocationCOO()`

Collective

Input Parameters:
- `x`     - vector being set
- `coo_v` - the value array
- `imode` - the insert mode

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecSetPreallocationCOO()`, `VecSetPreallocationCOOLocal()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecSetValuesCOO"))
"""
function VecSetValuesCOO(petsclib::PetscLibType, x::PetscVec, coo_v::Vector{PetscScalar}, imode::InsertMode) end

@for_petsc function VecSetValuesCOO(petsclib::$UnionPetscLib, x::PetscVec, coo_v::Vector{$PetscScalar}, imode::InsertMode )

    @chk ccall(
               (:VecSetValuesCOO, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}, InsertMode),
               x, coo_v, imode,
              )


	return nothing
end 

"""
	VecSetOption(petsclib::PetscLibType,x::PetscVec, op::VecOption, flag::PetscBool) 
Sets an option for controlling a vector's behavior.

Collective

Input Parameters:
- `x`    - the vector
- `op`   - the option
- `flag` - turn the option on or off

Supported Options:
- `VEC_IGNORE_OFF_PROC_ENTRIES` - which causes `VecSetValues()` to ignore
entries destined to be stored on a separate processor. This can be used
to eliminate the global reduction in the `VecAssemblyBegin()` if you know
that you have only used `VecSetValues()` to set local elements
- `VEC_IGNORE_NEGATIVE_INDICES` - which means you can pass negative indices
in ix in calls to `VecSetValues()` or `VecGetValues()`. These rows are simply
ignored.
- `VEC_SUBSET_OFF_PROC_ENTRIES` - which causes `VecAssemblyBegin()` to assume that the off-process
entries will always be a subset (possibly equal) of the off-process entries set on the
first assembly which had a true `VEC_SUBSET_OFF_PROC_ENTRIES` and the vector has not
changed this flag afterwards. If this assembly is not such first assembly, then this
assembly can reuse the communication pattern setup in that first assembly, thus avoiding
a global reduction. Subsequent assemblies setting off-process values should use the same
InsertMode as the first assembly.

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecSetOption"))
"""
function VecSetOption(petsclib::PetscLibType, x::PetscVec, op::VecOption, flag::PetscBool) end

@for_petsc function VecSetOption(petsclib::$UnionPetscLib, x::PetscVec, op::VecOption, flag::PetscBool )

    @chk ccall(
               (:VecSetOption, $petsc_library),
               PetscErrorCode,
               (CVec, VecOption, PetscBool),
               x, op, flag,
              )


	return nothing
end 

"""
	VecStashSetInitialSize(petsclib::PetscLibType,vec::PetscVec, size::PetscInt, bsize::PetscInt) 
sets the sizes of the vec
used during the assembly process to store values that belong to
other processors.

Not Collective, different processes can have different size stashes

Input Parameters:
- `vec`   - the vector
- `size`  - the initial size of the stash.
- `bsize` - the initial size of the block-stash(if used).

Options Database Keys:
- `-vecstash_initial_size <size> or <size0,size1,...sizep-1>`           - set initial size
- `-vecstash_block_initial_size <bsize> or <bsize0,bsize1,...bsizep-1>` - set initial block size

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecSetBlockSize()`, `VecSetValues()`, `VecSetValuesBlocked()`, `VecStashView()`

# External Links
$(_doc_external("Vec/VecStashSetInitialSize"))
"""
function VecStashSetInitialSize(petsclib::PetscLibType, vec::PetscVec, size::PetscInt, bsize::PetscInt) end

@for_petsc function VecStashSetInitialSize(petsclib::$UnionPetscLib, vec::PetscVec, size::$PetscInt, bsize::$PetscInt )

    @chk ccall(
               (:VecStashSetInitialSize, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt),
               vec, size, bsize,
              )


	return nothing
end 

"""
	VecSetRandom(petsclib::PetscLibType,x::PetscVec, rctx::PetscRandom) 
Sets all components of a vector to random numbers.

Logically Collective

Input Parameters:
- `x`    - the vector
- `rctx` - the random number context, formed by `PetscRandomCreate()`, or use `NULL` and it will create one internally.

Output Parameter:
- `x` - the vector

Example of Usage:
-seealso: [](ch_vectors), `Vec`, `VecSet()`, `VecSetValues()`, `PetscRandomCreate()`, `PetscRandomDestroy()`

# External Links
$(_doc_external("Vec/VecSetRandom"))
"""
function VecSetRandom(petsclib::PetscLibType, x::PetscVec, rctx::PetscRandom) end

@for_petsc function VecSetRandom(petsclib::$UnionPetscLib, x::PetscVec, rctx::PetscRandom )

    @chk ccall(
               (:VecSetRandom, $petsc_library),
               PetscErrorCode,
               (CVec, PetscRandom),
               x, rctx,
              )


	return nothing
end 

"""
	VecSetFromOptions(petsclib::PetscLibType,vec::PetscVec) 
Configures the vector from the options database.

Collective

Input Parameter:
- `vec` - The vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `VecSetOptionsPrefix()`

# External Links
$(_doc_external("Vec/VecSetFromOptions"))
"""
function VecSetFromOptions(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecSetFromOptions(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecSetSizes(petsclib::PetscLibType,v::PetscVec, n::PetscInt, N::PetscInt) 
Sets the local and global sizes, and checks to determine compatibility of the sizes

Collective

Input Parameters:
- `v` - the vector
- `n` - the local size (or `PETSC_DECIDE` to have it set)
- `N` - the global size (or `PETSC_DETERMINE` to have it set)

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `VecCreateSeq()`, `VecCreateMPI()`, `VecGetSize()`, `PetscSplitOwnership()`, `PetscLayout`,
`VecGetOwnershipRange()`, `VecGetOwnershipRanges()`, `MatSetSizes()`

# External Links
$(_doc_external("Vec/VecSetSizes"))
"""
function VecSetSizes(petsclib::PetscLibType, v::PetscVec, n::PetscInt, N::PetscInt) end

@for_petsc function VecSetSizes(petsclib::$UnionPetscLib, v::PetscVec, n::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:VecSetSizes, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt),
               v, n, N,
              )


	return nothing
end 

"""
	VecSetBlockSize(petsclib::PetscLibType,v::PetscVec, bs::PetscInt) 
Sets the block size for future calls to `VecSetValuesBlocked()`
and `VecSetValuesBlockedLocal()`.

Logically Collective

Input Parameters:
- `v`  - the vector
- `bs` - the blocksize

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSetValuesBlocked()`, `VecSetLocalToGlobalMapping()`, `VecGetBlockSize()`

# External Links
$(_doc_external("Vec/VecSetBlockSize"))
"""
function VecSetBlockSize(petsclib::PetscLibType, v::PetscVec, bs::PetscInt) end

@for_petsc function VecSetBlockSize(petsclib::$UnionPetscLib, v::PetscVec, bs::$PetscInt )

    @chk ccall(
               (:VecSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt),
               v, bs,
              )


	return nothing
end 

"""
	VecSetOptionsPrefix(petsclib::PetscLibType,v::PetscVec, prefix::Vector{Cchar}) 
Sets the prefix used for searching for all
`Vec` options in the database.

Logically Collective

Input Parameters:
- `v`      - the `Vec` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSetFromOptions()`

# External Links
$(_doc_external("Vec/VecSetOptionsPrefix"))
"""
function VecSetOptionsPrefix(petsclib::PetscLibType, v::PetscVec, prefix::Vector{Cchar}) end

@for_petsc function VecSetOptionsPrefix(petsclib::$UnionPetscLib, v::PetscVec, prefix::Vector{Cchar} )

    @chk ccall(
               (:VecSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Cchar}),
               v, prefix,
              )


	return nothing
end 

"""
	VecSetUp(petsclib::PetscLibType,v::PetscVec) 
Sets up the internal vector data structures for the later use.

Collective

Input Parameter:
- `v` - the `Vec` context

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `VecDestroy()`

# External Links
$(_doc_external("Vec/VecSetUp"))
"""
function VecSetUp(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecSetUp(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecSetUp, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecSwap(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
Swaps the values between two vectors, `x` and `y`.

Logically Collective

Input Parameters:
- `x` - the first vector
- `y` - the second vector

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSet()`

# External Links
$(_doc_external("Vec/VecSwap"))
"""
function VecSwap(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecSwap(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecSwap, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               x, y,
              )


	return nothing
end 

"""
	VecStashViewFromOptions(petsclib::PetscLibType,obj::PetscVec, bobj::PetscObject, optionname::Vector{Cchar}) 
Processes command line options to determine if/how a `VecStash` object is to be viewed.

Collective

Input Parameters:
- `obj`        - the `Vec` containing a stash
- `bobj`       - optional other object that provides the prefix
- `optionname` - option to activate viewing

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecStashSetInitialSize()`

# External Links
$(_doc_external("Vec/VecStashViewFromOptions"))
"""
function VecStashViewFromOptions(petsclib::PetscLibType, obj::PetscVec, bobj::PetscObject, optionname::Vector{Cchar}) end

@for_petsc function VecStashViewFromOptions(petsclib::$UnionPetscLib, obj::PetscVec, bobj::PetscObject, optionname::Vector{Cchar} )

    @chk ccall(
               (:VecStashViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CVec, PetscObject, Ptr{Cchar}),
               obj, bobj, optionname,
              )


	return nothing
end 

"""
	VecStashView(petsclib::PetscLibType,v::PetscVec, viewer::PetscViewer) 
Prints the entries in the vector stash and block stash.

Collective

Input Parameters:
- `v`      - the vector
- `viewer` - the viewer

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSetBlockSize()`, `VecSetValues()`, `VecSetValuesBlocked()`

# External Links
$(_doc_external("Vec/VecStashView"))
"""
function VecStashView(petsclib::PetscLibType, v::PetscVec, viewer::PetscViewer) end

@for_petsc function VecStashView(petsclib::$UnionPetscLib, v::PetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecStashView, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               v, viewer,
              )


	return nothing
end 

"""
	VecSetLayout(petsclib::PetscLibType,x::PetscVec, map::PetscLayout) 
set `PetscLayout` describing vector layout

Not Collective

Input Parameters:
- `x`   - the vector
- `map` - the layout

Level: developer

-seealso: [](ch_vectors), `Vec`, `PetscLayout`, `VecGetLayout()`, `VecGetSize()`, `VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecSetLayout"))
"""
function VecSetLayout(petsclib::PetscLibType, x::PetscVec, map::PetscLayout) end

@for_petsc function VecSetLayout(petsclib::$UnionPetscLib, x::PetscVec, map::PetscLayout )

    @chk ccall(
               (:VecSetLayout, $petsc_library),
               PetscErrorCode,
               (CVec, PetscLayout),
               x, map,
              )


	return nothing
end 

"""
	VecSetBindingPropagates(petsclib::PetscLibType,v::PetscVec, flg::PetscBool) 
Sets whether the state of being bound to the CPU for a GPU vector type propagates to child and some other associated objects

Input Parameters:
- `v`   - the vector
- `flg` - flag indicating whether the boundtocpu flag should be propagated

Level: developer

-seealso: [](ch_vectors), `Vec`, `MatSetBindingPropagates()`, `VecGetBindingPropagates()`

# External Links
$(_doc_external("Vec/VecSetBindingPropagates"))
"""
function VecSetBindingPropagates(petsclib::PetscLibType, v::PetscVec, flg::PetscBool) end

@for_petsc function VecSetBindingPropagates(petsclib::$UnionPetscLib, v::PetscVec, flg::PetscBool )

    @chk ccall(
               (:VecSetBindingPropagates, $petsc_library),
               PetscErrorCode,
               (CVec, PetscBool),
               v, flg,
              )


	return nothing
end 

"""
	VecSetPinnedMemoryMin(petsclib::PetscLibType,v::PetscVec, mbytes::Csize_t) 
Set the minimum data size for which pinned memory will be used for host (CPU) allocations.

Logically Collective

Input Parameters:
- `v`      - the vector
- `mbytes` - minimum data size in bytes

Options Database Key:
- `-vec_pinned_memory_min <size>` - minimum size (in bytes) for an allocation to use pinned memory on host.

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetPinnedMemoryMin()`

# External Links
$(_doc_external("Vec/VecSetPinnedMemoryMin"))
"""
function VecSetPinnedMemoryMin(petsclib::PetscLibType, v::PetscVec, mbytes::Csize_t) end

@for_petsc function VecSetPinnedMemoryMin(petsclib::$UnionPetscLib, v::PetscVec, mbytes::Csize_t )

    @chk ccall(
               (:VecSetPinnedMemoryMin, $petsc_library),
               PetscErrorCode,
               (CVec, Csize_t),
               v, mbytes,
              )


	return nothing
end 

"""
	VecSetType(petsclib::PetscLibType,vec::PetscVec, newType::VecType) 
Builds a vector, for a particular vector implementation.

Collective

Input Parameters:
- `vec`     - The vector object
- `newType` - The name of the vector type

Options Database Key:
- `-vec_type <type>` - Sets the vector type; use -help for a list
of available types

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecGetType()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecSetType"))
"""
function VecSetType(petsclib::PetscLibType, vec::PetscVec, newType::VecType) end

@for_petsc function VecSetType(petsclib::$UnionPetscLib, vec::PetscVec, newType::VecType )

    @chk ccall(
               (:VecSetType, $petsc_library),
               PetscErrorCode,
               (CVec, VecType),
               vec, newType,
              )


	return nothing
end 

"""
	VecStrideSet(petsclib::PetscLibType,v::PetscVec, start::PetscInt, s::PetscScalar) 
Sets a subvector of a vector defined
by a starting point and a stride with a given value

Logically Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)
- `s`     - value to set for each entry in that subvector

Level: advanced

-seealso: `Vec`, `VecNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideScale()`

# External Links
$(_doc_external("Vec/VecStrideSet"))
"""
function VecStrideSet(petsclib::PetscLibType, v::PetscVec, start::PetscInt, s::PetscScalar) end

@for_petsc function VecStrideSet(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt, s::$PetscScalar )

    @chk ccall(
               (:VecStrideSet, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscScalar),
               v, start, s,
              )


	return nothing
end 

"""
	VecStrideScale(petsclib::PetscLibType,v::PetscVec, start::PetscInt, scale::PetscScalar) 
Scales a subvector of a vector defined
by a starting point and a stride.

Logically Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)
- `scale` - value to multiply each subvector entry by

Level: advanced

-seealso: `Vec`, `VecNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideScale"))
"""
function VecStrideScale(petsclib::PetscLibType, v::PetscVec, start::PetscInt, scale::PetscScalar) end

@for_petsc function VecStrideScale(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt, scale::$PetscScalar )

    @chk ccall(
               (:VecStrideScale, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscScalar),
               v, start, scale,
              )


	return nothing
end 

"""
	nrm::PetscReal = VecStrideNorm(petsclib::PetscLibType,v::PetscVec, start::PetscInt, ntype::NormType) 
Computes the norm of subvector of a vector defined
by a starting point and a stride.

Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)
- `ntype` - type of norm, one of `NORM_1`, `NORM_2`, `NORM_INFINITY`

Output Parameter:
- `nrm` - the norm

Level: advanced

-seealso: `Vec`, `VecNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideNorm"))
"""
function VecStrideNorm(petsclib::PetscLibType, v::PetscVec, start::PetscInt, ntype::NormType) end

@for_petsc function VecStrideNorm(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt, ntype::NormType )
	nrm_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStrideNorm, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, NormType, Ptr{$PetscReal}),
               v, start, ntype, nrm_,
              )

	nrm = nrm_[]

	return nrm
end 

"""
	idex::PetscInt,nrm::PetscReal = VecStrideMax(petsclib::PetscLibType,v::PetscVec, start::PetscInt) 
Computes the maximum of subvector of a vector defined
by a starting point and a stride and optionally its location.

Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)

Output Parameters:
- `idex` - the location where the maximum occurred  (pass `NULL` if not required)
- `nrm`  - the maximum value in the subvector

Level: advanced

-seealso: `Vec`, `VecMax()`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`

# External Links
$(_doc_external("Vec/VecStrideMax"))
"""
function VecStrideMax(petsclib::PetscLibType, v::PetscVec, start::PetscInt) end

@for_petsc function VecStrideMax(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt )
	idex_ = Ref{$PetscInt}()
	nrm_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStrideMax, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}),
               v, start, idex_, nrm_,
              )

	idex = idex_[]
	nrm = nrm_[]

	return idex,nrm
end 

"""
	idex::PetscInt,nrm::PetscReal = VecStrideMin(petsclib::PetscLibType,v::PetscVec, start::PetscInt) 
Computes the minimum of subvector of a vector defined
by a starting point and a stride and optionally its location.

Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)

Output Parameters:
- `idex` - the location where the minimum occurred. (pass `NULL` if not required)
- `nrm`  - the minimum value in the subvector

Level: advanced

-seealso: `Vec`, `VecMin()`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideMin"))
"""
function VecStrideMin(petsclib::PetscLibType, v::PetscVec, start::PetscInt) end

@for_petsc function VecStrideMin(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt )
	idex_ = Ref{$PetscInt}()
	nrm_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStrideMin, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}),
               v, start, idex_, nrm_,
              )

	idex = idex_[]
	nrm = nrm_[]

	return idex,nrm
end 

"""
	sum::PetscScalar = VecStrideSum(petsclib::PetscLibType,v::PetscVec, start::PetscInt) 
Computes the sum of subvector of a vector defined
by a starting point and a stride.

Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)

Output Parameter:
- `sum` - the sum

Level: advanced

-seealso: `Vec`, `VecSum()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideSum"))
"""
function VecStrideSum(petsclib::PetscLibType, v::PetscVec, start::PetscInt) end

@for_petsc function VecStrideSum(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt )
	sum_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecStrideSum, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscScalar}),
               v, start, sum_,
              )

	sum = sum_[]

	return sum
end 

"""
	VecStrideScaleAll(petsclib::PetscLibType,v::PetscVec, scales::PetscScalar) 
Scales the subvectors of a vector defined
by a starting point and a stride.

Logically Collective

Input Parameters:
- `v`      - the vector
- `scales` - values to multiply each subvector entry by

Level: advanced

-seealso: `Vec`, `VecNorm()`, `VecStrideScale()`, `VecScale()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideScaleAll"))
"""
function VecStrideScaleAll(petsclib::PetscLibType, v::PetscVec, scales::PetscScalar) end

@for_petsc function VecStrideScaleAll(petsclib::$UnionPetscLib, v::PetscVec, scales::$PetscScalar )

    @chk ccall(
               (:VecStrideScaleAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               v, scales,
              )


	return nothing
end 

"""
	nrm::Vector{PetscReal} = VecStrideNormAll(petsclib::PetscLibType,v::PetscVec, ntype::NormType) 
Computes the norms of subvectors of a vector defined
by a starting point and a stride.

Collective

Input Parameters:
- `v`     - the vector
- `ntype` - type of norm, one of `NORM_1`, `NORM_2`, `NORM_INFINITY`

Output Parameter:
- `nrm` - the norms

Level: advanced

-seealso: `Vec`, `VecNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideNormAll"))
"""
function VecStrideNormAll(petsclib::PetscLibType, v::PetscVec, ntype::NormType) end

@for_petsc function VecStrideNormAll(petsclib::$UnionPetscLib, v::PetscVec, ntype::NormType )
	nrm = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecStrideNormAll, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               v, ntype, nrm,
              )


	return nrm
end 

"""
	idex::Vector{PetscInt},nrm::Vector{PetscReal} = VecStrideMaxAll(petsclib::PetscLibType,v::PetscVec) 
Computes the maximums of subvectors of a vector defined
by a starting point and a stride and optionally its location.

Collective

Input Parameter:
- `v` - the vector

Output Parameters:
- `idex` - the location where the maximum occurred (not supported, pass `NULL`,
if you need this, send mail to petsc-maint@mcs.anl.gov to request it)
- `nrm`  - the maximum values of each subvector

Level: advanced

-seealso: `Vec`, `VecMax()`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`

# External Links
$(_doc_external("Vec/VecStrideMaxAll"))
"""
function VecStrideMaxAll(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecStrideMaxAll(petsclib::$UnionPetscLib, v::PetscVec )
	idex = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!
	nrm = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecStrideMaxAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscReal}),
               v, idex, nrm,
              )


	return idex,nrm
end 

"""
	idex::Vector{PetscInt},nrm::Vector{PetscReal} = VecStrideMinAll(petsclib::PetscLibType,v::PetscVec) 
Computes the minimum of subvector of a vector defined
by a starting point and a stride and optionally its location.

Collective

Input Parameter:
- `v` - the vector

Output Parameters:
- `idex` - the location where the minimum occurred (not supported, pass `NULL`,
if you need this, send mail to petsc-maint@mcs.anl.gov to request it)
- `nrm`  - the minimums of each subvector

Level: advanced

-seealso: `Vec`, `VecMin()`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideMinAll"))
"""
function VecStrideMinAll(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecStrideMinAll(petsclib::$UnionPetscLib, v::PetscVec )
	idex = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!
	nrm = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecStrideMinAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscReal}),
               v, idex, nrm,
              )


	return idex,nrm
end 

"""
	sums::Vector{PetscScalar} = VecStrideSumAll(petsclib::PetscLibType,v::PetscVec) 
Computes the sums of subvectors of a vector defined by a stride.

Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `sums` - the sums

Level: advanced

-seealso: `Vec`, `VecSum()`, `VecStrideGather()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`

# External Links
$(_doc_external("Vec/VecStrideSumAll"))
"""
function VecStrideSumAll(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecStrideSumAll(petsclib::$UnionPetscLib, v::PetscVec )
	sums = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecStrideSumAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               v, sums,
              )


	return sums
end 

"""
	VecStrideGatherAll(petsclib::PetscLibType,v::PetscVec, s::Vector{PetscVec}, addv::InsertMode) 
Gathers all the single components from a multi
separate vectors.

Collective

Input Parameters:
- `v`    - the vector
- `addv` - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `s` - the location where the subvectors are stored

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGather()`,
`VecStrideScatterAll()`

# External Links
$(_doc_external("Vec/VecStrideGatherAll"))
"""
function VecStrideGatherAll(petsclib::PetscLibType, v::PetscVec, s::Vector{PetscVec}, addv::InsertMode) end

@for_petsc function VecStrideGatherAll(petsclib::$UnionPetscLib, v::PetscVec, s::Vector{PetscVec}, addv::InsertMode )

    @chk ccall(
               (:VecStrideGatherAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}, InsertMode),
               v, s, addv,
              )


	return nothing
end 

"""
	VecStrideScatterAll(petsclib::PetscLibType,s::Vector{PetscVec}, v::PetscVec, addv::InsertMode) 
Scatters all the single components from separate vectors into
a multi-component vector.

Collective

Input Parameters:
- `s`    - the location where the subvectors are stored
- `addv` - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `v` - the multicomponent vector

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGather()`,


# External Links
$(_doc_external("Vec/VecStrideScatterAll"))
"""
function VecStrideScatterAll(petsclib::PetscLibType, s::Vector{PetscVec}, v::PetscVec, addv::InsertMode) end

@for_petsc function VecStrideScatterAll(petsclib::$UnionPetscLib, s::Vector{PetscVec}, v::PetscVec, addv::InsertMode )

    @chk ccall(
               (:VecStrideScatterAll, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec}, CVec, InsertMode),
               s, v, addv,
              )


	return nothing
end 

"""
	VecStrideGather(petsclib::PetscLibType,v::PetscVec, start::PetscInt, s::PetscVec, addv::InsertMode) 
Gathers a single component from a multi
another vector.

Collective

Input Parameters:
- `v`     - the vector
- `start` - starting point of the subvector (defined by a stride)
- `addv`  - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `s` - the location where the subvector is stored

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideScatter()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGatherAll()`,
`VecStrideScatterAll()`

# External Links
$(_doc_external("Vec/VecStrideGather"))
"""
function VecStrideGather(petsclib::PetscLibType, v::PetscVec, start::PetscInt, s::PetscVec, addv::InsertMode) end

@for_petsc function VecStrideGather(petsclib::$UnionPetscLib, v::PetscVec, start::$PetscInt, s::PetscVec, addv::InsertMode )

    @chk ccall(
               (:VecStrideGather, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, CVec, InsertMode),
               v, start, s, addv,
              )


	return nothing
end 

"""
	VecStrideScatter(petsclib::PetscLibType,s::PetscVec, start::PetscInt, v::PetscVec, addv::InsertMode) 
Scatters a single component from a vector into a multi

Collective

Input Parameters:
- `s`     - the single-component vector
- `start` - starting point of the subvector (defined by a stride)
- `addv`  - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `v` - the location where the subvector is scattered (the multi-component vector)

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGatherAll()`,
`VecStrideScatterAll()`, `VecStrideSubSetScatter()`, `VecStrideSubSetGather()`

# External Links
$(_doc_external("Vec/VecStrideScatter"))
"""
function VecStrideScatter(petsclib::PetscLibType, s::PetscVec, start::PetscInt, v::PetscVec, addv::InsertMode) end

@for_petsc function VecStrideScatter(petsclib::$UnionPetscLib, s::PetscVec, start::$PetscInt, v::PetscVec, addv::InsertMode )

    @chk ccall(
               (:VecStrideScatter, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, CVec, InsertMode),
               s, start, v, addv,
              )


	return nothing
end 

"""
	VecStrideSubSetGather(petsclib::PetscLibType,v::PetscVec, nidx::PetscInt, idxv::Vector{PetscInt}, idxs::Vector{PetscInt}, s::PetscVec, addv::InsertMode) 
Gathers a subset of components from a multi
another vector.

Collective

Input Parameters:
- `v`    - the vector
- `nidx` - the number of indices
- `idxv` - the indices of the components 0 <= idxv[0] ...idxv[nidx-1] < bs(v), they need not be sorted
- `idxs` - the indices of the components 0 <= idxs[0] ...idxs[nidx-1] < bs(s), they need not be sorted, may be null if nidx == bs(s) or is `PETSC_DETERMINE`
- `addv` - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `s` - the location where the subvector is stored

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideScatter()`, `VecStrideGather()`, `VecStrideSubSetScatter()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGatherAll()`,
`VecStrideScatterAll()`

# External Links
$(_doc_external("Vec/VecStrideSubSetGather"))
"""
function VecStrideSubSetGather(petsclib::PetscLibType, v::PetscVec, nidx::PetscInt, idxv::Vector{PetscInt}, idxs::Vector{PetscInt}, s::PetscVec, addv::InsertMode) end

@for_petsc function VecStrideSubSetGather(petsclib::$UnionPetscLib, v::PetscVec, nidx::$PetscInt, idxv::Vector{$PetscInt}, idxs::Vector{$PetscInt}, s::PetscVec, addv::InsertMode )

    @chk ccall(
               (:VecStrideSubSetGather, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, CVec, InsertMode),
               v, nidx, idxv, idxs, s, addv,
              )


	return nothing
end 

"""
	VecStrideSubSetScatter(petsclib::PetscLibType,s::PetscVec, nidx::PetscInt, idxs::Vector{PetscInt}, idxv::Vector{PetscInt}, v::PetscVec, addv::InsertMode) 
Scatters components from a vector into a subset of components of a multi

Collective

Input Parameters:
- `s`    - the smaller-component vector
- `nidx` - the number of indices in idx
- `idxs` - the indices of the components in the smaller-component vector, 0 <= idxs[0] ...idxs[nidx-1] < bs(s) they need not be sorted, may be null if nidx == bs(s) or is `PETSC_DETERMINE`
- `idxv` - the indices of the components in the larger-component vector, 0 <= idx[0] ...idx[nidx-1] < bs(v) they need not be sorted
- `addv` - one of `ADD_VALUES`, `INSERT_VALUES`, `MAX_VALUES`

Output Parameter:
- `v` - the location where the subvector is into scattered (the multi-component vector)

Level: advanced

-seealso: `Vec`, `VecStrideNorm()`, `VecStrideGather()`, `VecStrideSubSetGather()`, `VecStrideMin()`, `VecStrideMax()`, `VecStrideGatherAll()`,
`VecStrideScatterAll()`

# External Links
$(_doc_external("Vec/VecStrideSubSetScatter"))
"""
function VecStrideSubSetScatter(petsclib::PetscLibType, s::PetscVec, nidx::PetscInt, idxs::Vector{PetscInt}, idxv::Vector{PetscInt}, v::PetscVec, addv::InsertMode) end

@for_petsc function VecStrideSubSetScatter(petsclib::$UnionPetscLib, s::PetscVec, nidx::$PetscInt, idxs::Vector{$PetscInt}, idxv::Vector{$PetscInt}, v::PetscVec, addv::InsertMode )

    @chk ccall(
               (:VecStrideSubSetScatter, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, CVec, InsertMode),
               s, nidx, idxs, idxv, v, addv,
              )


	return nothing
end 

"""
	VecSqrtAbs(petsclib::PetscLibType,v::PetscVec) 
Replaces each component of a vector by the square root of its magnitude.

Not Collective

Input Parameter:
- `v` - The vector

Level: beginner

-seealso: `Vec`, `VecLog()`, `VecExp()`, `VecReciprocal()`, `VecAbs()`


# External Links
$(_doc_external("Vec/VecSqrtAbs"))
"""
function VecSqrtAbs(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecSqrtAbs(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecSqrtAbs, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	sum::PetscScalar = VecSum(petsclib::PetscLibType,v::PetscVec) 
Computes the sum of all the components of a vector.

Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `sum` - the result

Level: beginner

-seealso: `Vec`, `VecMean()`, `VecNorm()`

# External Links
$(_doc_external("Vec/VecSum"))
"""
function VecSum(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecSum(petsclib::$UnionPetscLib, v::PetscVec )
	sum_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecSum, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               v, sum_,
              )

	sum = sum_[]

	return sum
end 

"""
	VecShift(petsclib::PetscLibType,v::PetscVec, shift::PetscScalar) 
Shifts all of the components of a vector by computing
`x[i] = x[i] + shift`.

Logically Collective

Input Parameters:
- `v`     - the vector
- `shift` - the shift

Level: intermediate

-seealso: `Vec`, `VecISShift()`

# External Links
$(_doc_external("Vec/VecShift"))
"""
function VecShift(petsclib::PetscLibType, v::PetscVec, shift::PetscScalar) end

@for_petsc function VecShift(petsclib::$UnionPetscLib, v::PetscVec, shift::$PetscScalar )

    @chk ccall(
               (:VecShift, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar),
               v, shift,
              )


	return nothing
end 

"""
	VecScatterSetUp(petsclib::PetscLibType,sf::VecScatter) 
Sets up the `VecScatter` to be able to actually scatter information between vectors

Collective

Input Parameter:
- `sf` - the scatter context

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterCopy()`

# External Links
$(_doc_external("Vec/VecScatterSetUp"))
"""
function VecScatterSetUp(petsclib::PetscLibType, sf::VecScatter) end

@for_petsc function VecScatterSetUp(petsclib::$UnionPetscLib, sf::VecScatter )

    @chk ccall(
               (:VecScatterSetUp, $petsc_library),
               PetscErrorCode,
               (VecScatter,),
               sf,
              )


	return nothing
end 

"""
	VecScatterSetType(petsclib::PetscLibType,sf::VecScatter, type::VecScatterType) 
Builds a vector scatter, for a particular vector scatter implementation.

Collective

Input Parameters:
- `sf`   - The `VecScatter` object
- `type` - The name of the vector scatter type

Options Database Key:
- `-sf_type <type>` - Sets the `VecScatterType`

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterType`, `VecScatterGetType()`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterSetType"))
"""
function VecScatterSetType(petsclib::PetscLibType, sf::VecScatter, type::VecScatterType) end

@for_petsc function VecScatterSetType(petsclib::$UnionPetscLib, sf::VecScatter, type::VecScatterType )

    @chk ccall(
               (:VecScatterSetType, $petsc_library),
               PetscErrorCode,
               (VecScatter, VecScatterType),
               sf, type,
              )


	return nothing
end 

"""
	type::VecScatterType = VecScatterGetType(petsclib::PetscLibType,sf::VecScatter) 
Gets the vector scatter type name (as a string) from the `VecScatter`.

Not Collective

Input Parameter:
- `sf` - The vector scatter

Output Parameter:
- `type` - The vector scatter type name

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterType`, `VecScatterSetType()`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterGetType"))
"""
function VecScatterGetType(petsclib::PetscLibType, sf::VecScatter) end

@for_petsc function VecScatterGetType(petsclib::$UnionPetscLib, sf::VecScatter )
	type_ = Ref{VecScatterType}()

    @chk ccall(
               (:VecScatterGetType, $petsc_library),
               PetscErrorCode,
               (VecScatter, Ptr{VecScatterType}),
               sf, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	VecScatterRegister(petsclib::PetscLibType,sname::Vector{Cchar}, fnc::external) 
Adds a new vector scatter component implementation

Not Collective

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

Level: advanced

-seealso: [](sec_scatter), `VecScatter`, `VecScatterType`, `VecRegister()`

# External Links
$(_doc_external("Vec/VecScatterRegister"))
"""
function VecScatterRegister(petsclib::PetscLibType, sname::Vector{Cchar}, fnc::external) end

@for_petsc function VecScatterRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, fnc::external )

    @chk ccall(
               (:VecScatterRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	flg::PetscBool = VecScatterGetMerged(petsclib::PetscLibType,sf::VecScatter) 
Returns true if the scatter is completed in the `VecScatterBegin()`
and the `VecScatterEnd()` does nothing

Not Collective

Input Parameter:
- `sf` - scatter context created with `VecScatterCreate()`

Output Parameter:
- `flg` - `PETSC_TRUE` if the `VecScatterBegin()`/`VecScatterEnd()` are all done during the `VecScatterBegin()`

Level: developer

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterEnd()`, `VecScatterBegin()`

# External Links
$(_doc_external("Vec/VecScatterGetMerged"))
"""
function VecScatterGetMerged(petsclib::PetscLibType, sf::VecScatter) end

@for_petsc function VecScatterGetMerged(petsclib::$UnionPetscLib, sf::VecScatter )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:VecScatterGetMerged, $petsc_library),
               PetscErrorCode,
               (VecScatter, Ptr{PetscBool}),
               sf, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	VecScatterDestroy(petsclib::PetscLibType,sf::VecScatter) 
Destroys a scatter context created by `VecScatterCreate()`

Collective

Input Parameter:
- `sf` - the scatter context

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterCopy()`

# External Links
$(_doc_external("Vec/VecScatterDestroy"))
"""
function VecScatterDestroy(petsclib::PetscLibType, sf::VecScatter) end

@for_petsc function VecScatterDestroy(petsclib::$UnionPetscLib, sf::VecScatter )

    sf_ = Ref{VecScatter}(sf)

    @chk ccall(
               (:VecScatterDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{VecScatter},),
               sf_,
              )


	return nothing
end 

"""
	VecScatterCopy(petsclib::PetscLibType,sf::VecScatter, newsf::VecScatter) 
Makes a copy of a scatter context.

Collective

Input Parameter:
- `sf` - the scatter context

Output Parameter:
- `newsf` - the context copy

Level: advanced

-seealso: [](sec_scatter), `VecScatter`, `VecScatterType`, `VecScatterCreate()`, `VecScatterDestroy()`

# External Links
$(_doc_external("Vec/VecScatterCopy"))
"""
function VecScatterCopy(petsclib::PetscLibType, sf::VecScatter, newsf::VecScatter) end

@for_petsc function VecScatterCopy(petsclib::$UnionPetscLib, sf::VecScatter, newsf::VecScatter )

    @chk ccall(
               (:VecScatterCopy, $petsc_library),
               PetscErrorCode,
               (VecScatter, Ptr{VecScatter}),
               sf, newsf,
              )


	return nothing
end 

"""
	VecScatterViewFromOptions(petsclib::PetscLibType,sf::VecScatter, obj::PetscObject, name::Vector{Cchar}) 
View a `VecScatter` object based on values in the options database

Collective

Input Parameters:
- `sf`   - the scatter context
- `obj`  - Optional object
- `name` - command line option

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterView()`, `PetscObjectViewFromOptions()`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterViewFromOptions"))
"""
function VecScatterViewFromOptions(petsclib::PetscLibType, sf::VecScatter, obj::PetscObject, name::Vector{Cchar}) end

@for_petsc function VecScatterViewFromOptions(petsclib::$UnionPetscLib, sf::VecScatter, obj::PetscObject, name::Vector{Cchar} )

    @chk ccall(
               (:VecScatterViewFromOptions, $petsc_library),
               PetscErrorCode,
               (VecScatter, PetscObject, Ptr{Cchar}),
               sf, obj, name,
              )


	return nothing
end 

"""
	VecScatterView(petsclib::PetscLibType,sf::VecScatter, viewer::PetscViewer) 
Views a vector scatter context.

Collective

Input Parameters:
- `sf`     - the scatter context
- `viewer` - the viewer for displaying the context

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `PetscViewer`, `VecScatterViewFromOptions()`, `PetscObjectViewFromOptions()`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterView"))
"""
function VecScatterView(petsclib::PetscLibType, sf::VecScatter, viewer::PetscViewer) end

@for_petsc function VecScatterView(petsclib::$UnionPetscLib, sf::VecScatter, viewer::PetscViewer )

    @chk ccall(
               (:VecScatterView, $petsc_library),
               PetscErrorCode,
               (VecScatter, PetscViewer),
               sf, viewer,
              )


	return nothing
end 

"""
	VecScatterRemap(petsclib::PetscLibType,sf::VecScatter, tomap::Vector{PetscInt}, frommap::Vector{PetscInt}) 
Remaps the "from" and "to" indices in a
vector scatter context.

Collective

Input Parameters:
- `sf`      - vector scatter context
- `tomap`   - remapping plan for "to" indices (may be `NULL`).
- `frommap` - remapping plan for "from" indices (may be `NULL`)

Level: developer

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterRemap"))
"""
function VecScatterRemap(petsclib::PetscLibType, sf::VecScatter, tomap::Vector{PetscInt}, frommap::Vector{PetscInt}) end

@for_petsc function VecScatterRemap(petsclib::$UnionPetscLib, sf::VecScatter, tomap::Vector{$PetscInt}, frommap::Vector{$PetscInt} )

    @chk ccall(
               (:VecScatterRemap, $petsc_library),
               PetscErrorCode,
               (VecScatter, Ptr{$PetscInt}, Ptr{$PetscInt}),
               sf, tomap, frommap,
              )


	return nothing
end 

"""
	VecScatterSetFromOptions(petsclib::PetscLibType,sf::VecScatter) 
Configures the vector scatter from values in the options database.

Collective

Input Parameter:
- `sf` - The vector scatter

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterDestroy()`, `VecScatterSetUp()`

# External Links
$(_doc_external("Vec/VecScatterSetFromOptions"))
"""
function VecScatterSetFromOptions(petsclib::PetscLibType, sf::VecScatter) end

@for_petsc function VecScatterSetFromOptions(petsclib::$UnionPetscLib, sf::VecScatter )

    @chk ccall(
               (:VecScatterSetFromOptions, $petsc_library),
               PetscErrorCode,
               (VecScatter,),
               sf,
              )


	return nothing
end 

"""
	newsf::VecScatter = VecScatterCreate(petsclib::PetscLibType,x::PetscVec, ix::IS, y::PetscVec, iy::IS) 
Creates a vector scatter `VecScatter` context that is used to communicate entries between two vectors `Vec`

Collective

Input Parameters:
- `x`  - a vector that defines the shape (parallel data layout of the vector) of vectors from which we scatter
- `y`  - a vector that defines the shape (parallel data layout of the vector) of vectors to which we scatter
- `ix` - the indices of `x` to scatter (if `NULL` scatters all values)
- `iy` - the indices of `y` to hold results (if `NULL` fills entire vector `yin` in order)

Output Parameter:
- `newsf` - location to store the new scatter context

Options Database Keys:
- `-vecscatter_view`              - Prints detail of communications
- `-vecscatter_view ::ascii_info` - Print less details about communication
- `-vecscatter_merge`             - `VecScatterBegin()` handles all of the communication, `VecScatterEnd()` is a nop
eliminates the chance for overlap of computation and communication

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterDestroy()`, `VecScatterCreateToAll()`, `VecScatterCreateToZero()`, `PetscSFCreate()`,
`VecScatterType`, `InsertMode`, `ScatterMode`, `VecScatterBegin()`, `VecScatterEnd()`

# External Links
$(_doc_external("Vec/VecScatterCreate"))
"""
function VecScatterCreate(petsclib::PetscLibType, x::PetscVec, ix::IS, y::PetscVec, iy::IS) end

@for_petsc function VecScatterCreate(petsclib::$UnionPetscLib, x::PetscVec, ix::IS, y::PetscVec, iy::IS )
	newsf_ = Ref{VecScatter}()

    @chk ccall(
               (:VecScatterCreate, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, CVec, CIS, Ptr{VecScatter}),
               x, ix, y, iy, newsf_,
              )

	newsf = newsf_[]

	return newsf
end 

"""
	ctx::VecScatter,vout::PetscVec = VecScatterCreateToAll(petsclib::PetscLibType,vin::PetscVec) 
Creates a vector and a scatter context that copies all
vector values to each processor

Collective

Input Parameter:
- `vin` - an `MPIVEC`

Output Parameters:
- `ctx`  - scatter context
- `vout` - output `SEQVEC` that is large enough to scatter into

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterCreateToZero()`, `VecScatterBegin()`, `VecScatterEnd()`

# External Links
$(_doc_external("Vec/VecScatterCreateToAll"))
"""
function VecScatterCreateToAll(petsclib::PetscLibType, vin::PetscVec) end

@for_petsc function VecScatterCreateToAll(petsclib::$UnionPetscLib, vin::PetscVec )
	ctx_ = Ref{VecScatter}()
	vout_ = Ref{CVec}()

    @chk ccall(
               (:VecScatterCreateToAll, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{VecScatter}, Ptr{CVec}),
               vin, ctx_, vout_,
              )

	ctx = ctx_[]
	vout = PetscVec(vout_[], petsclib)

	return ctx,vout
end 

"""
	ctx::VecScatter,vout::PetscVec = VecScatterCreateToZero(petsclib::PetscLibType,vin::PetscVec) 
Creates an output vector and a scatter context used to
copy all vector values into the output vector on the zeroth processor

Collective

Input Parameter:
- `vin` - `Vec` of type `MPIVEC`

Output Parameters:
- `ctx`  - scatter context
- `vout` - output `SEQVEC` that is large enough to scatter into on processor 0 and
of length zero on all other processors

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterCreateToAll()`, `VecScatterBegin()`, `VecScatterEnd()`

# External Links
$(_doc_external("Vec/VecScatterCreateToZero"))
"""
function VecScatterCreateToZero(petsclib::PetscLibType, vin::PetscVec) end

@for_petsc function VecScatterCreateToZero(petsclib::$UnionPetscLib, vin::PetscVec )
	ctx_ = Ref{VecScatter}()
	vout_ = Ref{CVec}()

    @chk ccall(
               (:VecScatterCreateToZero, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{VecScatter}, Ptr{CVec}),
               vin, ctx_, vout_,
              )

	ctx = ctx_[]
	vout = PetscVec(vout_[], petsclib)

	return ctx,vout
end 

"""
	VecScatterBegin(petsclib::PetscLibType,sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode) 
Begins a generalized scatter from one vector to
another. Complete the scattering phase with `VecScatterEnd()`.

Neighbor-wise Collective

Input Parameters:
- `sf`   - scatter context generated by `VecScatterCreate()`
- `x`    - the vector from which we scatter
- `y`    - the vector to which we scatter
- `addv` - either `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`, with `INSERT_VALUES` mode any location
not scattered to retains its old value; i.e. the vector is NOT first zeroed.
- `mode` - the scattering mode, usually `SCATTER_FORWARD`.  The available modes are: `SCATTER_FORWARD` or `SCATTER_REVERSE`

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterCreate()`, `VecScatterEnd()`, `InsertMode`, `ScatterMode`

# External Links
$(_doc_external("Vec/VecScatterBegin"))
"""
function VecScatterBegin(petsclib::PetscLibType, sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode) end

@for_petsc function VecScatterBegin(petsclib::$UnionPetscLib, sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode )

    @chk ccall(
               (:VecScatterBegin, $petsc_library),
               PetscErrorCode,
               (VecScatter, CVec, CVec, InsertMode, ScatterMode),
               sf, x, y, addv, mode,
              )


	return nothing
end 

"""
	VecScatterEnd(petsclib::PetscLibType,sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode) 
Ends a generalized scatter from one vector to another. Call
after first calling `VecScatterBegin()`.

Neighbor-wise Collective

Input Parameters:
- `sf`   - scatter context generated by `VecScatterCreate()`
- `x`    - the vector from which we scatter
- `y`    - the vector to which we scatter
- `addv` - one of `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`
- `mode` - the scattering mode, usually `SCATTER_FORWARD`.  The available modes are: `SCATTER_FORWARD`, `SCATTER_REVERSE`

Level: intermediate

-seealso: [](sec_scatter), `VecScatter`, `VecScatterBegin()`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecScatterEnd"))
"""
function VecScatterEnd(petsclib::PetscLibType, sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode) end

@for_petsc function VecScatterEnd(petsclib::$UnionPetscLib, sf::VecScatter, x::PetscVec, y::PetscVec, addv::InsertMode, mode::ScatterMode )

    @chk ccall(
               (:VecScatterEnd, $petsc_library),
               PetscErrorCode,
               (VecScatter, CVec, CVec, InsertMode, ScatterMode),
               sf, x, y, addv, mode,
              )


	return nothing
end 

"""
	stepmax::PetscReal = VecStepMaxBounded(petsclib::PetscLibType,X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec) 
See below

Collective

Input Parameters:
- `X`  - vector with no negative entries
- `XL` - lower bounds
- `XU` - upper bounds
- `DX` - step direction, can have negative, positive or zero entries

Output Parameter:
- `stepmax` - minimum value so that X[i] + stepmax*DX[i] <= XL[i]  or  XU[i] <= X[i] + stepmax*DX[i]

Level: intermediate

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecStepMaxBounded"))
"""
function VecStepMaxBounded(petsclib::PetscLibType, X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec) end

@for_petsc function VecStepMaxBounded(petsclib::$UnionPetscLib, X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec )
	stepmax_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStepMaxBounded, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, Ptr{$PetscReal}),
               X, DX, XL, XU, stepmax_,
              )

	stepmax = stepmax_[]

	return stepmax
end 

"""
	boundmin::PetscReal,wolfemin::PetscReal,boundmax::PetscReal = VecStepBoundInfo(petsclib::PetscLibType,X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec) 
See below

Collective

Input Parameters:
- `X`  - vector with no negative entries
- `XL` - lower bounds
- `XU` - upper bounds
- `DX` - step direction, can have negative, positive or zero entries

Output Parameters:
- `boundmin` - (may be `NULL` this it is not computed) maximum value so that   XL[i] <= X[i] + boundmax*DX[i] <= XU[i]
- `wolfemin` - (may be `NULL` this it is not computed)
- `boundmax` - (may be `NULL` this it is not computed) minimum value so that X[i] + boundmax*DX[i] <= XL[i]  or  XU[i] <= X[i] + boundmax*DX[i]

Level: advanced

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecStepBoundInfo"))
"""
function VecStepBoundInfo(petsclib::PetscLibType, X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec) end

@for_petsc function VecStepBoundInfo(petsclib::$UnionPetscLib, X::PetscVec, DX::PetscVec, XL::PetscVec, XU::PetscVec )
	boundmin_ = Ref{$PetscReal}()
	wolfemin_ = Ref{$PetscReal}()
	boundmax_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStepBoundInfo, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               X, DX, XL, XU, boundmin_, wolfemin_, boundmax_,
              )

	boundmin = boundmin_[]
	wolfemin = wolfemin_[]
	boundmax = boundmax_[]

	return boundmin,wolfemin,boundmax
end 

"""
	step::PetscReal = VecStepMax(petsclib::PetscLibType,X::PetscVec, DX::PetscVec) 
Returns the largest value so that x[i] + step*DX[i] >= 0 for all i

Collective

Input Parameters:
- `X`  - vector with no negative entries
- `DX` - a step direction, can have negative, positive or zero entries

Output Parameter:
- `step` - largest value such that x[i] + step*DX[i] >= 0 for all i

Level: advanced

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecStepMax"))
"""
function VecStepMax(petsclib::PetscLibType, X::PetscVec, DX::PetscVec) end

@for_petsc function VecStepMax(petsclib::$UnionPetscLib, X::PetscVec, DX::PetscVec )
	step_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecStepMax, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscReal}),
               X, DX, step_,
              )

	step = step_[]

	return step
end 

"""
	VecsDestroy(petsclib::PetscLibType,x::Vecs) 

# External Links
$(_doc_external("Vec/VecsDestroy"))
"""
function VecsDestroy(petsclib::PetscLibType, x::Vecs) end

@for_petsc function VecsDestroy(petsclib::$UnionPetscLib, x::Vecs )

    @chk ccall(
               (:VecsDestroy, $petsc_library),
               PetscErrorCode,
               (Vecs,),
               x,
              )


	return nothing
end 

"""
	x::Vecs = VecsCreateSeq(petsclib::PetscLibType,comm::MPI_Comm, p::PetscInt, m::PetscInt) 

# External Links
$(_doc_external("Vec/VecsCreateSeq"))
"""
function VecsCreateSeq(petsclib::PetscLibType, comm::MPI_Comm, p::PetscInt, m::PetscInt) end

@for_petsc function VecsCreateSeq(petsclib::$UnionPetscLib, comm::MPI_Comm, p::$PetscInt, m::$PetscInt )
	x_ = Ref{Vecs}()

    @chk ccall(
               (:VecsCreateSeq, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{Vecs}),
               comm, p, m, x_,
              )

	x = x_[]

	return x
end 

"""
	a::PetscScalar,x::Vecs = VecsCreateSeqWithArray(petsclib::PetscLibType,comm::MPI_Comm, p::PetscInt, m::PetscInt) 

# External Links
$(_doc_external("Vec/VecsCreateSeqWithArray"))
"""
function VecsCreateSeqWithArray(petsclib::PetscLibType, comm::MPI_Comm, p::PetscInt, m::PetscInt) end

@for_petsc function VecsCreateSeqWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, p::$PetscInt, m::$PetscInt )
	a_ = Ref{$PetscScalar}()
	x_ = Ref{Vecs}()

    @chk ccall(
               (:VecsCreateSeqWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{Vecs}),
               comm, p, m, a_, x_,
              )

	a = a_[]
	x = x_[]

	return a,x
end 

"""
	y::Vecs = VecsDuplicate(petsclib::PetscLibType,x::Vecs) 

# External Links
$(_doc_external("Vec/VecsDuplicate"))
"""
function VecsDuplicate(petsclib::PetscLibType, x::Vecs) end

@for_petsc function VecsDuplicate(petsclib::$UnionPetscLib, x::Vecs )
	y_ = Ref{Vecs}()

    @chk ccall(
               (:VecsDuplicate, $petsc_library),
               PetscErrorCode,
               (Vecs, Ptr{Vecs}),
               x, y_,
              )

	y = y_[]

	return y
end 

"""
	VecSetValuesSection(petsclib::PetscLibType,v::PetscVec, s::PetscSection, point::PetscInt, values::Vector{PetscScalar}, mode::InsertMode) 
Sets all the values associated with a given point, according to the section, in the given `Vec`

Not Collective

Input Parameters:
- `v`      - the `Vec`
- `s`      - the organizing `PetscSection`
- `point`  - the point
- `values` - the array of input values
- `mode`   - the insertion mode, either `ADD_VALUES` or `INSERT_VALUES`

Level: developer

-seealso: `PetscSection`, `PetscSectionCreate()`, `VecGetValuesSection()`

# External Links
$(_doc_external("Vec/VecSetValuesSection"))
"""
function VecSetValuesSection(petsclib::PetscLibType, v::PetscVec, s::PetscSection, point::PetscInt, values::Vector{PetscScalar}, mode::InsertMode) end

@for_petsc function VecSetValuesSection(petsclib::$UnionPetscLib, v::PetscVec, s::PetscSection, point::$PetscInt, values::Vector{$PetscScalar}, mode::InsertMode )

    @chk ccall(
               (:VecSetValuesSection, $petsc_library),
               PetscErrorCode,
               (CVec, PetscSection, $PetscInt, Ptr{$PetscScalar}, InsertMode),
               v, s, point, values, mode,
              )


	return nothing
end 

"""
	VecSetDM(petsclib::PetscLibType,v::PetscVec, dm::PetscDM) 
Sets the `DM` defining the data layout of the vector.

Not Collective

Input Parameters:
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
$(_doc_external("Dm/VecSetDM"))
"""
function VecSetDM(petsclib::PetscLibType, v::PetscVec, dm::PetscDM) end

@for_petsc function VecSetDM(petsclib::$UnionPetscLib, v::PetscVec, dm::PetscDM )

    @chk ccall(
               (:VecSetDM, $petsc_library),
               PetscErrorCode,
               (CVec, CDM),
               v, dm,
              )


	return nothing
end 

"""
	VecSFischer(petsclib::PetscLibType,X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, mu::PetscReal, FB::PetscVec) 
Evaluates the Smoothed Fischer
complementarity problems.

Logically Collective

Input Parameters:
- `X`  - current point
- `F`  - function evaluated at x
- `L`  - lower bounds
- `U`  - upper bounds
- `mu` - smoothing parameter

Output Parameter:
- `FB` - The Smoothed Fischer-Burmeister function vector

-seealso: `Vec`, `VecFischer()`, `MatDFischer()`, `MatDSFischer()`

# External Links
$(_doc_external("Tao/VecSFischer"))
"""
function VecSFischer(petsclib::PetscLibType, X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, mu::PetscReal, FB::PetscVec) end

@for_petsc function VecSFischer(petsclib::$UnionPetscLib, X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, mu::$PetscReal, FB::PetscVec )

    @chk ccall(
               (:VecSFischer, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, $PetscReal, CVec),
               X, F, L, U, mu, FB,
              )


	return nothing
end 

"""
	VecSetValue(petsclib::PetscLibType,v::PetscVec, i::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Vec/VecSetValue"))
"""
function VecSetValue(petsclib::PetscLibType, v::PetscVec, i::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function VecSetValue(petsclib::$UnionPetscLib, v::PetscVec, i::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:VecSetValue, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscScalar, InsertMode),
               v, i, va, mode,
              )


	return nothing
end 

"""
	VecSetValueLocal(petsclib::PetscLibType,v::PetscVec, i::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Vec/VecSetValueLocal"))
"""
function VecSetValueLocal(petsclib::PetscLibType, v::PetscVec, i::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function VecSetValueLocal(petsclib::$UnionPetscLib, v::PetscVec, i::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:VecSetValueLocal, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscScalar, InsertMode),
               v, i, va, mode,
              )


	return nothing
end 

"""
	VecSetErrorIfLocked(petsclib::PetscLibType,x::PetscVec, arg::PetscInt) 

# External Links
$(_doc_external("Vec/VecSetErrorIfLocked"))
"""
function VecSetErrorIfLocked(petsclib::PetscLibType, x::PetscVec, arg::PetscInt) end

@for_petsc function VecSetErrorIfLocked(petsclib::$UnionPetscLib, x::PetscVec, arg::$PetscInt )

    @chk ccall(
               (:VecSetErrorIfLocked, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt),
               x, arg,
              )


	return nothing
end 

