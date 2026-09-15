"""
	isout::IS = ISAllGather(petsclib::PetscLibType,is::AbstractIS) 
Given an index set `IS` on each processor, generates a large
index set (same on each processor) by concatenating together each
processors index set.

Collective

Input Parameter:
- `is` - the distributed index set

Output Parameter:
- `isout` - the concatenated index set (same on all processors)

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`

# External Links
$(_doc_external("IS/ISAllGather"))
"""
function ISAllGather(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISAllGather(petsclib::$UnionPetscLib, is::AbstractIS )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISAllGather, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{CIS}),
               is, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	outN::PetscInt,outindices::Ptr{ISColoringValue} = ISAllGatherColors(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, lindices::Vector{ISColoringValue}) 
Given a set of colors on each processor, generates a large
set (same on each processor) by concatenating together each processors colors

Collective

Input Parameters:
- `comm`     - communicator to share the indices
- `n`        - local size of set
- `lindices` - local colors

Output Parameters:
- `outN`       - total number of indices
- `outindices` - all of the colors

Level: intermediate

-seealso: `ISColoringValue`, `ISColoring()`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`, `ISAllGather()`

# External Links
$(_doc_external("IS/ISAllGatherColors"))
"""
function ISAllGatherColors(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, lindices::Vector{ISColoringValue}) end

@for_petsc function ISAllGatherColors(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, lindices::Vector{ISColoringValue} )
	outN_ = Ref{$PetscInt}()
	outindices_ = Ref{Ptr{ISColoringValue}}()

    @chk ccall(
               (:ISAllGatherColors, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{ISColoringValue}, Ptr{$PetscInt}, Ptr{Ptr{ISColoringValue}}),
               comm, n, lindices, outN_, outindices_,
              )

	outN = outN_[]
	outindices = outindices_[]

	return outN,outindices
end 

"""
	idx::Vector{PetscInt} = ISBlockGetIndices(petsclib::PetscLibType,is::AbstractIS) 
Gets the indices associated with each block in an `ISBLOCK`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `idx` - the integer indices, one for each block and count of block not indices

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISBLOCK`, `ISGetIndices()`, `ISBlockRestoreIndices()`, `ISBlockSetIndices()`, `ISCreateBlock()`

# External Links
$(_doc_external("IS/ISBlockGetIndices"))
"""
function ISBlockGetIndices(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISBlockGetIndices(petsclib::$UnionPetscLib, is::AbstractIS )
	idx_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:ISBlockGetIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, idx_,
              )

	n = div(ISGetLocalSize(petsclib, is), ISGetBlockSize(petsclib, is))
	idx = unsafe_wrap(Array, idx_[], n; own = false)

	return idx
end 

"""
	size::PetscInt = ISBlockGetLocalSize(petsclib::PetscLibType,is::AbstractIS) 
Returns the local number of blocks in the index set of `ISType` `ISBLOCK`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `size` - the local number of blocks

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGetBlockSize()`, `ISBlockGetSize()`, `ISGetSize()`, `ISCreateBlock()`, `ISBLOCK`

# External Links
$(_doc_external("IS/ISBlockGetLocalSize"))
"""
function ISBlockGetLocalSize(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISBlockGetLocalSize(petsclib::$UnionPetscLib, is::AbstractIS )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISBlockGetLocalSize, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}),
               is, size_,
              )

	size = size_[]

	return size
end 

"""
	size::PetscInt = ISBlockGetSize(petsclib::PetscLibType,is::AbstractIS) 
Returns the global number of blocks in parallel in the index set of `ISType` `ISBLOCK`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `size` - the global number of blocks

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGetBlockSize()`, `ISBlockGetLocalSize()`, `ISGetSize()`, `ISCreateBlock()`, `ISBLOCK`

# External Links
$(_doc_external("IS/ISBlockGetSize"))
"""
function ISBlockGetSize(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISBlockGetSize(petsclib::$UnionPetscLib, is::AbstractIS )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISBlockGetSize, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}),
               is, size_,
              )

	size = size_[]

	return size
end 

"""
	idx::Ptr{PetscInt} = ISBlockRestoreIndices(petsclib::PetscLibType,is::AbstractIS) 
Restores the indices associated with each block  in an `ISBLOCK` obtained with `ISBlockGetIndices()`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `idx` - the integer indices

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISBLOCK`, `ISRestoreIndices()`, `ISBlockGetIndices()`

# External Links
$(_doc_external("IS/ISBlockRestoreIndices"))
"""
function ISBlockRestoreIndices(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISBlockRestoreIndices(petsclib::$UnionPetscLib, is::AbstractIS )
	idx_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:ISBlockRestoreIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, idx_,
              )

	idx = idx_[]

	return idx
end 

"""
	ISBlockSetIndices(petsclib::PetscLibType,is::AbstractIS, bs::PetscInt, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) 
Set integers representing blocks of indices in an index set of `ISType` `ISBLOCK`

Collective

Input Parameters:
- `is`   - the index set
- `bs`   - number of elements in each block
- `n`    - the length of the index set (the number of blocks)
- `idx`  - the list of integers, one for each block, the integers contain the index of the first index of each block divided by the block size
- `mode` - see `PetscCopyMode`, only `PETSC_COPY_VALUES` and `PETSC_OWN_POINTER` are supported

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISCreateStride()`, `ISCreateGeneral()`, `ISAllGather()`, `ISCreateBlock()`, `ISBLOCK`, `ISGeneralSetIndices()`

# External Links
$(_doc_external("IS/ISBlockSetIndices"))
"""
function ISBlockSetIndices(petsclib::PetscLibType, is::AbstractIS, bs::PetscInt, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) end

@for_petsc function ISBlockSetIndices(petsclib::$UnionPetscLib, is::AbstractIS, bs::$PetscInt, n::$PetscInt, idx::Vector{$PetscInt}, mode::PetscCopyMode )

    @chk ccall(
               (:ISBlockSetIndices, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode),
               is, bs, n, idx, mode,
              )


	return nothing
end 

"""
	rows::IS = ISBuildTwoSided(petsclib::PetscLibType,ito::AbstractIS, toindx::AbstractIS) 
Takes an `IS` that describes where each element will be mapped globally over all ranks.
Generates an `IS` that contains new numbers from remote or local on the `IS`.

Collective

Input Parameters:
- `ito`    - an `IS` describes to which rank each entry will be mapped. Negative target rank will be ignored
- `toindx` - an `IS` describes what indices should send. `NULL` means sending natural numbering

Output Parameter:
- `rows` - contains new numbers from remote or local

Level: advanced

-seealso: [](sec_scatter), `IS`, `MatPartitioningCreate()`, `ISPartitioningToNumbering()`, `ISPartitioningCount()`

# External Links
$(_doc_external("IS/ISBuildTwoSided"))
"""
function ISBuildTwoSided(petsclib::PetscLibType, ito::AbstractIS, toindx::AbstractIS) end

@for_petsc function ISBuildTwoSided(petsclib::$UnionPetscLib, ito::AbstractIS, toindx::AbstractIS )
	rows_ = Ref{CIS}()

    @chk ccall(
               (:ISBuildTwoSided, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               ito, toindx, rows_,
              )

	rows = IS(rows_[], petsclib)

	return rows
end 

"""
	ISClearInfoCache(petsclib::PetscLibType,is::AbstractIS, clear_permanent_local::PetscBool) 
clear the cache of computed index set properties

Not Collective

Input Parameters:
- `is`                    - the index set
- `clear_permanent_local` - whether to remove the permanent status of local properties

Level: developer

-seealso: `IS`, `ISInfo`, `ISInfoType`, `ISSetInfo()`

# External Links
$(_doc_external("IS/ISClearInfoCache"))
"""
function ISClearInfoCache(petsclib::PetscLibType, is::AbstractIS, clear_permanent_local::PetscBool) end

@for_petsc function ISClearInfoCache(petsclib::$UnionPetscLib, is::AbstractIS, clear_permanent_local::PetscBool )

    @chk ccall(
               (:ISClearInfoCache, $petsc_library),
               PetscErrorCode,
               (CIS, PetscBool),
               is, clear_permanent_local,
              )


	return nothing
end 

"""
	isout::IS = ISComplement(petsclib::PetscLibType,is::AbstractIS, nmin::PetscInt, nmax::PetscInt) 
Given an index set `IS` generates the complement index set. That is
all indices that are NOT in the given set.

Collective

Input Parameters:
- `is`   - the index set
- `nmin` - the first index desired in the local part of the complement
- `nmax` - the largest index desired in the local part of the complement (note that all indices in `is` must be greater or equal to `nmin` and less than `nmax`)

Output Parameter:
- `isout` - the complement

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`, `ISAllGather()`

# External Links
$(_doc_external("IS/ISComplement"))
"""
function ISComplement(petsclib::PetscLibType, is::AbstractIS, nmin::PetscInt, nmax::PetscInt) end

@for_petsc function ISComplement(petsclib::$UnionPetscLib, is::AbstractIS, nmin::$PetscInt, nmax::$PetscInt )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISComplement, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, Ptr{CIS}),
               is, nmin, nmax, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	T::IS = ISComplementVec(petsclib::PetscLibType,S::AbstractIS, V::AbstractPetscVec) 
Creates the complement of the index set relative to a layout defined by a `Vec`

Collective

Input Parameters:
- `S` - a PETSc `IS`
- `V` - the reference vector space

Output Parameter:
- `T` - the complement of S

Level: advanced

-seealso: `IS`, `Vec`, `ISCreateGeneral()`

# External Links
$(_doc_external("Vec/ISComplementVec"))
"""
function ISComplementVec(petsclib::PetscLibType, S::AbstractIS, V::AbstractPetscVec) end

@for_petsc function ISComplementVec(petsclib::$UnionPetscLib, S::AbstractIS, V::AbstractPetscVec )
	T_ = Ref{CIS}()

    @chk ccall(
               (:ISComplementVec, $petsc_library),
               PetscErrorCode,
               (CIS, CVec, Ptr{CIS}),
               S, V, T_,
              )

	T = IS(T_[], petsclib)

	return T
end 

"""
	ISCompressIndicesGeneral(petsclib::PetscLibType,n::PetscInt, nkeys::PetscInt, bs::PetscInt, imax::PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS}) 
convert the indices of an array of `IS` into an array of `ISGENERAL` of block indices

Input Parameters:
- `n`     - maximum possible length of the index set
- `nkeys` - expected number of keys when using `PETSC_USE_CTABLE`
- `bs`    - the size of block
- `imax`  - the number of index sets
- `is_in` - the non-blocked array of index sets

Output Parameter:
- `is_out` - the blocked new index set, as `ISGENERAL`, not as `ISBLOCK`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGENERAL`, `ISExpandIndicesGeneral()`

# External Links
$(_doc_external("IS/ISCompressIndicesGeneral"))
"""
function ISCompressIndicesGeneral(petsclib::PetscLibType, n::PetscInt, nkeys::PetscInt, bs::PetscInt, imax::PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS}) end

@for_petsc function ISCompressIndicesGeneral(petsclib::$UnionPetscLib, n::$PetscInt, nkeys::$PetscInt, bs::$PetscInt, imax::$PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS} )

    @chk ccall(
               (:ISCompressIndicesGeneral, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CIS}, Ptr{CIS}),
               n, nkeys, bs, imax, is_in, is_out,
              )


	return nothing
end 

"""
	isout::IS = ISConcatenate(petsclib::PetscLibType,comm::MPI_Comm, len::PetscInt, islist::Vector{<:AbstractIS}) 
Forms a new `IS` by locally concatenating the indices from an `IS` list without reordering.

Collective

Input Parameters:
- `comm`   - communicator of the concatenated `IS`.
- `len`    - size of islist array (nonnegative)
- `islist` - array of index sets

Output Parameter:
- `isout` - The concatenated index set; empty, if `len` == 0.

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISDifference()`, `ISSum()`, `ISExpand()`, `ISIntersect()`

# External Links
$(_doc_external("IS/ISConcatenate"))
"""
function ISConcatenate(petsclib::PetscLibType, comm::MPI_Comm, len::PetscInt, islist::Vector{<:AbstractIS}) end

@for_petsc function ISConcatenate(petsclib::$UnionPetscLib, comm::MPI_Comm, len::$PetscInt, islist::Vector{<:AbstractIS} )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISConcatenate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CIS}, Ptr{CIS}),
               comm, len, islist, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	start::PetscInt,contig::PetscBool = ISContiguousLocal(petsclib::PetscLibType,is::AbstractIS, gstart::PetscInt, gend::PetscInt) 
Locates an index set with contiguous range within a global range, if possible

Not Collective

Input Parameters:
- `is`     - the index set
- `gstart` - global start
- `gend`   - global end

Output Parameters:
- `start`  - start of contiguous block, as an offset from `gstart`
- `contig` - `PETSC_TRUE` if the index set refers to contiguous entries on this process, else `PETSC_FALSE`

Level: developer

-seealso: `IS`, `ISGetLocalSize()`, `VecGetOwnershipRange()`

# External Links
$(_doc_external("IS/ISContiguousLocal"))
"""
function ISContiguousLocal(petsclib::PetscLibType, is::AbstractIS, gstart::PetscInt, gend::PetscInt) end

@for_petsc function ISContiguousLocal(petsclib::$UnionPetscLib, is::AbstractIS, gstart::$PetscInt, gend::$PetscInt )
	start_ = Ref{$PetscInt}()
	contig_ = Ref{PetscBool}()

    @chk ccall(
               (:ISContiguousLocal, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               is, gstart, gend, start_, contig_,
              )

	start = start_[]
	contig = contig_[]

	return start,contig
end 

"""
	ISCopy(petsclib::PetscLibType,is::AbstractIS, isy::AbstractIS) 
Copies an index set.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `isy` - the copy of the index set

Level: beginner

-seealso: `IS`, `ISDuplicate()`, `ISShift()`

# External Links
$(_doc_external("IS/ISCopy"))
"""
function ISCopy(petsclib::PetscLibType, is::AbstractIS, isy::AbstractIS) end

@for_petsc function ISCopy(petsclib::$UnionPetscLib, is::AbstractIS, isy::AbstractIS )

    @chk ccall(
               (:ISCopy, $petsc_library),
               PetscErrorCode,
               (CIS, CIS),
               is, isy,
              )


	return nothing
end 

"""
	is::IS = ISCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create an index set object. `IS`, index sets, are PETSc objects used to do efficient indexing into other data structures such as `Vec` and `Mat`

Collective

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `is` - the new index set

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISType()`, `ISSetType()`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`, `ISAllGather()`

# External Links
$(_doc_external("IS/ISCreate"))
"""
function ISCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function ISCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	is_ = Ref{CIS}()

    @chk ccall(
               (:ISCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CIS}),
               comm, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	is::IS = ISCreateBlock(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) 
Creates a data structure for an index set containing
a list of integers. Each integer represents a fixed block size set of indices.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `bs`   - number of elements in each block
- `n`    - the length of the index set (the number of blocks)
- `idx`  - the list of integers, one for each block, the integers contain the index of the first entry of each block divided by the block size
- `mode` - see `PetscCopyMode`, only `PETSC_COPY_VALUES` and `PETSC_OWN_POINTER` are supported in this routine

Output Parameter:
- `is` - the new index set

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISCreateStride()`, `ISCreateGeneral()`, `ISAllGather()`, `ISBlockSetIndices()`, `ISBLOCK`, `ISGENERAL`

# External Links
$(_doc_external("IS/ISCreateBlock"))
"""
function ISCreateBlock(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) end

@for_petsc function ISCreateBlock(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, idx::Vector{$PetscInt}, mode::PetscCopyMode )
	is_ = Ref{CIS}()

    @chk ccall(
               (:ISCreateBlock, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{CIS}),
               comm, bs, n, idx, mode, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	is::IS = ISCreateGeneral(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) 
Creates a data structure for an index set containing a list of integers.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `n`    - the length of the index set
- `idx`  - the list of integers
- `mode` - `PETSC_COPY_VALUES`, `PETSC_OWN_POINTER`, or `PETSC_USE_POINTER`; see `PetscCopyMode` for meaning of this flag.

Output Parameter:
- `is` - the new index set

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISGENERAL`, `ISCreateStride()`, `ISCreateBlock()`, `ISAllGather()`, `PETSC_COPY_VALUES`, `PETSC_OWN_POINTER`,
`PETSC_USE_POINTER`, `PetscCopyMode`, `ISGeneralSetIndicesFromMask()`

# External Links
$(_doc_external("IS/ISCreateGeneral"))
"""
function ISCreateGeneral(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) end

@for_petsc function ISCreateGeneral(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, idx::Vector{$PetscInt}, mode::PetscCopyMode )
	is_ = Ref{CIS}()

    @chk ccall(
               (:ISCreateGeneral, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{CIS}),
               comm, n, idx, mode, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	is::IS = ISCreateStride(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, first::PetscInt, step::PetscInt) 
Creates a data structure for an index set containing a list of evenly spaced integers.

Collective

Input Parameters:
- `comm`  - the MPI communicator
- `n`     - the length of the locally owned portion of the index set
- `first` - the first element of the locally owned portion of the index set
- `step`  - the change to the next index

Output Parameter:
- `is` - the new index set

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISStrideSetStride()`, `ISCreateGeneral()`, `ISCreateBlock()`, `ISAllGather()`, `ISSTRIDE`

# External Links
$(_doc_external("IS/ISCreateStride"))
"""
function ISCreateStride(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, first::PetscInt, step::PetscInt) end

@for_petsc function ISCreateStride(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, first::$PetscInt, step::$PetscInt )
	is_ = Ref{CIS}()

    @chk ccall(
               (:ISCreateStride, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{CIS}),
               comm, n, first, step, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	subis::IS = ISCreateSubIS(petsclib::PetscLibType,is::AbstractIS, comps::AbstractIS) 
Create a sub index set from a global index set selecting some components.

Collective

Input Parameters:
- `is`    - the index set
- `comps` - which components we will extract from `is`

Output Parameters:
- `subis` - the new sub index set

Example usage:
We have an index set `is` living on 3 processes with the following values:
| 4 9 0 | 2 6 7 | 10 11 1|
and another index set `comps` used to indicate which components of is  we want to take,
| 7 5  | 1 2 | 0 4|
The output index set `subis` should look like:
| 11 7 | 9 0 | 4 6|

Level: intermediate

-seealso: `IS`, `VecGetSubVector()`, `MatCreateSubMatrix()`

# External Links
$(_doc_external("IS/ISCreateSubIS"))
"""
function ISCreateSubIS(petsclib::PetscLibType, is::AbstractIS, comps::AbstractIS) end

@for_petsc function ISCreateSubIS(petsclib::$UnionPetscLib, is::AbstractIS, comps::AbstractIS )
	subis_ = Ref{CIS}()

    @chk ccall(
               (:ISCreateSubIS, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               is, comps, subis_,
              )

	subis = IS(subis_[], petsclib)

	return subis
end 

"""
	ISDestroy(petsclib::PetscLibType,is::AbstractIS) 
Destroys an index set.

Collective

Input Parameter:
- `is` - the index set

Level: beginner

-seealso: `IS`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`

# External Links
$(_doc_external("IS/ISDestroy"))
"""
function ISDestroy(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISDestroy(petsclib::$UnionPetscLib, is::AbstractIS )
	is_ = Ref(is.ptr)

    @chk ccall(
               (:ISDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CIS},),
               is_,
              )

	is.ptr = C_NULL

	return nothing
end 

"""
	isout::IS = ISDifference(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Computes the difference between two index sets.

Collective

Input Parameters:
- `is1` - first index, to have items removed from it
- `is2` - index values to be removed

Output Parameter:
- `isout` - is1 - is2

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISDestroy()`, `ISView()`, `ISSum()`, `ISExpand()`

# External Links
$(_doc_external("IS/ISDifference"))
"""
function ISDifference(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISDifference(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISDifference, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               is1, is2, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	newIS::IS = ISDuplicate(petsclib::PetscLibType,is::AbstractIS) 
Creates a duplicate copy of an index set.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `newIS` - the copy of the index set

Level: beginner

-seealso: `IS`, `ISCreateGeneral()`, `ISCopy()`

# External Links
$(_doc_external("IS/ISDuplicate"))
"""
function ISDuplicate(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISDuplicate(petsclib::$UnionPetscLib, is::AbstractIS )
	newIS_ = Ref{CIS}()

    @chk ccall(
               (:ISDuplicate, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{CIS}),
               is, newIS_,
              )

	newIS = IS(newIS_[], petsclib)

	return newIS
end 

"""
	c::IS = ISEmbed(petsclib::PetscLibType,a::AbstractIS, b::AbstractIS, drop::PetscBool) 
Embed `IS` `a` into `IS` `b` by finding the locations in `b` that have the same indices as in `a`.
If `c` is the `IS` of these locations, we have `a = b*c`, regarded as a composition of the
corresponding `ISLocalToGlobalMapping`.

Not Collective

Input Parameters:
- `a`    - `IS` to embed
- `b`    - `IS` to embed into
- `drop` - flag indicating whether to drop indices of `a` that are not in `b`.

Output Parameter:
- `c` - local embedding indices

Level: developer

-seealso: `IS`, `ISLocalToGlobalMapping`

# External Links
$(_doc_external("IS/ISEmbed"))
"""
function ISEmbed(petsclib::PetscLibType, a::AbstractIS, b::AbstractIS, drop::PetscBool) end

@for_petsc function ISEmbed(petsclib::$UnionPetscLib, a::AbstractIS, b::AbstractIS, drop::PetscBool )
	c_ = Ref{CIS}()

    @chk ccall(
               (:ISEmbed, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, PetscBool, Ptr{CIS}),
               a, b, drop, c_,
              )

	c = IS(c_[], petsclib)

	return c
end 

"""
	flg::PetscBool = ISEqual(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Compares if two index sets have the same set of indices.

Collective

Input Parameters:
- `is1` - first index set to compare
- `is2` - second index set to compare

Output Parameter:
- `flg` - output flag, either `PETSC_TRUE` (if both index sets have the
same indices), or `PETSC_FALSE` if the index sets differ by size
or by the set of indices)

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISEqualUnsorted()`

# External Links
$(_doc_external("IS/ISEqual"))
"""
function ISEqual(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISEqual(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:ISEqual, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{PetscBool}),
               is1, is2, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = ISEqualUnsorted(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Compares if two index sets have the same indices.

Collective

Input Parameters:
- `is1` - first index set to compare
- `is2` - second index set to compare

Output Parameter:
- `flg` - output flag, either `PETSC_TRUE` (if both index sets have the
same indices), or `PETSC_FALSE` if the index sets differ by size
or by the set of indices)

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISEqual()`

# External Links
$(_doc_external("IS/ISEqualUnsorted"))
"""
function ISEqualUnsorted(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISEqualUnsorted(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:ISEqualUnsorted, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{PetscBool}),
               is1, is2, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	isout::IS = ISExpand(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Computes the union of two index sets, by concatenating 2 lists and
removing duplicates.

Collective

Input Parameters:
- `is1` - first index set
- `is2` - index values to be added

Output Parameter:
- `isout` - `is1` + `is2` The index set `is2` is appended to `is1` removing duplicates

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISDestroy()`, `ISView()`, `ISDifference()`, `ISSum()`, `ISIntersect()`

# External Links
$(_doc_external("IS/ISExpand"))
"""
function ISExpand(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISExpand(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISExpand, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               is1, is2, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	ISExpandIndicesGeneral(petsclib::PetscLibType,n::PetscInt, nkeys::PetscInt, bs::PetscInt, imax::PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS}) 
convert the indices of an array `IS` into non

Input Parameters:
- `n`     - the length of the index set (not being used)
- `nkeys` - expected number of keys when `PETSC_USE_CTABLE` is used
- `bs`    - the size of block
- `imax`  - the number of index sets
- `is_in` - the blocked array of index sets, must be as large as `imax`

Output Parameter:
- `is_out` - the non-blocked new index set, as `ISGENERAL`, must be as large as `imax`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGENERAL`, `ISCompressIndicesGeneral()`

# External Links
$(_doc_external("IS/ISExpandIndicesGeneral"))
"""
function ISExpandIndicesGeneral(petsclib::PetscLibType, n::PetscInt, nkeys::PetscInt, bs::PetscInt, imax::PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS}) end

@for_petsc function ISExpandIndicesGeneral(petsclib::$UnionPetscLib, n::$PetscInt, nkeys::$PetscInt, bs::$PetscInt, imax::$PetscInt, is_in::Vector{<:AbstractIS}, is_out::Vector{<:AbstractIS} )

    @chk ccall(
               (:ISExpandIndicesGeneral, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CIS}, Ptr{CIS}),
               n, nkeys, bs, imax, is_in, is_out,
              )


	return nothing
end 

"""
	ISFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `IS` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: `PetscFinalize()`

# External Links
$(_doc_external("Vec/ISFinalizePackage"))
"""
function ISFinalizePackage(petsclib::PetscLibType) end

@for_petsc function ISFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:ISFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	ISGeneralFilter(petsclib::PetscLibType,is::AbstractIS, start::PetscInt, end_::PetscInt) 
Remove all indices outside of [start, end) from an `ISGENERAL`

Collective

Input Parameters:
- `is`    - the index set
- `start` - the lowest index kept
- `end`   - one more than the highest index kept, `start` \\le `end`

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISGENERAL`, `ISCreateGeneral()`, `ISGeneralSetIndices()`

# External Links
$(_doc_external("IS/ISGeneralFilter"))
"""
function ISGeneralFilter(petsclib::PetscLibType, is::AbstractIS, start::PetscInt, end_::PetscInt) end

@for_petsc function ISGeneralFilter(petsclib::$UnionPetscLib, is::AbstractIS, start::$PetscInt, end_::$PetscInt )

    @chk ccall(
               (:ISGeneralFilter, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt),
               is, start, end_,
              )


	return nothing
end 

"""
	ISGeneralSetIndices(petsclib::PetscLibType,is::AbstractIS, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) 
Sets the indices for an `ISGENERAL` index set

Logically Collective

Input Parameters:
- `is`   - the index set
- `n`    - the length of the index set
- `idx`  - the list of integers
- `mode` - see `PetscCopyMode` for meaning of this flag.

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISBLOCK`, `ISCreateGeneral()`, `ISGeneralSetIndicesFromMask()`, `ISBlockSetIndices()`, `ISGENERAL`, `PetscCopyMode`

# External Links
$(_doc_external("IS/ISGeneralSetIndices"))
"""
function ISGeneralSetIndices(petsclib::PetscLibType, is::AbstractIS, n::PetscInt, idx::Vector{PetscInt}, mode::PetscCopyMode) end

@for_petsc function ISGeneralSetIndices(petsclib::$UnionPetscLib, is::AbstractIS, n::$PetscInt, idx::Vector{$PetscInt}, mode::PetscCopyMode )

    @chk ccall(
               (:ISGeneralSetIndices, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, Ptr{$PetscInt}, PetscCopyMode),
               is, n, idx, mode,
              )


	return nothing
end 

"""
	ISGeneralSetIndicesFromMask(petsclib::PetscLibType,is::AbstractIS, rstart::PetscInt, rend::PetscInt, mask::Vector{PetscBool}) 
Sets the indices for an `ISGENERAL` index set using a boolean mask

Collective

Input Parameters:
- `is`     - the index set
- `rstart` - the range start index (inclusive)
- `rend`   - the range end index (exclusive)
- `mask`   - the boolean mask array of length rend-rstart, indices will be set for each `PETSC_TRUE` value in the array

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISCreateGeneral()`, `ISGeneralSetIndices()`, `ISGENERAL`

# External Links
$(_doc_external("IS/ISGeneralSetIndicesFromMask"))
"""
function ISGeneralSetIndicesFromMask(petsclib::PetscLibType, is::AbstractIS, rstart::PetscInt, rend::PetscInt, mask::Vector{PetscBool}) end

@for_petsc function ISGeneralSetIndicesFromMask(petsclib::$UnionPetscLib, is::AbstractIS, rstart::$PetscInt, rend::$PetscInt, mask::Vector{PetscBool} )

    @chk ccall(
               (:ISGeneralSetIndicesFromMask, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, Ptr{PetscBool}),
               is, rstart, rend, mask,
              )


	return nothing
end 

"""
	size::PetscInt = ISGetBlockSize(petsclib::PetscLibType,is::AbstractIS) 
Returns the number of elements in a block.

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `size` - the number of elements in a block

Level: intermediate

-seealso: `IS`, `ISBlockGetSize()`, `ISGetSize()`, `ISCreateBlock()`, `ISSetBlockSize()`

# External Links
$(_doc_external("IS/ISGetBlockSize"))
"""
function ISGetBlockSize(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetBlockSize(petsclib::$UnionPetscLib, is::AbstractIS )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGetBlockSize, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}),
               is, size_,
              )

	size = size_[]

	return size
end 

"""
	compress::PetscBool = ISGetCompressOutput(petsclib::PetscLibType,is::AbstractIS) 
Returns the flag for output compression

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `compress` - the flag to compress output

Level: intermediate

-seealso: `IS`, `ISSetCompressOutput()`, `ISView()`

# External Links
$(_doc_external("IS/ISGetCompressOutput"))
"""
function ISGetCompressOutput(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetCompressOutput(petsclib::$UnionPetscLib, is::AbstractIS )
	compress_ = Ref{PetscBool}()

    @chk ccall(
               (:ISGetCompressOutput, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{PetscBool}),
               is, compress_,
              )

	compress = compress_[]

	return compress
end 

"""
	ptr::Vector{PetscInt} = ISGetIndices(petsclib::PetscLibType,is::AbstractIS) 
Returns a pointer to the indices.  The user should call
`ISRestoreIndices()` after having looked at the indices.  The user should
NOT change the indices.

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `ptr` - the location to put the pointer to the indices

Level: intermediate

-seealso: `IS`, `ISRestoreIndices()`

# External Links
$(_doc_external("IS/ISGetIndices"))
"""
function ISGetIndices(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetIndices(petsclib::$UnionPetscLib, is::AbstractIS )
	ptr_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:ISGetIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, ptr_,
              )

	ptr = unsafe_wrap(Array, ptr_[], ISGetLocalSize(petsclib, is); own = false)

	return ptr
end 

"""
	flg::PetscBool = ISGetInfo(petsclib::PetscLibType,is::AbstractIS, info::ISInfo, type::ISInfoType, compute::PetscBool) 
Determine whether an index set satisfies a given property

Collective or Logically Collective if the type is `IS_GLOBAL` (logically collective if the value of the property has been permanently set with `ISSetInfo()`)

Input Parameters:
- `is`      - the index set
- `info`    - describing a property of the index set, one of those listed in the documentation of `ISSetInfo()`
- `compute` - if `PETSC_FALSE`, the property will not be computed if it is not already known and the property will be assumed to be false
- `type`    - whether the property is local (`IS_LOCAL`) or global (`IS_GLOBAL`)

Output Parameter:
- `flg` - whether the property is true (`PETSC_TRUE`) or false (`PETSC_FALSE`)

Level: advanced

-seealso: `IS`, `ISInfo`, `ISInfoType`, `ISSetInfo()`, `ISClearInfoCache()`

# External Links
$(_doc_external("IS/ISGetInfo"))
"""
function ISGetInfo(petsclib::PetscLibType, is::AbstractIS, info::ISInfo, type::ISInfoType, compute::PetscBool) end

@for_petsc function ISGetInfo(petsclib::$UnionPetscLib, is::AbstractIS, info::ISInfo, type::ISInfoType, compute::PetscBool )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:ISGetInfo, $petsc_library),
               PetscErrorCode,
               (CIS, ISInfo, ISInfoType, PetscBool, Ptr{PetscBool}),
               is, info, type, compute, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	map::PetscLayout = ISGetLayout(petsclib::PetscLibType,is::AbstractIS) 
get `PetscLayout` describing index set layout

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `map` - the layout

Level: developer

-seealso: `IS`, `PetscLayout`, `ISSetLayout()`, `ISGetSize()`, `ISGetLocalSize()`

# External Links
$(_doc_external("IS/ISGetLayout"))
"""
function ISGetLayout(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetLayout(petsclib::$UnionPetscLib, is::AbstractIS )
	map_ = Ref{PetscLayout}()

    @chk ccall(
               (:ISGetLayout, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{PetscLayout}),
               is, map_,
              )

	map = map_[]

	return map
end 

"""
	size::PetscInt = ISGetLocalSize(petsclib::PetscLibType,is::AbstractIS) 
Returns the local (processor) length of an index set.

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `size` - the local size

Level: beginner

-seealso: `IS`, `ISGetSize()`

# External Links
$(_doc_external("IS/ISGetLocalSize"))
"""
function ISGetLocalSize(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetLocalSize(petsclib::$UnionPetscLib, is::AbstractIS )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGetLocalSize, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}),
               is, size_,
              )

	size = size_[]

	return size
end 

"""
	min::PetscInt,max::PetscInt = ISGetMinMax(petsclib::PetscLibType,is::AbstractIS) 
Gets the minimum and maximum values in an `IS`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameters:
- `min` - the minimum value, you may pass `NULL`
- `max` - the maximum value, you may pass `NULL`

Level: intermediate

-seealso: `IS`, `ISGetIndices()`, `ISRestoreIndices()`

# External Links
$(_doc_external("IS/ISGetMinMax"))
"""
function ISGetMinMax(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetMinMax(petsclib::$UnionPetscLib, is::AbstractIS )
	min_ = Ref{$PetscInt}()
	max_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGetMinMax, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}, Ptr{$PetscInt}),
               is, min_, max_,
              )

	min = min_[]
	max = max_[]

	return min,max
end 

"""
	complement::IS = ISGetNonlocalIS(petsclib::PetscLibType,is::AbstractIS) 
Gather all nonlocal indices for this `IS` and present
them as another sequential index set.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `complement` - sequential `IS` with indices identical to the result of
`ISGetNonlocalIndices()`

Level: intermediate

-seealso: `IS`, `ISGetNonlocalIndices()`, `ISRestoreNonlocalIndices()`, `ISAllGather()`, `ISGetSize()`

# External Links
$(_doc_external("IS/ISGetNonlocalIS"))
"""
function ISGetNonlocalIS(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetNonlocalIS(petsclib::$UnionPetscLib, is::AbstractIS )
	complement_ = Ref{CIS}()

    @chk ccall(
               (:ISGetNonlocalIS, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{CIS}),
               is, complement_,
              )

	complement = IS(complement_[], petsclib)

	return complement
end 

"""
	indices::Vector{PetscInt} = ISGetNonlocalIndices(petsclib::PetscLibType,is::AbstractIS) 
Retrieve an array of indices from remote processors
in this communicator.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `indices` - indices with rank 0 indices first, and so on,  omitting
the current rank.  Total number of indices is the difference
total and local, obtained with `ISGetSize()` and `ISGetLocalSize()`,
respectively.

Level: intermediate

-seealso: `IS`, `ISGetTotalIndices()`, `ISRestoreNonlocalIndices()`, `ISGetSize()`, `ISGetLocalSize().`

# External Links
$(_doc_external("IS/ISGetNonlocalIndices"))
"""
function ISGetNonlocalIndices(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetNonlocalIndices(petsclib::$UnionPetscLib, is::AbstractIS )
	indices_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:ISGetNonlocalIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, indices_,
              )

	n = ISGetSize(petsclib, is) - ISGetLocalSize(petsclib, is)
	indices = unsafe_wrap(Array, indices_[], n; own = false)

	return indices
end 

"""
	pStart::PetscInt,pEnd::PetscInt,points::Ptr{PetscInt} = ISGetPointRange(petsclib::PetscLibType,pointIS::AbstractIS) 
Returns a description of the points in an `IS` suitable for traversal

Not Collective

Input Parameter:
- `pointIS` - The `IS` object

Output Parameters:
- `pStart` - The first index, see notes
- `pEnd`   - One past the last index, see notes
- `points` - The indices, see notes

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISRestorePointRange()`, `ISGetPointSubrange()`, `ISGetIndices()`, `ISCreateStride()`

# External Links
$(_doc_external("IS/ISGetPointRange"))
"""
function ISGetPointRange(petsclib::PetscLibType, pointIS::AbstractIS) end

@for_petsc function ISGetPointRange(petsclib::$UnionPetscLib, pointIS::AbstractIS )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()
	points_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:ISGetPointRange, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               pointIS, pStart_, pEnd_, points_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]
	points = points_[]

	return pStart,pEnd,points
end 

"""
	ISGetPointSubrange(petsclib::PetscLibType,subpointIS::AbstractIS, pStart::PetscInt, pEnd::PetscInt, points::Vector{PetscInt}) 
Configures the input `IS` to be a subrange for the traversal information given

Not Collective

Input Parameters:
- `subpointIS` - The `IS` object to be configured
- `pStart`     - The first index of the subrange
- `pEnd`       - One past the last index for the subrange
- `points`     - The indices for the entire range, from `ISGetPointRange()`

Output Parameters:
- `subpointIS` - The `IS` object now configured to be a subrange

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGetPointRange()`, `ISRestorePointRange()`, `ISGetIndices()`, `ISCreateStride()`

# External Links
$(_doc_external("IS/ISGetPointSubrange"))
"""
function ISGetPointSubrange(petsclib::PetscLibType, subpointIS::AbstractIS, pStart::PetscInt, pEnd::PetscInt, points::Vector{PetscInt}) end

@for_petsc function ISGetPointSubrange(petsclib::$UnionPetscLib, subpointIS::AbstractIS, pStart::$PetscInt, pEnd::$PetscInt, points::Vector{$PetscInt} )

    @chk ccall(
               (:ISGetPointSubrange, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               subpointIS, pStart, pEnd, points,
              )


	return nothing
end 

"""
	size::PetscInt = ISGetSize(petsclib::PetscLibType,is::AbstractIS) 
Returns the global length of an index set.

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `size` - the global size

Level: beginner

-seealso: `IS`

# External Links
$(_doc_external("IS/ISGetSize"))
"""
function ISGetSize(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetSize(petsclib::$UnionPetscLib, is::AbstractIS )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGetSize, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}),
               is, size_,
              )

	size = size_[]

	return size
end 

"""
	indices::Vector{PetscInt} = ISGetTotalIndices(petsclib::PetscLibType,is::AbstractIS) 
Retrieve an array containing all indices across the communicator.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `indices` - total indices with rank 0 indices first, and so on; total array size is
the same as returned with `ISGetSize()`.

Level: intermediate

-seealso: `IS`, `ISRestoreTotalIndices()`, `ISGetNonlocalIndices()`, `ISGetSize()`

# External Links
$(_doc_external("IS/ISGetTotalIndices"))
"""
function ISGetTotalIndices(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetTotalIndices(petsclib::$UnionPetscLib, is::AbstractIS )
	indices_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:ISGetTotalIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, indices_,
              )

	n = ISGetSize(petsclib, is)
	indices = unsafe_wrap(Array, indices_[], n; own = false)

	return indices
end 

"""
	type::ISType = ISGetType(petsclib::PetscLibType,is::AbstractIS) 
Gets the index set type name, `ISType`, (as a string) from the `IS`.

Not Collective

Input Parameter:
- `is` - The index set

Output Parameter:
- `type` - The index set type name

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISType`, `ISSetType()`, `ISCreate()`

# External Links
$(_doc_external("IS/ISGetType"))
"""
function ISGetType(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISGetType(petsclib::$UnionPetscLib, is::AbstractIS )
	type_ = Ref{ISType}()

    @chk ccall(
               (:ISGetType, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{ISType}),
               is, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	nout::PetscInt = ISGlobalToLocalMappingApply(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::PetscInt, idx::Vector{PetscInt}, idxout::Vector{PetscInt}) 
Provides the local numbering for a list of integers
specified with a global numbering.

Not Collective

Input Parameters:
- `mapping` - mapping between local and global numbering
- `type`    - `IS_GTOLM_MASK` - maps global indices with no local value to -1 in the output list (i.e., mask them)
`IS_GTOLM_DROP` - drops the indices with no local value from the output list
- `n`       - number of global indices to map
- `idx`     - global indices to map

Output Parameters:
- `nout`   - number of indices in output array (if type == `IS_GTOLM_MASK` then nout = n)
- `idxout` - local index of each global index, one must pass in an array long enough
to hold all the indices. You can call `ISGlobalToLocalMappingApply()` with
idxout == NULL to determine the required length (returned in nout)
and then allocate the required space and call `ISGlobalToLocalMappingApply()`
a second time to set the values.

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingApply()`, `ISGlobalToLocalMappingApplyBlock()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingDestroy()`

# External Links
$(_doc_external("IS/ISGlobalToLocalMappingApply"))
"""
function ISGlobalToLocalMappingApply(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::PetscInt, idx::Vector{PetscInt}, idxout::Vector{PetscInt}) end

@for_petsc function ISGlobalToLocalMappingApply(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::$PetscInt, idx::Vector{$PetscInt}, idxout::Vector{$PetscInt} )
	nout_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGlobalToLocalMappingApply, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, ISGlobalToLocalMappingMode, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mapping, type, n, idx, nout_, idxout,
              )

	nout = nout_[]

	return nout
end 

"""
	nout::PetscInt = ISGlobalToLocalMappingApplyBlock(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::PetscInt, idx::Vector{PetscInt}, idxout::Vector{PetscInt}) 
Provides the local block numbering for a list of integers
specified with a block global numbering.

Not Collective

Input Parameters:
- `mapping` - mapping between local and global numbering
- `type`    - `IS_GTOLM_MASK` - maps global indices with no local value to -1 in the output list (i.e., mask them)
`IS_GTOLM_DROP` - drops the indices with no local value from the output list
- `n`       - number of global indices to map
- `idx`     - global indices to map

Output Parameters:
- `nout`   - number of indices in output array (if type == `IS_GTOLM_MASK` then nout = n)
- `idxout` - local index of each global index, one must pass in an array long enough
to hold all the indices. You can call `ISGlobalToLocalMappingApplyBlock()` with
idxout == NULL to determine the required length (returned in nout)
and then allocate the required space and call `ISGlobalToLocalMappingApplyBlock()`
a second time to set the values.

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingApply()`, `ISGlobalToLocalMappingApply()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingDestroy()`

# External Links
$(_doc_external("IS/ISGlobalToLocalMappingApplyBlock"))
"""
function ISGlobalToLocalMappingApplyBlock(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::PetscInt, idx::Vector{PetscInt}, idxout::Vector{PetscInt}) end

@for_petsc function ISGlobalToLocalMappingApplyBlock(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, n::$PetscInt, idx::Vector{$PetscInt}, idxout::Vector{$PetscInt} )
	nout_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISGlobalToLocalMappingApplyBlock, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, ISGlobalToLocalMappingMode, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mapping, type, n, idx, nout_, idxout,
              )

	nout = nout_[]

	return nout
end 

"""
	newis::IS = ISGlobalToLocalMappingApplyIS(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, is::AbstractIS) 
Creates from an `IS` in the global numbering
a new index set using the local numbering defined in an `ISLocalToGlobalMapping`
context.

Not Collective

Input Parameters:
- `mapping` - mapping between local and global numbering
- `type`    - `IS_GTOLM_MASK` - maps global indices with no local value to -1 in the output list (i.e., mask them)
`IS_GTOLM_DROP` - drops the indices with no local value from the output list
- `is`      - index set in global numbering

Output Parameter:
- `newis` - index set in local numbering

Level: advanced

-seealso: [](sec_scatter), `ISGlobalToLocalMapping`, `ISGlobalToLocalMappingApply()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingDestroy()`

# External Links
$(_doc_external("IS/ISGlobalToLocalMappingApplyIS"))
"""
function ISGlobalToLocalMappingApplyIS(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, is::AbstractIS) end

@for_petsc function ISGlobalToLocalMappingApplyIS(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, type::ISGlobalToLocalMappingMode, is::AbstractIS )
	newis_ = Ref{CIS}()

    @chk ccall(
               (:ISGlobalToLocalMappingApplyIS, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, ISGlobalToLocalMappingMode, CIS, Ptr{CIS}),
               mapping, type, is, newis_,
              )

	newis = IS(newis_[], petsclib)

	return newis
end 

"""
	ident::PetscBool = ISIdentity(petsclib::PetscLibType,is::AbstractIS) 
Determines whether index set is the identity mapping.

Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `ident` - `PETSC_TRUE` if an identity, else `PETSC_FALSE`

Level: intermediate

-seealso: `IS`, `ISSetIdentity()`, `ISGetInfo()`

# External Links
$(_doc_external("IS/ISIdentity"))
"""
function ISIdentity(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISIdentity(petsclib::$UnionPetscLib, is::AbstractIS )
	ident_ = Ref{PetscBool}()

    @chk ccall(
               (:ISIdentity, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{PetscBool}),
               is, ident_,
              )

	ident = ident_[]

	return ident
end 

"""
	ISInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `IS` package. It is called
from PetscDLLibraryRegister_petscvec() when using dynamic libraries, and on the first call to ISCreateXXXX()
when using shared or static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/ISInitializePackage"))
"""
function ISInitializePackage(petsclib::PetscLibType) end

@for_petsc function ISInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:ISInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	isout::IS = ISIntersect(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Computes the intersection of two index sets, by sorting and comparing.

Collective

Input Parameters:
- `is1` - first index set
- `is2` - second index set

Output Parameter:
- `isout` - the sorted intersection of `is1` and `is2`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISDestroy()`, `ISView()`, `ISDifference()`, `ISSum()`, `ISExpand()`, `ISConcatenate()`

# External Links
$(_doc_external("IS/ISIntersect"))
"""
function ISIntersect(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISIntersect(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISIntersect, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               is1, is2, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	isout::IS = ISInvertPermutation(petsclib::PetscLibType,is::AbstractIS, nlocal::PetscInt) 
Creates a new permutation that is the inverse of
a given permutation.

Collective

Input Parameters:
- `is`     - the index set
- `nlocal` - number of indices on this processor in result (ignored for 1 processor) or
use `PETSC_DECIDE`

Output Parameter:
- `isout` - the inverse permutation

Level: intermediate

-seealso: `IS`, `ISGetInfo()`, `ISSetPermutation()`, `ISGetPermutation()`

# External Links
$(_doc_external("IS/ISInvertPermutation"))
"""
function ISInvertPermutation(petsclib::PetscLibType, is::AbstractIS, nlocal::PetscInt) end

@for_petsc function ISInvertPermutation(petsclib::$UnionPetscLib, is::AbstractIS, nlocal::$PetscInt )
	isout_ = Ref{CIS}()

    @chk ccall(
               (:ISInvertPermutation, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, Ptr{CIS}),
               is, nlocal, isout_,
              )

	isout = IS(isout_[], petsclib)

	return isout
end 

"""
	xis::IS,yis::IS = ISListToPair(petsclib::PetscLibType,comm::MPI_Comm, listlen::PetscInt, islist::Vector{<:AbstractIS}) 
Convert an `IS` list to a pair of `IS` of equal length defining an equivalent integer multimap.
Each `IS` in `islist` is assigned an integer j so that all of the indices of that `IS` are
mapped to j.

Collective

Input Parameters:
- `comm`    - `MPI_Comm`
- `listlen` - `IS` list length
- `islist`  - `IS` list

Output Parameters:
- `xis` - domain `IS`
- `yis` - range  `IS`

Level: developer

-seealso: `IS`, `ISPairToList()`

# External Links
$(_doc_external("IS/ISListToPair"))
"""
function ISListToPair(petsclib::PetscLibType, comm::MPI_Comm, listlen::PetscInt, islist::Vector{<:AbstractIS}) end

@for_petsc function ISListToPair(petsclib::$UnionPetscLib, comm::MPI_Comm, listlen::$PetscInt, islist::Vector{<:AbstractIS} )
	xis_ = Ref{CIS}()
	yis_ = Ref{CIS}()

    @chk ccall(
               (:ISListToPair, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CIS}, Ptr{CIS}, Ptr{CIS}),
               comm, listlen, islist, xis_, yis_,
              )

	xis = IS(xis_[], petsclib)
	yis = IS(yis_[], petsclib)

	return xis,yis
end 

"""
	ISLoad(petsclib::PetscLibType,is::AbstractIS, viewer::PetscViewer) 
Loads an index set that has been stored in binary or HDF5 format with `ISView()`.

Collective

Input Parameters:
- `is`     - the newly loaded index set, this needs to have been created with `ISCreate()` or some related function before a call to `ISLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()` or HDF5 file viewer, obtained from `PetscViewerHDF5Open()`

Level: intermediate

-seealso: `IS`, `PetscViewerBinaryOpen()`, `ISView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("IS/ISLoad"))
"""
function ISLoad(petsclib::PetscLibType, is::AbstractIS, viewer::PetscViewer) end

@for_petsc function ISLoad(petsclib::$UnionPetscLib, is::AbstractIS, viewer::PetscViewer )

    @chk ccall(
               (:ISLoad, $petsc_library),
               PetscErrorCode,
               (CIS, PetscViewer),
               is, viewer,
              )


	return nothing
end 

"""
	location::PetscInt = ISLocate(petsclib::PetscLibType,is::AbstractIS, key::PetscInt) 
determine the location of an index within the local component of an index set

Not Collective

Input Parameters:
- `is`  - the index set
- `key` - the search key

Output Parameter:
- `location` - if >= 0, a location within the index set that is equal to the key, otherwise the key is not in the index set

Level: intermediate

-seealso: `IS`

# External Links
$(_doc_external("IS/ISLocate"))
"""
function ISLocate(petsclib::PetscLibType, is::AbstractIS, key::PetscInt) end

@for_petsc function ISLocate(petsclib::$UnionPetscLib, is::AbstractIS, key::$PetscInt )
	location_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISLocate, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, Ptr{$PetscInt}),
               is, key, location_,
              )

	location = location_[]

	return location
end 

"""
	newis::IS = ISOnComm(petsclib::PetscLibType,is::AbstractIS, comm::MPI_Comm, mode::PetscCopyMode) 
Split a parallel `IS` on subcomms (usually self) or concatenate index sets on subcomms into a parallel index set

Collective

Input Parameters:
- `is`   - index set
- `comm` - communicator for new index set
- `mode` - copy semantics, `PETSC_USE_POINTER` for no-copy if possible, otherwise `PETSC_COPY_VALUES`

Output Parameter:
- `newis` - new `IS` on `comm`

Level: advanced

-seealso: `IS`

# External Links
$(_doc_external("IS/ISOnComm"))
"""
function ISOnComm(petsclib::PetscLibType, is::AbstractIS, comm::MPI_Comm, mode::PetscCopyMode) end

@for_petsc function ISOnComm(petsclib::$UnionPetscLib, is::AbstractIS, comm::MPI_Comm, mode::PetscCopyMode )
	newis_ = Ref{CIS}()

    @chk ccall(
               (:ISOnComm, $petsc_library),
               PetscErrorCode,
               (CIS, MPI_Comm, PetscCopyMode, Ptr{CIS}),
               is, comm, mode, newis_,
              )

	newis = IS(newis_[], petsclib)

	return newis
end 

"""
	listlen::PetscInt,islist::IS = ISPairToList(petsclib::PetscLibType,xis::AbstractIS, yis::AbstractIS) 
Convert an `IS` pair encoding an integer map to a list of `IS`.

Collective

Input Parameters:
- `xis` - domain `IS`
- `yis` - range `IS`, the maximum value must be less than `PETSC_MPI_INT_MAX`

Output Parameters:
- `listlen` - length of `islist`
- `islist`  - list of `IS`s breaking up indices by color

Level: developer

-seealso: `IS`, `ISListToPair()`

# External Links
$(_doc_external("IS/ISPairToList"))
"""
function ISPairToList(petsclib::PetscLibType, xis::AbstractIS, yis::AbstractIS) end

@for_petsc function ISPairToList(petsclib::$UnionPetscLib, xis::AbstractIS, yis::AbstractIS )
	listlen_ = Ref{$PetscInt}()
	islist_ = Ref{CIS}()

    @chk ccall(
               (:ISPairToList, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
               xis, yis, listlen_, islist_,
              )

	listlen = listlen_[]
	islist = IS(islist_[], petsclib)

	return listlen,islist
end 

"""
	ISPartitioningCount(petsclib::PetscLibType,part::AbstractIS, len::PetscInt, count::Vector{PetscInt}) 
Takes a `IS` that represents a partitioning (the MPI rank that each local entry belongs to) and determines the number of
resulting elements on each (partition) rank

Collective

Input Parameters:
- `part` - a partitioning as generated by `MatPartitioningApply()` or `MatPartitioningApplyND()`
- `len`  - length of the array count, this is the total number of partitions

Output Parameter:
- `count` - array of length size, to contain the number of elements assigned
to each partition, where size is the number of partitions generated
(see notes below).

Level: advanced

-seealso: [](sec_scatter), `IS`, `MatPartitioningCreate()`, `AOCreateBasic()`, `ISPartitioningToNumbering()`,
`MatPartitioningSetNParts()`, `MatPartitioningApply()`, `MatPartitioningApplyND()`

# External Links
$(_doc_external("IS/ISPartitioningCount"))
"""
function ISPartitioningCount(petsclib::PetscLibType, part::AbstractIS, len::PetscInt, count::Vector{PetscInt}) end

@for_petsc function ISPartitioningCount(petsclib::$UnionPetscLib, part::AbstractIS, len::$PetscInt, count::Vector{$PetscInt} )

    @chk ccall(
               (:ISPartitioningCount, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, Ptr{$PetscInt}),
               part, len, count,
              )


	return nothing
end 

"""
	is::IS = ISPartitioningToNumbering(petsclib::PetscLibType,part::AbstractIS) 
Takes an `IS' that represents a partitioning (the MPI rank that each local entry belongs to) and on each MPI process
generates an `IS` that contains a new global node number in the new ordering for each entry

Collective

Input Parameter:
- `part` - a partitioning as generated by `MatPartitioningApply()` or `MatPartitioningApplyND()`

Output Parameter:
- `is` - on each processor the index set that defines the global numbers
(in the new numbering) for all the nodes currently (before the partitioning)
on that processor

Level: advanced

-seealso: [](sec_scatter), `IS`, `MatPartitioningCreate()`, `AOCreateBasic()`, `ISPartitioningCount()`

# External Links
$(_doc_external("IS/ISPartitioningToNumbering"))
"""
function ISPartitioningToNumbering(petsclib::PetscLibType, part::AbstractIS) end

@for_petsc function ISPartitioningToNumbering(petsclib::$UnionPetscLib, part::AbstractIS )
	is_ = Ref{CIS}()

    @chk ccall(
               (:ISPartitioningToNumbering, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{CIS}),
               part, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	perm::PetscBool = ISPermutation(petsclib::PetscLibType,is::AbstractIS) 
`PETSC_TRUE` or `PETSC_FALSE` depending on whether the
index set has been declared to be a permutation.

Logically Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `perm` - `PETSC_TRUE` if a permutation, else `PETSC_FALSE`

Level: intermediate

-seealso: `IS`, `ISSetPermutation()`, `ISGetInfo()`

# External Links
$(_doc_external("IS/ISPermutation"))
"""
function ISPermutation(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISPermutation(petsclib::$UnionPetscLib, is::AbstractIS )
	perm_ = Ref{PetscBool}()

    @chk ccall(
               (:ISPermutation, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{PetscBool}),
               is, perm_,
              )

	perm = perm_[]

	return perm
end 

"""
	ISRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new index set implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine itself

-seealso: [](sec_scatter), `IS`, `ISType`, `ISSetType()`, `ISRegisterAll()`, `ISRegisterDestroy()`

# External Links
$(_doc_external("IS/ISRegister"))
"""
function ISRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function ISRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:ISRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	ISRegisterAll(petsclib::PetscLibType) 
Registers all of the index set components in the `IS` package.

Not Collective

Level: advanced

-seealso: [](sec_scatter), `IS`, `ISType`, `ISRegister()`

# External Links
$(_doc_external("IS/ISRegisterAll"))
"""
function ISRegisterAll(petsclib::PetscLibType) end

@for_petsc function ISRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:ISRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	N::PetscInt,subset_n::IS = ISRenumber(petsclib::PetscLibType,subset::AbstractIS, subset_mult::AbstractIS) 
Renumbers the non

Collective

Input Parameters:
- `subset`      - the index set
- `subset_mult` - the multiplicity of each entry in subset (optional, can be `NULL`)

Output Parameters:
- `N`        - one past the largest entry of the new `IS`
- `subset_n` - the new `IS`

Level: intermediate

-seealso: `IS`

# External Links
$(_doc_external("IS/ISRenumber"))
"""
function ISRenumber(petsclib::PetscLibType, subset::AbstractIS, subset_mult::AbstractIS) end

@for_petsc function ISRenumber(petsclib::$UnionPetscLib, subset::AbstractIS, subset_mult::AbstractIS )
	N_ = Ref{$PetscInt}()
	subset_n_ = Ref{CIS}()

    @chk ccall(
               (:ISRenumber, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{$PetscInt}, Ptr{CIS}),
               subset, subset_mult, N_, subset_n_,
              )

	N = N_[]
	subset_n = IS(subset_n_[], petsclib)

	return N,subset_n
end 

"""
	ISRestoreIndices(petsclib::PetscLibType,is::AbstractIS, ptr::Vector{PetscInt}) 
Restores an index set to a usable state after a call to `ISGetIndices()`.

Not Collective

Input Parameters:
- `is`  - the index set
- `ptr` - the pointer obtained by `ISGetIndices()`

Level: intermediate

-seealso: `IS`, `ISGetIndices()`

# External Links
$(_doc_external("IS/ISRestoreIndices"))
"""
function ISRestoreIndices(petsclib::PetscLibType, is::AbstractIS, ptr::Vector{PetscInt}) end

@for_petsc function ISRestoreIndices(petsclib::$UnionPetscLib, is::AbstractIS, ptr::Vector{$PetscInt} )
	ptr_ = Ref(pointer(ptr))

    @chk ccall(
               (:ISRestoreIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, ptr_,
              )


	return nothing
end 

"""
	ISRestoreNonlocalIS(petsclib::PetscLibType,is::AbstractIS, complement::AbstractIS) 
Restore the `IS` obtained with `ISGetNonlocalIS()`.

Not collective.

Input Parameters:
- `is`         - the index set
- `complement` - index set of `is`'s nonlocal indices

Level: intermediate

-seealso: `IS`, `ISGetNonlocalIS()`, `ISGetNonlocalIndices()`, `ISRestoreNonlocalIndices()`

# External Links
$(_doc_external("IS/ISRestoreNonlocalIS"))
"""
function ISRestoreNonlocalIS(petsclib::PetscLibType, is::AbstractIS, complement::AbstractIS) end

@for_petsc function ISRestoreNonlocalIS(petsclib::$UnionPetscLib, is::AbstractIS, complement::AbstractIS )
	complement_ = Ref(complement.ptr)

    @chk ccall(
               (:ISRestoreNonlocalIS, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{CIS}),
               is, complement_,
              )

	complement.ptr = complement_[]

	return nothing
end 

"""
	ISRestoreNonlocalIndices(petsclib::PetscLibType,is::AbstractIS, indices::Vector{PetscInt}) 
Restore the index array obtained with `ISGetNonlocalIndices()`.

Not Collective.

Input Parameters:
- `is`      - the index set
- `indices` - index array; must be the array obtained with `ISGetNonlocalIndices()`

Level: intermediate

-seealso: `IS`, `ISGetTotalIndices()`, `ISGetNonlocalIndices()`, `ISRestoreTotalIndices()`

# External Links
$(_doc_external("IS/ISRestoreNonlocalIndices"))
"""
function ISRestoreNonlocalIndices(petsclib::PetscLibType, is::AbstractIS, indices::Vector{PetscInt}) end

@for_petsc function ISRestoreNonlocalIndices(petsclib::$UnionPetscLib, is::AbstractIS, indices::Vector{$PetscInt} )
	indices_ = Ref(pointer(indices))

    @chk ccall(
               (:ISRestoreNonlocalIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, indices_,
              )


	return nothing
end 

"""
	ISRestorePointRange(petsclib::PetscLibType,pointIS::AbstractIS, pStart::PetscInt, pEnd::PetscInt, points::Vector{PetscInt}) 
Destroys the traversal description created with `ISGetPointRange()`

Not Collective

Input Parameters:
- `pointIS` - The `IS` object
- `pStart`  - The first index, from `ISGetPointRange()`
- `pEnd`    - One past the last index, from `ISGetPointRange()`
- `points`  - The indices, from `ISGetPointRange()`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGetPointRange()`, `ISGetPointSubrange()`, `ISGetIndices()`, `ISCreateStride()`

# External Links
$(_doc_external("IS/ISRestorePointRange"))
"""
function ISRestorePointRange(petsclib::PetscLibType, pointIS::AbstractIS, pStart::PetscInt, pEnd::PetscInt, points::Vector{PetscInt}) end

@for_petsc function ISRestorePointRange(petsclib::$UnionPetscLib, pointIS::AbstractIS, pStart::$PetscInt, pEnd::$PetscInt, points::Vector{$PetscInt} )
	pStart_ = Ref{$PetscInt}(pStart)
	pEnd_ = Ref{$PetscInt}(pEnd)
	points_ = Ref(pointer(points))

    @chk ccall(
               (:ISRestorePointRange, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               pointIS, pStart_, pEnd_, points_,
              )


	return nothing
end 

"""
	ISRestoreTotalIndices(petsclib::PetscLibType,is::AbstractIS, indices::Vector{PetscInt}) 
Restore the index array obtained with `ISGetTotalIndices()`.

Not Collective.

Input Parameters:
- `is`      - the index set
- `indices` - index array; must be the array obtained with `ISGetTotalIndices()`

Level: intermediate

-seealso: `IS`, `ISGetNonlocalIndices()`

# External Links
$(_doc_external("IS/ISRestoreTotalIndices"))
"""
function ISRestoreTotalIndices(petsclib::PetscLibType, is::AbstractIS, indices::Vector{PetscInt}) end

@for_petsc function ISRestoreTotalIndices(petsclib::$UnionPetscLib, is::AbstractIS, indices::Vector{$PetscInt} )
	indices_ = Ref(pointer(indices))

    @chk ccall(
               (:ISRestoreTotalIndices, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{Ptr{$PetscInt}}),
               is, indices_,
              )


	return nothing
end 

"""
	ISSetBlockSize(petsclib::PetscLibType,is::AbstractIS, bs::PetscInt) 
informs an index set that it has a given block size

Logicall Collective

Input Parameters:
- `is` - index set
- `bs` - block size

Level: intermediate

-seealso: `IS`, `ISGetBlockSize()`, `ISCreateBlock()`, `ISBlockGetIndices()`,

# External Links
$(_doc_external("IS/ISSetBlockSize"))
"""
function ISSetBlockSize(petsclib::PetscLibType, is::AbstractIS, bs::PetscInt) end

@for_petsc function ISSetBlockSize(petsclib::$UnionPetscLib, is::AbstractIS, bs::$PetscInt )

    @chk ccall(
               (:ISSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt),
               is, bs,
              )


	return nothing
end 

"""
	ISSetCompressOutput(petsclib::PetscLibType,is::AbstractIS, compress::PetscBool) 
set the flag for output compression

Logicall Collective

Input Parameters:
- `is`       - index set
- `compress` - flag for output compression

Level: intermediate

-seealso: `IS`, `ISGetCompressOutput()`, `ISView()`

# External Links
$(_doc_external("IS/ISSetCompressOutput"))
"""
function ISSetCompressOutput(petsclib::PetscLibType, is::AbstractIS, compress::PetscBool) end

@for_petsc function ISSetCompressOutput(petsclib::$UnionPetscLib, is::AbstractIS, compress::PetscBool )

    @chk ccall(
               (:ISSetCompressOutput, $petsc_library),
               PetscErrorCode,
               (CIS, PetscBool),
               is, compress,
              )


	return nothing
end 

"""
	ISSetIdentity(petsclib::PetscLibType,is::AbstractIS) 
Informs the index set that it is an identity.

Logically Collective

Input Parameter:
- `is` - the index set

Level: intermediate

-seealso: `IS`, `ISIdentity()`, `ISSetInfo()`, `ISClearInfoCache()`

# External Links
$(_doc_external("IS/ISSetIdentity"))
"""
function ISSetIdentity(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISSetIdentity(petsclib::$UnionPetscLib, is::AbstractIS )

    @chk ccall(
               (:ISSetIdentity, $petsc_library),
               PetscErrorCode,
               (CIS,),
               is,
              )


	return nothing
end 

"""
	ISSetInfo(petsclib::PetscLibType,is::AbstractIS, info::ISInfo, type::ISInfoType, permanent::PetscBool, flg::PetscBool) 
Set known information about an index set.

Logically Collective if `ISInfoType` is `IS_GLOBAL`

Input Parameters:
- `is`        - the index set
- `info`      - describing a property of the index set, one of those listed below,
- `type`      - `IS_LOCAL` if the information describes the local portion of the index set,
`IS_GLOBAL` if it describes the whole index set
- `permanent` - `PETSC_TRUE` if it is known that the property will persist through changes to the index set, `PETSC_FALSE` otherwise
If the user sets a property as permanently known, it will bypass computation of that property
- `flg`       - set the described property as true (`PETSC_TRUE`) or false (`PETSC_FALSE`)

Values of `info` Describing `IS` Structure:
- `IS_SORTED`      - the [local part of the] index set is sorted in ascending order
- `IS_UNIQUE`      - each entry in the [local part of the] index set is unique
- `IS_PERMUTATION` - the [local part of the] index set is a permutation of the integers {0, 1, ..., N-1}, where N is the size of the [local part of the] index set
- `IS_INTERVAL`    - the [local part of the] index set is equal to a contiguous range of integers {f, f + 1, ..., f + N-1}
- `IS_IDENTITY`    - the [local part of the] index set is equal to the integers {0, 1, ..., N-1}

Level: advanced

-seealso: `ISInfo`, `ISInfoType`, `IS`

# External Links
$(_doc_external("IS/ISSetInfo"))
"""
function ISSetInfo(petsclib::PetscLibType, is::AbstractIS, info::ISInfo, type::ISInfoType, permanent::PetscBool, flg::PetscBool) end

@for_petsc function ISSetInfo(petsclib::$UnionPetscLib, is::AbstractIS, info::ISInfo, type::ISInfoType, permanent::PetscBool, flg::PetscBool )

    @chk ccall(
               (:ISSetInfo, $petsc_library),
               PetscErrorCode,
               (CIS, ISInfo, ISInfoType, PetscBool, PetscBool),
               is, info, type, permanent, flg,
              )


	return nothing
end 

"""
	ISSetLayout(petsclib::PetscLibType,is::AbstractIS, map::PetscLayout) 
set `PetscLayout` describing index set layout

Collective

Input Parameters:
- `is`  - the index set
- `map` - the layout

Level: developer

-seealso: `IS`, `PetscLayout`, `ISCreate()`, `ISGetLayout()`, `ISGetSize()`, `ISGetLocalSize()`

# External Links
$(_doc_external("IS/ISSetLayout"))
"""
function ISSetLayout(petsclib::PetscLibType, is::AbstractIS, map::PetscLayout) end

@for_petsc function ISSetLayout(petsclib::$UnionPetscLib, is::AbstractIS, map::PetscLayout )

    @chk ccall(
               (:ISSetLayout, $petsc_library),
               PetscErrorCode,
               (CIS, PetscLayout),
               is, map,
              )


	return nothing
end 

"""
	ISSetPermutation(petsclib::PetscLibType,is::AbstractIS) 
Informs the index set that it is a permutation.

Logically Collective

Input Parameter:
- `is` - the index set

Level: intermediate

-seealso: `IS`, `ISPermutation()`, `ISSetInfo()`, `ISClearInfoCache().`

# External Links
$(_doc_external("IS/ISSetPermutation"))
"""
function ISSetPermutation(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISSetPermutation(petsclib::$UnionPetscLib, is::AbstractIS )

    @chk ccall(
               (:ISSetPermutation, $petsc_library),
               PetscErrorCode,
               (CIS,),
               is,
              )


	return nothing
end 

"""
	ISSetType(petsclib::PetscLibType,is::AbstractIS, method::ISType) 
Builds a index set, for a particular `ISType`

Collective

Input Parameters:
- `is`     - The index set object
- `method` - The name of the index set type

Options Database Key:
- `-is_type <type>` - Sets the index set type; use `-help` for a list of available types

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISGENERAL`, `ISBLOCK`, `ISGetType()`, `ISCreate()`, `ISCreateGeneral()`, `ISCreateStride()`, `ISCreateBlock()`

# External Links
$(_doc_external("IS/ISSetType"))
"""
function ISSetType(petsclib::PetscLibType, is::AbstractIS, method::ISType) end

@for_petsc function ISSetType(petsclib::$UnionPetscLib, is::AbstractIS, method::ISType )

    @chk ccall(
               (:ISSetType, $petsc_library),
               PetscErrorCode,
               (CIS, ISType),
               is, method,
              )


	return nothing
end 

"""
	ISShift(petsclib::PetscLibType,is::AbstractIS, offset::PetscInt, isy::AbstractIS) 
Shift all indices by given offset

Collective

Input Parameters:
- `is`     - the index set
- `offset` - the offset

Output Parameter:
- `isy` - the shifted copy of the input index set

Level: beginner

-seealso: `ISDuplicate()`, `ISCopy()`

# External Links
$(_doc_external("IS/ISShift"))
"""
function ISShift(petsclib::PetscLibType, is::AbstractIS, offset::PetscInt, isy::AbstractIS) end

@for_petsc function ISShift(petsclib::$UnionPetscLib, is::AbstractIS, offset::$PetscInt, isy::AbstractIS )

    @chk ccall(
               (:ISShift, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, CIS),
               is, offset, isy,
              )


	return nothing
end 

"""
	ISSort(petsclib::PetscLibType,is::AbstractIS) 
Sorts the indices of an index set.

Collective

Input Parameter:
- `is` - the index set

Level: intermediate

-seealso: `IS`, `ISSortRemoveDups()`, `ISSorted()`

# External Links
$(_doc_external("IS/ISSort"))
"""
function ISSort(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISSort(petsclib::$UnionPetscLib, is::AbstractIS )

    @chk ccall(
               (:ISSort, $petsc_library),
               PetscErrorCode,
               (CIS,),
               is,
              )


	return nothing
end 

"""
	h::IS = ISSortPermutation(petsclib::PetscLibType,f::AbstractIS, always::PetscBool) 
calculate the permutation of the indices into a nondecreasing order.

Not Collective

Input Parameters:
- `f`      - `IS` to sort
- `always` - build the permutation even when `f`'s indices are nondecreasing.

Output Parameter:
- `h` - permutation or `NULL`, if `f` is nondecreasing and `always` == `PETSC_FALSE`.

Level: advanced

-seealso: `IS`, `ISLocalToGlobalMapping`, `ISSort()`

# External Links
$(_doc_external("IS/ISSortPermutation"))
"""
function ISSortPermutation(petsclib::PetscLibType, f::AbstractIS, always::PetscBool) end

@for_petsc function ISSortPermutation(petsclib::$UnionPetscLib, f::AbstractIS, always::PetscBool )
	h_ = Ref{CIS}()

    @chk ccall(
               (:ISSortPermutation, $petsc_library),
               PetscErrorCode,
               (CIS, PetscBool, Ptr{CIS}),
               f, always, h_,
              )

	h = IS(h_[], petsclib)

	return h
end 

"""
	ISSortRemoveDups(petsclib::PetscLibType,is::AbstractIS) 
Sorts the indices of an index set, removing duplicates.

Collective

Input Parameter:
- `is` - the index set

Level: intermediate

-seealso: `IS`, `ISSort()`, `ISSorted()`

# External Links
$(_doc_external("IS/ISSortRemoveDups"))
"""
function ISSortRemoveDups(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISSortRemoveDups(petsclib::$UnionPetscLib, is::AbstractIS )

    @chk ccall(
               (:ISSortRemoveDups, $petsc_library),
               PetscErrorCode,
               (CIS,),
               is,
              )


	return nothing
end 

"""
	flg::PetscBool = ISSorted(petsclib::PetscLibType,is::AbstractIS) 
Checks the indices to determine whether they have been sorted.

Not Collective

Input Parameter:
- `is` - the index set

Output Parameter:
- `flg` - output flag, either `PETSC_TRUE` if the index set is sorted,
or `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `ISSort()`, `ISSortRemoveDups()`

# External Links
$(_doc_external("IS/ISSorted"))
"""
function ISSorted(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISSorted(petsclib::$UnionPetscLib, is::AbstractIS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:ISSorted, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{PetscBool}),
               is, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	first::PetscInt,step::PetscInt = ISStrideGetInfo(petsclib::PetscLibType,is::AbstractIS) 
Returns the first index in a stride index set and the stride width from an `IS` of `ISType` `ISSTRIDE`

Not Collective

Input Parameter:
- `is` - the index set

Output Parameters:
- `first` - the first index
- `step`  - the stride width

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISCreateStride()`, `ISGetSize()`, `ISSTRIDE`

# External Links
$(_doc_external("IS/ISStrideGetInfo"))
"""
function ISStrideGetInfo(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISStrideGetInfo(petsclib::$UnionPetscLib, is::AbstractIS )
	first_ = Ref{$PetscInt}()
	step_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISStrideGetInfo, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{$PetscInt}, Ptr{$PetscInt}),
               is, first_, step_,
              )

	first = first_[]
	step = step_[]

	return first,step
end 

"""
	ISStrideSetStride(petsclib::PetscLibType,is::AbstractIS, n::PetscInt, first::PetscInt, step::PetscInt) 
Sets the stride information for a stride index set.

Logically Collective

Input Parameters:
- `is`    - the index set
- `n`     - the length of the locally owned portion of the index set
- `first` - the first element of the locally owned portion of the index set
- `step`  - the change to the next index

Level: beginner

-seealso: [](sec_scatter), `IS`, `ISCreateGeneral()`, `ISCreateBlock()`, `ISAllGather()`, `ISSTRIDE`, `ISCreateStride()`, `ISStrideGetInfo()`

# External Links
$(_doc_external("IS/ISStrideSetStride"))
"""
function ISStrideSetStride(petsclib::PetscLibType, is::AbstractIS, n::PetscInt, first::PetscInt, step::PetscInt) end

@for_petsc function ISStrideSetStride(petsclib::$UnionPetscLib, is::AbstractIS, n::$PetscInt, first::$PetscInt, step::$PetscInt )

    @chk ccall(
               (:ISStrideSetStride, $petsc_library),
               PetscErrorCode,
               (CIS, $PetscInt, $PetscInt, $PetscInt),
               is, n, first, step,
              )


	return nothing
end 

"""
	is3::IS = ISSum(petsclib::PetscLibType,is1::AbstractIS, is2::AbstractIS) 
Computes the sum (union) of two index sets.

Only sequential version (at the moment)

Input Parameters:
- `is1` - index set to be extended
- `is2` - index values to be added

Output Parameter:
- `is3` - the sum; this can not be `is1` or `is2`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `ISDestroy()`, `ISView()`, `ISDifference()`, `ISExpand()`

# External Links
$(_doc_external("IS/ISSum"))
"""
function ISSum(petsclib::PetscLibType, is1::AbstractIS, is2::AbstractIS) end

@for_petsc function ISSum(petsclib::$UnionPetscLib, is1::AbstractIS, is2::AbstractIS )
	is3_ = Ref{CIS}()

    @chk ccall(
               (:ISSum, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CIS}),
               is1, is2, is3_,
              )

	is3 = IS(is3_[], petsclib)

	return is3
end 

"""
	ISToGeneral(petsclib::PetscLibType,is::AbstractIS) 
Converts an IS object of any type to `ISGENERAL` type

Collective

Input Parameter:
- `is` - the index set

Level: intermediate

-seealso: `IS`, `ISSorted()`

# External Links
$(_doc_external("IS/ISToGeneral"))
"""
function ISToGeneral(petsclib::PetscLibType, is::AbstractIS) end

@for_petsc function ISToGeneral(petsclib::$UnionPetscLib, is::AbstractIS )

    @chk ccall(
               (:ISToGeneral, $petsc_library),
               PetscErrorCode,
               (CIS,),
               is,
              )


	return nothing
end 

"""
	ISView(petsclib::PetscLibType,is::AbstractIS, viewer::PetscViewer) 
Displays an index set.

Collective

Input Parameters:
- `is`     - the index set
- `viewer` - viewer used to display the set, for example `PETSC_VIEWER_STDOUT_SELF`.

Level: intermediate

-seealso: `IS`, `PetscViewer`, `PetscViewerASCIIOpen()`, `ISViewFromOptions()`

# External Links
$(_doc_external("IS/ISView"))
"""
function ISView(petsclib::PetscLibType, is::AbstractIS, viewer::PetscViewer) end

@for_petsc function ISView(petsclib::$UnionPetscLib, is::AbstractIS, viewer::PetscViewer )

    @chk ccall(
               (:ISView, $petsc_library),
               PetscErrorCode,
               (CIS, PetscViewer),
               is, viewer,
              )


	return nothing
end 

"""
	ISViewFromOptions(petsclib::PetscLibType,A::AbstractIS, obj::PetscObject, name::String) 
View an `IS` based on options in the options database

Collective

Input Parameters:
- `A`    - the index set
- `obj`  - Optional object that provides the prefix for the options database
- `name` - command line option

Level: intermediate

-seealso: `IS`, `ISView()`, `PetscObjectViewFromOptions()`, `ISCreate()`

# External Links
$(_doc_external("IS/ISViewFromOptions"))
"""
function ISViewFromOptions(petsclib::PetscLibType, A::AbstractIS, obj::PetscObject, name::String) end

@for_petsc function ISViewFromOptions(petsclib::$UnionPetscLib, A::AbstractIS, obj::PetscObject, name::String )

    @chk ccall(
               (:ISViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CIS, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

