# autodefined type arguments for class ------
# -------------------------------------------------------
"""
	map::PetscLayout = PetscLayoutCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Allocates `PetscLayout` object

Collective

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `map` - the new `PetscLayout`

Level: advanced

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutSetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutGetLocalSize()`,
`PetscLayout`, `PetscLayoutDestroy()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`, `PetscLayoutSetUp()`,
`PetscLayoutCreateFromSizes()`

# External Links
$(_doc_external("Vec/PetscLayoutCreate"))
"""
function PetscLayoutCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscLayoutCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	map_ = Ref{PetscLayout}()

    @chk ccall(
               (:PetscLayoutCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscLayout}),
               comm, map_,
              )

	map = map_[]

	return map
end 

"""
	map::PetscLayout = PetscLayoutCreateFromSizes(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt, bs::PetscInt) 
Allocates `PetscLayout` object and sets the layout sizes, and sets the layout up.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `n`    - the local size (or `PETSC_DECIDE`)
- `N`    - the global size (or `PETSC_DECIDE`)
- `bs`   - the block size (or `PETSC_DECIDE`)

Output Parameter:
- `map` - the new `PetscLayout`

Level: advanced

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutGetLocalSize()`, `PetscLayout`, `PetscLayoutDestroy()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`, `PetscLayoutSetUp()`, `PetscLayoutCreateFromRanges()`

# External Links
$(_doc_external("Vec/PetscLayoutCreateFromSizes"))
"""
function PetscLayoutCreateFromSizes(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt, bs::PetscInt) end

@for_petsc function PetscLayoutCreateFromSizes(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt, bs::$PetscInt )
	map_ = Ref{PetscLayout}()

    @chk ccall(
               (:PetscLayoutCreateFromSizes, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{PetscLayout}),
               comm, n, N, bs, map_,
              )

	map = map_[]

	return map
end 

"""
	PetscLayoutDestroy(petsclib::PetscLibType,map::PetscLayout) 
Frees a `PetscLayout` object and frees its range if that exists.

Collective

Input Parameter:
- `map` - the `PetscLayout`

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutSetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutGetLocalSize()`,
`PetscLayout`, `PetscLayoutCreate()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`, `PetscLayoutSetUp()`

# External Links
$(_doc_external("Vec/PetscLayoutDestroy"))
"""
function PetscLayoutDestroy(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutDestroy(petsclib::$UnionPetscLib, map::PetscLayout )

    @chk ccall(
               (:PetscLayoutDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLayout},),
               map,
              )


	return nothing
end 

"""
	newmap::PetscLayout = PetscLayoutCreateFromRanges(petsclib::PetscLibType,comm::MPI_Comm, range::Vector{PetscInt}, mode::PetscCopyMode, bs::PetscInt) 
Creates a new `PetscLayout` with the given ownership ranges and sets it up.

Collective

Input Parameters:
- `comm`  - the MPI communicator
- `range` - the array of ownership ranges for each rank with length commsize+1
- `mode`  - the copy mode for range
- `bs`    - the block size (or `PETSC_DECIDE`)

Output Parameter:
- `newmap` - the new `PetscLayout`

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`,
`PetscLayoutGetLocalSize()`, `PetscLayout`, `PetscLayoutDestroy()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`, `PetscLayoutSetUp()`, `PetscLayoutCreateFromSizes()`

# External Links
$(_doc_external("Vec/PetscLayoutCreateFromRanges"))
"""
function PetscLayoutCreateFromRanges(petsclib::PetscLibType, comm::MPI_Comm, range::Vector{PetscInt}, mode::PetscCopyMode, bs::PetscInt) end

@for_petsc function PetscLayoutCreateFromRanges(petsclib::$UnionPetscLib, comm::MPI_Comm, range::Vector{$PetscInt}, mode::PetscCopyMode, bs::$PetscInt )
	newmap_ = Ref{PetscLayout}()

    @chk ccall(
               (:PetscLayoutCreateFromRanges, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscInt}, PetscCopyMode, $PetscInt, Ptr{PetscLayout}),
               comm, range, mode, bs, newmap_,
              )

	newmap = newmap_[]

	return newmap
end 

"""
	PetscLayoutSetUp(petsclib::PetscLibType,map::PetscLayout) 
given a map where you have set either the global or local
size sets up the map so that it may be used.

Collective

Input Parameter:
- `map` - pointer to the map

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutSetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutGetLocalSize()`,
`PetscLayout`, `PetscLayoutDestroy()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`, `PetscLayoutCreate()`, `PetscSplitOwnership()`

# External Links
$(_doc_external("Vec/PetscLayoutSetUp"))
"""
function PetscLayoutSetUp(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutSetUp(petsclib::$UnionPetscLib, map::PetscLayout )

    @chk ccall(
               (:PetscLayoutSetUp, $petsc_library),
               PetscErrorCode,
               (PetscLayout,),
               map,
              )


	return nothing
end 

"""
	out::PetscLayout = PetscLayoutDuplicate(petsclib::PetscLibType,in::PetscLayout) 
creates a new `PetscLayout` with the same information as a given one. If the `PetscLayout` already exists it is destroyed first.

Collective

Input Parameter:
- `in` - input `PetscLayout` to be duplicated

Output Parameter:
- `out` - the copy

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutDestroy()`, `PetscLayoutSetUp()`, `PetscLayoutReference()`

# External Links
$(_doc_external("Vec/PetscLayoutDuplicate"))
"""
function PetscLayoutDuplicate(petsclib::PetscLibType, in::PetscLayout) end

@for_petsc function PetscLayoutDuplicate(petsclib::$UnionPetscLib, in::PetscLayout )
	out_ = Ref{PetscLayout}()

    @chk ccall(
               (:PetscLayoutDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{PetscLayout}),
               in, out_,
              )

	out = out_[]

	return out
end 

"""
	PetscLayoutReference(petsclib::PetscLibType,in::PetscLayout, out::PetscLayout) 
Causes a PETSc `Vec` or `Mat` to share a `PetscLayout` with one that already exists.

Collective

Input Parameter:
- `in` - input `PetscLayout` to be copied

Output Parameter:
- `out` - the reference location

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutDestroy()`, `PetscLayoutSetUp()`, `PetscLayoutDuplicate()`

# External Links
$(_doc_external("Vec/PetscLayoutReference"))
"""
function PetscLayoutReference(petsclib::PetscLibType, in::PetscLayout, out::PetscLayout) end

@for_petsc function PetscLayoutReference(petsclib::$UnionPetscLib, in::PetscLayout, out::PetscLayout )

    @chk ccall(
               (:PetscLayoutReference, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{PetscLayout}),
               in, out,
              )


	return nothing
end 

"""
	PetscLayoutSetISLocalToGlobalMapping(petsclib::PetscLibType,in::PetscLayout, ltog::ISLocalToGlobalMapping) 
sets a `ISLocalGlobalMapping` into a `PetscLayout`

Collective

Input Parameters:
- `in`   - input `PetscLayout`
- `ltog` - the local to global mapping

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutDestroy()`, `PetscLayoutSetUp()`, `PetscLayoutDuplicate()`

# External Links
$(_doc_external("Vec/PetscLayoutSetISLocalToGlobalMapping"))
"""
function PetscLayoutSetISLocalToGlobalMapping(petsclib::PetscLibType, in::PetscLayout, ltog::ISLocalToGlobalMapping) end

@for_petsc function PetscLayoutSetISLocalToGlobalMapping(petsclib::$UnionPetscLib, in::PetscLayout, ltog::ISLocalToGlobalMapping )

    @chk ccall(
               (:PetscLayoutSetISLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (PetscLayout, ISLocalToGlobalMapping),
               in, ltog,
              )


	return nothing
end 

"""
	PetscLayoutSetLocalSize(petsclib::PetscLibType,map::PetscLayout, n::PetscInt) 
Sets the local size for a `PetscLayout` object.

Collective

Input Parameters:
- `map` - pointer to the map
- `n`   - the local size, pass `PETSC_DECIDE` (the default) to have this value determined by the global size set with `PetscLayoutSetSize()`

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutSetUp()`
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`

# External Links
$(_doc_external("Vec/PetscLayoutSetLocalSize"))
"""
function PetscLayoutSetLocalSize(petsclib::PetscLibType, map::PetscLayout, n::PetscInt) end

@for_petsc function PetscLayoutSetLocalSize(petsclib::$UnionPetscLib, map::PetscLayout, n::$PetscInt )

    @chk ccall(
               (:PetscLayoutSetLocalSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt),
               map, n,
              )


	return nothing
end 

"""
	n::PetscInt = PetscLayoutGetLocalSize(petsclib::PetscLibType,map::PetscLayout) 
Gets the local size for a `PetscLayout` object.

Not Collective

Input Parameter:
- `map` - pointer to the map

Output Parameter:
- `n` - the local size

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutSetUp()`
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`

# External Links
$(_doc_external("Vec/PetscLayoutGetLocalSize"))
"""
function PetscLayoutGetLocalSize(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutGetLocalSize(petsclib::$UnionPetscLib, map::PetscLayout )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLayoutGetLocalSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{$PetscInt}),
               map, n_,
              )

	n = n_[]

	return n
end 

"""
	PetscLayoutSetSize(petsclib::PetscLibType,map::PetscLayout, n::PetscInt) 
Sets the global size for a `PetscLayout` object.

Logically Collective

Input Parameters:
- `map` - pointer to the map
- `n`   - the global size, use `PETSC_DETERMINE` (the default) to have this value computed as the sum of the local sizes set with `PetscLayoutSetLocalSize()`

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutGetSize()`, `PetscLayoutSetUp()`
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`

# External Links
$(_doc_external("Vec/PetscLayoutSetSize"))
"""
function PetscLayoutSetSize(petsclib::PetscLibType, map::PetscLayout, n::PetscInt) end

@for_petsc function PetscLayoutSetSize(petsclib::$UnionPetscLib, map::PetscLayout, n::$PetscInt )

    @chk ccall(
               (:PetscLayoutSetSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt),
               map, n,
              )


	return nothing
end 

"""
	n::PetscInt = PetscLayoutGetSize(petsclib::PetscLibType,map::PetscLayout) 
Gets the global size for a `PetscLayout` object.

Not Collective

Input Parameter:
- `map` - pointer to the map

Output Parameter:
- `n` - the global size

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutSetUp()`
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetBlockSize()`

# External Links
$(_doc_external("Vec/PetscLayoutGetSize"))
"""
function PetscLayoutGetSize(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutGetSize(petsclib::$UnionPetscLib, map::PetscLayout )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLayoutGetSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{$PetscInt}),
               map, n_,
              )

	n = n_[]

	return n
end 

"""
	PetscLayoutSetBlockSize(petsclib::PetscLibType,map::PetscLayout, bs::PetscInt) 
Sets the block size for a `PetscLayout` object.

Logically Collective

Input Parameters:
- `map` - pointer to the map
- `bs`  - the size

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutGetBlockSize()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutSetUp()`

# External Links
$(_doc_external("Vec/PetscLayoutSetBlockSize"))
"""
function PetscLayoutSetBlockSize(petsclib::PetscLibType, map::PetscLayout, bs::PetscInt) end

@for_petsc function PetscLayoutSetBlockSize(petsclib::$UnionPetscLib, map::PetscLayout, bs::$PetscInt )

    @chk ccall(
               (:PetscLayoutSetBlockSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt),
               map, bs,
              )


	return nothing
end 

"""
	bs::PetscInt = PetscLayoutGetBlockSize(petsclib::PetscLibType,map::PetscLayout) 
Gets the block size for a `PetscLayout` object.

Not Collective

Input Parameter:
- `map` - pointer to the map

Output Parameter:
- `bs` - the size

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutSetSize()`, `PetscLayoutSetUp()`
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutGetSize()`

# External Links
$(_doc_external("Vec/PetscLayoutGetBlockSize"))
"""
function PetscLayoutGetBlockSize(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutGetBlockSize(petsclib::$UnionPetscLib, map::PetscLayout )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLayoutGetBlockSize, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{$PetscInt}),
               map, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	rstart::PetscInt,rend_::PetscInt = PetscLayoutGetRange(petsclib::PetscLibType,map::PetscLayout) 
gets the range of values owned by this process

Not Collective

Input Parameter:
- `map` - pointer to the map

Output Parameters:
- `rstart` - first index owned by this process
- `rend`   - one more than the last index owned by this process

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutSetSize()`,
`PetscLayoutGetSize()`, `PetscLayoutGetRanges()`, `PetscLayoutSetBlockSize()`, `PetscLayoutSetUp()`

# External Links
$(_doc_external("Vec/PetscLayoutGetRange"))
"""
function PetscLayoutGetRange(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutGetRange(petsclib::$UnionPetscLib, map::PetscLayout )
	rstart_ = Ref{$PetscInt}()
	rend__ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLayoutGetRange, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{$PetscInt}, Ptr{$PetscInt}),
               map, rstart_, rend__,
              )

	rstart = rstart_[]
	rend_ = rend__[]

	return rstart,rend_
end 

"""
	range::Vector{PetscInt} = PetscLayoutGetRanges(petsclib::PetscLibType,map::PetscLayout) 
gets the ranges of values owned by all processes

Not Collective

Input Parameter:
- `map` - pointer to the map

Output Parameter:
- `range` - start of each processors range of indices (the final entry is one more than the
last index on the last process). The length of the array is one more than the number of processes in the MPI
communicator owned by `map`

Level: developer

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutSetSize()`,
`PetscLayoutGetSize()`, `PetscLayoutGetRange()`, `PetscLayoutSetBlockSize()`, `PetscLayoutSetUp()`

# External Links
$(_doc_external("Vec/PetscLayoutGetRanges"))
"""
function PetscLayoutGetRanges(petsclib::PetscLibType, map::PetscLayout) end

@for_petsc function PetscLayoutGetRanges(petsclib::$UnionPetscLib, map::PetscLayout )
	range_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscLayoutGetRanges, $petsc_library),
               PetscErrorCode,
               (PetscLayout, Ptr{Ptr{$PetscInt}}),
               map, range_,
              )

	range = unsafe_wrap(Array, range_[], VecGetLocalSize(petsclib, x); own = false)

	return range
end 

"""
	congruent::PetscBool = PetscLayoutCompare(petsclib::PetscLibType,mapa::PetscLayout, mapb::PetscLayout) 
Compares two layouts

Not Collective

Input Parameters:
- `mapa` - pointer to the first map
- `mapb` - pointer to the second map

Output Parameter:
- `congruent` - `PETSC_TRUE` if the two layouts are congruent, `PETSC_FALSE` otherwise

Level: beginner

-seealso: [PetscLayout](sec_matlayout), `PetscLayoutCreate()`, `PetscLayoutSetLocalSize()`, `PetscLayoutGetLocalSize()`, `PetscLayoutGetBlockSize()`,
`PetscLayoutGetRange()`, `PetscLayoutGetRanges()`, `PetscLayoutSetSize()`, `PetscLayoutGetSize()`, `PetscLayoutSetUp()`

# External Links
$(_doc_external("Vec/PetscLayoutCompare"))
"""
function PetscLayoutCompare(petsclib::PetscLibType, mapa::PetscLayout, mapb::PetscLayout) end

@for_petsc function PetscLayoutCompare(petsclib::$UnionPetscLib, mapa::PetscLayout, mapb::PetscLayout )
	congruent_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLayoutCompare, $petsc_library),
               PetscErrorCode,
               (PetscLayout, PetscLayout, Ptr{PetscBool}),
               mapa, mapb, congruent_,
              )

	congruent = congruent_[]

	return congruent
end 

"""
	PetscLayoutFindOwner(petsclib::PetscLibType,map::PetscLayout, idx::PetscInt, owner::PetscMPIInt) 
Find the owning MPI process for a global index

Not Collective; No Fortran Support

Input Parameters:
- `map` - the layout
- `idx` - global index to find the owner of

Output Parameter:
- `owner` - the owning rank

Level: developer

-seealso: `PetscLayout`, `PetscLayoutFindOwnerIndex()`

# External Links
$(_doc_external("Vec/PetscLayoutFindOwner"))
"""
function PetscLayoutFindOwner(petsclib::PetscLibType, map::PetscLayout, idx::PetscInt, owner::PetscMPIInt) end

@for_petsc function PetscLayoutFindOwner(petsclib::$UnionPetscLib, map::PetscLayout, idx::$PetscInt, owner::PetscMPIInt )

    @chk ccall(
               (:PetscLayoutFindOwner, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt, Ptr{PetscMPIInt}),
               map, idx, owner,
              )


	return nothing
end 

"""
	lidx::PetscInt = PetscLayoutFindOwnerIndex(petsclib::PetscLibType,map::PetscLayout, idx::PetscInt, owner::PetscMPIInt) 
Find the owning MPI process and the local index on that process for a global index

Not Collective; No Fortran Support

Input Parameters:
- `map` - the layout
- `idx` - global index to find the owner of

Output Parameters:
- `owner` - the owning rank
- `lidx`  - local index used by the owner for `idx`

Level: developer

-seealso: `PetscLayout`, `PetscLayoutFindOwner()`

# External Links
$(_doc_external("Vec/PetscLayoutFindOwnerIndex"))
"""
function PetscLayoutFindOwnerIndex(petsclib::PetscLibType, map::PetscLayout, idx::PetscInt, owner::PetscMPIInt) end

@for_petsc function PetscLayoutFindOwnerIndex(petsclib::$UnionPetscLib, map::PetscLayout, idx::$PetscInt, owner::PetscMPIInt )
	lidx_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLayoutFindOwnerIndex, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
               map, idx, owner, lidx_,
              )

	lidx = lidx_[]

	return lidx
end 

"""
	on::PetscInt,oidxs::Vector{PetscInt},ogidxs::Vector{PetscInt} = PetscLayoutMapLocal(petsclib::PetscLibType,map::PetscLayout, N::PetscInt, idxs::Vector{PetscInt}) 

# External Links
$(_doc_external("Vec/PetscLayoutMapLocal"))
"""
function PetscLayoutMapLocal(petsclib::PetscLibType, map::PetscLayout, N::PetscInt, idxs::Vector{PetscInt}) end

@for_petsc function PetscLayoutMapLocal(petsclib::$UnionPetscLib, map::PetscLayout, N::$PetscInt, idxs::Vector{$PetscInt} )
	on_ = Ref{$PetscInt}()
	oidxs_ = Ref{Ptr{$PetscInt}}()
	ogidxs_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscLayoutMapLocal, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               map, N, idxs, on_, oidxs_, ogidxs_,
              )

	on = on_[]
	oidxs = unsafe_wrap(Array, oidxs_[], VecGetLocalSize(petsclib, x); own = false)
	ogidxs = unsafe_wrap(Array, ogidxs_[], VecGetLocalSize(petsclib, x); own = false)

	return on,oidxs,ogidxs
end 

