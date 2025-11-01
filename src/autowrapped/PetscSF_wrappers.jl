# autodefined type arguments for class ------
mutable struct _n_MPI_Group end
const MPI_Group = Ptr{_n_MPI_Group}



"""
	sf::PetscSF = PetscSFCreate(petsclib::PetscLibType,comm::MPI_Comm) 
create a star forest communication context

Collective

Input Parameter:
- `comm` - communicator on which the star forest will operate

Output Parameter:
- `sf` - new star forest context

Options Database Key:
- `-sf_type basic`                 - Use MPI persistent Isend/Irecv for communication (Default)
- `-sf_type window`                - Use MPI-3 one-sided window for communication
- `-sf_type neighbor`              - Use MPI-3 neighborhood collectives for communication
- `-sf_neighbor_persistent <bool>` - If true, use MPI-4 persistent neighborhood collectives for communication (used along with -sf_type neighbor)

Level: intermediate

-seealso: `PetscSF`, `PetscSFSetType`, `PetscSFSetGraph()`, `PetscSFSetGraphWithPattern()`, `PetscSFDestroy()`

# External Links
$(_doc_external("Vec/PetscSFCreate"))
"""
function PetscSFCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSFCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	sf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSF}),
               comm, sf_,
              )

	sf = sf_[]

	return sf
end 

"""
	PetscSFReset(petsclib::PetscLibType,sf::PetscSF) 
Reset a star forest so that different sizes or neighbors can be used

Collective

Input Parameter:
- `sf` - star forest

Level: advanced

-seealso: `PetscSF`, `PetscSFCreate()`, `PetscSFSetGraph()`, `PetscSFDestroy()`

# External Links
$(_doc_external("Vec/PetscSFReset"))
"""
function PetscSFReset(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFReset(petsclib::$UnionPetscLib, sf::PetscSF )

    @chk ccall(
               (:PetscSFReset, $petsc_library),
               PetscErrorCode,
               (PetscSF,),
               sf,
              )


	return nothing
end 

"""
	PetscSFSetType(petsclib::PetscLibType,sf::PetscSF, type::PetscSFType) 
Set the `PetscSF` communication implementation

Collective

Input Parameters:
- `sf`   - the `PetscSF` context
- `type` - a known method
-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFSetType"))
"""
function PetscSFSetType(petsclib::PetscLibType, sf::PetscSF, type::PetscSFType) end

@for_petsc function PetscSFSetType(petsclib::$UnionPetscLib, sf::PetscSF, type::PetscSFType )

    @chk ccall(
               (:PetscSFSetType, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFType),
               sf, type,
              )


	return nothing
end 

"""
	type::PetscSFType = PetscSFGetType(petsclib::PetscLibType,sf::PetscSF) 
Get the `PetscSF` communication implementation

Not Collective

Input Parameter:
- `sf` - the `PetscSF` context

Output Parameter:
- `type` - the `PetscSF` type name

Level: intermediate

-seealso: `PetscSF`, `PetscSFType`, `PetscSFSetType()`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFGetType"))
"""
function PetscSFGetType(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFGetType(petsclib::$UnionPetscLib, sf::PetscSF )
	type_ = Ref{PetscSFType}()

    @chk ccall(
               (:PetscSFGetType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFType}),
               sf, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscSFDestroy(petsclib::PetscLibType,sf::PetscSF) 
destroy a star forest

Collective

Input Parameter:
- `sf` - address of star forest

Level: intermediate

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFReset()`

# External Links
$(_doc_external("Vec/PetscSFDestroy"))
"""
function PetscSFDestroy(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFDestroy(petsclib::$UnionPetscLib, sf::PetscSF )

    @chk ccall(
               (:PetscSFDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSF},),
               sf,
              )


	return nothing
end 

"""
	PetscSFSetUp(petsclib::PetscLibType,sf::PetscSF) 
set up communication structures for a `PetscSF`, after this is done it may be used to perform communication

Collective

Input Parameter:
- `sf` - star forest communication object

Level: beginner

-seealso: `PetscSF`, `PetscSFType`, `PetscSFSetFromOptions()`, `PetscSFSetType()`

# External Links
$(_doc_external("Vec/PetscSFSetUp"))
"""
function PetscSFSetUp(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFSetUp(petsclib::$UnionPetscLib, sf::PetscSF )

    @chk ccall(
               (:PetscSFSetUp, $petsc_library),
               PetscErrorCode,
               (PetscSF,),
               sf,
              )


	return nothing
end 

"""
	PetscSFSetFromOptions(petsclib::PetscLibType,sf::PetscSF) 
set `PetscSF` options using the options database

Logically Collective

Input Parameter:
- `sf` - star forest

Options Database Keys:
- `-sf_type`                      - implementation type, see `PetscSFSetType()`
- `-sf_rank_order`                - sort composite points for gathers and scatters in rank order, gathers are non-deterministic otherwise
- `-sf_use_default_stream`        - Assume callers of `PetscSF` computed the input root/leafdata with the default CUDA stream. `PetscSF` will also
use the default stream to process data. Therefore, no stream synchronization is needed between `PetscSF` and its caller (default: true).
If true, this option only works with `-use_gpu_aware_mpi 1`.
- `-sf_use_stream_aware_mpi`      - Assume the underlying MPI is CUDA-stream aware and `PetscSF` won't sync streams for send/recv buffers passed to MPI (default: false).
If true, this option only works with `-use_gpu_aware_mpi 1`.

- `-sf_backend <cuda,hip,kokkos>` - Select the device backend`PetscSF` uses. Currently `PetscSF` has these backends: cuda - hip and Kokkos.
On CUDA (HIP) devices, one can choose cuda (hip) or kokkos with the default being kokkos. On other devices,
the only available is kokkos.

Level: intermediate

-seealso: `PetscSF`, `PetscSFCreate()`, `PetscSFSetType()`

# External Links
$(_doc_external("Vec/PetscSFSetFromOptions"))
"""
function PetscSFSetFromOptions(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFSetFromOptions(petsclib::$UnionPetscLib, sf::PetscSF )

    @chk ccall(
               (:PetscSFSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSF,),
               sf,
              )


	return nothing
end 

"""
	PetscSFSetRankOrder(petsclib::PetscLibType,sf::PetscSF, flg::PetscBool) 
sort multi

Logically Collective

Input Parameters:
- `sf`  - star forest
- `flg` - `PETSC_TRUE` to sort, `PETSC_FALSE` to skip sorting (lower setup cost, but non-deterministic)

Level: advanced

-seealso: `PetscSF`, `PetscSFType`, `PetscSFGatherBegin()`, `PetscSFScatterBegin()`

# External Links
$(_doc_external("Vec/PetscSFSetRankOrder"))
"""
function PetscSFSetRankOrder(petsclib::PetscLibType, sf::PetscSF, flg::PetscBool) end

@for_petsc function PetscSFSetRankOrder(petsclib::$UnionPetscLib, sf::PetscSF, flg::PetscBool )

    @chk ccall(
               (:PetscSFSetRankOrder, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscBool),
               sf, flg,
              )


	return nothing
end 

"""
	PetscSFSetGraph(petsclib::PetscLibType,sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, iloc::Vector{PetscInt}, locmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode) 
Set a parallel star forest

Collective

Input Parameters:
- `sf`         - star forest
- `nroots`     - number of root vertices on the current process (these are possible targets for other process to attach leaves)
- `nleaves`    - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`     - locations of leaves in leafdata buffers, pass `NULL` for contiguous storage (locations must be >= 0, enforced
during setup in debug mode)
- `localmode`  - copy mode for `ilocal`
- `iremote`    - remote locations of root vertices for each leaf on the current process, length is 2 `nleaves'
(locations must be >= 0, enforced during setup in debug mode)
- `remotemode` - copy mode for `iremote`

Level: intermediate

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFSetGraph"))
"""
function PetscSFSetGraph(petsclib::PetscLibType, sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, iloc::Vector{PetscInt}, locmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode) end

@for_petsc function PetscSFSetGraph(petsclib::$UnionPetscLib, sf::PetscSF, nroots::$PetscInt, nleaves::$PetscInt, iloc::Vector{$PetscInt}, locmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode )

    @chk ccall(
               (:PetscSFSetGraph, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{PetscSFNode}, PetscCopyMode),
               sf, nroots, nleaves, iloc, locmode, iremote, remotemode,
              )


	return nothing
end 

"""
	PetscSFSetGraphWithPattern(petsclib::PetscLibType,sf::PetscSF, map::PetscLayout, pattern::PetscSFPattern) 
Sets the graph of a `PetscSF` with a specific pattern

Collective

Input Parameters:
- `sf`      - The `PetscSF`
- `map`     - Layout of roots over all processes (insignificant when pattern is `PETSCSF_PATTERN_ALLTOALL`)
- `pattern` - One of `PETSCSF_PATTERN_ALLGATHER`, `PETSCSF_PATTERN_GATHER`, `PETSCSF_PATTERN_ALLTOALL`

Level: intermediate

-seealso: `PetscSF`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFSetGraphWithPattern"))
"""
function PetscSFSetGraphWithPattern(petsclib::PetscLibType, sf::PetscSF, map::PetscLayout, pattern::PetscSFPattern) end

@for_petsc function PetscSFSetGraphWithPattern(petsclib::$UnionPetscLib, sf::PetscSF, map::PetscLayout, pattern::PetscSFPattern )

    @chk ccall(
               (:PetscSFSetGraphWithPattern, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscLayout, PetscSFPattern),
               sf, map, pattern,
              )


	return nothing
end 

"""
	isf::PetscSF = PetscSFCreateInverseSF(petsclib::PetscLibType,sf::PetscSF) 
given a `PetscSF` in which all vertices have degree 1, creates the inverse map

Collective

Input Parameter:
- `sf` - star forest to invert

Output Parameter:
- `isf` - inverse of `sf`

Level: advanced

-seealso: `PetscSF`, `PetscSFType`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFCreateInverseSF"))
"""
function PetscSFCreateInverseSF(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFCreateInverseSF(petsclib::$UnionPetscLib, sf::PetscSF )
	isf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateInverseSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSF}),
               sf, isf_,
              )

	isf = isf_[]

	return isf
end 

"""
	newsf::PetscSF = PetscSFDuplicate(petsclib::PetscLibType,sf::PetscSF, opt::PetscSFDuplicateOption) 
duplicate a `PetscSF`, optionally preserving rank connectivity and graph

Collective

Input Parameters:
- `sf`  - communication object to duplicate
- `opt` - `PETSCSF_DUPLICATE_CONFONLY`, `PETSCSF_DUPLICATE_RANKS`, or `PETSCSF_DUPLICATE_GRAPH` (see `PetscSFDuplicateOption`)

Output Parameter:
- `newsf` - new communication object

Level: beginner

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFSetType()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFDuplicate"))
"""
function PetscSFDuplicate(petsclib::PetscLibType, sf::PetscSF, opt::PetscSFDuplicateOption) end

@for_petsc function PetscSFDuplicate(petsclib::$UnionPetscLib, sf::PetscSF, opt::PetscSFDuplicateOption )
	newsf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFDuplicateOption, Ptr{PetscSF}),
               sf, opt, newsf_,
              )

	newsf = newsf_[]

	return newsf
end 

"""
	nroots::PetscInt,nleaves::PetscInt,iloc::Vector{PetscInt} = PetscSFGetGraph(petsclib::PetscLibType,sf::PetscSF, iremote::Vector{PetscSFNode}) 
Get the graph specifying a parallel star forest

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `nroots`  - number of root vertices on the current process (these are possible targets for other process to attach leaves)
- `nleaves` - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`  - locations of leaves in leafdata buffers (if returned value is `NULL`, it means leaves are in contiguous storage)
- `iremote` - remote locations of root vertices for each leaf on the current process

Level: intermediate

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraph"))
"""
function PetscSFGetGraph(petsclib::PetscLibType, sf::PetscSF, iremote::Vector{PetscSFNode}) end

@for_petsc function PetscSFGetGraph(petsclib::$UnionPetscLib, sf::PetscSF, iremote::Vector{PetscSFNode} )
	nroots_ = Ref{$PetscInt}()
	nleaves_ = Ref{$PetscInt}()
	iloc_ = Ref{Ptr{$PetscInt}}()
	iremote_ = Ref(pointer(iremote))

    @chk ccall(
               (:PetscSFGetGraph, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{PetscSFNode}}),
               sf, nroots_, nleaves_, iloc_, iremote_,
              )

	nroots = nroots_[]
	nleaves = nleaves_[]
	iloc = unsafe_wrap(Array, iloc_[], VecGetLocalSize(petsclib, x); own = false)

	return nroots,nleaves,iloc
end 

"""
	minleaf::PetscInt,maxleaf::PetscInt = PetscSFGetLeafRange(petsclib::PetscLibType,sf::PetscSF) 
Get the active leaf ranges

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `minleaf` - minimum active leaf on this process. Returns 0 if there are no leaves.
- `maxleaf` - maximum active leaf on this process. Returns -1 if there are no leaves.

Level: developer

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetLeafRange"))
"""
function PetscSFGetLeafRange(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFGetLeafRange(petsclib::$UnionPetscLib, sf::PetscSF )
	minleaf_ = Ref{$PetscInt}()
	maxleaf_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSFGetLeafRange, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}),
               sf, minleaf_, maxleaf_,
              )

	minleaf = minleaf_[]
	maxleaf = maxleaf_[]

	return minleaf,maxleaf
end 

"""
	PetscSFViewFromOptions(petsclib::PetscLibType,A::PetscSF, obj::PetscObject, name::String) 
View a `PetscSF` based on arguments in the options database

Collective

Input Parameters:
- `A`    - the star forest
- `obj`  - Optional object that provides the prefix for the option names
- `name` - command line option

Level: intermediate

-seealso: `PetscSF`, `PetscSFView`, `PetscObjectViewFromOptions()`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFViewFromOptions"))
"""
function PetscSFViewFromOptions(petsclib::PetscLibType, A::PetscSF, obj::PetscObject, name::String) end

@for_petsc function PetscSFViewFromOptions(petsclib::$UnionPetscLib, A::PetscSF, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscSFViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscSFView(petsclib::PetscLibType,sf::PetscSF, viewer::PetscViewer) 
view a star forest

Collective

Input Parameters:
- `sf`     - star forest
- `viewer` - viewer to display graph, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: beginner

-seealso: `PetscSF`, `PetscViewer`, `PetscSFCreate()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFView"))
"""
function PetscSFView(petsclib::PetscLibType, sf::PetscSF, viewer::PetscViewer) end

@for_petsc function PetscSFView(petsclib::$UnionPetscLib, sf::PetscSF, viewer::PetscViewer )

    @chk ccall(
               (:PetscSFView, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscViewer),
               sf, viewer,
              )


	return nothing
end 

"""
	roffset::Vector{PetscInt},rmine::Vector{PetscInt},rremote::Vector{PetscInt} = PetscSFGetRootRanks(petsclib::PetscLibType,sf::PetscSF, nranks::PetscMPIInt, ranks::Vector{PetscMPIInt}) 
Get root ranks and number of vertices referenced by leaves on this process

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `nranks`  - number of ranks referenced by local part
- `ranks`   - [`nranks`] array of ranks
- `roffset` - [`nranks`+1] offset in `rmine`/`rremote` for each rank
- `rmine`   - [`roffset`[`nranks`]] concatenated array holding local indices referencing each remote rank, or `NULL`
- `rremote` - [`roffset`[`nranks`]] concatenated array holding remote indices referenced for each remote rank, or `NULL`

Level: developer

-seealso: `PetscSF`, `PetscSFGetLeafRanks()`

# External Links
$(_doc_external("Vec/PetscSFGetRootRanks"))
"""
function PetscSFGetRootRanks(petsclib::PetscLibType, sf::PetscSF, nranks::PetscMPIInt, ranks::Vector{PetscMPIInt}) end

@for_petsc function PetscSFGetRootRanks(petsclib::$UnionPetscLib, sf::PetscSF, nranks::PetscMPIInt, ranks::Vector{PetscMPIInt} )
	ranks_ = Ref(pointer(ranks))
	roffset_ = Ref{Ptr{$PetscInt}}()
	rmine_ = Ref{Ptr{$PetscInt}}()
	rremote_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetRootRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscMPIInt}, Ptr{Ptr{PetscMPIInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, nranks, ranks_, roffset_, rmine_, rremote_,
              )

	roffset = unsafe_wrap(Array, roffset_[], VecGetLocalSize(petsclib, x); own = false)
	rmine = unsafe_wrap(Array, rmine_[], VecGetLocalSize(petsclib, x); own = false)
	rremote = unsafe_wrap(Array, rremote_[], VecGetLocalSize(petsclib, x); own = false)

	return roffset,rmine,rremote
end 

"""
	ioffset::Vector{PetscInt},irootloc::Vector{PetscInt} = PetscSFGetLeafRanks(petsclib::PetscLibType,sf::PetscSF, niranks::PetscMPIInt, iranks::Vector{PetscMPIInt}) 
Get leaf ranks referencing roots on this process

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `niranks`  - number of leaf ranks referencing roots on this process
- `iranks`   - [`niranks`] array of ranks
- `ioffset`  - [`niranks`+1] offset in `irootloc` for each rank
- `irootloc` - [`ioffset`[`niranks`]] concatenated array holding local indices of roots referenced by each leaf rank

Level: developer

-seealso: `PetscSF`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("Vec/PetscSFGetLeafRanks"))
"""
function PetscSFGetLeafRanks(petsclib::PetscLibType, sf::PetscSF, niranks::PetscMPIInt, iranks::Vector{PetscMPIInt}) end

@for_petsc function PetscSFGetLeafRanks(petsclib::$UnionPetscLib, sf::PetscSF, niranks::PetscMPIInt, iranks::Vector{PetscMPIInt} )
	iranks_ = Ref(pointer(iranks))
	ioffset_ = Ref{Ptr{$PetscInt}}()
	irootloc_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetLeafRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscMPIInt}, Ptr{Ptr{PetscMPIInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, niranks, iranks_, ioffset_, irootloc_,
              )

	ioffset = unsafe_wrap(Array, ioffset_[], VecGetLocalSize(petsclib, x); own = false)
	irootloc = unsafe_wrap(Array, irootloc_[], VecGetLocalSize(petsclib, x); own = false)

	return ioffset,irootloc
end 

"""
	PetscSFSetUpRanks(petsclib::PetscLibType,sf::PetscSF, dgroup::MPI_Group) 
Set up data structures associated with ranks; this is for internal use by `PetscSF` implementations.

Collective

Input Parameters:
- `sf`     - `PetscSF` to set up; `PetscSFSetGraph()` must have been called
- `dgroup` - `MPI_Group` of ranks to be distinguished (e.g., for self or shared memory exchange)

Level: developer

-seealso: `PetscSF`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("Vec/PetscSFSetUpRanks"))
"""
function PetscSFSetUpRanks(petsclib::PetscLibType, sf::PetscSF, dgroup::MPI_Group) end

@for_petsc function PetscSFSetUpRanks(petsclib::$UnionPetscLib, sf::PetscSF, dgroup::MPI_Group )

    @chk ccall(
               (:PetscSFSetUpRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, MPI_Group),
               sf, dgroup,
              )


	return nothing
end 

"""
	PetscSFGetGroups(petsclib::PetscLibType,sf::PetscSF, incoming::MPI_Group, outgoing::MPI_Group) 
gets incoming and outgoing process groups

Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `incoming` - group of origin processes for incoming edges (leaves that reference my roots)
- `outgoing` - group of destination processes for outgoing edges (roots that I reference)

Level: developer

-seealso: `PetscSF`, `PetscSFGetWindow()`, `PetscSFRestoreWindow()`

# External Links
$(_doc_external("Vec/PetscSFGetGroups"))
"""
function PetscSFGetGroups(petsclib::PetscLibType, sf::PetscSF, incoming::MPI_Group, outgoing::MPI_Group) end

@for_petsc function PetscSFGetGroups(petsclib::$UnionPetscLib, sf::PetscSF, incoming::MPI_Group, outgoing::MPI_Group )

    @chk ccall(
               (:PetscSFGetGroups, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{MPI_Group}, Ptr{MPI_Group}),
               sf, incoming, outgoing,
              )


	return nothing
end 

"""
	PetscSFGetRanksSF(petsclib::PetscLibType,sf::PetscSF, rsf::PetscSF) 
gets the `PetscSF` to perform communications with root ranks

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `rsf` - the star forest with a single root per process to perform communications

Level: developer

-seealso: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("Vec/PetscSFGetRanksSF"))
"""
function PetscSFGetRanksSF(petsclib::PetscLibType, sf::PetscSF, rsf::PetscSF) end

@for_petsc function PetscSFGetRanksSF(petsclib::$UnionPetscLib, sf::PetscSF, rsf::PetscSF )

    @chk ccall(
               (:PetscSFGetRanksSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSF}),
               sf, rsf,
              )


	return nothing
end 

"""
	PetscSFGetMultiSF(petsclib::PetscLibType,sf::PetscSF, multi::PetscSF) 
gets the inner `PetscSF` implementing gathers and scatters

Collective

Input Parameter:
- `sf` - star forest that may contain roots with 0 or with more than 1 vertex

Output Parameter:
- `multi` - star forest with split roots, such that each root has degree exactly 1

Level: developer

-seealso: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGatherBegin()`, `PetscSFScatterBegin()`, `PetscSFComputeMultiRootOriginalNumbering()`

# External Links
$(_doc_external("Vec/PetscSFGetMultiSF"))
"""
function PetscSFGetMultiSF(petsclib::PetscLibType, sf::PetscSF, multi::PetscSF) end

@for_petsc function PetscSFGetMultiSF(petsclib::$UnionPetscLib, sf::PetscSF, multi::PetscSF )

    @chk ccall(
               (:PetscSFGetMultiSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSF}),
               sf, multi,
              )


	return nothing
end 

"""
	esf::PetscSF = PetscSFCreateEmbeddedRootSF(petsclib::PetscLibType,sf::PetscSF, nselected::PetscInt, selected::PetscInt) 
removes edges from all but the selected roots of a `PetscSF`, does not remap indices

Collective

Input Parameters:
- `sf`        - original star forest
- `nselected` - number of selected roots on this process
- `selected`  - indices of the selected roots on this process

Output Parameter:
- `esf` - new star forest

Level: advanced

-seealso: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFCreateEmbeddedRootSF"))
"""
function PetscSFCreateEmbeddedRootSF(petsclib::PetscLibType, sf::PetscSF, nselected::PetscInt, selected::PetscInt) end

@for_petsc function PetscSFCreateEmbeddedRootSF(petsclib::$UnionPetscLib, sf::PetscSF, nselected::$PetscInt, selected::$PetscInt )
	esf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateEmbeddedRootSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
               sf, nselected, selected, esf_,
              )

	esf = esf_[]

	return esf
end 

"""
	newsf::PetscSF = PetscSFCreateEmbeddedLeafSF(petsclib::PetscLibType,sf::PetscSF, nselected::PetscInt, selected::PetscInt) 
removes edges from all but the selected leaves of a `PetscSF`, does not remap indices

Collective

Input Parameters:
- `sf`        - original star forest
- `nselected` - number of selected leaves on this process
- `selected`  - indices of the selected leaves on this process

Output Parameter:
- `newsf` - new star forest

Level: advanced

-seealso: `PetscSF`, `PetscSFCreateEmbeddedRootSF()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFCreateEmbeddedLeafSF"))
"""
function PetscSFCreateEmbeddedLeafSF(petsclib::PetscLibType, sf::PetscSF, nselected::PetscInt, selected::PetscInt) end

@for_petsc function PetscSFCreateEmbeddedLeafSF(petsclib::$UnionPetscLib, sf::PetscSF, nselected::$PetscInt, selected::$PetscInt )
	newsf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateEmbeddedLeafSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
               sf, nselected, selected, newsf_,
              )

	newsf = newsf_[]

	return newsf
end 

"""
	degree::Vector{PetscInt} = PetscSFComputeDegreeBegin(petsclib::PetscLibType,sf::PetscSF) 
begin computation of degree for each root vertex, to be completed with `PetscSFComputeDegreeEnd()`

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `degree` - degree of each root vertex

Level: advanced

-seealso: `PetscSF`, `PetscSFGatherBegin()`, `PetscSFComputeDegreeEnd()`

# External Links
$(_doc_external("Vec/PetscSFComputeDegreeBegin"))
"""
function PetscSFComputeDegreeBegin(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFComputeDegreeBegin(petsclib::$UnionPetscLib, sf::PetscSF )
	degree_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFComputeDegreeBegin, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{Ptr{$PetscInt}}),
               sf, degree_,
              )

	degree = unsafe_wrap(Array, degree_[], VecGetLocalSize(petsclib, x); own = false)

	return degree
end 

"""
	degree::Vector{PetscInt} = PetscSFComputeDegreeEnd(petsclib::PetscLibType,sf::PetscSF) 
complete computation of degree for each root vertex, started with `PetscSFComputeDegreeBegin()`

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `degree` - degree of each root vertex

Level: developer

-seealso: `PetscSF`, `PetscSFGatherBegin()`, `PetscSFComputeDegreeBegin()`

# External Links
$(_doc_external("Vec/PetscSFComputeDegreeEnd"))
"""
function PetscSFComputeDegreeEnd(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFComputeDegreeEnd(petsclib::$UnionPetscLib, sf::PetscSF )
	degree_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFComputeDegreeEnd, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{Ptr{$PetscInt}}),
               sf, degree_,
              )

	degree = unsafe_wrap(Array, degree_[], VecGetLocalSize(petsclib, x); own = false)

	return degree
end 

"""
	nMultiRoots::PetscInt,multiRootsOrigNumbering::Vector{PetscInt} = PetscSFComputeMultiRootOriginalNumbering(petsclib::PetscLibType,sf::PetscSF, degree::Vector{PetscInt}) 
Returns original numbering of multi
Each multi-root is assigned index of the corresponding original root.

Collective

Input Parameters:
- `sf`     - star forest
- `degree` - degree of each root vertex, computed with `PetscSFComputeDegreeBegin()`/`PetscSFComputeDegreeEnd()`

Output Parameters:
- `nMultiRoots`             - (optional) number of multi-roots (roots of multi-`PetscSF`)
- `multiRootsOrigNumbering` - original indices of multi-roots; length of this array is `nMultiRoots`

Level: developer

-seealso: `PetscSF`, `PetscSFComputeDegreeBegin()`, `PetscSFComputeDegreeEnd()`, `PetscSFGetMultiSF()`

# External Links
$(_doc_external("Vec/PetscSFComputeMultiRootOriginalNumbering"))
"""
function PetscSFComputeMultiRootOriginalNumbering(petsclib::PetscLibType, sf::PetscSF, degree::Vector{PetscInt}) end

@for_petsc function PetscSFComputeMultiRootOriginalNumbering(petsclib::$UnionPetscLib, sf::PetscSF, degree::Vector{$PetscInt} )
	nMultiRoots_ = Ref{$PetscInt}()
	multiRootsOrigNumbering_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFComputeMultiRootOriginalNumbering, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               sf, degree, nMultiRoots_, multiRootsOrigNumbering_,
              )

	nMultiRoots = nMultiRoots_[]
	multiRootsOrigNumbering = unsafe_wrap(Array, multiRootsOrigNumbering_[], VecGetLocalSize(petsclib, x); own = false)

	return nMultiRoots,multiRootsOrigNumbering
end 

"""
	PetscSFCompose(petsclib::PetscLibType,sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF) 
Compose a new `PetscSF` by putting the second `PetscSF` under the first one in a top (roots) down (leaves) view

Input Parameters:
- `sfA` - The first `PetscSF`
- `sfB` - The second `PetscSF`

Output Parameter:
- `sfBA` - The composite `PetscSF`

Level: developer

-seealso: `PetscSF`, `PetscSFComposeInverse()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFCompose"))
"""
function PetscSFCompose(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF) end

@for_petsc function PetscSFCompose(petsclib::$UnionPetscLib, sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF )

    @chk ccall(
               (:PetscSFCompose, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfA, sfB, sfBA,
              )


	return nothing
end 

"""
	PetscSFComposeInverse(petsclib::PetscLibType,sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF) 
Compose a new `PetscSF` by putting the inverse of the second `PetscSF` under the first one

Input Parameters:
- `sfA` - The first `PetscSF`
- `sfB` - The second `PetscSF`

Output Parameter:
- `sfBA` - The composite `PetscSF`.

Level: developer

-seealso: `PetscSF`, `PetscSFCompose()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`, `PetscSFCreateInverseSF()`

# External Links
$(_doc_external("Vec/PetscSFComposeInverse"))
"""
function PetscSFComposeInverse(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF) end

@for_petsc function PetscSFComposeInverse(petsclib::$UnionPetscLib, sfA::PetscSF, sfB::PetscSF, sfBA::PetscSF )

    @chk ccall(
               (:PetscSFComposeInverse, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfA, sfB, sfBA,
              )


	return nothing
end 

"""
	PetscSFConcatenate(petsclib::PetscLibType,comm::MPI_Comm, nsfs::PetscInt, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::Vector{PetscInt}, newsf::PetscSF) 
concatenate multiple `PetscSF` into one

Input Parameters:
- `comm`        - the communicator
- `nsfs`        - the number of input `PetscSF`
- `sfs`         - the array of input `PetscSF`
- `rootMode`    - the root mode specifying how roots are handled
- `leafOffsets` - the array of local leaf offsets, one for each input `PetscSF`, or `NULL` for contiguous storage

Output Parameter:
- `newsf` - The resulting `PetscSF`

Level: advanced

-seealso: `PetscSF`, `PetscSFCompose()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`, `PetscSFConcatenateRootMode`

# External Links
$(_doc_external("Vec/PetscSFConcatenate"))
"""
function PetscSFConcatenate(petsclib::PetscLibType, comm::MPI_Comm, nsfs::PetscInt, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::Vector{PetscInt}, newsf::PetscSF) end

@for_petsc function PetscSFConcatenate(petsclib::$UnionPetscLib, comm::MPI_Comm, nsfs::$PetscInt, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::Vector{$PetscInt}, newsf::PetscSF )

    @chk ccall(
               (:PetscSFConcatenate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{PetscSF}, PetscSFConcatenateRootMode, Ptr{$PetscInt}, Ptr{PetscSF}),
               comm, nsfs, sfs, rootMode, leafOffsets, newsf,
              )


	return nothing
end 

"""
	PetscSFRegister(petsclib::PetscLibType,name::String, create::external) 
Adds an implementation of the `PetscSF` communication protocol.

Not Collective, No Fortran Support

Input Parameters:
- `name`   - name of a new user-defined implementation
- `create` - routine to create method context

-seealso: `PetscSF`, `PetscSFType`, `PetscSFRegisterAll()`, `PetscSFInitializePackage()`

# External Links
$(_doc_external("Vec/PetscSFRegister"))
"""
function PetscSFRegister(petsclib::PetscLibType, name::String, create::external) end

@for_petsc function PetscSFRegister(petsclib::$UnionPetscLib, name::String, create::external )

    @chk ccall(
               (:PetscSFRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               name, create,
              )


	return nothing
end 

"""
	PetscSFInitializePackage(petsclib::PetscLibType) 
Initialize `PetscSF` package

Logically Collective

Level: developer

-seealso: `PetscSF`, `PetscSFFinalizePackage()`

# External Links
$(_doc_external("Vec/PetscSFInitializePackage"))
"""
function PetscSFInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscSFInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSFInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSFFinalizePackage(petsclib::PetscLibType) 
Finalize `PetscSF` package, it is called from `PetscFinalize()`

Logically Collective

Level: developer

-seealso: `PetscSF`, `PetscSFInitializePackage()`

# External Links
$(_doc_external("Vec/PetscSFFinalizePackage"))
"""
function PetscSFFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscSFFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSFFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSFWindowSetFlavorType(petsclib::PetscLibType,sf::PetscSF, flavor::PetscSFWindowFlavorType) 
Set flavor type for `MPI_Win` creation

Logically Collective

Input Parameters:
- `sf`     - star forest for communication of type `PETSCSFWINDOW`
- `flavor` - flavor type

Options Database Key:
- `-sf_window_flavor <flavor>` - sets the flavor type CREATE, DYNAMIC, ALLOCATE or SHARED (see `PetscSFWindowFlavorType`)

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetFlavorType()`

# External Links
$(_doc_external("Vec/PetscSFWindowSetFlavorType"))
"""
function PetscSFWindowSetFlavorType(petsclib::PetscLibType, sf::PetscSF, flavor::PetscSFWindowFlavorType) end

@for_petsc function PetscSFWindowSetFlavorType(petsclib::$UnionPetscLib, sf::PetscSF, flavor::PetscSFWindowFlavorType )

    @chk ccall(
               (:PetscSFWindowSetFlavorType, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFWindowFlavorType),
               sf, flavor,
              )


	return nothing
end 

"""
	flavor::PetscSFWindowFlavorType = PetscSFWindowGetFlavorType(petsclib::PetscLibType,sf::PetscSF) 
Get  `PETSCSFWINDOW` flavor type for `PetscSF` communication

Logically Collective

Input Parameter:
- `sf` - star forest for communication of type `PETSCSFWINDOW`

Output Parameter:
- `flavor` - flavor type

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetFlavorType()`

# External Links
$(_doc_external("Vec/PetscSFWindowGetFlavorType"))
"""
function PetscSFWindowGetFlavorType(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFWindowGetFlavorType(petsclib::$UnionPetscLib, sf::PetscSF )
	flavor_ = Ref{PetscSFWindowFlavorType}()

    @chk ccall(
               (:PetscSFWindowGetFlavorType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFWindowFlavorType}),
               sf, flavor_,
              )

	flavor = unsafe_string(flavor_[])

	return flavor
end 

"""
	PetscSFWindowSetSyncType(petsclib::PetscLibType,sf::PetscSF, sync::PetscSFWindowSyncType) 
Set synchronization type for `PetscSF` communication of type  `PETSCSFWINDOW`

Logically Collective

Input Parameters:
- `sf`   - star forest for communication
- `sync` - synchronization type

Options Database Key:
- `-sf_window_sync <sync>` - sets the synchronization type FENCE, LOCK, or ACTIVE (see `PetscSFWindowSyncType`)

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetSyncType()`, `PetscSFWindowSyncType`

# External Links
$(_doc_external("Vec/PetscSFWindowSetSyncType"))
"""
function PetscSFWindowSetSyncType(petsclib::PetscLibType, sf::PetscSF, sync::PetscSFWindowSyncType) end

@for_petsc function PetscSFWindowSetSyncType(petsclib::$UnionPetscLib, sf::PetscSF, sync::PetscSFWindowSyncType )

    @chk ccall(
               (:PetscSFWindowSetSyncType, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFWindowSyncType),
               sf, sync,
              )


	return nothing
end 

"""
	sync::PetscSFWindowSyncType = PetscSFWindowGetSyncType(petsclib::PetscLibType,sf::PetscSF) 
Get synchronization type for `PetscSF` communication of type `PETSCSFWINDOW`

Logically Collective

Input Parameter:
- `sf` - star forest for communication

Output Parameter:
- `sync` - synchronization type

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetSyncType()`, `PetscSFWindowSyncType`

# External Links
$(_doc_external("Vec/PetscSFWindowGetSyncType"))
"""
function PetscSFWindowGetSyncType(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFWindowGetSyncType(petsclib::$UnionPetscLib, sf::PetscSF )
	sync_ = Ref{PetscSFWindowSyncType}()

    @chk ccall(
               (:PetscSFWindowGetSyncType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFWindowSyncType}),
               sf, sync_,
              )

	sync = unsafe_string(sync_[])

	return sync
end 

"""
	PetscSFWindowSetInfo(petsclib::PetscLibType,sf::PetscSF, info::MPI_Info) 
Set the `MPI_Info` handle that will be used for subsequent windows allocation

Logically Collective

Input Parameters:
- `sf`   - star forest for communication
- `info` - `MPI_Info` handle

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetInfo()`

# External Links
$(_doc_external("Vec/PetscSFWindowSetInfo"))
"""
function PetscSFWindowSetInfo(petsclib::PetscLibType, sf::PetscSF, info::MPI_Info) end

@for_petsc function PetscSFWindowSetInfo(petsclib::$UnionPetscLib, sf::PetscSF, info::MPI_Info )

    @chk ccall(
               (:PetscSFWindowSetInfo, $petsc_library),
               PetscErrorCode,
               (PetscSF, MPI_Info),
               sf, info,
              )


	return nothing
end 

"""
	PetscSFWindowGetInfo(petsclib::PetscLibType,sf::PetscSF, info::MPI_Info) 
Get the `MPI_Info` handle used for windows allocation

Logically Collective

Input Parameter:
- `sf` - star forest for communication

Output Parameter:
- `info` - `MPI_Info` handle

Level: advanced

-seealso: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetInfo()`

# External Links
$(_doc_external("Vec/PetscSFWindowGetInfo"))
"""
function PetscSFWindowGetInfo(petsclib::PetscLibType, sf::PetscSF, info::MPI_Info) end

@for_petsc function PetscSFWindowGetInfo(petsclib::$UnionPetscLib, sf::PetscSF, info::MPI_Info )

    @chk ccall(
               (:PetscSFWindowGetInfo, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{MPI_Info}),
               sf, info,
              )


	return nothing
end 

"""
	PetscSFSetGraphFromCoordinates(petsclib::PetscLibType,sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, dim::PetscInt, tol::PetscReal, rootcoords::PetscReal, leafcoords::PetscReal) 
Create SF by fuzzy matching leaf coordinates to root coordinates

Collective

Input Parameters:
- `sf`         - PetscSF to set graph on
- `nroots`     - number of root coordinates
- `nleaves`    - number of leaf coordinates
- `dim`        - spatial dimension of coordinates
- `tol`        - positive tolerance for matching
- `rootcoords` - array of root coordinates in which root i component d is [i*dim+d]
- `leafcoords` - array of root coordinates in which leaf i component d is [i*dim+d]

-seealso: `PetscSFCreate()`, `PetscSFSetGraph()`, `PetscSFCreateByMatchingIndices()`

# External Links
$(_doc_external("Vec/PetscSFSetGraphFromCoordinates"))
"""
function PetscSFSetGraphFromCoordinates(petsclib::PetscLibType, sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, dim::PetscInt, tol::PetscReal, rootcoords::PetscReal, leafcoords::PetscReal) end

@for_petsc function PetscSFSetGraphFromCoordinates(petsclib::$UnionPetscLib, sf::PetscSF, nroots::$PetscInt, nleaves::$PetscInt, dim::$PetscInt, tol::$PetscReal, rootcoords::$PetscReal, leafcoords::$PetscReal )

    @chk ccall(
               (:PetscSFSetGraphFromCoordinates, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, $PetscInt, $PetscInt, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sf, nroots, nleaves, dim, tol, rootcoords, leafcoords,
              )


	return nothing
end 

"""
	iloc::PetscInt = PetscSFSetGraphLayout(petsclib::PetscLibType,sf::PetscSF, layout::PetscLayout, nleaves::PetscInt, locmode::PetscCopyMode, gremote::PetscInt) 
Set a parallel star forest via global indices and a `PetscLayout`

Collective

Input Parameters:
- `sf`        - star forest
- `layout`    - `PetscLayout` defining the global space for roots
- `nleaves`   - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`    - locations of leaves in leafdata buffers, pass NULL for contiguous storage
- `localmode` - copy mode for ilocal
- `gremote`   - root vertices in global numbering corresponding to leaves in ilocal

Level: intermediate

-seealso: `PetscSF`, `PetscSFGetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFSetGraphLayout"))
"""
function PetscSFSetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout, nleaves::PetscInt, locmode::PetscCopyMode, gremote::PetscInt) end

@for_petsc function PetscSFSetGraphLayout(petsclib::$UnionPetscLib, sf::PetscSF, layout::PetscLayout, nleaves::$PetscInt, locmode::PetscCopyMode, gremote::$PetscInt )
	iloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSFSetGraphLayout, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscLayout, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{$PetscInt}),
               sf, layout, nleaves, iloc_, locmode, gremote,
              )

	iloc = iloc_[]

	return iloc
end 

"""
	nleaves::PetscInt,iloc::Vector{PetscInt},gremote::Vector{PetscInt} = PetscSFGetGraphLayout(petsclib::PetscLibType,sf::PetscSF, layout::PetscLayout) 
Get the global indices and `PetscLayout` that describe this star forest

Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `layout`  - `PetscLayout` defining the global space for roots
- `nleaves` - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`  - locations of leaves in leafdata buffers, or `NULL` for contiguous storage
- `gremote` - root vertices in global numbering corresponding to leaves in ilocal

Level: intermediate

-seealso: `PetscSF`, `PetscSFSetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraphLayout"))
"""
function PetscSFGetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout) end

@for_petsc function PetscSFGetGraphLayout(petsclib::$UnionPetscLib, sf::PetscSF, layout::PetscLayout )
	nleaves_ = Ref{$PetscInt}()
	iloc_ = Ref{Ptr{$PetscInt}}()
	gremote_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetGraphLayout, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscLayout}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, layout, nleaves_, iloc_, gremote_,
              )

	nleaves = nleaves_[]
	iloc = unsafe_wrap(Array, iloc_[], VecGetLocalSize(petsclib, x); own = false)
	gremote = unsafe_wrap(Array, gremote_[], VecGetLocalSize(petsclib, x); own = false)

	return nleaves,iloc,gremote
end 

"""
	PetscSFSetGraphSection(petsclib::PetscLibType,sf::PetscSF, locSection::PetscSection, globSection::PetscSection) 
Sets the `PetscSF` graph encoding the parallel dof overlap based upon the `PetscSection` describing the data layout.

Input Parameters:
- `sf`            - The `PetscSF`
- `localSection`  - `PetscSection` describing the local data layout
- `globalSection` - `PetscSection` describing the global data layout

Level: developer

-seealso: `PetscSF`, `PetscSFSetGraph()`, `PetscSFSetGraphLayout()`

# External Links
$(_doc_external("Vec/PetscSFSetGraphSection"))
"""
function PetscSFSetGraphSection(petsclib::PetscLibType, sf::PetscSF, locSection::PetscSection, globSection::PetscSection) end

@for_petsc function PetscSFSetGraphSection(petsclib::$UnionPetscLib, sf::PetscSF, locSection::PetscSection, globSection::PetscSection )

    @chk ccall(
               (:PetscSFSetGraphSection, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, PetscSection),
               sf, locSection, globSection,
              )


	return nothing
end 

"""
	remoteOffsets::Vector{PetscInt} = PetscSFDistributeSection(petsclib::PetscLibType,sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) 
Create a new `PetscSection` reorganized, moving from the root to the leaves of the `PetscSF`

Collective

Input Parameters:
- `sf`          - The `PetscSF`
- `rootSection` - Section defined on root space

Output Parameters:
- `remoteOffsets` - root offsets in leaf storage, or `NULL`, its length will be the size of the chart of `leafSection`
- `leafSection`   - Section defined on the leaf space

Level: advanced

-seealso: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFDistributeSection"))
"""
function PetscSFDistributeSection(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) end

@for_petsc function PetscSFDistributeSection(petsclib::$UnionPetscLib, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection )
	remoteOffsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFDistributeSection, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, Ptr{Ptr{$PetscInt}}, PetscSection),
               sf, rootSection, remoteOffsets_, leafSection,
              )

	remoteOffsets = unsafe_wrap(Array, remoteOffsets_[], VecGetLocalSize(petsclib, x); own = false)

	return remoteOffsets
end 

"""
	remoteOffsets::Vector{PetscInt} = PetscSFCreateRemoteOffsets(petsclib::PetscLibType,sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) 
Create offsets for point data on remote processes

Collective

Input Parameters:
- `sf`          - The `PetscSF`
- `rootSection` - Data layout of remote points for outgoing data (this is layout for roots)
- `leafSection` - Data layout of local points for incoming data  (this is layout for leaves)

Output Parameter:
- `remoteOffsets` - Offsets for point data on remote processes (these are offsets from the root section), or `NULL`

Level: developer

-seealso: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFCreateRemoteOffsets"))
"""
function PetscSFCreateRemoteOffsets(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) end

@for_petsc function PetscSFCreateRemoteOffsets(petsclib::$UnionPetscLib, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection )
	remoteOffsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFCreateRemoteOffsets, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, PetscSection, Ptr{Ptr{$PetscInt}}),
               sf, rootSection, leafSection, remoteOffsets_,
              )

	remoteOffsets = unsafe_wrap(Array, remoteOffsets_[], VecGetLocalSize(petsclib, x); own = false)

	return remoteOffsets
end 

"""
	sectionSF::PetscSF = PetscSFCreateSectionSF(petsclib::PetscLibType,sf::PetscSF, rootSection::PetscSection, remoteOffsets::Vector{PetscInt}, leafSection::PetscSection) 
Create an expanded `PetscSF` of dofs, assuming the input `PetscSF` relates points

Collective

Input Parameters:
- `sf`            - The `PetscSF`
- `rootSection`   - Data layout of remote points for outgoing data (this is usually the serial section)
- `remoteOffsets` - Offsets for point data on remote processes (these are offsets from the root section), or NULL
- `leafSection`   - Data layout of local points for incoming data  (this is the distributed section)

Output Parameter:
- `sectionSF` - The new `PetscSF`

Level: advanced

-seealso: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFCreateSectionSF"))
"""
function PetscSFCreateSectionSF(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, remoteOffsets::Vector{PetscInt}, leafSection::PetscSection) end

@for_petsc function PetscSFCreateSectionSF(petsclib::$UnionPetscLib, sf::PetscSF, rootSection::PetscSection, remoteOffsets::Vector{$PetscInt}, leafSection::PetscSection )
	sectionSF_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateSectionSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, Ptr{$PetscInt}, PetscSection, Ptr{PetscSF}),
               sf, rootSection, remoteOffsets, leafSection, sectionSF_,
              )

	sectionSF = sectionSF_[]

	return sectionSF
end 

"""
	sf::PetscSF = PetscSFCreateFromLayouts(petsclib::PetscLibType,rmap::PetscLayout, lmap::PetscLayout) 
Creates a parallel star forest mapping two `PetscLayout` objects

Collective

Input Parameters:
- `rmap` - `PetscLayout` defining the global root space
- `lmap` - `PetscLayout` defining the global leaf space

Output Parameter:
- `sf` - The parallel star forest

Level: intermediate

-seealso: `PetscSF`, `PetscSFCreate()`, `PetscLayoutCreate()`, `PetscSFSetGraphLayout()`

# External Links
$(_doc_external("Vec/PetscSFCreateFromLayouts"))
"""
function PetscSFCreateFromLayouts(petsclib::PetscLibType, rmap::PetscLayout, lmap::PetscLayout) end

@for_petsc function PetscSFCreateFromLayouts(petsclib::$UnionPetscLib, rmap::PetscLayout, lmap::PetscLayout )
	sf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateFromLayouts, $petsc_library),
               PetscErrorCode,
               (PetscLayout, PetscLayout, Ptr{PetscSF}),
               rmap, lmap, sf_,
              )

	sf = sf_[]

	return sf
end 

"""
	sfA::PetscSF,sf::PetscSF = PetscSFCreateByMatchingIndices(petsclib::PetscLibType,layout::PetscLayout, numRootIndices::PetscInt, rootIndices::PetscInt, rootLocalIndices::PetscInt, rootLocalOffset::PetscInt, numLeafIndices::PetscInt, leafIndices::PetscInt, leafLocalIndices::PetscInt, leafLocalOffset::PetscInt) 
Create `PetscSF` by matching root and leaf indices

Collective

Input Parameters:
- `layout`           - `PetscLayout` defining the global index space and the rank that brokers each index
- `numRootIndices`   - size of rootIndices
- `rootIndices`      - `PetscInt` array of global indices of which this process requests ownership
- `rootLocalIndices` - root local index permutation (NULL if no permutation)
- `rootLocalOffset`  - offset to be added to root local indices
- `numLeafIndices`   - size of leafIndices
- `leafIndices`      - `PetscInt` array of global indices with which this process requires data associated
- `leafLocalIndices` - leaf local index permutation (NULL if no permutation)
- `leafLocalOffset`  - offset to be added to leaf local indices

Output Parameters:
- `sfA` - star forest representing the communication pattern from the layout space to the leaf space (NULL if not needed)
- `sf`  - star forest representing the communication pattern from the root space to the leaf space

Level: advanced

Example 1:
-seealso: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("Vec/PetscSFCreateByMatchingIndices"))
"""
function PetscSFCreateByMatchingIndices(petsclib::PetscLibType, layout::PetscLayout, numRootIndices::PetscInt, rootIndices::PetscInt, rootLocalIndices::PetscInt, rootLocalOffset::PetscInt, numLeafIndices::PetscInt, leafIndices::PetscInt, leafLocalIndices::PetscInt, leafLocalOffset::PetscInt) end

@for_petsc function PetscSFCreateByMatchingIndices(petsclib::$UnionPetscLib, layout::PetscLayout, numRootIndices::$PetscInt, rootIndices::$PetscInt, rootLocalIndices::$PetscInt, rootLocalOffset::$PetscInt, numLeafIndices::$PetscInt, leafIndices::$PetscInt, leafLocalIndices::$PetscInt, leafLocalOffset::$PetscInt )
	sfA_ = Ref{PetscSF}()
	sf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateByMatchingIndices, $petsc_library),
               PetscErrorCode,
               (PetscLayout, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, Ptr{PetscSF}, Ptr{PetscSF}),
               layout, numRootIndices, rootIndices, rootLocalIndices, rootLocalOffset, numLeafIndices, leafIndices, leafLocalIndices, leafLocalOffset, sfA_, sf_,
              )

	sfA = sfA_[]
	sf = sf_[]

	return sfA,sf
end 

"""
	PetscSFMerge(petsclib::PetscLibType,sfa::PetscSF, sfb::PetscSF, merged::PetscSF) 
append/merge indices of `sfb` into `sfa`, with preference for `sfb`

Collective

Input Parameters:
- `sfa` - default `PetscSF`
- `sfb` - additional edges to add/replace edges in sfa

Output Parameter:
- `merged` - new `PetscSF` with combined edges

Level: intermediate

-seealso: `PetscSFCompose()`

# External Links
$(_doc_external("Vec/PetscSFMerge"))
"""
function PetscSFMerge(petsclib::PetscLibType, sfa::PetscSF, sfb::PetscSF, merged::PetscSF) end

@for_petsc function PetscSFMerge(petsclib::$UnionPetscLib, sfa::PetscSF, sfb::PetscSF, merged::PetscSF )

    @chk ccall(
               (:PetscSFMerge, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfa, sfb, merged,
              )


	return nothing
end 

"""
	vsf::PetscSF = PetscSFCreateStridedSF(petsclib::PetscLibType,sf::PetscSF, bs::PetscInt, ldr::PetscInt, ldl::PetscInt) 
Create an `PetscSF` to communicate interleaved blocks of data

Collective

Input Parameters:
- `sf`  - star forest
- `bs`  - stride
- `ldr` - leading dimension of root space
- `ldl` - leading dimension of leaf space

Output Parameter:
- `vsf` - the new `PetscSF`

Level: intermediate

-seealso: `PetscSF`, `PetscSFCreate()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFCreateStridedSF"))
"""
function PetscSFCreateStridedSF(petsclib::PetscLibType, sf::PetscSF, bs::PetscInt, ldr::PetscInt, ldl::PetscInt) end

@for_petsc function PetscSFCreateStridedSF(petsclib::$UnionPetscLib, sf::PetscSF, bs::$PetscInt, ldr::$PetscInt, ldl::$PetscInt )
	vsf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateStridedSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, $PetscInt, $PetscInt, Ptr{PetscSF}),
               sf, bs, ldr, ldl, vsf_,
              )

	vsf = vsf_[]

	return vsf
end 

"""
	PetscSFGetSubSF(petsclib::PetscLibType,mainsf::PetscSF, map::ISLocalToGlobalMapping, subSF::PetscSF) 
Returns an `PetscSF` for a specific subset of points. Leaves are re

Collective

Input Parameters:
- `mainsf` - `PetscSF` structure
- `map`    - a `ISLocalToGlobalMapping` that contains the subset of points

Output Parameter:
- `subSF` - a subset of the `mainSF` for the desired subset.

Level: intermediate

-seealso: `PetscSF`

# External Links
$(_doc_external("Dm/PetscSFGetSubSF"))
"""
function PetscSFGetSubSF(petsclib::PetscLibType, mainsf::PetscSF, map::ISLocalToGlobalMapping, subSF::PetscSF) end

@for_petsc function PetscSFGetSubSF(petsclib::$UnionPetscLib, mainsf::PetscSF, map::ISLocalToGlobalMapping, subSF::PetscSF )

    @chk ccall(
               (:PetscSFGetSubSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, ISLocalToGlobalMapping, Ptr{PetscSF}),
               mainsf, map, subSF,
              )


	return nothing
end 

"""
	PetscSFGetRanks(petsclib::PetscLibType,sf::PetscSF, nranks::PetscMPIInt, ranks::PetscMPIInt, roffset::PetscInt, rmine::PetscInt, rremote::PetscInt) 

# External Links
$(_doc_external("Vec/PetscSFGetRanks"))
"""
function PetscSFGetRanks(petsclib::PetscLibType, sf::PetscSF, nranks::PetscMPIInt, ranks::PetscMPIInt, roffset::PetscInt, rmine::PetscInt, rremote::PetscInt) end

@for_petsc function PetscSFGetRanks(petsclib::$UnionPetscLib, sf::PetscSF, nranks::PetscMPIInt, ranks::PetscMPIInt, roffset::$PetscInt, rmine::$PetscInt, rremote::$PetscInt )

    @chk ccall(
               (:PetscSFGetRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscMPIInt}, PetscMPIInt, $PetscInt, $PetscInt, $PetscInt),
               sf, nranks, ranks, roffset, rmine, rremote,
              )


	return nothing
end 

"""
	selected::PetscInt,esf::PetscSF = PetscSFCreateEmbeddedSF(petsclib::PetscLibType,sf::PetscSF, nselected::PetscInt) 

# External Links
$(_doc_external("Vec/PetscSFCreateEmbeddedSF"))
"""
function PetscSFCreateEmbeddedSF(petsclib::PetscLibType, sf::PetscSF, nselected::PetscInt) end

@for_petsc function PetscSFCreateEmbeddedSF(petsclib::$UnionPetscLib, sf::PetscSF, nselected::$PetscInt )
	selected_ = Ref{$PetscInt}()
	esf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCreateEmbeddedSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSF}),
               sf, nselected, selected_, esf_,
              )

	selected = selected_[]
	esf = esf_[]

	return selected,esf
end 

