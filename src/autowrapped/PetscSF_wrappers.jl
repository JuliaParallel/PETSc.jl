"""
	sfBA::PetscSF = PetscSFCompose(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF) 
Compose a new `PetscSF` by putting the second `PetscSF` under the first one in a top (roots) down (leaves) view

Input Parameters:
- `sfA` - The first `PetscSF`
- `sfB` - The second `PetscSF`

Output Parameter:
- `sfBA` - The composite `PetscSF`

Level: developer

See also: `PetscSF`, `PetscSFComposeInverse()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFCompose"))
"""
function PetscSFCompose(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF)
    error("PetscSFCompose: no generated method for these argument types")
end

@for_petsc function PetscSFCompose(petsclib::$UnionPetscLib, sfA::PetscSF, sfB::PetscSF )
	sfBA_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFCompose, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfA, sfB, sfBA_,
              )

	sfBA = sfBA_[]

	return sfBA
end 

"""
	sfBA::PetscSF = PetscSFComposeInverse(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF) 
Compose a new `PetscSF` by putting the inverse of the second `PetscSF` under the first one

Input Parameters:
- `sfA` - The first `PetscSF`
- `sfB` - The second `PetscSF`

Output Parameter:
- `sfBA` - The composite `PetscSF`.

Level: developer

See also: `PetscSF`, `PetscSFCompose()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`, `PetscSFCreateInverseSF()`

# External Links
$(_doc_external("PetscSF/PetscSFComposeInverse"))
"""
function PetscSFComposeInverse(petsclib::PetscLibType, sfA::PetscSF, sfB::PetscSF)
    error("PetscSFComposeInverse: no generated method for these argument types")
end

@for_petsc function PetscSFComposeInverse(petsclib::$UnionPetscLib, sfA::PetscSF, sfB::PetscSF )
	sfBA_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFComposeInverse, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfA, sfB, sfBA_,
              )

	sfBA = sfBA_[]

	return sfBA
end 

"""
	degree::Ptr{PetscInt} = PetscSFComputeDegreeBegin(petsclib::PetscLibType, sf::PetscSF) 
begin computation of the degree of each root vertex, to be completed with `PetscSFComputeDegreeEnd()`

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `degree` - degree (the number of leaves) of each root vertex

Level: advanced

See also: `PetscSF`, `PetscSFGatherBegin()`, `PetscSFComputeDegreeEnd()`

# External Links
$(_doc_external("PetscSF/PetscSFComputeDegreeBegin"))
"""
function PetscSFComputeDegreeBegin(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFComputeDegreeBegin: no generated method for these argument types")
end

@for_petsc function PetscSFComputeDegreeBegin(petsclib::$UnionPetscLib, sf::PetscSF )
	degree_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFComputeDegreeBegin, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{Ptr{$PetscInt}}),
               sf, degree_,
              )

	degree = degree_[]

	return degree
end 

"""
	degree::Ptr{PetscInt} = PetscSFComputeDegreeEnd(petsclib::PetscLibType, sf::PetscSF) 
complete computation of degree for each root vertex, started with `PetscSFComputeDegreeBegin()`

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `degree` - degree of each root vertex

Level: developer

See also: `PetscSF`, `PetscSFGatherBegin()`, `PetscSFComputeDegreeBegin()`

# External Links
$(_doc_external("PetscSF/PetscSFComputeDegreeEnd"))
"""
function PetscSFComputeDegreeEnd(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFComputeDegreeEnd: no generated method for these argument types")
end

@for_petsc function PetscSFComputeDegreeEnd(petsclib::$UnionPetscLib, sf::PetscSF )
	degree_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFComputeDegreeEnd, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{Ptr{$PetscInt}}),
               sf, degree_,
              )

	degree = degree_[]

	return degree
end 

"""
	nMultiRoots::PetscInt,multiRootsOrigNumbering::Vector{PetscInt} = PetscSFComputeMultiRootOriginalNumbering(petsclib::PetscLibType, sf::PetscSF, degree::Vector{PetscInt}) 
Returns original numbering of multi-roots (roots of multi-`PetscSF` returned by `PetscSFGetMultiSF()`).
Each multi-root is assigned index of the corresponding original root.

Collective

Input Parameters:
- `sf`     - star forest
- `degree` - degree of each root vertex, computed with `PetscSFComputeDegreeBegin()` and `PetscSFComputeDegreeEnd()`

Output Parameters:
- `nMultiRoots`             - (optional) number of multi-roots (roots of multi-`PetscSF`)
- `multiRootsOrigNumbering` - original indices of multi-roots; length of this array is `nMultiRoots`

Level: developer

See also: `PetscSF`, `PetscSFComputeDegreeBegin()`, `PetscSFComputeDegreeEnd()`, `PetscSFGetMultiSF()`

# External Links
$(_doc_external("PetscSF/PetscSFComputeMultiRootOriginalNumbering"))
"""
function PetscSFComputeMultiRootOriginalNumbering(petsclib::PetscLibType, sf::PetscSF, degree::AbstractVector{<:Number})
    error("PetscSFComputeMultiRootOriginalNumbering: no generated method for these argument types")
end

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
	multiRootsOrigNumbering = multiRootsOrigNumbering_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, multiRootsOrigNumbering_[], nMultiRoots; own = false)

	return nMultiRoots,multiRootsOrigNumbering
end 

"""
	newsf::PetscSF = PetscSFConcatenate(petsclib::PetscLibType, comm::MPI_Comm, nsfs::PetscInt, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::Vector{PetscInt}) 
concatenate multiple `PetscSF` into a new `PetscSF`

Input Parameters:
- `comm`        - the communicator
- `nsfs`        - the number of input `PetscSF`
- `sfs`         - the array of input `PetscSF`
- `rootMode`    - the root mode specifying how roots are handled
- `leafOffsets` - the array of local leaf offsets, one for each input `PetscSF`, or `NULL` for contiguous storage

Output Parameter:
- `newsf` - The resulting `PetscSF`

Level: advanced

See also: `PetscSF`, `PetscSFCompose()`, `PetscSFGetGraph()`, `PetscSFSetGraph()`, `PetscSFConcatenateRootMode`

# External Links
$(_doc_external("PetscSF/PetscSFConcatenate"))
"""
function PetscSFConcatenate(petsclib::PetscLibType, comm::MPI_Comm, nsfs::Integer, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::AbstractVector{<:Number})
    error("PetscSFConcatenate: no generated method for these argument types")
end

@for_petsc function PetscSFConcatenate(petsclib::$UnionPetscLib, comm::MPI_Comm, nsfs::$PetscInt, sfs::Vector{PetscSF}, rootMode::PetscSFConcatenateRootMode, leafOffsets::Vector{$PetscInt} )
	newsf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFConcatenate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{PetscSF}, PetscSFConcatenateRootMode, Ptr{$PetscInt}, Ptr{PetscSF}),
               comm, nsfs, sfs, rootMode, leafOffsets, newsf_,
              )

	newsf = newsf_[]

	return newsf
end 

"""
	sf::PetscSF = PetscSFCreate(petsclib::PetscLibType, comm::MPI_Comm) 
create a star forest communication context

Collective

Input Parameter:
- `comm` - communicator on which the star forest will operate

Output Parameter:
- `sf` - new star forest context

Options Database Key:
- `-sf_type (basic|window|neighbor)`     - Use MPI persistent Isend/Irecv, or MPI-3 one-sided window, or MPI-3 neighborhood collectives for communication
- `-sf_neighbor_persistent (true|false)` - Use MPI-4 persistent neighborhood collectives for communication (used along with `-sf_type neighbor`)

Level: intermediate

See also: `PetscSF`, `PetscSFSetType`, `PetscSFSetGraph()`, `PetscSFSetGraphWithPattern()`, `PetscSFDestroy()`

# External Links
$(_doc_external("PetscSF/PetscSFCreate"))
"""
function PetscSFCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscSFCreate: no generated method for these argument types")
end

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
	sfA::PetscSF,sf::PetscSF = PetscSFCreateByMatchingIndices(petsclib::PetscLibType, layout::PetscLayout, numRootIndices::PetscInt, rootIndices::Vector{PetscInt}, rootLocalIndices::Vector{PetscInt}, rootLocalOffset::PetscInt, numLeafIndices::PetscInt, leafIndices::Vector{PetscInt}, leafLocalIndices::Vector{PetscInt}, leafLocalOffset::PetscInt) 
Create `PetscSF` by matching root and leaf indices

Collective

Input Parameters:
- `layout`           - `PetscLayout` defining the global index space and the MPI rank that brokers each index
- `numRootIndices`   - size of `rootIndices`
- `rootIndices`      - array of global indices of which this process requests ownership
- `rootLocalIndices` - root local index permutation (`NULL` if no permutation)
- `rootLocalOffset`  - offset to be added to `rootLocalIndices`
- `numLeafIndices`   - size of `leafIndices`
- `leafIndices`      - array of global indices with which this process requires data associated
- `leafLocalIndices` - leaf local index permutation (`NULL` if no permutation)
- `leafLocalOffset`  - offset to be added to `leafLocalIndices`

Output Parameters:
- `sfA` - star forest representing the communication pattern from the layout space to the leaf space (`NULL` if not needed)
- `sf`  - star forest representing the communication pattern from the root space to the leaf space

Level: advanced

Example 1:
``
rank             : 0            1            2
rootIndices      : [1 0 2]      [3]          [3]
rootLocalOffset  : 100          200          300
layout           : [0 1]        [2]          [3]
leafIndices      : [0]          [2]          [0 3]
leafLocalOffset  : 400          500          600

would build the following PetscSF

[0] 400 <- (0,101)
[1] 500 <- (0,102)
[2] 600 <- (0,101)
[2] 601 <- (2,300)
``

Example 2:
``
rank             : 0               1               2
rootIndices      : [1 0 2]         [3]             [3]
rootLocalOffset  : 100             200             300
layout           : [0 1]           [2]             [3]
leafIndices      : rootIndices     rootIndices     rootIndices
leafLocalOffset  : rootLocalOffset rootLocalOffset rootLocalOffset

would build the following PetscSF

[1] 200 <- (2,300)
``

Example 3:
``
No process requests ownership of global index 1, but no process needs it.

rank             : 0            1            2
numRootIndices   : 2            1            1
rootIndices      : [0 2]        [3]          [3]
rootLocalOffset  : 100          200          300
layout           : [0 1]        [2]          [3]
numLeafIndices   : 1            1            2
leafIndices      : [0]          [2]          [0 3]
leafLocalOffset  : 400          500          600

would build the following PetscSF

[0] 400 <- (0,100)
[1] 500 <- (0,101)
[2] 600 <- (0,100)
[2] 601 <- (2,300)
``

See also: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateByMatchingIndices"))
"""
function PetscSFCreateByMatchingIndices(petsclib::PetscLibType, layout::PetscLayout, numRootIndices::Integer, rootIndices::AbstractVector{<:Number}, rootLocalIndices::AbstractVector{<:Number}, rootLocalOffset::Integer, numLeafIndices::Integer, leafIndices::AbstractVector{<:Number}, leafLocalIndices::AbstractVector{<:Number}, leafLocalOffset::Integer)
    error("PetscSFCreateByMatchingIndices: no generated method for these argument types")
end

@for_petsc function PetscSFCreateByMatchingIndices(petsclib::$UnionPetscLib, layout::PetscLayout, numRootIndices::$PetscInt, rootIndices::Vector{$PetscInt}, rootLocalIndices::Vector{$PetscInt}, rootLocalOffset::$PetscInt, numLeafIndices::$PetscInt, leafIndices::Vector{$PetscInt}, leafLocalIndices::Vector{$PetscInt}, leafLocalOffset::$PetscInt )
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
	newsf::PetscSF = PetscSFCreateEmbeddedLeafSF(petsclib::PetscLibType, sf::PetscSF, nselected::PetscInt, selected::Vector{PetscInt}) 
removes edges from all but the selected leaves of a `PetscSF`, does not remap indices

Collective

Input Parameters:
- `sf`        - original star forest
- `nselected` - number of selected leaves on this MPI process
- `selected`  - indices of the selected leaves on this MPI process

Output Parameter:
- `newsf` - new star forest

Level: advanced

See also: `PetscSF`, `PetscSFCreateEmbeddedRootSF()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateEmbeddedLeafSF"))
"""
function PetscSFCreateEmbeddedLeafSF(petsclib::PetscLibType, sf::PetscSF, nselected::Integer, selected::AbstractVector{<:Number})
    error("PetscSFCreateEmbeddedLeafSF: no generated method for these argument types")
end

@for_petsc function PetscSFCreateEmbeddedLeafSF(petsclib::$UnionPetscLib, sf::PetscSF, nselected::$PetscInt, selected::Vector{$PetscInt} )
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
	esf::PetscSF = PetscSFCreateEmbeddedRootSF(petsclib::PetscLibType, sf::PetscSF, nselected::PetscInt, selected::Vector{PetscInt}) 
removes edges from all but the selected roots of a `PetscSF`, does not remap indices

Collective

Input Parameters:
- `sf`        - original star forest
- `nselected` - number of selected roots on this MPI process
- `selected`  - indices of the selected roots on this MPI process

Output Parameter:
- `esf` - new star forest

Level: advanced

See also: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateEmbeddedRootSF"))
"""
function PetscSFCreateEmbeddedRootSF(petsclib::PetscLibType, sf::PetscSF, nselected::Integer, selected::AbstractVector{<:Number})
    error("PetscSFCreateEmbeddedRootSF: no generated method for these argument types")
end

@for_petsc function PetscSFCreateEmbeddedRootSF(petsclib::$UnionPetscLib, sf::PetscSF, nselected::$PetscInt, selected::Vector{$PetscInt} )
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
	sf::PetscSF = PetscSFCreateFromLayouts(petsclib::PetscLibType, rmap::PetscLayout, lmap::PetscLayout) 
Creates a parallel star forest mapping between two `PetscLayout` objects

Collective

Input Parameters:
- `rmap` - `PetscLayout` defining the global root space
- `lmap` - `PetscLayout` defining the global leaf space

Output Parameter:
- `sf` - The parallel star forest

Level: intermediate

See also: `PetscSF`, `PetscLayout`, `PetscSFCreate()`, `PetscSFSetGraph()`, `PetscLayoutCreate()`, `PetscSFSetGraphLayout()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateFromLayouts"))
"""
function PetscSFCreateFromLayouts(petsclib::PetscLibType, rmap::PetscLayout, lmap::PetscLayout)
    error("PetscSFCreateFromLayouts: no generated method for these argument types")
end

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
	isf::PetscSF = PetscSFCreateInverseSF(petsclib::PetscLibType, sf::PetscSF) 
given a `PetscSF` in which all roots have degree 1 (exactly one leaf), creates the inverse map

Collective

Input Parameter:
- `sf` - star forest to invert

Output Parameter:
- `isf` - inverse of `sf`

Level: advanced

See also: `PetscSF`, `PetscSFType`, `PetscSFSetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateInverseSF"))
"""
function PetscSFCreateInverseSF(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFCreateInverseSF: no generated method for these argument types")
end

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
	remoteOffsets::Ptr{PetscInt} = PetscSFCreateRemoteOffsets(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) 
Create offsets for point data on remote processes

Collective

Input Parameters:
- `sf`          - The `PetscSF`
- `rootSection` - Data layout of remote points for outgoing data (this is layout for roots)
- `leafSection` - Data layout of local points for incoming data  (this is layout for leaves)

Output Parameter:
- `remoteOffsets` - Offsets for point data on remote processes (these are offsets from the root section), or `NULL`

Level: developer

See also: `PetscSF`, `PetscSFCreate()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateRemoteOffsets"))
"""
function PetscSFCreateRemoteOffsets(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection)
    error("PetscSFCreateRemoteOffsets: no generated method for these argument types")
end

@for_petsc function PetscSFCreateRemoteOffsets(petsclib::$UnionPetscLib, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection )
	remoteOffsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFCreateRemoteOffsets, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, PetscSection, Ptr{Ptr{$PetscInt}}),
               sf, rootSection, leafSection, remoteOffsets_,
              )

	remoteOffsets = remoteOffsets_[]

	return remoteOffsets
end 

"""
	sectionSF::PetscSF = PetscSFCreateSectionSF(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, remoteOffsets::Vector{PetscInt}, leafSection::PetscSection) 
Create an expanded `PetscSF` of dofs, assuming the input `PetscSF` relates points

Collective

Input Parameters:
- `sf`            - The `PetscSF`
- `rootSection`   - Data layout of remote points for outgoing data (this is usually the serial section)
- `remoteOffsets` - Offsets for point data on remote processes (these are offsets from the root section), or `NULL`
- `leafSection`   - Data layout of local points for incoming data  (this is the distributed section)

Output Parameter:
- `sectionSF` - The new `PetscSF`

Level: advanced

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFDistributeSection()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateSectionSF"))
"""
function PetscSFCreateSectionSF(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, remoteOffsets::AbstractVector{<:Number}, leafSection::PetscSection)
    error("PetscSFCreateSectionSF: no generated method for these argument types")
end

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
	vsf::PetscSF = PetscSFCreateStridedSF(petsclib::PetscLibType, sf::PetscSF, bs::PetscInt, ldr::PetscInt, ldl::PetscInt) 
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

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFCreateStridedSF"))
"""
function PetscSFCreateStridedSF(petsclib::PetscLibType, sf::PetscSF, bs::Integer, ldr::Integer, ldl::Integer)
    error("PetscSFCreateStridedSF: no generated method for these argument types")
end

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
	PetscSFDestroy(petsclib::PetscLibType, sf::Union{PetscSF, Ref{PetscSF}}) 
destroy a star forest

Collective

Input Parameter:
- `sf` - address of star forest

Level: intermediate

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFReset()`

# External Links
$(_doc_external("PetscSF/PetscSFDestroy"))
"""
function PetscSFDestroy(petsclib::PetscLibType, sf::Union{PetscSF, Ref{PetscSF}})
    error("PetscSFDestroy: no generated method for these argument types")
end

@for_petsc function PetscSFDestroy(petsclib::$UnionPetscLib, sf::Union{PetscSF, Ref{PetscSF}} )
	sf_ = sf isa Base.RefValue ? sf : Ref{PetscSF}(sf)

    @chk ccall(
               (:PetscSFDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSF},),
               sf_,
              )


	return nothing
end 

"""
	remoteOffsets::Ptr{PetscInt} = PetscSFDistributeSection(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection) 
Create a new `PetscSection` reorganized, moving from the root to the leaves of the `PetscSF`

Collective

Input Parameters:
- `sf`          - The `PetscSF`
- `rootSection` - Section defined on root space

Output Parameters:
- `remoteOffsets` - root offsets in leaf storage, or `NULL`, its length will be the size of the chart of `leafSection`
- `leafSection`   - Section defined on the leaf space

Level: advanced

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFCreateSectionSF()`

# External Links
$(_doc_external("PetscSF/PetscSFDistributeSection"))
"""
function PetscSFDistributeSection(petsclib::PetscLibType, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection)
    error("PetscSFDistributeSection: no generated method for these argument types")
end

@for_petsc function PetscSFDistributeSection(petsclib::$UnionPetscLib, sf::PetscSF, rootSection::PetscSection, leafSection::PetscSection )
	remoteOffsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFDistributeSection, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, Ptr{Ptr{$PetscInt}}, PetscSection),
               sf, rootSection, remoteOffsets_, leafSection,
              )

	remoteOffsets = remoteOffsets_[]

	return remoteOffsets
end 

"""
	newsf::PetscSF = PetscSFDuplicate(petsclib::PetscLibType, sf::PetscSF, opt::PetscSFDuplicateOption) 
duplicate a `PetscSF`, optionally preserving rank connectivity and graph

Collective

Input Parameters:
- `sf`  - communication object to duplicate
- `opt` - `PETSCSF_DUPLICATE_CONFONLY`, `PETSCSF_DUPLICATE_RANKS`, or `PETSCSF_DUPLICATE_GRAPH` (see `PetscSFDuplicateOption`)

Output Parameter:
- `newsf` - new communication object

Level: beginner

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFSetType()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFDuplicate"))
"""
function PetscSFDuplicate(petsclib::PetscLibType, sf::PetscSF, opt::PetscSFDuplicateOption)
    error("PetscSFDuplicate: no generated method for these argument types")
end

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
	PetscSFFinalizePackage(petsclib::PetscLibType) 
Finalize `PetscSF` package, it is called from `PetscFinalize()`

Logically Collective

Level: developer

See also: `PetscSF`, `PetscSFInitializePackage()`

# External Links
$(_doc_external("PetscSF/PetscSFFinalizePackage"))
"""
function PetscSFFinalizePackage(petsclib::PetscLibType)
    error("PetscSFFinalizePackage: no generated method for these argument types")
end

@for_petsc function PetscSFFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSFFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

# override for PetscSFGetGraph; C signature: PetscSFGetGraph(PetscSF sf, PetscInt* nroots, PetscInt* nleaves, PetscInt* ilocal[], PetscSFNode* iremote[])
"""
	nroots::PetscInt,nleaves::PetscInt,iloc::Vector{PetscInt} = PetscSFGetGraph(petsclib::PetscLibType, sf::PetscSF, iremote::Vector{PetscSFNode}) 
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

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraph"))
"""
function PetscSFGetGraph(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFGetGraph(petsclib::$UnionPetscLib, sf::PetscSF )
	nroots_ = Ref{$PetscInt}()
	nleaves_ = Ref{$PetscInt}()
	ilocal_ = Ref{Ptr{$PetscInt}}(C_NULL)
	iremote_ = Ref{Ptr{PetscSFNode}}(C_NULL)
    @chk ccall(
               (:PetscSFGetGraph, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{PetscSFNode}}),
               sf, nroots_, nleaves_, ilocal_, iremote_,
              )
	nroots = nroots_[]
	nleaves = nleaves_[]
	# ilocal == NULL means the leaves are contiguous [0, nleaves)
	ilocal = ilocal_[] == C_NULL ? nothing : unsafe_wrap(Array, ilocal_[], max(nleaves, 0); own = false)
	iremote = iremote_[] == C_NULL ? nothing : unsafe_wrap(Array, iremote_[], max(nleaves, 0); own = false)
	return nroots,nleaves,ilocal,iremote
end

# override for PetscSFGetGraphLayout; C signature: PetscSFGetGraphLayout(PetscSF sf, PetscLayout* layout, PetscInt* nleaves, PetscInt* ilocal[], PetscInt* gremote[])
"""
	nleaves::PetscInt,iloc::Vector{PetscInt},gremote::Vector{PetscInt} = PetscSFGetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout) 
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

See also: `PetscSF`, `PetscSFSetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

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
	incoming::MPI_Group,outgoing::MPI_Group = PetscSFGetGroups(petsclib::PetscLibType, sf::PetscSF) 
gets incoming and outgoing process groups

Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `incoming` - group of origin processes for incoming edges (leaves that reference my roots)
- `outgoing` - group of destination processes for outgoing edges (roots that I reference)

Level: developer

See also: `PetscSF`, `PetscSFGetWindow()`, `PetscSFRestoreWindow()`

# External Links
$(_doc_external("PetscSF/PetscSFGetGroups"))
"""
function PetscSFGetGroups(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetGroups: no generated method for these argument types")
end

@for_petsc function PetscSFGetGroups(petsclib::$UnionPetscLib, sf::PetscSF )
	incoming_ = Ref{MPI_Group}()
	outgoing_ = Ref{MPI_Group}()

    @chk ccall(
               (:PetscSFGetGroups, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{MPI_Group}, Ptr{MPI_Group}),
               sf, incoming_, outgoing_,
              )

	incoming = incoming_[]
	outgoing = outgoing_[]

	return incoming,outgoing
end 

"""
	minleaf::PetscInt,maxleaf::PetscInt = PetscSFGetLeafRange(petsclib::PetscLibType, sf::PetscSF) 
Get the active leaf ranges

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `minleaf` - minimum active leaf on this MPI process. Returns 0 if there are no leaves.
- `maxleaf` - maximum active leaf on this MPI process. Returns -1 if there are no leaves.

Level: developer

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFGetLeafRange"))
"""
function PetscSFGetLeafRange(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetLeafRange: no generated method for these argument types")
end

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
	niranks::PetscMPIInt,iranks::Vector{PetscMPIInt},ioffset::Vector{PetscInt},irootloc::Vector{PetscInt} = PetscSFGetLeafRanks(petsclib::PetscLibType, sf::PetscSF) 
Get leaf MPI ranks referencing roots on this process

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `niranks`  - number of leaf MPI processes referencing roots on this process
- `iranks`   - [`niranks`] array of MPI ranks
- `ioffset`  - [`niranks`+1] offset in `irootloc` for each MPI process
- `irootloc` - [`ioffset`[`niranks`]] concatenated array holding local indices of roots referenced by each leaf MPI process

Level: developer

See also: `PetscSF`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("PetscSF/PetscSFGetLeafRanks"))
"""
function PetscSFGetLeafRanks(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetLeafRanks: no generated method for these argument types")
end

@for_petsc function PetscSFGetLeafRanks(petsclib::$UnionPetscLib, sf::PetscSF )
	niranks_ = Ref{PetscMPIInt}()
	iranks_ = Ref{Ptr{PetscMPIInt}}()
	ioffset_ = Ref{Ptr{$PetscInt}}()
	irootloc_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetLeafRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscMPIInt}, Ptr{Ptr{PetscMPIInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, niranks_, iranks_, ioffset_, irootloc_,
              )

	niranks = niranks_[]
	iranks = iranks_[] == C_NULL ? PetscMPIInt[] : unsafe_wrap(Array, iranks_[], niranks; own = false)
	ioffset = ioffset_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, ioffset_[], niranks + 1; own = false)
	irootloc = irootloc_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, irootloc_[], ioffset[end]; own = false)

	return niranks,iranks,ioffset,irootloc
end 

"""
	multi::PetscSF = PetscSFGetMultiSF(petsclib::PetscLibType, sf::PetscSF) 
gets the inner `PetscSF` implementing gathers and scatters

Collective

Input Parameter:
- `sf` - star forest that may contain roots with 0 or with more than 1 vertex

Output Parameter:
- `multi` - star forest with split roots, such that each root has degree exactly 1 (has one leaf)

Level: developer

See also: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGatherBegin()`, `PetscSFScatterBegin()`, `PetscSFComputeMultiRootOriginalNumbering()`

# External Links
$(_doc_external("PetscSF/PetscSFGetMultiSF"))
"""
function PetscSFGetMultiSF(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetMultiSF: no generated method for these argument types")
end

@for_petsc function PetscSFGetMultiSF(petsclib::$UnionPetscLib, sf::PetscSF )
	multi_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFGetMultiSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSF}),
               sf, multi_,
              )

	multi = multi_[]

	return multi
end 

"""
	rsf::PetscSF = PetscSFGetRanksSF(petsclib::PetscLibType, sf::PetscSF) 
gets the `PetscSF` to perform communications with root ranks

Collective

Input Parameter:
- `sf` - star forest

Output Parameter:
- `rsf` - the star forest with a single root per MPI process to perform communications

Level: developer

See also: `PetscSF`, `PetscSFSetGraph()`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("PetscSF/PetscSFGetRanksSF"))
"""
function PetscSFGetRanksSF(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetRanksSF: no generated method for these argument types")
end

@for_petsc function PetscSFGetRanksSF(petsclib::$UnionPetscLib, sf::PetscSF )
	rsf_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFGetRanksSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSF}),
               sf, rsf_,
              )

	rsf = rsf_[]

	return rsf
end 

"""
	nranks::PetscMPIInt,ranks::Vector{PetscMPIInt},roffset::Vector{PetscInt},rmine::Vector{PetscInt},rremote::Vector{PetscInt} = PetscSFGetRootRanks(petsclib::PetscLibType, sf::PetscSF) 
Get the root MPI ranks and number of vertices referenced by leaves on this process

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `nranks`  - number of MPI processes referenced by local part
- `ranks`   - [`nranks`] array of MPI ranks
- `roffset` - [`nranks`+1] offset in `rmine` and `rremote` for each MPI process
- `rmine`   - [`roffset`[`nranks`]] concatenated array holding local indices referencing each remote MPI process, or `NULL`
- `rremote` - [`roffset`[`nranks`]] concatenated array holding remote indices referenced for each remote MPI process, or `NULL`

Level: developer

See also: `PetscSF`, `PetscSFGetLeafRanks()`

# External Links
$(_doc_external("PetscSF/PetscSFGetRootRanks"))
"""
function PetscSFGetRootRanks(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetRootRanks: no generated method for these argument types")
end

@for_petsc function PetscSFGetRootRanks(petsclib::$UnionPetscLib, sf::PetscSF )
	nranks_ = Ref{PetscMPIInt}()
	ranks_ = Ref{Ptr{PetscMPIInt}}()
	roffset_ = Ref{Ptr{$PetscInt}}()
	rmine_ = Ref{Ptr{$PetscInt}}()
	rremote_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetRootRanks, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscMPIInt}, Ptr{Ptr{PetscMPIInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, nranks_, ranks_, roffset_, rmine_, rremote_,
              )

	nranks = nranks_[]
	ranks = ranks_[] == C_NULL ? PetscMPIInt[] : unsafe_wrap(Array, ranks_[], nranks; own = false)
	roffset = roffset_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, roffset_[], nranks + 1; own = false)
	rmine = rmine_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, rmine_[], roffset[end]; own = false)
	rremote = rremote_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, rremote_[], roffset[end]; own = false)

	return nranks,ranks,roffset,rmine,rremote
end 

"""
	subSF::PetscSF = PetscSFGetSubSF(petsclib::PetscLibType, mainsf::PetscSF, map::ISLocalToGlobalMapping) 
Returns an `PetscSF` for a specific subset of points. Leaves are re-numbered to reflect the new ordering

Collective

Input Parameters:
- `mainsf` - `PetscSF` structure
- `map`    - a `ISLocalToGlobalMapping` that contains the subset of points

Output Parameter:
- `subSF` - a subset of the `mainSF` for the desired subset.

Level: intermediate

See also: `PetscSF`

# External Links
$(_doc_external("DMNetwork/PetscSFGetSubSF"))
"""
function PetscSFGetSubSF(petsclib::PetscLibType, mainsf::PetscSF, map::ISLocalToGlobalMapping)
    error("PetscSFGetSubSF: no generated method for these argument types")
end

@for_petsc function PetscSFGetSubSF(petsclib::$UnionPetscLib, mainsf::PetscSF, map::ISLocalToGlobalMapping )
	subSF_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFGetSubSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, ISLocalToGlobalMapping, Ptr{PetscSF}),
               mainsf, map, subSF_,
              )

	subSF = subSF_[]

	return subSF
end 

"""
	type::String = PetscSFGetType(petsclib::PetscLibType, sf::PetscSF) 
Get the `PetscSF` communication implementation

Not Collective

Input Parameter:
- `sf` - the `PetscSF` context

Output Parameter:
- `type` - the `PetscSF` type name

Level: intermediate

See also: `PetscSF`, `PetscSFType`, `PetscSFSetType()`, `PetscSFCreate()`

# External Links
$(_doc_external("PetscSF/PetscSFGetType"))
"""
function PetscSFGetType(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFGetType: no generated method for these argument types")
end

@for_petsc function PetscSFGetType(petsclib::$UnionPetscLib, sf::PetscSF )
	type_ = Ref{PetscSFType}()

    @chk ccall(
               (:PetscSFGetType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFType}),
               sf, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	PetscSFInitializePackage(petsclib::PetscLibType) 
Initialize `PetscSF` package

Logically Collective

Level: developer

See also: `PetscSF`, `PetscSFFinalizePackage()`

# External Links
$(_doc_external("Sys/PetscSFInitializePackage"))
"""
function PetscSFInitializePackage(petsclib::PetscLibType)
    error("PetscSFInitializePackage: no generated method for these argument types")
end

@for_petsc function PetscSFInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSFInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	merged::PetscSF = PetscSFMerge(petsclib::PetscLibType, sfa::PetscSF, sfb::PetscSF) 
append/merge indices of `sfb` into `sfa`, with preference for `sfb`

Collective

Input Parameters:
- `sfa` - default `PetscSF`
- `sfb` - additional edges to add/replace edges in `sfa`

Output Parameter:
- `merged` - new `PetscSF` with combined edges

Level: intermediate

See also: `PetscSF`, `PetscSFCompose()`

# External Links
$(_doc_external("PetscSF/PetscSFMerge"))
"""
function PetscSFMerge(petsclib::PetscLibType, sfa::PetscSF, sfb::PetscSF)
    error("PetscSFMerge: no generated method for these argument types")
end

@for_petsc function PetscSFMerge(petsclib::$UnionPetscLib, sfa::PetscSF, sfb::PetscSF )
	merged_ = Ref{PetscSF}()

    @chk ccall(
               (:PetscSFMerge, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSF, Ptr{PetscSF}),
               sfa, sfb, merged_,
              )

	merged = merged_[]

	return merged
end 

"""
	PetscSFRegister(petsclib::PetscLibType, name::String, create::external) 
Adds an implementation of the `PetscSF` communication protocol.

Not Collective, No Fortran Support

Input Parameters:
- `name`   - name of a new user-defined implementation
- `create` - routine to create method context

See also: `PetscSF`, `PetscSFType`, `PetscSFRegisterAll()`, `PetscSFInitializePackage()`

# External Links
$(_doc_external("PetscSF/PetscSFRegister"))
"""
function PetscSFRegister(petsclib::PetscLibType, name::String, create::external)
    error("PetscSFRegister: no generated method for these argument types")
end

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
	PetscSFReset(petsclib::PetscLibType, sf::PetscSF) 
Reset a star forest so that different sizes or neighbors can be used

Collective

Input Parameter:
- `sf` - star forest

Level: advanced

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFSetGraph()`, `PetscSFDestroy()`

# External Links
$(_doc_external("PetscSF/PetscSFReset"))
"""
function PetscSFReset(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFReset: no generated method for these argument types")
end

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
	PetscSFSetFromOptions(petsclib::PetscLibType, sf::PetscSF) 
set `PetscSF` options using the options database

Logically Collective

Input Parameter:
- `sf` - star forest

Options Database Keys:
- `-sf_type (basic|window|neighbor)` - implementation type, see `PetscSFSetType()`
- `-sf_rank_order (true|false)`      - sort composite points for gathers and scatters in MPI rank order, gathers are non-deterministic otherwise
- `-sf_use_default_stream`           - Assume callers of `PetscSF` computed the input root/leafdata with the default CUDA stream. `PetscSF` will also
use the default stream to process data. Therefore, no stream synchronization is needed between `PetscSF` and its caller (default: true).
If true, this option only works with `-use_gpu_aware_mpi 1`.
- `-sf_use_stream_aware_mpi`         - Assume the underlying MPI is CUDA-stream aware and `PetscSF` won't sync streams for send/recv buffers passed to MPI (default: false).
If true, this option only works with `-use_gpu_aware_mpi 1`.
- `-sf_backend (cuda|hip|kokkos)`    - Select the device backend `PetscSF` uses. On CUDA (HIP) devices, one can choose `cuda` (`hip`) or `kokkos` with the default being `kokkos`.
On other devices, the only available is `kokkos`.

Level: intermediate

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFSetType()`

# External Links
$(_doc_external("PetscSF/PetscSFSetFromOptions"))
"""
function PetscSFSetFromOptions(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFSetFromOptions: no generated method for these argument types")
end

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
	PetscSFSetGraph(petsclib::PetscLibType, sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, ilocal::Vector{PetscInt}, localmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode) 
Set a parallel star forest

Collective

Input Parameters:
- `sf`         - star forest
- `nroots`     - number of root vertices on the current MPI process (these are possible targets for other process to attach leaves)
- `nleaves`    - number of leaf vertices on the current MPI process, each of these references a root on any process
- `ilocal`     - locations of leaves in leafdata buffers (locations must be >= 0, enforced during setup in debug mode), pass `NULL` for contiguous storage (same as passing (0, 1, 2, ..., nleaves-1))
- `localmode`  - copy mode for `ilocal`
- `iremote`    - remote locations of root vertices for each leaf on the current process, length is `nleaves' (locations must be >= 0, enforced during setup in debug mode)
- `remotemode` - copy mode for `iremote`

Level: intermediate

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFGetGraph()`, `PetscSFSetGraphWithPattern()`

# External Links
$(_doc_external("PetscSF/PetscSFSetGraph"))
"""
function PetscSFSetGraph(petsclib::PetscLibType, sf::PetscSF, nroots::Integer, nleaves::Integer, ilocal::AbstractVector{<:Number}, localmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode)
    error("PetscSFSetGraph: no generated method for these argument types")
end

@for_petsc function PetscSFSetGraph(petsclib::$UnionPetscLib, sf::PetscSF, nroots::$PetscInt, nleaves::$PetscInt, ilocal::Vector{$PetscInt}, localmode::PetscCopyMode, iremote::Vector{PetscSFNode}, remotemode::PetscCopyMode )

    @chk ccall(
               (:PetscSFSetGraph, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{PetscSFNode}, PetscCopyMode),
               sf, nroots, nleaves, ilocal, localmode, iremote, remotemode,
              )


	return nothing
end 

"""
	PetscSFSetGraphFromCoordinates(petsclib::PetscLibType, sf::PetscSF, nroots::PetscInt, nleaves::PetscInt, dim::PetscInt, tol::PetscReal, rootcoords::Vector{PetscReal}, leafcoords::Vector{PetscReal}) 
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

See also: `PetscSFCreate()`, `PetscSFSetGraph()`, `PetscSFCreateByMatchingIndices()`

# External Links
$(_doc_external("PetscSF/PetscSFSetGraphFromCoordinates"))
"""
function PetscSFSetGraphFromCoordinates(petsclib::PetscLibType, sf::PetscSF, nroots::Integer, nleaves::Integer, dim::Integer, tol::Real, rootcoords::AbstractVector{<:Number}, leafcoords::AbstractVector{<:Number})
    error("PetscSFSetGraphFromCoordinates: no generated method for these argument types")
end

@for_petsc function PetscSFSetGraphFromCoordinates(petsclib::$UnionPetscLib, sf::PetscSF, nroots::$PetscInt, nleaves::$PetscInt, dim::$PetscInt, tol::$PetscReal, rootcoords::Vector{$PetscReal}, leafcoords::Vector{$PetscReal} )

    @chk ccall(
               (:PetscSFSetGraphFromCoordinates, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, $PetscInt, $PetscInt, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sf, nroots, nleaves, dim, tol, rootcoords, leafcoords,
              )


	return nothing
end 

"""
	PetscSFSetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout, nleaves::PetscInt, ilocal::Vector{PetscInt}, localmode::PetscCopyMode, gremote::Vector{PetscInt}) 
Set a `PetscSF` communication pattern using global indices and a `PetscLayout`

Collective

Input Parameters:
- `sf`        - star forest
- `layout`    - `PetscLayout` defining the global space for roots, i.e. which roots are owned by each MPI process
- `nleaves`   - number of leaf vertices on the current process, each of these references a root on any MPI process
- `ilocal`    - locations of leaves in leafdata buffers, pass `NULL` for contiguous storage, that is the locations are in [0,`nleaves`)
- `localmode` - copy mode for `ilocal`
- `gremote`   - root vertices in global numbering corresponding to the leaves

Level: intermediate

See also: `PetscSF`, `PetscSFGetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFSetGraphLayout"))
"""
function PetscSFSetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout, nleaves::Integer, ilocal::AbstractVector{<:Number}, localmode::PetscCopyMode, gremote::AbstractVector{<:Number})
    error("PetscSFSetGraphLayout: no generated method for these argument types")
end

@for_petsc function PetscSFSetGraphLayout(petsclib::$UnionPetscLib, sf::PetscSF, layout::PetscLayout, nleaves::$PetscInt, ilocal::Vector{$PetscInt}, localmode::PetscCopyMode, gremote::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSFSetGraphLayout, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscLayout, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{$PetscInt}),
               sf, layout, nleaves, ilocal, localmode, gremote,
              )


	return nothing
end 

"""
	PetscSFSetGraphSection(petsclib::PetscLibType, sf::PetscSF, localSection::PetscSection, globalSection::PetscSection) 
Sets the `PetscSF` graph (communication pattern) encoding the parallel dof overlap based upon the `PetscSection` describing the data layout.

Input Parameters:
- `sf`            - The `PetscSF`
- `localSection`  - `PetscSection` describing the local data layout
- `globalSection` - `PetscSection` describing the global data layout

Level: developer

See also: `PetscSF`, `PetscSFSetGraph()`, `PetscSFSetGraphLayout()`

# External Links
$(_doc_external("PetscSF/PetscSFSetGraphSection"))
"""
function PetscSFSetGraphSection(petsclib::PetscLibType, sf::PetscSF, localSection::PetscSection, globalSection::PetscSection)
    error("PetscSFSetGraphSection: no generated method for these argument types")
end

@for_petsc function PetscSFSetGraphSection(petsclib::$UnionPetscLib, sf::PetscSF, localSection::PetscSection, globalSection::PetscSection )

    @chk ccall(
               (:PetscSFSetGraphSection, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSection, PetscSection),
               sf, localSection, globalSection,
              )


	return nothing
end 

"""
	PetscSFSetGraphWithPattern(petsclib::PetscLibType, sf::PetscSF, map::PetscLayout, pattern::PetscSFPattern) 
Sets the graph of a `PetscSF` with a specific pattern

Collective

Input Parameters:
- `sf`      - The `PetscSF`
- `map`     - Layout of roots over all processes (not used when pattern is `PETSCSF_PATTERN_ALLTOALL`)
- `pattern` - One of `PETSCSF_PATTERN_ALLGATHER`, `PETSCSF_PATTERN_GATHER`, `PETSCSF_PATTERN_ALLTOALL`

Level: intermediate

See also: `PetscSF`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFSetGraphWithPattern"))
"""
function PetscSFSetGraphWithPattern(petsclib::PetscLibType, sf::PetscSF, map::PetscLayout, pattern::PetscSFPattern)
    error("PetscSFSetGraphWithPattern: no generated method for these argument types")
end

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
	PetscSFSetRankOrder(petsclib::PetscLibType, sf::PetscSF, flg::PetscBool) 
sort multi-points for gathers and scatters by MPI rank order

Logically Collective

Input Parameters:
- `sf`  - star forest
- `flg` - `PETSC_TRUE` to sort, `PETSC_FALSE` to skip sorting (false has a lower setup cost, but is non-deterministic)

Level: advanced

See also: `PetscSF`, `PetscSFType`, `PetscSFGatherBegin()`, `PetscSFScatterBegin()`

# External Links
$(_doc_external("PetscSF/PetscSFSetRankOrder"))
"""
function PetscSFSetRankOrder(petsclib::PetscLibType, sf::PetscSF, flg::PetscBool)
    error("PetscSFSetRankOrder: no generated method for these argument types")
end

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
	PetscSFSetType(petsclib::PetscLibType, sf::PetscSF, type::String) 
Set the `PetscSF` communication implementation

Collective

Input Parameters:
- `sf`   - the `PetscSF` context
- `type` - a known method
``
PETSCSFWINDOW - MPI-2/3 one-sided
PETSCSFBASIC - basic implementation using MPI-1 two-sided
``

Options Database Key:
- `-sf_type (basic|window|neighbor)` - Sets the method; see `PetscSFType`

Level: intermediate

See also: `PetscSF`, `PetscSFType`, `PetscSFCreate()`

# External Links
$(_doc_external("PetscSF/PetscSFSetType"))
"""
function PetscSFSetType(petsclib::PetscLibType, sf::PetscSF, type::String)
    error("PetscSFSetType: no generated method for these argument types")
end

@for_petsc function PetscSFSetType(petsclib::$UnionPetscLib, sf::PetscSF, type::String )

    @chk ccall(
               (:PetscSFSetType, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFType),
               sf, type,
              )


	return nothing
end 

"""
	PetscSFSetUp(petsclib::PetscLibType, sf::PetscSF) 
set up communication structures for a `PetscSF`, after this is done it may be used to perform communication

Collective

Input Parameter:
- `sf` - star forest communication object

Level: beginner

See also: `PetscSF`, `PetscSFType`, `PetscSFSetFromOptions()`, `PetscSFSetType()`

# External Links
$(_doc_external("PetscSF/PetscSFSetUp"))
"""
function PetscSFSetUp(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFSetUp: no generated method for these argument types")
end

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
	PetscSFSetUpRanks(petsclib::PetscLibType, sf::PetscSF, dgroup::MPI_Group) 
Set up data structures associated with MPI ranks; this is for internal use by `PetscSF` implementations.

Collective

Input Parameters:
- `sf`     - `PetscSF` to set up; `PetscSFSetGraph()` must have been called
- `dgroup` - `MPI_Group` of ranks to be distinguished (e.g., for self or shared memory exchange)

Level: developer

See also: `PetscSF`, `PetscSFGetRootRanks()`

# External Links
$(_doc_external("PetscSF/PetscSFSetUpRanks"))
"""
function PetscSFSetUpRanks(petsclib::PetscLibType, sf::PetscSF, dgroup::MPI_Group)
    error("PetscSFSetUpRanks: no generated method for these argument types")
end

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
	PetscSFView(petsclib::PetscLibType, sf::PetscSF, viewer::PetscViewer) 
view a star forest

Collective

Input Parameters:
- `sf`     - star forest
- `viewer` - viewer to display graph, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: beginner

See also: `PetscSF`, `PetscViewer`, `PetscSFCreate()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("PetscSF/PetscSFView"))
"""
function PetscSFView(petsclib::PetscLibType, sf::PetscSF, viewer::PetscViewer)
    error("PetscSFView: no generated method for these argument types")
end

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
	PetscSFViewFromOptions(petsclib::PetscLibType, A::PetscSF, obj, name::String) 
View a `PetscSF` based on arguments in the options database

Collective

Input Parameters:
- `A`    - the star forest
- `obj`  - Optional object that provides the prefix for the option names
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `PetscSF`, `PetscSFView`, `PetscObjectViewFromOptions()`, `PetscSFCreate()`

# External Links
$(_doc_external("PetscSF/PetscSFViewFromOptions"))
"""
function PetscSFViewFromOptions(petsclib::PetscLibType, A::PetscSF, obj, name::String)
    error("PetscSFViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscSFViewFromOptions(petsclib::$UnionPetscLib, A::PetscSF, obj, name::String )

    @chk ccall(
               (:PetscSFViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	flavor::PetscSFWindowFlavorType = PetscSFWindowGetFlavorType(petsclib::PetscLibType, sf::PetscSF) 
Get  `PETSCSFWINDOW` flavor type for `PetscSF` communication

Logically Collective

Input Parameter:
- `sf` - star forest for communication of type `PETSCSFWINDOW`

Output Parameter:
- `flavor` - flavor type

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetFlavorType()`

# External Links
$(_doc_external("PetscSF/PetscSFWindowGetFlavorType"))
"""
function PetscSFWindowGetFlavorType(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFWindowGetFlavorType: no generated method for these argument types")
end

@for_petsc function PetscSFWindowGetFlavorType(petsclib::$UnionPetscLib, sf::PetscSF )
	flavor_ = Ref{PetscSFWindowFlavorType}()

    @chk ccall(
               (:PetscSFWindowGetFlavorType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFWindowFlavorType}),
               sf, flavor_,
              )

	flavor = flavor_[]

	return flavor
end 

"""
	info::MPI_Info = PetscSFWindowGetInfo(petsclib::PetscLibType, sf::PetscSF) 
Get the `MPI_Info` handle used for windows allocation

Logically Collective

Input Parameter:
- `sf` - star forest for communication

Output Parameter:
- `info` - `MPI_Info` handle

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetInfo()`

# External Links
$(_doc_external("PetscSF/PetscSFWindowGetInfo"))
"""
function PetscSFWindowGetInfo(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFWindowGetInfo: no generated method for these argument types")
end

@for_petsc function PetscSFWindowGetInfo(petsclib::$UnionPetscLib, sf::PetscSF )
	info_ = Ref{MPI_Info}()

    @chk ccall(
               (:PetscSFWindowGetInfo, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{MPI_Info}),
               sf, info_,
              )

	info = info_[]

	return info
end 

"""
	sync::PetscSFWindowSyncType = PetscSFWindowGetSyncType(petsclib::PetscLibType, sf::PetscSF) 
Get synchronization type for `PetscSF` communication of type `PETSCSFWINDOW`

Logically Collective

Input Parameter:
- `sf` - star forest for communication

Output Parameter:
- `sync` - synchronization type

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowSetSyncType()`, `PetscSFWindowSyncType`

# External Links
$(_doc_external("PetscSF/PetscSFWindowGetSyncType"))
"""
function PetscSFWindowGetSyncType(petsclib::PetscLibType, sf::PetscSF)
    error("PetscSFWindowGetSyncType: no generated method for these argument types")
end

@for_petsc function PetscSFWindowGetSyncType(petsclib::$UnionPetscLib, sf::PetscSF )
	sync_ = Ref{PetscSFWindowSyncType}()

    @chk ccall(
               (:PetscSFWindowGetSyncType, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscSFWindowSyncType}),
               sf, sync_,
              )

	sync = sync_[]

	return sync
end 

"""
	PetscSFWindowSetFlavorType(petsclib::PetscLibType, sf::PetscSF, flavor::PetscSFWindowFlavorType) 
Set flavor type for `MPI_Win` creation

Logically Collective

Input Parameters:
- `sf`     - star forest for communication of type `PETSCSFWINDOW`
- `flavor` - flavor type

Options Database Key:
- `-sf_window_flavor flavor` - sets the flavor type CREATE, DYNAMIC, ALLOCATE or SHARED (see `PetscSFWindowFlavorType`)

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetFlavorType()`

# External Links
$(_doc_external("PetscSF/PetscSFWindowSetFlavorType"))
"""
function PetscSFWindowSetFlavorType(petsclib::PetscLibType, sf::PetscSF, flavor::PetscSFWindowFlavorType)
    error("PetscSFWindowSetFlavorType: no generated method for these argument types")
end

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
	PetscSFWindowSetInfo(petsclib::PetscLibType, sf::PetscSF, info::MPI_Info) 
Set the `MPI_Info` handle that will be used for subsequent windows allocation

Logically Collective

Input Parameters:
- `sf`   - star forest for communication
- `info` - `MPI_Info` handle

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetInfo()`

# External Links
$(_doc_external("PetscSF/PetscSFWindowSetInfo"))
"""
function PetscSFWindowSetInfo(petsclib::PetscLibType, sf::PetscSF, info::MPI_Info)
    error("PetscSFWindowSetInfo: no generated method for these argument types")
end

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
	PetscSFWindowSetSyncType(petsclib::PetscLibType, sf::PetscSF, sync::PetscSFWindowSyncType) 
Set synchronization type for `PetscSF` communication of type  `PETSCSFWINDOW`

Logically Collective

Input Parameters:
- `sf`   - star forest for communication
- `sync` - synchronization type

Options Database Key:
- `-sf_window_sync sync` - sets the synchronization type FENCE, LOCK, or ACTIVE (see `PetscSFWindowSyncType`)

Level: advanced

See also: `PetscSF`, `PETSCSFWINDOW`, `PetscSFSetFromOptions()`, `PetscSFWindowGetSyncType()`, `PetscSFWindowSyncType`

# External Links
$(_doc_external("PetscSF/PetscSFWindowSetSyncType"))
"""
function PetscSFWindowSetSyncType(petsclib::PetscLibType, sf::PetscSF, sync::PetscSFWindowSyncType)
    error("PetscSFWindowSetSyncType: no generated method for these argument types")
end

@for_petsc function PetscSFWindowSetSyncType(petsclib::$UnionPetscLib, sf::PetscSF, sync::PetscSFWindowSyncType )

    @chk ccall(
               (:PetscSFWindowSetSyncType, $petsc_library),
               PetscErrorCode,
               (PetscSF, PetscSFWindowSyncType),
               sf, sync,
              )


	return nothing
end 

