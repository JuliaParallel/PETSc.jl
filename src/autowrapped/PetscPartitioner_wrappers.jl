mutable struct _n_PetscPartitioner end
const PetscPartitioner = Ptr{_n_PetscPartitioner}

"""
	PetscPartitionerRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new PetscPartitioner implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: `PetscPartitionerRegisterAll()`


# External Links
$(_doc_external("Mat/PetscPartitionerRegister"))
"""
function PetscPartitionerRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscPartitionerRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscPartitionerRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscPartitionerFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the PetscPartitioner package.
It is called from PetscFinalize().

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Mat/PetscPartitionerFinalizePackage"))
"""
function PetscPartitionerFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscPartitionerFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPartitionerFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscPartitionerInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the PetscPartitioner package.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Mat/PetscPartitionerInitializePackage"))
"""
function PetscPartitionerInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscPartitionerInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPartitionerInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscPartitionerSetType(petsclib::PetscLibType,part::PetscPartitioner, name::PetscPartitionerType) 
Builds a particular `PetscPartitioner`

Collective

Input Parameters:
- `part` - The `PetscPartitioner` object
- `name` - The kind of partitioner

Options Database Key:
- `-petscpartitioner_type <type>` - Sets the `PetscPartitioner` type

Level: intermediate

-seealso: `PetscPartitionerGetType()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerSetType"))
"""
function PetscPartitionerSetType(petsclib::PetscLibType, part::PetscPartitioner, name::PetscPartitionerType) end

@for_petsc function PetscPartitionerSetType(petsclib::$UnionPetscLib, part::PetscPartitioner, name::PetscPartitionerType )

    @chk ccall(
               (:PetscPartitionerSetType, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, PetscPartitionerType),
               part, name,
              )


	return nothing
end 

"""
	name::PetscPartitionerType = PetscPartitionerGetType(petsclib::PetscLibType,part::PetscPartitioner) 
Gets the PetscPartitioner type name (as a string) from the object.

Not Collective

Input Parameter:
- `part` - The PetscPartitioner

Output Parameter:
- `name` - The PetscPartitioner type name

Level: intermediate

-seealso: `PetscPartitionerSetType()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerGetType"))
"""
function PetscPartitionerGetType(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerGetType(petsclib::$UnionPetscLib, part::PetscPartitioner )
	name_ = Ref{PetscPartitionerType}()

    @chk ccall(
               (:PetscPartitionerGetType, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, Ptr{PetscPartitionerType}),
               part, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscPartitionerViewFromOptions(petsclib::PetscLibType,A::PetscPartitioner, obj::PetscObject, name::String) 
View a `PetscPartitioner` object based on options in the options database

Collective

Input Parameters:
- `A`    - the `PetscPartitioner` object
- `obj`  - Optional `PetscObject` that provides the options prefix
- `name` - command line option

Level: intermediate

-seealso: `PetscPartitionerView()`, `PetscObjectViewFromOptions()`

# External Links
$(_doc_external("Mat/PetscPartitionerViewFromOptions"))
"""
function PetscPartitionerViewFromOptions(petsclib::PetscLibType, A::PetscPartitioner, obj::PetscObject, name::String) end

@for_petsc function PetscPartitionerViewFromOptions(petsclib::$UnionPetscLib, A::PetscPartitioner, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscPartitionerViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscPartitionerView(petsclib::PetscLibType,part::PetscPartitioner, v::PetscViewer) 
Views a `PetscPartitioner`

Collective

Input Parameters:
- `part` - the `PetscPartitioner` object to view
- `v`    - the viewer

Level: developer

-seealso: `PetscPartitionerDestroy()`

# External Links
$(_doc_external("Mat/PetscPartitionerView"))
"""
function PetscPartitionerView(petsclib::PetscLibType, part::PetscPartitioner, v::PetscViewer) end

@for_petsc function PetscPartitionerView(petsclib::$UnionPetscLib, part::PetscPartitioner, v::PetscViewer )

    @chk ccall(
               (:PetscPartitionerView, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, PetscViewer),
               part, v,
              )


	return nothing
end 

"""
	PetscPartitionerSetFromOptions(petsclib::PetscLibType,part::PetscPartitioner) 
sets parameters in a `PetscPartitioner` from the options database

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to set options for

Options Database Keys:
- `-petscpartitioner_type <type>`        - Sets the `PetscPartitioner` type; use -help for a list of available types
- `-petscpartitioner_use_vertex_weights` - Uses weights associated with the graph vertices
- `-petscpartitioner_view_graph`         - View the graph each time PetscPartitionerPartition is called. Viewer can be customized, see `PetscOptionsCreateViewer()`

Level: developer

-seealso: `PetscPartitionerView()`, `PetscPartitionerSetType()`, `PetscPartitionerPartition()`

# External Links
$(_doc_external("Mat/PetscPartitionerSetFromOptions"))
"""
function PetscPartitionerSetFromOptions(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerSetFromOptions(petsclib::$UnionPetscLib, part::PetscPartitioner )

    @chk ccall(
               (:PetscPartitionerSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner,),
               part,
              )


	return nothing
end 

"""
	PetscPartitionerSetUp(petsclib::PetscLibType,part::PetscPartitioner) 
Construct data structures for the `PetscPartitioner`

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to setup

Level: developer

-seealso: `PetscPartitionerView()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("Mat/PetscPartitionerSetUp"))
"""
function PetscPartitionerSetUp(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerSetUp(petsclib::$UnionPetscLib, part::PetscPartitioner )

    @chk ccall(
               (:PetscPartitionerSetUp, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner,),
               part,
              )


	return nothing
end 

"""
	PetscPartitionerReset(petsclib::PetscLibType,part::PetscPartitioner) 
Resets data structures for the `PetscPartitioner`

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to reset

Level: developer

-seealso: `PetscPartitionerSetUp()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("Mat/PetscPartitionerReset"))
"""
function PetscPartitionerReset(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerReset(petsclib::$UnionPetscLib, part::PetscPartitioner )

    @chk ccall(
               (:PetscPartitionerReset, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner,),
               part,
              )


	return nothing
end 

"""
	PetscPartitionerDestroy(petsclib::PetscLibType,part::PetscPartitioner) 
Destroys a `PetscPartitioner` object

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to destroy

Level: developer

-seealso: `PetscPartitionerView()`

# External Links
$(_doc_external("Mat/PetscPartitionerDestroy"))
"""
function PetscPartitionerDestroy(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerDestroy(petsclib::$UnionPetscLib, part::PetscPartitioner )

    @chk ccall(
               (:PetscPartitionerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscPartitioner},),
               part,
              )


	return nothing
end 

"""
	PetscPartitionerPartition(petsclib::PetscLibType,part::PetscPartitioner, nparts::PetscInt, numVertices::PetscInt, start::Vector{PetscInt}, adjacency::Vector{PetscInt}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection, partition::IS) 
Partition a graph

Collective

Input Parameters:
- `part`          - The `PetscPartitioner`
- `nparts`        - Number of partitions
- `numVertices`   - Number of vertices in the local part of the graph
- `start`         - row pointers for the local part of the graph (CSR style)
- `adjacency`     - adjacency list (CSR style)
- `vertexSection` - PetscSection describing the absolute weight of each local vertex (can be `NULL`)
- `edgeSection`   - PetscSection describing the absolute weight of each local edge (can be `NULL`)
- `targetSection` - PetscSection describing the absolute weight of each partition (can be `NULL`)

Output Parameters:
- `partSection` - The `PetscSection` giving the division of points by partition
- `partition`   - The list of points by partition

Options Database Keys:
- `-petscpartitioner_view`       - View the partitioner information
- `-petscpartitioner_view_graph` - View the graph we are partitioning

Level: developer

-seealso: `PetscPartitionerCreate()`, `PetscPartitionerSetType()`, `PetscSectionCreate()`, `PetscSectionSetChart()`, `PetscSectionSetDof()`

# External Links
$(_doc_external("Mat/PetscPartitionerPartition"))
"""
function PetscPartitionerPartition(petsclib::PetscLibType, part::PetscPartitioner, nparts::PetscInt, numVertices::PetscInt, start::Vector{PetscInt}, adjacency::Vector{PetscInt}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection, partition::IS) end

@for_petsc function PetscPartitionerPartition(petsclib::$UnionPetscLib, part::PetscPartitioner, nparts::$PetscInt, numVertices::$PetscInt, start::Vector{$PetscInt}, adjacency::Vector{$PetscInt}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection, partition::IS )
	partition_ = Ref(partition.ptr)

    @chk ccall(
               (:PetscPartitionerPartition, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, PetscSection, PetscSection, PetscSection, PetscSection, Ptr{CIS}),
               part, nparts, numVertices, start, adjacency, vertexSection, edgeSection, targetSection, partSection, partition_,
              )

	partition.ptr = C_NULL

	return nothing
end 

"""
	part::PetscPartitioner = PetscPartitionerCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscPartitioner` object. The type can then be set with `PetscPartitionerSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscPartitioner` object

Output Parameter:
- `part` - The `PetscPartitioner` object

Level: beginner

-seealso: `PetscPartitionerSetType()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("Mat/PetscPartitionerCreate"))
"""
function PetscPartitionerCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscPartitionerCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	part_ = Ref{PetscPartitioner}()

    @chk ccall(
               (:PetscPartitionerCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscPartitioner}),
               comm, part_,
              )

	part = part_[]

	return part
end 

"""
	PetscPartitionerShellSetPartition(petsclib::PetscLibType,part::PetscPartitioner, size::PetscInt, sizes::Vector{PetscInt}, points::Vector{PetscInt}) 
Set an artificial partition for a mesh

Collective

Input Parameters:
- `part`   - The `PetscPartitioner`
- `size`   - The number of partitions
- `sizes`  - array of length size (or `NULL`) providing the number of points in each partition
- `points` - array of length sum(sizes) (may be `NULL` iff sizes is `NULL`), a permutation of the points that groups those assigned to each partition in order (i.e., partition 0 first, partition 1 next, etc.)

Level: developer

-seealso: `DMPlexDistribute()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerShellSetPartition"))
"""
function PetscPartitionerShellSetPartition(petsclib::PetscLibType, part::PetscPartitioner, size::PetscInt, sizes::Vector{PetscInt}, points::Vector{PetscInt}) end

@for_petsc function PetscPartitionerShellSetPartition(petsclib::$UnionPetscLib, part::PetscPartitioner, size::$PetscInt, sizes::Vector{$PetscInt}, points::Vector{$PetscInt} )

    @chk ccall(
               (:PetscPartitionerShellSetPartition, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               part, size, sizes, points,
              )


	return nothing
end 

"""
	PetscPartitionerShellSetRandom(petsclib::PetscLibType,part::PetscPartitioner, random::PetscBool) 
Set the flag to use a random partition

Collective

Input Parameters:
- `part`   - The `PetscPartitioner`
- `random` - The flag to use a random partition

Level: intermediate

-seealso: `PetscPartitionerShellGetRandom()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerShellSetRandom"))
"""
function PetscPartitionerShellSetRandom(petsclib::PetscLibType, part::PetscPartitioner, random::PetscBool) end

@for_petsc function PetscPartitionerShellSetRandom(petsclib::$UnionPetscLib, part::PetscPartitioner, random::PetscBool )

    @chk ccall(
               (:PetscPartitionerShellSetRandom, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, PetscBool),
               part, random,
              )


	return nothing
end 

"""
	random::PetscBool = PetscPartitionerShellGetRandom(petsclib::PetscLibType,part::PetscPartitioner) 
get the flag to use a random partition

Collective

Input Parameter:
- `part` - The `PetscPartitioner`

Output Parameter:
- `random` - The flag to use a random partition

Level: intermediate

-seealso: `PetscPartitionerShellSetRandom()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerShellGetRandom"))
"""
function PetscPartitionerShellGetRandom(petsclib::PetscLibType, part::PetscPartitioner) end

@for_petsc function PetscPartitionerShellGetRandom(petsclib::$UnionPetscLib, part::PetscPartitioner )
	random_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscPartitionerShellGetRandom, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, Ptr{PetscBool}),
               part, random_,
              )

	random = random_[]

	return random
end 

"""
	PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::PetscLibType,part::PetscPartitioner, mp::MatPartitioning) 
Get a `MatPartitioning` instance wrapped by this `PetscPartitioner`.

Not Collective

Input Parameter:
- `part` - The `PetscPartitioner`

Output Parameter:
- `mp` - The `MatPartitioning`

Level: developer

-seealso: `DMPlexDistribute()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("Mat/PetscPartitionerMatPartitioningGetMatPartitioning"))
"""
function PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::PetscLibType, part::PetscPartitioner, mp::MatPartitioning) end

@for_petsc function PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::$UnionPetscLib, part::PetscPartitioner, mp::MatPartitioning )

    @chk ccall(
               (:PetscPartitionerMatPartitioningGetMatPartitioning, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, Ptr{MatPartitioning}),
               part, mp,
              )


	return nothing
end 

"""
	PetscPartitionerDMPlexPartition(petsclib::PetscLibType,part::PetscPartitioner, dm::PetscDM, targetSection::PetscSection, partSection::PetscSection, partition::IS) 
Create a non

Collective

Input Parameters:
- `part`          - The `PetscPartitioner`
- `targetSection` - The `PetscSection` describing the absolute weight of each partition (can be `NULL`)
- `dm`            - The mesh `DM`

Output Parameters:
- `partSection` - The `PetscSection` giving the division of points by partition
- `partition`   - The list of points by partition

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `PetscPartitioner`, `PetscSection`, `DMPlexDistribute()`, `PetscPartitionerCreate()`, `PetscSectionCreate()`,
`PetscSectionSetChart()`, `PetscPartitionerPartition()`

# External Links
$(_doc_external("Dm/PetscPartitionerDMPlexPartition"))
"""
function PetscPartitionerDMPlexPartition(petsclib::PetscLibType, part::PetscPartitioner, dm::PetscDM, targetSection::PetscSection, partSection::PetscSection, partition::IS) end

@for_petsc function PetscPartitionerDMPlexPartition(petsclib::$UnionPetscLib, part::PetscPartitioner, dm::PetscDM, targetSection::PetscSection, partSection::PetscSection, partition::IS )
	partition_ = Ref(partition.ptr)

    @chk ccall(
               (:PetscPartitionerDMPlexPartition, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, CDM, PetscSection, PetscSection, Ptr{CIS}),
               part, dm, targetSection, partSection, partition_,
              )

	partition.ptr = C_NULL

	return nothing
end 

