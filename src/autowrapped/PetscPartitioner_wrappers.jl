"""
	part::PetscPartitioner = PetscPartitionerCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates an empty `PetscPartitioner` object. The type can then be set with `PetscPartitionerSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscPartitioner` object

Output Parameter:
- `part` - The `PetscPartitioner` object

Level: beginner

See also: `PetscPartitionerSetType()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerCreate"))
"""
function PetscPartitionerCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscPartitionerCreate: no generated method for these argument types")
end

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
	partition::IS = PetscPartitionerDMPlexPartition(petsclib::PetscLibType, part::PetscPartitioner, dm::AbstractPetscDM, targetSection::PetscSection, partSection::PetscSection) 
Create a non-overlapping partition of the cells in the mesh

Collective

Input Parameters:
- `part`          - The `PetscPartitioner`
- `targetSection` - The `PetscSection` describing the absolute weight of each partition (can be `NULL`)
- `dm`            - The mesh `DM`

Output Parameters:
- `partSection` - The `PetscSection` giving the division of points by partition
- `partition`   - The list of points by partition

Level: developer

See also: `DM`, `DMPLEX`, `PetscPartitioner`, `PetscSection`, `DMPlexDistribute()`, `PetscPartitionerCreate()`, `PetscSectionCreate()`,
`PetscSectionSetChart()`, `PetscPartitionerPartition()`

# External Links
$(_doc_external("DMPlex/PetscPartitionerDMPlexPartition"))
"""
function PetscPartitionerDMPlexPartition(petsclib::PetscLibType, part::PetscPartitioner, dm::AbstractPetscDM, targetSection::PetscSection, partSection::PetscSection)
    error("PetscPartitionerDMPlexPartition: no generated method for these argument types")
end

@for_petsc function PetscPartitionerDMPlexPartition(petsclib::$UnionPetscLib, part::PetscPartitioner, dm::AbstractPetscDM, targetSection::PetscSection, partSection::PetscSection )
	partition_ = Ref{CIS}()

    @chk ccall(
               (:PetscPartitionerDMPlexPartition, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, CDM, PetscSection, PetscSection, Ptr{CIS}),
               part, dm, targetSection, partSection, partition_,
              )

	partition = IS(partition_[], petsclib)

	return partition
end 

"""
	PetscPartitionerDestroy(petsclib::PetscLibType, part::Union{PetscPartitioner, Ref{PetscPartitioner}}) 
Destroys a `PetscPartitioner` object

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to destroy

Level: developer

See also: `PetscPartitionerView()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerDestroy"))
"""
function PetscPartitionerDestroy(petsclib::PetscLibType, part::Union{PetscPartitioner, Ref{PetscPartitioner}})
    error("PetscPartitionerDestroy: no generated method for these argument types")
end

@for_petsc function PetscPartitionerDestroy(petsclib::$UnionPetscLib, part::Union{PetscPartitioner, Ref{PetscPartitioner}} )
	part_ = part isa Base.RefValue ? part : Ref{PetscPartitioner}(part)

    @chk ccall(
               (:PetscPartitionerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscPartitioner},),
               part_,
              )


	return nothing
end 

"""
	PetscPartitionerFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the PetscPartitioner package.
It is called from PetscFinalize().

Level: developer

See also: `PetscInitialize()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerFinalizePackage"))
"""
function PetscPartitionerFinalizePackage(petsclib::PetscLibType)
    error("PetscPartitionerFinalizePackage: no generated method for these argument types")
end

@for_petsc function PetscPartitionerFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPartitionerFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	name::PetscPartitionerType = PetscPartitionerGetType(petsclib::PetscLibType, part::PetscPartitioner) 
Gets the PetscPartitioner type name (as a string) from the object.

Not Collective

Input Parameter:
- `part` - The PetscPartitioner

Output Parameter:
- `name` - The PetscPartitioner type name

Level: intermediate

See also: `PetscPartitionerSetType()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerGetType"))
"""
function PetscPartitionerGetType(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerGetType: no generated method for these argument types")
end

@for_petsc function PetscPartitionerGetType(petsclib::$UnionPetscLib, part::PetscPartitioner )
	name_ = Ref{PetscPartitionerType}()

    @chk ccall(
               (:PetscPartitionerGetType, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, Ptr{PetscPartitionerType}),
               part, name_,
              )

	name = name_[] == C_NULL ? "" : unsafe_string(name_[])

	return name
end 

"""
	PetscPartitionerInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the PetscPartitioner package.

Level: developer

See also: `PetscInitialize()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerInitializePackage"))
"""
function PetscPartitionerInitializePackage(petsclib::PetscLibType)
    error("PetscPartitionerInitializePackage: no generated method for these argument types")
end

@for_petsc function PetscPartitionerInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPartitionerInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	mp::MatPartitioning = PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::PetscLibType, part::PetscPartitioner) 
Get a `MatPartitioning` instance wrapped by this `PetscPartitioner`.

Not Collective

Input Parameter:
- `part` - The `PetscPartitioner`

Output Parameter:
- `mp` - The `MatPartitioning`

Level: developer

See also: `DMPlexDistribute()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerMatPartitioningGetMatPartitioning"))
"""
function PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerMatPartitioningGetMatPartitioning: no generated method for these argument types")
end

@for_petsc function PetscPartitionerMatPartitioningGetMatPartitioning(petsclib::$UnionPetscLib, part::PetscPartitioner )
	mp_ = Ref{MatPartitioning}()

    @chk ccall(
               (:PetscPartitionerMatPartitioningGetMatPartitioning, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, Ptr{MatPartitioning}),
               part, mp_,
              )

	mp = mp_[]

	return mp
end 

"""
	partition::IS = PetscPartitionerPartition(petsclib::PetscLibType, part::PetscPartitioner, nparts::PetscInt, numVertices::PetscInt, start::Vector{PetscInt}, adjacency::Vector{PetscInt}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection) 
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

See also: `PetscPartitionerCreate()`, `PetscPartitionerSetType()`, `PetscSectionCreate()`, `PetscSectionSetChart()`, `PetscSectionSetDof()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerPartition"))
"""
function PetscPartitionerPartition(petsclib::PetscLibType, part::PetscPartitioner, nparts::Integer, numVertices::Integer, start::AbstractVector{<:Number}, adjacency::AbstractVector{<:Number}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection)
    error("PetscPartitionerPartition: no generated method for these argument types")
end

@for_petsc function PetscPartitionerPartition(petsclib::$UnionPetscLib, part::PetscPartitioner, nparts::$PetscInt, numVertices::$PetscInt, start::Vector{$PetscInt}, adjacency::Vector{$PetscInt}, vertexSection::PetscSection, edgeSection::PetscSection, targetSection::PetscSection, partSection::PetscSection )
	partition_ = Ref{CIS}()

    @chk ccall(
               (:PetscPartitionerPartition, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, PetscSection, PetscSection, PetscSection, PetscSection, Ptr{CIS}),
               part, nparts, numVertices, start, adjacency, vertexSection, edgeSection, targetSection, partSection, partition_,
              )

	partition = IS(partition_[], petsclib)

	return partition
end 

"""
	PetscPartitionerRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a new PetscPartitioner implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

See also: `PetscPartitionerRegisterAll()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerRegister"))
"""
function PetscPartitionerRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscPartitionerRegister: no generated method for these argument types")
end

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
	PetscPartitionerReset(petsclib::PetscLibType, part::PetscPartitioner) 
Resets data structures for the `PetscPartitioner`

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to reset

Level: developer

See also: `PetscPartitionerSetUp()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerReset"))
"""
function PetscPartitionerReset(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerReset: no generated method for these argument types")
end

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
	PetscPartitionerSetFromOptions(petsclib::PetscLibType, part::PetscPartitioner) 
sets parameters in a `PetscPartitioner` from the options database

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to set options for

Options Database Keys:
- `-petscpartitioner_type type`          - Sets the `PetscPartitionerType`
- `-petscpartitioner_use_vertex_weights` - Uses weights associated with the graph vertices
- `-petscpartitioner_view_graph`         - View the graph each time PetscPartitionerPartition is called. Viewer can be customized, see `PetscOptionsCreateViewer()`

Level: developer

See also: `PetscPartitionerView()`, `PetscPartitionerSetType()`, `PetscPartitionerPartition()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerSetFromOptions"))
"""
function PetscPartitionerSetFromOptions(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerSetFromOptions: no generated method for these argument types")
end

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
	PetscPartitionerSetType(petsclib::PetscLibType, part::PetscPartitioner, name::PetscPartitionerType) 
Builds a particular `PetscPartitioner`

Collective

Input Parameters:
- `part` - The `PetscPartitioner` object
- `name` - The kind of partitioner

Options Database Key:
- `-petscpartitioner_type type` - Sets the `PetscPartitionerType`

Level: intermediate

See also: `PetscPartitionerGetType()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerSetType"))
"""
function PetscPartitionerSetType(petsclib::PetscLibType, part::PetscPartitioner, name::PetscPartitionerType)
    error("PetscPartitionerSetType: no generated method for these argument types")
end

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
	PetscPartitionerSetUp(petsclib::PetscLibType, part::PetscPartitioner) 
Construct data structures for the `PetscPartitioner`

Collective

Input Parameter:
- `part` - the `PetscPartitioner` object to setup

Level: developer

See also: `PetscPartitionerView()`, `PetscPartitionerDestroy()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerSetUp"))
"""
function PetscPartitionerSetUp(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerSetUp: no generated method for these argument types")
end

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
	random::PetscBool = PetscPartitionerShellGetRandom(petsclib::PetscLibType, part::PetscPartitioner) 
get the flag to use a random partition

Collective

Input Parameter:
- `part` - The `PetscPartitioner`

Output Parameter:
- `random` - The flag to use a random partition

Level: intermediate

See also: `PetscPartitionerShellSetRandom()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerShellGetRandom"))
"""
function PetscPartitionerShellGetRandom(petsclib::PetscLibType, part::PetscPartitioner)
    error("PetscPartitionerShellGetRandom: no generated method for these argument types")
end

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
	PetscPartitionerShellSetPartition(petsclib::PetscLibType, part::PetscPartitioner, size::PetscInt, sizes::Vector{PetscInt}, points::Vector{PetscInt}) 
Set an artificial partition for a mesh

Collective

Input Parameters:
- `part`   - The `PetscPartitioner`
- `size`   - The number of partitions
- `sizes`  - array of length size (or `NULL`) providing the number of points in each partition
- `points` - array of length sum(sizes) (may be `NULL` iff sizes is `NULL`), a permutation of the points that groups those assigned to each partition in order (i.e., partition 0 first, partition 1 next, etc.)

Level: developer

See also: `DMPlexDistribute()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerShellSetPartition"))
"""
function PetscPartitionerShellSetPartition(petsclib::PetscLibType, part::PetscPartitioner, size::Integer, sizes::AbstractVector{<:Number}, points::AbstractVector{<:Number})
    error("PetscPartitionerShellSetPartition: no generated method for these argument types")
end

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
	PetscPartitionerShellSetRandom(petsclib::PetscLibType, part::PetscPartitioner, random::PetscBool) 
Set the flag to use a random partition

Collective

Input Parameters:
- `part`   - The `PetscPartitioner`
- `random` - The flag to use a random partition

Level: intermediate

See also: `PetscPartitionerShellGetRandom()`, `PetscPartitionerCreate()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerShellSetRandom"))
"""
function PetscPartitionerShellSetRandom(petsclib::PetscLibType, part::PetscPartitioner, random::PetscBool)
    error("PetscPartitionerShellSetRandom: no generated method for these argument types")
end

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
	PetscPartitionerView(petsclib::PetscLibType, part::PetscPartitioner, v::PetscViewer) 
Views a `PetscPartitioner`

Collective

Input Parameters:
- `part` - the `PetscPartitioner` object to view
- `v`    - the viewer

Level: developer

See also: `PetscPartitionerDestroy()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerView"))
"""
function PetscPartitionerView(petsclib::PetscLibType, part::PetscPartitioner, v::PetscViewer)
    error("PetscPartitionerView: no generated method for these argument types")
end

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
	PetscPartitionerViewFromOptions(petsclib::PetscLibType, A::PetscPartitioner, obj, name::String) 
View a `PetscPartitioner` object based on options in the options database

Collective

Input Parameters:
- `A`    - the `PetscPartitioner` object
- `obj`  - Optional `PetscObject` that provides the options prefix
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `PetscPartitionerView()`, `PetscObjectViewFromOptions()`

# External Links
$(_doc_external("MatGraphOperations/PetscPartitionerViewFromOptions"))
"""
function PetscPartitionerViewFromOptions(petsclib::PetscLibType, A::PetscPartitioner, obj, name::String)
    error("PetscPartitionerViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscPartitionerViewFromOptions(petsclib::$UnionPetscLib, A::PetscPartitioner, obj, name::String )

    @chk ccall(
               (:PetscPartitionerViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscPartitioner, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

