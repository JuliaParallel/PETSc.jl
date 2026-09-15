"""
	iscoloring::ISColoring = ISColoringCreate(petsclib::PetscLibType,comm::MPI_Comm, ncolors::PetscInt, n::PetscInt, colors::Vector{ISColoringValue}, mode::PetscCopyMode) 
Generates an `ISColoring` context from lists (provided by each MPI process) of colors for each node.

Collective

Input Parameters:
- `comm`    - communicator for the processors creating the coloring
- `ncolors` - max color value
- `n`       - number of nodes on this processor
- `colors`  - array containing the colors for this MPI rank, color numbers begin at 0, for each local node
- `mode`    - see `PetscCopyMode` for meaning of this flag.

Output Parameter:
- `iscoloring` - the resulting coloring data structure

Options Database Key:
- `-is_coloring_view` - Activates `ISColoringView()`

Level: advanced

-seealso: `ISColoring`, `ISColoringValue`, `MatColoringCreate()`, `ISColoringView()`, `ISColoringDestroy()`, `ISColoringSetType()`

# External Links
$(_doc_external("IS/ISColoringCreate"))
"""
function ISColoringCreate(petsclib::PetscLibType, comm::MPI_Comm, ncolors::Integer, n::Integer, colors::Vector{ISColoringValue}, mode::PetscCopyMode)
    error("ISColoringCreate: no generated method for these argument types")
end

@for_petsc function ISColoringCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, ncolors::$PetscInt, n::$PetscInt, colors::Vector{ISColoringValue}, mode::PetscCopyMode )
	iscoloring_ = Ref{ISColoring}()

    @chk ccall(
               (:ISColoringCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{ISColoringValue}, PetscCopyMode, Ptr{ISColoring}),
               comm, ncolors, n, colors, mode, iscoloring_,
              )

	iscoloring = iscoloring_[]

	return iscoloring
end 

"""
	ISColoringDestroy(petsclib::PetscLibType,iscoloring::Union{ISColoring, Ref{ISColoring}}) 
Destroys an `ISColoring` coloring context.

Collective

Input Parameter:
- `iscoloring` - the coloring context

Level: advanced

-seealso: `ISColoring`, `ISColoringView()`, `MatColoring`

# External Links
$(_doc_external("IS/ISColoringDestroy"))
"""
function ISColoringDestroy(petsclib::PetscLibType, iscoloring::Union{ISColoring, Ref{ISColoring}})
    error("ISColoringDestroy: no generated method for these argument types")
end

@for_petsc function ISColoringDestroy(petsclib::$UnionPetscLib, iscoloring::Union{ISColoring, Ref{ISColoring}} )
	iscoloring_ = iscoloring isa Base.RefValue ? iscoloring : Ref{ISColoring}(iscoloring)

    @chk ccall(
               (:ISColoringDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{ISColoring},),
               iscoloring_,
              )


	return nothing
end 

# override for ISColoringGetColors; C signature: ISColoringGetColors(ISColoring iscoloring, PetscInt* n, PetscInt* nc, ISColoringValue** colors)
"""
	n::PetscInt, nc::PetscInt, colors::Vector{UInt16} = ISColoringGetColors(petsclib::PetscLibType, iscoloring::ISColoring)
Returns an array with the color for each local node

Not Collective

Input Parameter:
- `iscoloring` - the coloring context

Output Parameters:
- `n`      - number of nodes (DOFs)
- `nc`     - number of colors
- `colors` - copy of the color array (one `UInt16` entry per DOF)

Level: advanced

-seealso: `ISColoring`, `ISColoringValue`, `ISColoringRestoreIS()`, `ISColoringView()`, `ISColoringGetIS()`

# External Links
$(_doc_external("Vec/ISColoringGetColors"))
"""
function ISColoringGetColors(petsclib::PetscLibType, iscoloring::ISColoring) end

@for_petsc function ISColoringGetColors(petsclib::$UnionPetscLib, iscoloring::ISColoring)
    n_  = Ref{$PetscInt}()
    nc_ = Ref{$PetscInt}()
    # PETSc returns a pointer into its own storage; we must not free it.
    colors_ptr_ = Ref{ISColoringValue}(C_NULL)

    @chk ccall(
               (:ISColoringGetColors, $petsc_library),
               PetscErrorCode,
               (ISColoring, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{ISColoringValue}),
               iscoloring, n_, nc_, colors_ptr_,
              )

    n  = n_[]
    nc = nc_[]
    # ISColoringValue is a PETSc opaque pointer alias for `unsigned short *`.
    # Reinterpret as Ptr{UInt16} and copy to a Julia-owned Vector so the
    # caller can safely use the data after ISColoringDestroy.
    colors_raw = Ptr{UInt16}(UInt(colors_ptr_[]))
    colors = copy(unsafe_wrap(Vector{UInt16}, colors_raw, Int(n); own = false))

    return n, nc, colors
end

"""
	nn::PetscInt,isis::Ptr{IS} = ISColoringGetIS(petsclib::PetscLibType,iscoloring::ISColoring, mode::PetscCopyMode) 
Extracts index sets from the coloring context. Each is contains the nodes of one color

Collective

Input Parameters:
- `iscoloring` - the coloring context
- `mode`       - if this value is `PETSC_OWN_POINTER` then the caller owns the pointer and must free the array of `IS` and each `IS` in the array

Output Parameters:
- `nn`   - number of index sets in the coloring context
- `isis` - array of index sets

Level: advanced

-seealso: `ISColoring`, `IS`, `ISColoringRestoreIS()`, `ISColoringView()`, `ISColoringGetColoring()`, `ISColoringGetColors()`

# External Links
$(_doc_external("IS/ISColoringGetIS"))
"""
function ISColoringGetIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode)
    error("ISColoringGetIS: no generated method for these argument types")
end

@for_petsc function ISColoringGetIS(petsclib::$UnionPetscLib, iscoloring::ISColoring, mode::PetscCopyMode )
	nn_ = Ref{$PetscInt}()
	isis_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:ISColoringGetIS, $petsc_library),
               PetscErrorCode,
               (ISColoring, PetscCopyMode, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
               iscoloring, mode, nn_, isis_,
              )

	nn = nn_[]
	isis = isis_[]

	return nn,isis
end 

"""
	type::ISColoringType = ISColoringGetType(petsclib::PetscLibType,coloring::ISColoring) 
gets if the coloring is for the local representation (including ghost points) or the global representation

Collective

Input Parameter:
- `coloring` - the coloring object

Output Parameter:
- `type` - either `IS_COLORING_LOCAL` or `IS_COLORING_GLOBAL`

Level: intermediate

-seealso: `MatFDColoringCreate()`, `ISColoring`, `ISColoringType`, `ISColoringCreate()`, `IS_COLORING_LOCAL`, `IS_COLORING_GLOBAL`, `ISColoringSetType()`

# External Links
$(_doc_external("IS/ISColoringGetType"))
"""
function ISColoringGetType(petsclib::PetscLibType, coloring::ISColoring)
    error("ISColoringGetType: no generated method for these argument types")
end

@for_petsc function ISColoringGetType(petsclib::$UnionPetscLib, coloring::ISColoring )
	type_ = Ref{ISColoringType}()

    @chk ccall(
               (:ISColoringGetType, $petsc_library),
               PetscErrorCode,
               (ISColoring, Ptr{ISColoringType}),
               coloring, type_,
              )

	type = type_[]

	return type
end 

"""
	ISColoringReference(petsclib::PetscLibType,coloring::ISColoring) 
Increases the reference count of an `ISColoring` object by one

Logically collective

Input Parameter:
- `coloring` - the `ISColoring` object

Level: developer

-seealso: `ISColoring`, `ISColoringCreate()`, `ISColoringDestroy()`

# External Links
$(_doc_external("IS/ISColoringReference"))
"""
function ISColoringReference(petsclib::PetscLibType, coloring::ISColoring)
    error("ISColoringReference: no generated method for these argument types")
end

@for_petsc function ISColoringReference(petsclib::$UnionPetscLib, coloring::ISColoring )

    @chk ccall(
               (:ISColoringReference, $petsc_library),
               PetscErrorCode,
               (ISColoring,),
               coloring,
              )


	return nothing
end 

"""
	ISColoringRestoreIS(petsclib::PetscLibType,iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractArray{IS}}) 
Restores the index sets extracted from the coloring context with `ISColoringGetIS()` using `PETSC_USE_POINTER`

Collective

Input Parameters:
- `iscoloring` - the coloring context
- `mode`       - who retains ownership of the is
- `is`         - array of index sets

Level: advanced

-seealso: `ISColoring()`, `IS`, `ISColoringGetIS()`, `ISColoringView()`, `PetscCopyMode`

# External Links
$(_doc_external("IS/ISColoringRestoreIS"))
"""
function ISColoringRestoreIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractArray{IS}})
    error("ISColoringRestoreIS: no generated method for these argument types")
end

@for_petsc function ISColoringRestoreIS(petsclib::$UnionPetscLib, iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractArray{IS}} )
	is_ = Ref{Ptr{CIS}}(is isa Ptr ? is : pointer(is))

    @chk ccall(
               (:ISColoringRestoreIS, $petsc_library),
               PetscErrorCode,
               (ISColoring, PetscCopyMode, Ptr{Ptr{CIS}}),
               iscoloring, mode, is_,
              )


	return nothing
end 

"""
	ISColoringSetType(petsclib::PetscLibType,coloring::ISColoring, type::ISColoringType) 
indicates if the coloring is for the local representation (including ghost points) or the global representation of a `Mat`

Collective

Input Parameters:
- `coloring` - the coloring object
- `type`     - either `IS_COLORING_LOCAL` or `IS_COLORING_GLOBAL`

Level: intermediate

-seealso: `MatFDColoringCreate()`, `ISColoring`, `ISColoringType`, `ISColoringCreate()`, `IS_COLORING_LOCAL`, `IS_COLORING_GLOBAL`, `ISColoringGetType()`

# External Links
$(_doc_external("IS/ISColoringSetType"))
"""
function ISColoringSetType(petsclib::PetscLibType, coloring::ISColoring, type::ISColoringType)
    error("ISColoringSetType: no generated method for these argument types")
end

@for_petsc function ISColoringSetType(petsclib::$UnionPetscLib, coloring::ISColoring, type::ISColoringType )

    @chk ccall(
               (:ISColoringSetType, $petsc_library),
               PetscErrorCode,
               (ISColoring, ISColoringType),
               coloring, type,
              )


	return nothing
end 

"""
	ISColoringValueCast(petsclib::PetscLibType,a::PetscCount, b::ISColoringValue) 

# External Links
$(_doc_external("Vec/ISColoringValueCast"))
"""
function ISColoringValueCast(petsclib::PetscLibType, a::PetscCount, b::ISColoringValue)
    error("ISColoringValueCast: no generated method for these argument types")
end

@for_petsc function ISColoringValueCast(petsclib::$UnionPetscLib, a::PetscCount, b::ISColoringValue )

    @chk ccall(
               (:ISColoringValueCast, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{ISColoringValue}),
               a, b,
              )


	return nothing
end 

"""
	ISColoringView(petsclib::PetscLibType,iscoloring::ISColoring, viewer::PetscViewer) 
Views an `ISColoring` coloring context.

Collective

Input Parameters:
- `iscoloring` - the coloring context
- `viewer`     - the viewer

Level: advanced

-seealso: `ISColoring()`, `ISColoringViewFromOptions()`, `ISColoringDestroy()`, `ISColoringGetIS()`, `MatColoring`

# External Links
$(_doc_external("IS/ISColoringView"))
"""
function ISColoringView(petsclib::PetscLibType, iscoloring::ISColoring, viewer::PetscViewer)
    error("ISColoringView: no generated method for these argument types")
end

@for_petsc function ISColoringView(petsclib::$UnionPetscLib, iscoloring::ISColoring, viewer::PetscViewer )

    @chk ccall(
               (:ISColoringView, $petsc_library),
               PetscErrorCode,
               (ISColoring, PetscViewer),
               iscoloring, viewer,
              )


	return nothing
end 

"""
	ISColoringViewFromOptions(petsclib::PetscLibType,obj::ISColoring, bobj, name::String) 
Processes command line options to determine if/how an `ISColoring` object is to be viewed.

Collective

Input Parameters:
- `obj`  - the `ISColoring` object
- `bobj` - prefix to use for viewing, or `NULL` to use prefix of `mat`
- `name` - option to activate viewing

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

-seealso: `ISColoring`, `ISColoringView()`, `PetscObjectViewFromOptions()`

# External Links
$(_doc_external("IS/ISColoringViewFromOptions"))
"""
function ISColoringViewFromOptions(petsclib::PetscLibType, obj::ISColoring, bobj, name::String)
    error("ISColoringViewFromOptions: no generated method for these argument types")
end

@for_petsc function ISColoringViewFromOptions(petsclib::$UnionPetscLib, obj::ISColoring, bobj, name::String )

    @chk ccall(
               (:ISColoringViewFromOptions, $petsc_library),
               PetscErrorCode,
               (ISColoring, PetscObject, Ptr{Cchar}),
               obj, bobj, name,
              )


	return nothing
end 

"""
	out::Vector{PetscInt} = ISLocalToGlobalMappingApply(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, N::PetscInt, in::Vector{PetscInt}) 
Takes a list of integers in a local numbering
and converts them to the global numbering.

Not Collective

Input Parameters:
- `mapping` - the local to global mapping context
- `N`       - number of integers
- `in`      - input indices in local numbering

Output Parameter:
- `out` - indices in global numbering

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingApplyBlock()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingDestroy()`,
`ISLocalToGlobalMappingApplyIS()`, `AOCreateBasic()`, `AOApplicationToPetsc()`,
`AOPetscToApplication()`, `ISGlobalToLocalMappingApply()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingApply"))
"""
function ISLocalToGlobalMappingApply(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, N::Integer, in::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingApply: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingApply(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, N::$PetscInt, in::Vector{$PetscInt} )
	out = Vector{$PetscInt}(undef, Int(N))

    @chk ccall(
               (:ISLocalToGlobalMappingApply, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mapping, N, in, out,
              )


	return out
end 

"""
	out::Vector{PetscInt} = ISLocalToGlobalMappingApplyBlock(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, N::PetscInt, in::Vector{PetscInt}) 
Takes a list of integers in a local block numbering and converts them to the global block numbering

Not Collective

Input Parameters:
- `mapping` - the local to global mapping context
- `N`       - number of integers
- `in`      - input indices in local block numbering

Output Parameter:
- `out` - indices in global block numbering

Example:
If the index values are {0,1,6,7} set with a call to `ISLocalToGlobalMappingCreate`(`PETSC_COMM_SELF`,2,2,{0,3}) then the mapping applied to 0
(the first block) would produce 0 and the mapping applied to 1 (the second block) would produce 3.

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingApply()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingDestroy()`,
`ISLocalToGlobalMappingApplyIS()`, `AOCreateBasic()`, `AOApplicationToPetsc()`,
`AOPetscToApplication()`, `ISGlobalToLocalMappingApply()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingApplyBlock"))
"""
function ISLocalToGlobalMappingApplyBlock(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, N::Integer, in::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingApplyBlock: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingApplyBlock(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, N::$PetscInt, in::Vector{$PetscInt} )
	out = Vector{$PetscInt}(undef, Int(N))

    @chk ccall(
               (:ISLocalToGlobalMappingApplyBlock, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mapping, N, in, out,
              )


	return out
end 

"""
	newis::IS = ISLocalToGlobalMappingApplyIS(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, is::AbstractIS) 
Creates from an `IS` in the local numbering
a new index set using the global numbering defined in an `ISLocalToGlobalMapping`
context.

Collective

Input Parameters:
- `mapping` - mapping between local and global numbering
- `is`      - index set in local numbering

Output Parameter:
- `newis` - index set in global numbering

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingApply()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingDestroy()`, `ISGlobalToLocalMappingApply()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingApplyIS"))
"""
function ISLocalToGlobalMappingApplyIS(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, is::AbstractIS)
    error("ISLocalToGlobalMappingApplyIS: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingApplyIS(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, is::AbstractIS )
	newis_ = Ref{CIS}()

    @chk ccall(
               (:ISLocalToGlobalMappingApplyIS, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, CIS, Ptr{CIS}),
               mapping, is, newis_,
              )

	newis = IS(newis_[], petsclib)

	return newis
end 

"""
	ltogcat::ISLocalToGlobalMapping = ISLocalToGlobalMappingConcatenate(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, ltogs::Vector{ISLocalToGlobalMapping}) 
Create a new mapping that concatenates a list of mappings

Not Collective

Input Parameters:
- `comm`  - communicator for the new mapping, must contain the communicator of every mapping to concatenate
- `n`     - number of mappings to concatenate
- `ltogs` - local to global mappings

Output Parameter:
- `ltogcat` - new mapping

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingConcatenate"))
"""
function ISLocalToGlobalMappingConcatenate(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, ltogs::Vector{ISLocalToGlobalMapping})
    error("ISLocalToGlobalMappingConcatenate: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingConcatenate(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, ltogs::Vector{ISLocalToGlobalMapping} )
	ltogcat_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:ISLocalToGlobalMappingConcatenate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
               comm, n, ltogs, ltogcat_,
              )

	ltogcat = ltogcat_[]

	return ltogcat
end 

"""
	mapping::ISLocalToGlobalMapping = ISLocalToGlobalMappingCreate(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, indices::Vector{PetscInt}, mode::PetscCopyMode) 
Creates a mapping between a local (0 to n)
ordering and a global parallel ordering.

Not Collective, but communicator may have more than one process

Input Parameters:
- `comm`    - MPI communicator
- `bs`      - the block size
- `n`       - the number of local elements divided by the block size, or equivalently the number of block indices
- `indices` - the global index for each local element, these do not need to be in increasing order (sorted), these values should not be scaled (i.e. multiplied) by the blocksize bs
- `mode`    - see PetscCopyMode

Output Parameter:
- `mapping` - new mapping data structure

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingSetFromOptions()`,
`ISLOCALTOGLOBALMAPPINGBASIC`, `ISLOCALTOGLOBALMAPPINGHASH`,
`ISLocalToGlobalMappingSetType()`, `ISLocalToGlobalMappingType`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingCreate"))
"""
function ISLocalToGlobalMappingCreate(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, indices::AbstractVector{<:Number}, mode::PetscCopyMode)
    error("ISLocalToGlobalMappingCreate: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, indices::Vector{$PetscInt}, mode::PetscCopyMode )
	mapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:ISLocalToGlobalMappingCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscInt}, PetscCopyMode, Ptr{ISLocalToGlobalMapping}),
               comm, bs, n, indices, mode, mapping_,
              )

	mapping = mapping_[]

	return mapping
end 

"""
	mapping::ISLocalToGlobalMapping = ISLocalToGlobalMappingCreateIS(petsclib::PetscLibType,is::AbstractIS) 
Creates a mapping between a local (0 to n)
ordering and a global parallel ordering.

Not Collective

Input Parameter:
- `is` - index set containing the global numbers for each local number

Output Parameter:
- `mapping` - new mapping data structure

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingSetFromOptions()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingCreateIS"))
"""
function ISLocalToGlobalMappingCreateIS(petsclib::PetscLibType, is::AbstractIS)
    error("ISLocalToGlobalMappingCreateIS: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingCreateIS(petsclib::$UnionPetscLib, is::AbstractIS )
	mapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:ISLocalToGlobalMappingCreateIS, $petsc_library),
               PetscErrorCode,
               (CIS, Ptr{ISLocalToGlobalMapping}),
               is, mapping_,
              )

	mapping = mapping_[]

	return mapping
end 

"""
	mapping::ISLocalToGlobalMapping = ISLocalToGlobalMappingCreateSF(petsclib::PetscLibType,sf::PetscSF, start::PetscInt) 
Creates a mapping between a local (0 to n) ordering and a global parallel ordering induced by a star forest.

Collective

Input Parameters:
- `sf`    - star forest mapping contiguous local indices to (rank, offset)
- `start` - first global index on this process, or `PETSC_DECIDE` to compute contiguous global numbering automatically

Output Parameter:
- `mapping` - new mapping data structure

Level: advanced

-seealso: [](sec_scatter), `PetscSF`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingSetFromOptions()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingCreateSF"))
"""
function ISLocalToGlobalMappingCreateSF(petsclib::PetscLibType, sf::PetscSF, start::Integer)
    error("ISLocalToGlobalMappingCreateSF: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingCreateSF(petsclib::$UnionPetscLib, sf::PetscSF, start::$PetscInt )
	mapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:ISLocalToGlobalMappingCreateSF, $petsc_library),
               PetscErrorCode,
               (PetscSF, $PetscInt, Ptr{ISLocalToGlobalMapping}),
               sf, start, mapping_,
              )

	mapping = mapping_[]

	return mapping
end 

"""
	ISLocalToGlobalMappingDestroy(petsclib::PetscLibType,mapping::Union{ISLocalToGlobalMapping, Ref{ISLocalToGlobalMapping}}) 
Destroys a mapping between a local (0 to n)
ordering and a global parallel ordering.

Not Collective

Input Parameter:
- `mapping` - mapping data structure

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingDestroy"))
"""
function ISLocalToGlobalMappingDestroy(petsclib::PetscLibType, mapping::Union{ISLocalToGlobalMapping, Ref{ISLocalToGlobalMapping}})
    error("ISLocalToGlobalMappingDestroy: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingDestroy(petsclib::$UnionPetscLib, mapping::Union{ISLocalToGlobalMapping, Ref{ISLocalToGlobalMapping}} )
	mapping_ = mapping isa Base.RefValue ? mapping : Ref{ISLocalToGlobalMapping}(mapping)

    @chk ccall(
               (:ISLocalToGlobalMappingDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{ISLocalToGlobalMapping},),
               mapping_,
              )


	return nothing
end 

"""
	nltog::ISLocalToGlobalMapping = ISLocalToGlobalMappingDuplicate(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping) 
Duplicates the local to global mapping object

Not Collective

Input Parameter:
- `ltog` - local to global mapping

Output Parameter:
- `nltog` - the duplicated local to global mapping

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingDuplicate"))
"""
function ISLocalToGlobalMappingDuplicate(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingDuplicate: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingDuplicate(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping )
	nltog_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:ISLocalToGlobalMappingDuplicate, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{ISLocalToGlobalMapping}),
               ltog, nltog_,
              )

	nltog = nltog_[]

	return nltog
end 

"""
	array::Vector{PetscInt} = ISLocalToGlobalMappingGetBlockIndices(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping) 
Get global indices for every local block in a `ISLocalToGlobalMapping`

Not Collective

Input Parameter:
- `ltog` - local to global mapping

Output Parameter:
- `array` - array of indices

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingApply()`,
`ISLocalToGlobalMappingRestoreBlockIndices()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetBlockIndices"))
"""
function ISLocalToGlobalMappingGetBlockIndices(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetBlockIndices: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetBlockIndices(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping )
	array_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:ISLocalToGlobalMappingGetBlockIndices, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
               ltog, array_,
              )

	n = div(ISLocalToGlobalMappingGetSize(petsclib, ltog), ISLocalToGlobalMappingGetBlockSize(petsclib, ltog))
	array = unsafe_wrap(Array, array_[], n; own = false)

	return array
end 

"""
	nproc::PetscInt,procs::Ptr{PetscInt},numprocs::Ptr{PetscInt},indices::Ptr{Ptr{PetscInt}} = ISLocalToGlobalMappingGetBlockInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the neighbor information

Collective the first time it is called

Input Parameter:
- `mapping` - the mapping from local to global indexing

Output Parameters:
- `nproc`    - number of processes that are connected to the calling process
- `procs`    - neighboring processes
- `numprocs` - number of block indices for each process
- `indices`  - block indices (in local numbering) shared with neighbors (sorted by global numbering)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingRestoreBlockInfo()`, `ISLocalToGlobalMappingGetBlockMultiLeavesSF()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetBlockInfo"))
"""
function ISLocalToGlobalMappingGetBlockInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetBlockInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetBlockInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	nproc_ = Ref{$PetscInt}()
	procs_ = Ref{Ptr{$PetscInt}}()
	numprocs_ = Ref{Ptr{$PetscInt}}()
	indices_ = Ref{Ptr{Ptr{$PetscInt}}}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetBlockInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, nproc_, procs_, numprocs_, indices_,
              )

	nproc = nproc_[]
	procs = procs_[]
	numprocs = numprocs_[]
	indices = indices_[]

	return nproc,procs,numprocs,indices
end 

"""
	mlsf::PetscSF = ISLocalToGlobalMappingGetBlockMultiLeavesSF(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Get the star

Collective the first time it is called

Input Parameter:
- `mapping` - the mapping from local to global indexing

Output Parameter:
- `mlsf` - the `PetscSF`

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingGetBlockNodeInfo()`, `PetscSF`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetBlockMultiLeavesSF"))
"""
function ISLocalToGlobalMappingGetBlockMultiLeavesSF(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetBlockMultiLeavesSF: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetBlockMultiLeavesSF(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	mlsf_ = Ref{PetscSF}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetBlockMultiLeavesSF, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{PetscSF}),
               mapping, mlsf_,
              )

	mlsf = mlsf_[]

	return mlsf
end 

"""
	n::PetscInt,n_procs::Ptr{PetscInt},procs::Ptr{Ptr{PetscInt}} = ISLocalToGlobalMappingGetBlockNodeInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the neighbor information for each local block index

Collective the first time it is called

Input Parameter:
- `mapping` - the mapping from local to global indexing

Output Parameters:
- `n`       - number of local block nodes
- `n_procs` - an array storing the number of processes for each local block node (including self)
- `procs`   - the processes' rank for each local block node (sorted, self is first)

Level: advanced

-seealso: `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetBlockInfo()`, `ISLocalToGlobalMappingRestoreBlockNodeInfo()`, `ISLocalToGlobalMappingGetNodeInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetBlockNodeInfo"))
"""
function ISLocalToGlobalMappingGetBlockNodeInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetBlockNodeInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetBlockNodeInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	n_ = Ref{$PetscInt}()
	n_procs_ = Ref{Ptr{$PetscInt}}()
	procs_ = Ref{Ptr{Ptr{$PetscInt}}}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetBlockNodeInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, n_, n_procs_, procs_,
              )

	n = n_[]
	n_procs = n_procs_[]
	procs = procs_[]

	return n,n_procs,procs
end 

"""
	bs::PetscInt = ISLocalToGlobalMappingGetBlockSize(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the blocksize of the mapping
ordering and a global parallel ordering.

Not Collective

Input Parameter:
- `mapping` - mapping data structure

Output Parameter:
- `bs` - the blocksize

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetBlockSize"))
"""
function ISLocalToGlobalMappingGetBlockSize(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetBlockSize: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetBlockSize(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetBlockSize, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}),
               mapping, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	array::Vector{PetscInt} = ISLocalToGlobalMappingGetIndices(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping) 
Get global indices for every local point that is mapped

Not Collective

Input Parameter:
- `ltog` - local to global mapping

Output Parameter:
- `array` - array of indices, the length of this array may be obtained with `ISLocalToGlobalMappingGetSize()`

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingApply()`, `ISLocalToGlobalMappingRestoreIndices()`,
`ISLocalToGlobalMappingGetBlockIndices()`, `ISLocalToGlobalMappingRestoreBlockIndices()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetIndices"))
"""
function ISLocalToGlobalMappingGetIndices(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetIndices: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetIndices(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping )
	array_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:ISLocalToGlobalMappingGetIndices, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
               ltog, array_,
              )

	n = ISLocalToGlobalMappingGetSize(petsclib, ltog)
	array = unsafe_wrap(Array, array_[], n; own = false)

	return array
end 

"""
	nproc::PetscInt,procs::Ptr{PetscInt},numprocs::Ptr{PetscInt},indices::Ptr{Ptr{PetscInt}} = ISLocalToGlobalMappingGetInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the neighbor information for each process

Collective the first time it is called

Input Parameter:
- `mapping` - the mapping from local to global indexing

Output Parameters:
- `nproc`    - number of processes that are connected to the calling process
- `procs`    - neighboring processes
- `numprocs` - number of indices for each process
- `indices`  - indices (in local numbering) shared with neighbors (sorted by global numbering)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingRestoreInfo()`, `ISLocalToGlobalMappingGetNodeInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetInfo"))
"""
function ISLocalToGlobalMappingGetInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	nproc_ = Ref{$PetscInt}()
	procs_ = Ref{Ptr{$PetscInt}}()
	numprocs_ = Ref{Ptr{$PetscInt}}()
	indices_ = Ref{Ptr{Ptr{$PetscInt}}}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, nproc_, procs_, numprocs_, indices_,
              )

	nproc = nproc_[]
	procs = procs_[]
	numprocs = numprocs_[]
	indices = indices_[]

	return nproc,procs,numprocs,indices
end 

"""
	n::PetscInt,n_procs::Ptr{PetscInt},procs::Ptr{Ptr{PetscInt}} = ISLocalToGlobalMappingGetNodeInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the neighbor information of local nodes

Collective the first time it is called

Input Parameter:
- `mapping` - the mapping from local to global indexing

Output Parameters:
- `n`       - number of local nodes
- `n_procs` - an array storing the number of processes for each local node (including self)
- `procs`   - the processes' rank for each local node (sorted, self is first)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetInfo()`, `ISLocalToGlobalMappingRestoreNodeInfo()`, `ISLocalToGlobalMappingGetBlockNodeInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetNodeInfo"))
"""
function ISLocalToGlobalMappingGetNodeInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetNodeInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetNodeInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	n_ = Ref{$PetscInt}()
	n_procs_ = Ref{Ptr{$PetscInt}}()
	procs_ = Ref{Ptr{Ptr{$PetscInt}}}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetNodeInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, n_, n_procs_, procs_,
              )

	n = n_[]
	n_procs = n_procs_[]
	procs = procs_[]

	return n,n_procs,procs
end 

"""
	n::PetscInt = ISLocalToGlobalMappingGetSize(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Gets the local size of a local to global mapping

Not Collective

Input Parameter:
- `mapping` - local to global mapping

Output Parameter:
- `n` - the number of entries in the local mapping, `ISLocalToGlobalMappingGetIndices()` returns an array of this length

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetSize"))
"""
function ISLocalToGlobalMappingGetSize(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetSize: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetSize(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetSize, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}),
               mapping, n_,
              )

	n = n_[]

	return n
end 

"""
	type::ISLocalToGlobalMappingType = ISLocalToGlobalMappingGetType(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping) 
Get the type of the `ISLocalToGlobalMapping`

Not Collective

Input Parameter:
- `ltog` - the `ISLocalToGlobalMapping` object

Output Parameter:
- `type` - the type

Level: intermediate

-seealso: [](sec_scatter), `ISLocalToGlobalMappingType`, `ISLocalToGlobalMappingRegister()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingSetType()`,
`PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingGetType"))
"""
function ISLocalToGlobalMappingGetType(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingGetType: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingGetType(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping )
	type_ = Ref{ISLocalToGlobalMappingType}()

    @chk ccall(
               (:ISLocalToGlobalMappingGetType, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{ISLocalToGlobalMappingType}),
               ltog, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	ISLocalToGlobalMappingLoad(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, viewer::PetscViewer) 
Loads a local

Collective on viewer

Input Parameters:
- `mapping` - the newly loaded map, this needs to have been created with `ISLocalToGlobalMappingCreate()` or some related function before a call to `ISLocalToGlobalMappingLoad()`
- `viewer`  - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

-seealso: [](sec_scatter), `PetscViewer`, `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingView()`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingLoad"))
"""
function ISLocalToGlobalMappingLoad(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, viewer::PetscViewer)
    error("ISLocalToGlobalMappingLoad: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingLoad(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, viewer::PetscViewer )

    @chk ccall(
               (:ISLocalToGlobalMappingLoad, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, PetscViewer),
               mapping, viewer,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Registers a method for applying a global to local mapping with an `ISLocalToGlobalMapping`

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new method
- `function` - routine to create method context

-seealso: [](sec_scatter), `ISLocalToGlobalMappingRegisterAll()`, `ISLocalToGlobalMappingRegisterDestroy()`, `ISLOCALTOGLOBALMAPPINGBASIC`,
`ISLOCALTOGLOBALMAPPINGHASH`, `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingApply()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRegister"))
"""
function ISLocalToGlobalMappingRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("ISLocalToGlobalMappingRegister: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:ISLocalToGlobalMappingRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRegisterAll(petsclib::PetscLibType) 
Registers all of the local to global mapping components in the `IS` package.

Not Collective

Level: advanced

-seealso: [](sec_scatter), `ISRegister()`, `ISLocalToGlobalRegister()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRegisterAll"))
"""
function ISLocalToGlobalMappingRegisterAll(petsclib::PetscLibType)
    error("ISLocalToGlobalMappingRegisterAll: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:ISLocalToGlobalMappingRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreBlockIndices(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{PetscInt}}) 
Restore indices obtained with `ISLocalToGlobalMappingGetBlockIndices()`

Not Collective

Input Parameters:
- `ltog`  - local to global mapping
- `array` - array of indices

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingApply()`, `ISLocalToGlobalMappingGetIndices()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreBlockIndices"))
"""
function ISLocalToGlobalMappingRestoreBlockIndices(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{<:Number}})
    error("ISLocalToGlobalMappingRestoreBlockIndices: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreBlockIndices(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{$PetscInt}} )
	array_ = Ref{Ptr{$PetscInt}}(array isa Ptr ? array : pointer(array))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreBlockIndices, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
               ltog, array_,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreBlockInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, nproc::PetscInt, procs::Union{Ptr, AbstractArray{PetscInt}}, numprocs::Union{Ptr, AbstractArray{PetscInt}}, indices::Vector{PetscInt}) 
Frees the memory allocated by `ISLocalToGlobalMappingGetBlockInfo()`

Not Collective

Input Parameters:
- `mapping`  - the mapping from local to global indexing
- `nproc`    - number of processes that are connected to the calling process
- `procs`    - neighboring processes
- `numprocs` - number of block indices for each process
- `indices`  - block indices (in local numbering) shared with neighbors (sorted by global numbering)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreBlockInfo"))
"""
function ISLocalToGlobalMappingRestoreBlockInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, nproc::Integer, procs::Union{Ptr, AbstractArray{<:Number}}, numprocs::Union{Ptr, AbstractArray{<:Number}}, indices::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingRestoreBlockInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreBlockInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, nproc::$PetscInt, procs::Union{Ptr, AbstractArray{$PetscInt}}, numprocs::Union{Ptr, AbstractArray{$PetscInt}}, indices::Vector{$PetscInt} )
	nproc_ = Ref{$PetscInt}(nproc)
	procs_ = Ref{Ptr{$PetscInt}}(procs isa Ptr ? procs : pointer(procs))
	numprocs_ = Ref{Ptr{$PetscInt}}(numprocs isa Ptr ? numprocs : pointer(numprocs))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreBlockInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, nproc_, procs_, numprocs_, indices,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreBlockNodeInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, n::PetscInt, n_procs::Union{Ptr, AbstractArray{PetscInt}}, procs::Vector{PetscInt}) 
Frees the memory allocated by `ISLocalToGlobalMappingGetBlockNodeInfo()`

Not Collective

Input Parameters:
- `mapping` - the mapping from local to global indexing
- `n`       - number of local block nodes
- `n_procs` - an array storing the number of processes for each local block nodes (including self)
- `procs`   - the processes' rank for each local block node (sorted, self is first)

Level: advanced

-seealso: `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetBlockNodeInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreBlockNodeInfo"))
"""
function ISLocalToGlobalMappingRestoreBlockNodeInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, n::Integer, n_procs::Union{Ptr, AbstractArray{<:Number}}, procs::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingRestoreBlockNodeInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreBlockNodeInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, n::$PetscInt, n_procs::Union{Ptr, AbstractArray{$PetscInt}}, procs::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}(n)
	n_procs_ = Ref{Ptr{$PetscInt}}(n_procs isa Ptr ? n_procs : pointer(n_procs))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreBlockNodeInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, n_, n_procs_, procs,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreIndices(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{PetscInt}}) 
Restore indices obtained with `ISLocalToGlobalMappingGetIndices()`

Not Collective

Input Parameters:
- `ltog`  - local to global mapping
- `array` - array of indices

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingApply()`, `ISLocalToGlobalMappingGetIndices()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreIndices"))
"""
function ISLocalToGlobalMappingRestoreIndices(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{<:Number}})
    error("ISLocalToGlobalMappingRestoreIndices: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreIndices(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping, array::Union{Ptr, AbstractArray{$PetscInt}} )
	array_ = Ref{Ptr{$PetscInt}}(array isa Ptr ? array : pointer(array))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreIndices, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{Ptr{$PetscInt}}),
               ltog, array_,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, nproc::PetscInt, procs::Union{Ptr, AbstractArray{PetscInt}}, numprocs::Union{Ptr, AbstractArray{PetscInt}}, indices::Vector{PetscInt}) 
Frees the memory allocated by `ISLocalToGlobalMappingGetInfo()`

Not Collective

Input Parameters:
- `mapping`  - the mapping from local to global indexing
- `nproc`    - number of processes that are connected to the calling process
- `procs`    - neighboring processes
- `numprocs` - number of indices for each process
- `indices`  - indices (in local numbering) shared with neighbors (sorted by global numbering)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreInfo"))
"""
function ISLocalToGlobalMappingRestoreInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, nproc::Integer, procs::Union{Ptr, AbstractArray{<:Number}}, numprocs::Union{Ptr, AbstractArray{<:Number}}, indices::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingRestoreInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, nproc::$PetscInt, procs::Union{Ptr, AbstractArray{$PetscInt}}, numprocs::Union{Ptr, AbstractArray{$PetscInt}}, indices::Vector{$PetscInt} )
	nproc_ = Ref{$PetscInt}(nproc)
	procs_ = Ref{Ptr{$PetscInt}}(procs isa Ptr ? procs : pointer(procs))
	numprocs_ = Ref{Ptr{$PetscInt}}(numprocs isa Ptr ? numprocs : pointer(numprocs))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, nproc_, procs_, numprocs_, indices,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingRestoreNodeInfo(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, n::PetscInt, n_procs::Union{Ptr, AbstractArray{PetscInt}}, procs::Vector{PetscInt}) 
Frees the memory allocated by `ISLocalToGlobalMappingGetNodeInfo()`

Not Collective

Input Parameters:
- `mapping` - the mapping from local to global indexing
- `n`       - number of local nodes
- `n_procs` - an array storing the number of processes for each local node (including self)
- `procs`   - the processes' rank for each local node (sorted, self is first)

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`, `ISLocalToGlobalMappingCreate()`,
`ISLocalToGlobalMappingGetInfo()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingRestoreNodeInfo"))
"""
function ISLocalToGlobalMappingRestoreNodeInfo(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, n::Integer, n_procs::Union{Ptr, AbstractArray{<:Number}}, procs::AbstractVector{<:Number})
    error("ISLocalToGlobalMappingRestoreNodeInfo: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingRestoreNodeInfo(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, n::$PetscInt, n_procs::Union{Ptr, AbstractArray{$PetscInt}}, procs::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}(n)
	n_procs_ = Ref{Ptr{$PetscInt}}(n_procs isa Ptr ? n_procs : pointer(n_procs))

    @chk ccall(
               (:ISLocalToGlobalMappingRestoreNodeInfo, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{Ptr{$PetscInt}}}),
               mapping, n_, n_procs_, procs,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingSetBlockSize(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, bs::PetscInt) 
Sets the blocksize of the mapping

Not Collective

Input Parameters:
- `mapping` - mapping data structure
- `bs`      - the blocksize

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreateIS()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingSetBlockSize"))
"""
function ISLocalToGlobalMappingSetBlockSize(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, bs::Integer)
    error("ISLocalToGlobalMappingSetBlockSize: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingSetBlockSize(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, bs::$PetscInt )

    @chk ccall(
               (:ISLocalToGlobalMappingSetBlockSize, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, $PetscInt),
               mapping, bs,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingSetFromOptions(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping) 
Set mapping options from the options database.

Not Collective

Input Parameter:
- `mapping` - mapping data structure

Options Database Key:
- `-islocaltoglobalmapping_type (basic|hash)` - nonscalable and scalable versions

Level: advanced

-seealso: [](sec_scatter), `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`,
`ISLocalToGlobalMappingCreateIS()`, `ISLOCALTOGLOBALMAPPINGBASIC`,
`ISLOCALTOGLOBALMAPPINGHASH`, `ISLocalToGlobalMappingSetType()`, `ISLocalToGlobalMappingType`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingSetFromOptions"))
"""
function ISLocalToGlobalMappingSetFromOptions(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping)
    error("ISLocalToGlobalMappingSetFromOptions: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingSetFromOptions(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:ISLocalToGlobalMappingSetFromOptions, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping,),
               mapping,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingSetType(petsclib::PetscLibType,ltog::ISLocalToGlobalMapping, type::ISLocalToGlobalMappingType) 
Sets the implementation type `ISLocalToGlobalMapping` will use

Logically Collective

Input Parameters:
- `ltog` - the `ISLocalToGlobalMapping` object
- `type` - a known method

Options Database Key:
- `-islocaltoglobalmapping_type (basic|hash)` - Sets the method used for applying the mapping, see `ISLocalToGlobalMappingType`

Level: intermediate

-seealso: [](sec_scatter), `ISLocalToGlobalMappingType`, `ISLocalToGlobalMappingRegister()`, `ISLocalToGlobalMappingCreate()`, `ISLocalToGlobalMappingGetType()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingSetType"))
"""
function ISLocalToGlobalMappingSetType(petsclib::PetscLibType, ltog::ISLocalToGlobalMapping, type::ISLocalToGlobalMappingType)
    error("ISLocalToGlobalMappingSetType: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingSetType(petsclib::$UnionPetscLib, ltog::ISLocalToGlobalMapping, type::ISLocalToGlobalMappingType )

    @chk ccall(
               (:ISLocalToGlobalMappingSetType, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, ISLocalToGlobalMappingType),
               ltog, type,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingView(petsclib::PetscLibType,mapping::ISLocalToGlobalMapping, viewer::PetscViewer) 
View a local to global mapping

Collective on viewer

Input Parameters:
- `mapping` - local to global mapping
- `viewer`  - viewer

Level: intermediate

-seealso: [](sec_scatter), `PetscViewer`, `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingDestroy()`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingView"))
"""
function ISLocalToGlobalMappingView(petsclib::PetscLibType, mapping::ISLocalToGlobalMapping, viewer::PetscViewer)
    error("ISLocalToGlobalMappingView: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingView(petsclib::$UnionPetscLib, mapping::ISLocalToGlobalMapping, viewer::PetscViewer )

    @chk ccall(
               (:ISLocalToGlobalMappingView, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, PetscViewer),
               mapping, viewer,
              )


	return nothing
end 

"""
	ISLocalToGlobalMappingViewFromOptions(petsclib::PetscLibType,A::ISLocalToGlobalMapping, obj, name::String) 
View an `ISLocalToGlobalMapping` based on values in the options database

Collective

Input Parameters:
- `A`    - the local to global mapping object
- `obj`  - Optional object that provides the options prefix used for the options database query
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

-seealso: [](sec_scatter), `PetscViewer`, `ISLocalToGlobalMapping`, `ISLocalToGlobalMappingView`, `PetscObjectViewFromOptions()`, `ISLocalToGlobalMappingCreate()`

# External Links
$(_doc_external("IS/ISLocalToGlobalMappingViewFromOptions"))
"""
function ISLocalToGlobalMappingViewFromOptions(petsclib::PetscLibType, A::ISLocalToGlobalMapping, obj, name::String)
    error("ISLocalToGlobalMappingViewFromOptions: no generated method for these argument types")
end

@for_petsc function ISLocalToGlobalMappingViewFromOptions(petsclib::$UnionPetscLib, A::ISLocalToGlobalMapping, obj, name::String )

    @chk ccall(
               (:ISLocalToGlobalMappingViewFromOptions, $petsc_library),
               PetscErrorCode,
               (ISLocalToGlobalMapping, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

