

# Helper to convert points_per_proc tuples to PetscInt
to_petscint_tuple(t::Tuple, PetscInt) = map(arr -> PetscInt.(arr), t)

"""
    da = DMStag(
        petsclib::PetscLib
        comm::MPI.Comm,
        boundary_type::NTuple{D, DMBoundaryType},
        global_dim::NTuple{D, Integer},
        dof_per_node::NTuple{1 + D, Integer},
        stencil_width::Integer,
        stencil_type;
        points_per_proc::Tuple,
        processors::Tuple,
        setfromoptions = true,
        dmsetup = true,
        prefix = "",
        options...
    )

Creates a `D`-dimensional distributed staggered array with the options specified
using keyword arguments.

The Tuple `dof_per_node` specifies how many degrees of freedom are at all the
staggerings in the order:
 - 1D: `(vertex, element)`
 - 2D: `(vertex, edge, element)`
 - 3D: `(vertex, edge, face, element)`

If keyword argument `points_per_proc[k] isa Vector{petsclib.PetscInt}` then this
specifies the points per processor in dimension `k`.

If keyword argument `processors[k] isa Integer` then this specifies the number of
processors used in dimension `k`; ignored when `D == 1`.

If keyword argument `setfromoptions == true` then `setfromoptions!` called.

If keyword argument `dmsetup == true` then `setup!` is called.

When `D == 1` the `stencil_type` argument is not required and ignored if
specified.

# External Links
$(_doc_external("DMStag/DMStagCreate1d"))
$(_doc_external("DMStag/DMStagCreate2d"))
$(_doc_external("DMStag/DMStagCreate3d"))
"""
function DMStag(
    petsclib::PetscLib,
    comm::MPI.Comm,
    boundary_type::NTuple{N, DMBoundaryType},
    global_dim::NTuple{N, Integer},
    dof_per_node::NTuple{N1, Integer},
    stencil_width::Integer,
    stencil_type = DMSTAG_STENCIL_BOX;
    points_per_proc::Union{Tuple, Nothing} = nothing,
    processors = nothing,
    setfromoptions = true,
    dmsetup = true,
    prefix = "",
    options...,
) where {PetscLib, N, N1}
    @assert N1 == N + 1 "dof_per_node should have length N + 1"
    PetscInt = inttype(PetscLib)
    opts = Options(petsclib; options...)

    if isnothing(points_per_proc)
        points_per_proc = ntuple(_ -> nothing, N)
    end
    if isnothing(processors)
        processors = ntuple(_ -> PETSC_DECIDE, N)
    end

    ref_points_per_proc = ntuple(N) do d
        if isnothing(points_per_proc[d]) || points_per_proc[d] == PETSC_DECIDE
            C_NULL
        else
            @assert points_per_proc[d] isa Array
            @assert length(points_per_proc[d]) == MPI.Comm_size(comm)
            points_per_proc[d]
        end
    end
   # ref_points_per_proc = to_petscint_tuple(ref_points_per_proc, PetscInt)  

    if N==1  
        da = LibPETSc.DMStagCreate1d(petsclib,
                                   comm, 
                                   boundary_type[1], 
                                   PetscInt(global_dim[1]), 
                                   PetscInt(dof_per_node[1]), 
                                   PetscInt(dof_per_node[2]), 
                                   stencil_type,
                                   PetscInt(stencil_width), 
                                   ref_points_per_proc[1]
                                   )
    elseif N==2
        da =   LibPETSc.DMStagCreate2d(
                                    petsclib,
                                    comm,
                                    boundary_type[1], boundary_type[2],
                                    PetscInt(global_dim[1]), PetscInt(global_dim[2]),
                                    PetscInt(processors[1]), PetscInt(processors[2]),
                                    PetscInt(dof_per_node[1]), 
                                    PetscInt(dof_per_node[2]),
                                    PetscInt(dof_per_node[3]), 
                                    stencil_type,
                                    PetscInt(stencil_width),
                                    ref_points_per_proc[1], ref_points_per_proc[2]
                                )                                    
     elseif N==3
        da =   LibPETSc.DMStagCreate3d(
                                    petsclib,
                                    comm,
                                    boundary_type[1], boundary_type[2], boundary_type[3],
                                    PetscInt(global_dim[1]), PetscInt(global_dim[2]), PetscInt(global_dim[3]),
                                    PetscInt(processors[1]), PetscInt(processors[2]), PetscInt(processors[3]),
                                    PetscInt(dof_per_node[1]), 
                                    PetscInt(dof_per_node[2]),
                                    PetscInt(dof_per_node[3]),
                                    PetscInt(dof_per_node[4]), 
                                    stencil_type,
                                    PetscInt(stencil_width),
                                    ref_points_per_proc[1], ref_points_per_proc[2], ref_points_per_proc[3]
                                )                    
    end

    if !isempty(prefix)
        # options prefix
        LibPETSc.DMSetOptionsPrefix(petsclib, da, prefix)
    end

    if setfromoptions
        # set options (if any)
        opts = PETSc.Options(petsclib; options...);
        push!(opts)
        LibPETSc.DMSetFromOptions(PetscLib, da)
        pop!(opts)
    end
    
    if dmsetup
        # initialize dmda
        setup!(da)            
    end

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, da)
    end
    return da
end


function DMStag(
    dm::AbstractPetscDM{PetscLib},
    dof_per_node::Union{NTuple{2,Int},NTuple{3,Int},NTuple{4,Int}},
    dmsetfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib}
    @assert  PETSc.gettype(dm) == "stag" "DM must be of type DMStag"
    petsclib = getlib(PetscLib)
    PetscInt = petsclib.PetscInt
    dmnew = PetscDM{PetscLib}(C_NULL, petsclib.age)

    s = size(dof_per_node,1)

    dof_per_node_C = [0,0,0,0]

    for (i, value) in enumerate(dof_per_node)
        dof_per_node_C[i] = value
    end


    #with(dm.opts) do
    dmnew =  LibPETSc.DMStagCreateCompatibleDMStag(
                    PetscLib,
                    dm,
                    PetscInt(dof_per_node_C[1]),
                    PetscInt(dof_per_node_C[2]),
                    PetscInt(dof_per_node_C[3]),
                    PetscInt(dof_per_node_C[4]),
                    )
    #end

    #=
    dmsetfromoptions && setfromoptions!(dmnew)
    dmsetup && setup!(dmnew)

    =#
    comm  = getcomm(dm);

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dmnew)
    end    

    return dmnew
end


#Base.size(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm)
#globalsize(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm::AbstractDMStag)
#boundarytypes(dm::AbstractDMStag)  = DMStagGetBoundaryTypes(dm::AbstractDMStag) 

"""
    setuniformcoordinates_stag!(
        dm::AbstractDMStag,
        xyzmin::Union{NTuple{1,Int},NTuple{2,Int},NTuple{3,Int}},
        xyzmax::Union{NTuple{1,Int},NTuple{2,Int},NTuple{3,Int}},
    )

Sets uniform coordinates for the DMStag `dm` in the range specified by `xyzmin` and `xyzmax`.
"""
function setuniformcoordinates_stag!(
    dm::AbstractPetscDM{PetscLib},
    xyzmin::NTuple,
    xyzmax::NTuple,
    ) where {PetscLib}
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag"
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar

    xmin = PetscScalar(xyzmin[1])
    xmax = PetscScalar(xyzmax[1])

    s = size(xyzmin,1)

    ymin = (s > 1) ? PetscScalar(xyzmin[2]) : PetscScalar(0)
    ymax = (s > 1) ? PetscScalar(xyzmax[2]) : PetscScalar(0)

    zmin = (s > 2) ? PetscScalar(xyzmin[3]) : PetscScalar(0)
    zmax = (s > 2) ? PetscScalar(xyzmax[3]) : PetscScalar(0)
    
    #=
    LibPETSc.DMStagSetUniformCoordinatesProduct(
        getlib(PetscLib),
        dm,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    )
    =#
    petsclib=getlib(PetscLib)
    LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm, xmin, xmax, ymin, ymax, zmin, zmax)

    return nothing
end

"""
    corners = getcorners_dmstag(dm::AbstractPetscDM)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. Also included is `nextra` of
the number of extra partial elements in each direction.

# External Links
$(_doc_external("DMDA/DMStagGetCorners"))
"""
function getcorners_dmstag(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag"
    x, y, z, m, n, p, nExtrax, nExtray, nExtraz = LibPETSc.DMStagGetCorners(PetscLib, dm)
    return (
        lower  = CartesianIndex(x + 1, y + 1, z + 1),
        upper  = CartesianIndex(x + m, y + n, z + p),
        size   = (m, n, p),
        nextra = (nExtrax, nExtray, nExtraz),
    )
end


"""
    corners = getghostcorners_dmstag(dm::AbstractPetscDM)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`. Also included is `nextra` of
the number of extra partial elements in each direction.

# External Links
$(_doc_external("DMDA/DMStagGetCorners"))
"""
function getghostcorners_dmstag(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag"
    x, y, z, m, n, p = LibPETSc.DMStagGetGhostCorners(PetscLib, dm)
    return (
        lower = CartesianIndex(x + 1, y + 1, z + 1),
        upper = CartesianIndex(x + m, y + n, z + p),
        size  = (m, n, p),
    )
end


"""
    local_indices_dmstag(dm::DMStag)

Return indices for the central/vertex nodes of a local (ghosted) array built from the
input `dm`. This takes ghost points into account and provides index ranges for
accessing staggered data, so that e.g. `array[local_indices_dmstag(dm).center.x]`
correctly skips the ghost region on the low side.

# Returns

A `NamedTuple` with:
- `center`: Tuple of ranges `(x, y, z)` for cell-centered indices
- `vertex`: Tuple of ranges `(x, y, z)` for vertex indices

# Note

In Julia, array indices start at 1, whereas PETSc uses 0-based indexing with
possibly negative ghost indices. This function handles the conversion automatically.

# See also

[`global_indices_dmstag`](@ref) for the equivalent indices into a non-ghosted, global array.
"""
function local_indices_dmstag end

function local_indices_dmstag(dm::PetscDM{PetscLib}) where {PetscLib}
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag" 
    # In Julia, indices in arrays start @ 1, whereas they can go negative in C
    x, y, z, m, n, p, nx, ny, nz = LibPETSc.DMStagGetCorners(PetscLib, dm)
    gx, gy, gz, _, _, _ = LibPETSc.DMStagGetGhostCorners(PetscLib, dm)

    c  = CartesianIndex(x + 1, y + 1, z + 1)
    gc = CartesianIndex(gx + 1, gy + 1, gz + 1)

    lo = c + (c - gc)
    hi = lo + CartesianIndex(m - 1, n - 1, p - 1)

    return (
        center = (
            x = lo[1]:hi[1],
            y = lo[2]:hi[2],
            z = lo[3]:hi[3],
        ),
        vertex = (
            x = lo[1]:(hi[1] + nx),
            y = lo[2]:(hi[2] + ny),
            z = lo[3]:(hi[3] + nz),
        ),
    )

end



Base.@deprecate DMStagGetIndices(dm) local_indices_dmstag(dm)

"""
    global_indices_dmstag(dm::DMStag)

Return indices for the central/vertex nodes of the global (non-ghosted) array built
from the input `dm`, i.e. the process-local interior region only, excluding ghost
points.

# Returns

A `NamedTuple` with:
- `center`: Tuple of ranges `(x, y, z)` for cell-centered indices
- `vertex`: Tuple of ranges `(x, y, z)` for vertex indices

# Note

In Julia, array indices start at 1, whereas PETSc uses 0-based indexing. This function
handles the conversion automatically.

# See also

[`local_indices_dmstag`](@ref) for the equivalent indices into a ghosted, local array.
"""
function global_indices_dmstag(dm::PetscDM{PetscLib}) where {PetscLib}
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag"
    x, y, z, m, n, p, nx, ny, nz = LibPETSc.DMStagGetCorners(PetscLib, dm)

    return (
        center = (
            x = (x + 1):(x + m),
            y = (y + 1):(y + n),
            z = (z + 1):(z + p),
        ),
        vertex = (
            x = (x + 1):(x + m + nx),
            y = (y + 1):(y + n + ny),
            z = (z + 1):(z + p + nz),
        ),
    )

end

"""
    slot::Int = DMStagDOF_Slot(dm::PetscDM{PetscLib}, loc::LibPETSc.DMStagStencilLocation, dof::Int) 

Returns the location `slot` for a degree of freedom `dof` at a given stencil location `loc` in the DMStag `dm`.
Note that the returned `slot` is 1-based for Julia compatibility.    
"""
function DMStagDOF_Slot(dm::PetscDM{PetscLib}, loc::LibPETSc.DMStagStencilLocation, dof::Int) where {PetscLib} 
    @assert PETSc.gettype(dm) == "stag" "DM must be of type DMStag" 

    slot = LibPETSc.DMStagGetLocationSlot(getlib(PetscLib), dm, loc, PetscLib.PetscInt(dof))
    return slot+1
end
