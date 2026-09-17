

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

If keyword argument `setfromoptions == true` then `set_from_options!` called.

If keyword argument `dmsetup == true` then `setup!` is called.

When `D == 1` the `stencil_type` argument is not required and ignored if
specified.

# External Links
$(doc_external("DMStag/DMStagCreate1d"))
$(doc_external("DMStag/DMStagCreate2d"))
$(doc_external("DMStag/DMStagCreate3d"))
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
    N1 == N + 1 || throw(
        DimensionMismatch(
            "dof_per_node has length $N1, " *
            "but a $N-dimensional DMStag needs $(N + 1)",
        ),
    )
    PetscInt = inttype(PetscLib)
    stencil_type = stencil_type_enum(DMStagStencilType, stencil_type)

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
            points_per_proc[d] isa Array || throw(
                ArgumentError(
                    "points_per_proc[$d] must be an Array, " *
                    "got $(typeof(points_per_proc[d]))",
                ),
            )
            length(points_per_proc[d]) == MPI.Comm_size(comm) || throw(
                DimensionMismatch(
                    "points_per_proc[$d] has $(length(points_per_proc[d])) entries, " *
                    "but the communicator has $(MPI.Comm_size(comm)) ranks",
                ),
            )
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

    # Take ownership of the handle the creator returned: from here on `da` is a
    # `DMStag{PetscLib, N}` and every DMStag method dispatches on it.
    da = DMStag{PetscLib, N}(da.ptr, petsclib.age, true)

    if !isempty(prefix)
        # options prefix
        LibPETSc.DMSetOptionsPrefix(petsclib, da, prefix)
    end

    if setfromoptions
        # set options (if any)
        opts = PetscOptions(petsclib; options...)
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
    return own_dm!(da, comm)
end


"""
    DMStag(dm::DMStag, dof_per_node; setfromoptions = true, dmsetup = true, options...)

A `DMStag` compatible with `dm` — same dimension, communicator and layout — but
with the degrees of freedom per stratum given by `dof_per_node`.

The v0.4 spelling took `dmsetfromoptions` and slurped anything else into an
`options` collection it never used. Both are honoured here, and `options...` is
the same options database keyword set the other constructors take.

# External Links
$(doc_external("DMSTAG/DMStagCreateCompatibleDMStag"))
"""
function DMStag(
    dm::DMStag{PetscLib, N},
    dof_per_node::Union{NTuple{2, Integer}, NTuple{3, Integer}, NTuple{4, Integer}};
    setfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib, N}
    petsclib = getlib(PetscLib)
    PetscInt = petsclib.PetscInt

    dof_per_node_C = ntuple(i -> i <= length(dof_per_node) ? PetscInt(dof_per_node[i]) : PetscInt(0), 4)

    dmnew = LibPETSc.DMStagCreateCompatibleDMStag(
        PetscLib,
        dm,
        dof_per_node_C[1],
        dof_per_node_C[2],
        dof_per_node_C[3],
        dof_per_node_C[4],
    )
    dmnew = DMStag{PetscLib, N}(dmnew.ptr, petsclib.age, true)

    if setfromoptions
        opts = PetscOptions(petsclib; options...)
        push!(opts)
        LibPETSc.DMSetFromOptions(PetscLib, dmnew)
        pop!(opts)
    end

    dmsetup && setup!(dmnew)

    return own_dm!(dmnew, comm(dm))
end


#Base.size(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm)
#globalsize(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm::AbstractDMStag)
#boundarytypes(dm::AbstractDMStag)  = DMStagGetBoundaryTypes(dm::AbstractDMStag) 

"""
    set_uniform_coordinates!(dm::DMStag, xyzmin, xyzmax)

The [`DMStag`](@ref) method of [`set_uniform_coordinates!`](@ref): sets uniform
coordinates on `dm` in the range specified by `xyzmin` and `xyzmax`.

# External Links
$(doc_external("DMSTAG/DMStagSetUniformCoordinatesProduct"))
"""
function set_uniform_coordinates!(
    dm::DMStag{PetscLib},
    xyzmin::NTuple,
    xyzmax::NTuple,
    ) where {PetscLib}
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
    corners(dm::DMStag{PetscLib, N})

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. Also included is `nextra`,
the number of extra partial elements in each direction.

The result is dimension-correct (§12): `lower` and `upper` are
`CartesianIndex{N}`, `size` and `nextra` are `NTuple{N,Int}`, with no padding
to three entries. This is a break with no shim (§16).

# External Links
$(doc_external("DMSTAG/DMStagGetCorners"))
"""
function corners(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    x, y, z, m, n, p, nex, ney, nez = LibPETSc.DMStagGetCorners(PetscLib, dm)
    lo = (Int(x), Int(y), Int(z))
    sz = (Int(m), Int(n), Int(p))
    ne = (Int(nex), Int(ney), Int(nez))
    return (
        lower  = CartesianIndex(ntuple(i -> lo[i] + 1, Val(N))),
        upper  = CartesianIndex(ntuple(i -> lo[i] + sz[i], Val(N))),
        size   = ntuple(i -> sz[i], Val(N)),
        nextra = ntuple(i -> ne[i], Val(N)),
    )
end


"""
    ghost_corners(dm::DMStag{PetscLib, N})

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`.

There is no `nextra` field: `DMStagGetGhostCorners` does not report the extra
partial elements, and v0.4's docstring promised a field the function never
returned. Ask [`corners`](@ref) for `nextra`.

Dimension-correct like [`corners`](@ref): `CartesianIndex{N}` and `NTuple{N,Int}`.

# External Links
$(doc_external("DMSTAG/DMStagGetGhostCorners"))
"""
function ghost_corners(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    x, y, z, m, n, p = LibPETSc.DMStagGetGhostCorners(PetscLib, dm)
    lo = (Int(x), Int(y), Int(z))
    sz = (Int(m), Int(n), Int(p))
    return (
        lower = CartesianIndex(ntuple(i -> lo[i] + 1, Val(N))),
        upper = CartesianIndex(ntuple(i -> lo[i] + sz[i], Val(N))),
        size  = ntuple(i -> sz[i], Val(N)),
    )
end


"""
    local_indices(dm::DMStag)

Return indices for the central/vertex nodes of a local (ghosted) array built from the
input `dm`. This takes ghost points into account and provides index ranges for
accessing staggered data, so that e.g. `array[local_indices(dm).center.x]`
correctly skips the ghost region on the low side.

# Returns

A `NamedTuple` with:
- `center`: `NamedTuple` of ranges keyed `x`, `y`, `z` for cell-centered indices
- `vertex`: `NamedTuple` of ranges keyed `x`, `y`, `z` for vertex indices

Both are dimension-correct (§12): a 2D `DMStag` yields `(x = …, y = …)` with no
`z`. This is a break with no shim (§16).

# Note

In Julia, array indices start at 1, whereas PETSc uses 0-based indexing with
possibly negative ghost indices. This function handles the conversion automatically.

# See also

[`global_indices`](@ref) for the equivalent indices into a non-ghosted, global array.
"""
function local_indices(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    # In Julia, indices in arrays start @ 1, whereas they can go negative in C
    x, y, z, m, n, p, nex, ney, nez = LibPETSc.DMStagGetCorners(PetscLib, dm)
    gx, gy, gz, _, _, _ = LibPETSc.DMStagGetGhostCorners(PetscLib, dm)

    c  = (Int(x) + 1, Int(y) + 1, Int(z) + 1)
    gc = (Int(gx) + 1, Int(gy) + 1, Int(gz) + 1)
    sz = (Int(m), Int(n), Int(p))
    ne = (Int(nex), Int(ney), Int(nez))

    # The low-side ghost band is skipped by shifting the owned range by the
    # distance between the owned and the ghosted lower corner.
    lo = ntuple(i -> 2c[i] - gc[i], Val(N))
    hi = ntuple(i -> lo[i] + sz[i] - 1, Val(N))

    return (
        center = axis_ranges(i -> lo[i]:hi[i], Val(N)),
        vertex = axis_ranges(i -> lo[i]:(hi[i] + ne[i]), Val(N)),
    )
end

"""
    global_indices(dm::DMStag)

Return indices for the central/vertex nodes of the global (non-ghosted) array built
from the input `dm`, i.e. the process-local interior region only, excluding ghost
points.

# Returns

A `NamedTuple` with:
- `center`: `NamedTuple` of ranges keyed `x`, `y`, `z` for cell-centered indices
- `vertex`: `NamedTuple` of ranges keyed `x`, `y`, `z` for vertex indices

Both are dimension-correct (§12): a 2D `DMStag` yields `(x = …, y = …)` with no
`z`. This is a break with no shim (§16).

# Note

In Julia, array indices start at 1, whereas PETSc uses 0-based indexing. This function
handles the conversion automatically.

# See also

[`local_indices`](@ref) for the equivalent indices into a ghosted, local array.
"""
function global_indices(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    x, y, z, m, n, p, nex, ney, nez = LibPETSc.DMStagGetCorners(PetscLib, dm)
    lo = (Int(x), Int(y), Int(z))
    sz = (Int(m), Int(n), Int(p))
    ne = (Int(nex), Int(ney), Int(nez))

    return (
        center = axis_ranges(i -> (lo[i] + 1):(lo[i] + sz[i]), Val(N)),
        vertex = axis_ranges(i -> (lo[i] + 1):(lo[i] + sz[i] + ne[i]), Val(N)),
    )
end

"""
    slot::Int = dof_slot(dm::DMStag, loc::LibPETSc.DMStagStencilLocation, dof::Int) 

Returns the location `slot` for a degree of freedom `dof` at a given stencil location `loc` in the DMStag `dm`.
Note that the returned `slot` is 1-based for Julia compatibility.    
"""
function dof_slot(dm::DMStag{PetscLib}, loc::LibPETSc.DMStagStencilLocation, dof::Int) where {PetscLib} 
    slot = LibPETSc.DMStagGetLocationSlot(getlib(PetscLib), dm, loc, PetscLib.PetscInt(dof))
    return slot+1
end
