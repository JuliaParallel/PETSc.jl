

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
    check_initialized(getlib(PetscLib))
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
    check_initialized(petsclib)
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

    return dm
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
`dof` is PETSc's 0-based component number at that location, as in [`stencil`](@ref).
"""
function dof_slot(dm::DMStag{PetscLib}, loc::LibPETSc.DMStagStencilLocation, dof::Int) where {PetscLib} 
    slot = LibPETSc.DMStagGetLocationSlot(getlib(PetscLib), dm, loc, PetscLib.PetscInt(dof))
    return slot+1
end


# ============================================================================
#   Locations and stencils
# ============================================================================

"""
    vertex_location(dm::DMStag)

The location of the vertex DMStag stores with each element, its lower corner:
`DMSTAG_LEFT` in 1D, `DMSTAG_DOWN_LEFT` in 2D and `DMSTAG_BACK_DOWN_LEFT` in 3D.

See also [`face_location`](@ref), [`edge_location`](@ref), [`element_location`](@ref)
and [`stencil`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencilLocation"))
"""
vertex_location(::DMStag{PetscLib, 1}) where {PetscLib} = LibPETSc.DMSTAG_LEFT
vertex_location(::DMStag{PetscLib, 2}) where {PetscLib} = LibPETSc.DMSTAG_DOWN_LEFT
vertex_location(::DMStag{PetscLib, 3}) where {PetscLib} = LibPETSc.DMSTAG_BACK_DOWN_LEFT

"""
    edge_location(dm::DMStag, a::Integer, b::Integer)

The location of the edge touching the lower faces of axes `a` and `b`, in either
order: the edge a shear component τ_ab lives on. In 3D, `(1, 2)` is
`DMSTAG_DOWN_LEFT`, `(1, 3)` is `DMSTAG_BACK_LEFT` and `(2, 3)` is
`DMSTAG_BACK_DOWN`. In 2D the vertex plays that role, so `(1, 2)` returns
`DMSTAG_DOWN_LEFT`. A 1D DM has no edges, and an axis outside `1:ndims(dm)`, or
`a == b`, throws an `ArgumentError`.

See also [`vertex_location`](@ref), [`face_location`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencilLocation"))
"""
function edge_location(dm::DMStag{PetscLib, N}, a::Integer, b::Integer) where {PetscLib, N}
    N == 1 && throw(ArgumentError("a 1D DMStag has no edges"))
    check_axis(dm, a)
    check_axis(dm, b)
    a == b && throw(ArgumentError("an edge needs two different axes, got ($a, $b)"))
    lo, hi = minmax(a, b)
    lo == 1 && hi == 2 && return LibPETSc.DMSTAG_DOWN_LEFT
    lo == 1 && return LibPETSc.DMSTAG_BACK_LEFT
    return LibPETSc.DMSTAG_BACK_DOWN
end

"""
    face_location(dm::DMStag, axis::Integer)

The location of the face normal to `axis` on the lower side of each element:
`DMSTAG_LEFT` for axis 1, `DMSTAG_DOWN` for axis 2 and `DMSTAG_BACK` for axis 3.
Axes count the DM's own axes, so `face_location(dm, ndims(dm))` is the last axis in
any dimension. An axis outside `1:ndims(dm)` throws an `ArgumentError`.

See also [`vertex_location`](@ref), [`edge_location`](@ref), [`element_location`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencilLocation"))
"""
function face_location(dm::DMStag, axis::Integer)
    check_axis(dm, axis)
    axis == 1 && return LibPETSc.DMSTAG_LEFT
    axis == 2 && return LibPETSc.DMSTAG_DOWN
    return LibPETSc.DMSTAG_BACK
end

"""
    element_location(dm::DMStag)

The location of the element interior, `DMSTAG_ELEMENT` in every dimension.

See also [`vertex_location`](@ref), [`face_location`](@ref), [`edge_location`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencilLocation"))
"""
element_location(::DMStag) = LibPETSc.DMSTAG_ELEMENT

function check_axis(::DMStag{PetscLib, N}, axis::Integer) where {PetscLib, N}
    1 <= axis <= N || throw(ArgumentError("axis $axis is not an axis of a $(N)D DMStag"))
    return nothing
end

"""
    stencil(dm::DMStag, loc::LibPETSc.DMStagStencilLocation, I; dof = 0)

The `LibPETSc.DMStagStencil` addressing component `dof` at location `loc` of
element `I`. `I` is a `CartesianIndex{N}` or an `NTuple{N, Integer}` with the
1-based element indices [`corners`](@ref) and [`ghost_corners`](@ref) use, and the
stencil holds them 0-based, as PETSc reads them. Indices are not checked, so a ghost
element (0 or `N + 1` on a periodic axis) is valid. `dof` is PETSc's 0-based
component number, as in [`dof_slot`](@ref).

Allocation free, for use inside assembly loops.

```julia
lower = corners(dm).lower
row = stencil(dm, face_location(dm, 1), lower)             # the first owned x face
col = stencil(dm, element_location(dm), lower; dof = 1)    # second element component
```

See also [`set_values!`](@ref), [`zero_rows_local!`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencil"))
"""
@inline function stencil(
    dm::DMStag{PetscLib, N},
    loc::LibPETSc.DMStagStencilLocation,
    I::NTuple{N, Integer};
    dof::Integer = 0,
) where {PetscLib, N}
    T = inttype(PetscLib)
    i = ntuple(d -> d <= N ? T(I[d] - 1) : zero(T), Val(3))
    return LibPETSc.DMStagStencil(loc, i[1], i[2], i[3], T(dof))
end

@inline stencil(
    dm::DMStag{PetscLib, N},
    loc::LibPETSc.DMStagStencilLocation,
    I::CartesianIndex{N};
    dof::Integer = 0,
) where {PetscLib, N} = stencil(dm, loc, Tuple(I); dof)

# ============================================================================
#   Operations on stencils
# ============================================================================

"""
    set_values!(J::AbstractPetscMat, dm::DMStag, rows, cols, vals, mode = INSERT_VALUES)

Write the dense block `vals` into `J` at the stencils `rows` × `cols` and return
`J`. `vals` is row-major, with `length(rows) * length(cols)` entries. Under
`ADD_VALUES`, entries whose row and column repeat are summed. `rows`, `cols` and
`vals` can be any `AbstractVector`; one that is not a `Vector` is copied first.

# External Links
$(doc_external("DMStag/DMStagMatSetValuesStencil"))
"""
function set_values!(
    J::AbstractPetscMat{PetscLib},
    dm::DMStag{PetscLib},
    rows::AbstractVector{LibPETSc.DMStagStencil},
    cols::AbstractVector{LibPETSc.DMStagStencil},
    vals::AbstractVector,
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    length(vals) == length(rows) * length(cols) || throw(
        DimensionMismatch(
            "a $(length(rows))x$(length(cols)) block needs " *
            "$(length(rows) * length(cols)) values, got $(length(vals))",
        ),
    )
    T = inttype(PetscLib)
    LibPETSc.DMStagMatSetValuesStencil(
        getlib(PetscLib), dm, J,
        T(length(rows)), stencil_vector(rows),
        T(length(cols)), stencil_vector(cols),
        value_vector(PetscLib, vals), mode,
    )
    return J
end

"""
    set_values!(v::AbstractPetscVec, dm::DMStag, positions, vals, mode = INSERT_VALUES)

Write `vals` into `v` at the stencils `positions` and return `v`. Both can be any
`AbstractVector` of the same length; one that is not a `Vector` is copied first.

# External Links
$(doc_external("DMStag/DMStagVecSetValuesStencil"))
"""
function set_values!(
    v::AbstractPetscVec{PetscLib},
    dm::DMStag{PetscLib},
    positions::AbstractVector{LibPETSc.DMStagStencil},
    vals::AbstractVector,
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    length(vals) == length(positions) || throw(
        DimensionMismatch("$(length(positions)) positions, but $(length(vals)) values"),
    )
    LibPETSc.DMStagVecSetValuesStencil(
        getlib(PetscLib), dm, v, inttype(PetscLib)(length(positions)),
        stencil_vector(positions), value_vector(PetscLib, vals), mode,
    )
    return v
end

"""
    zero_rows_local!(J::AbstractPetscMat, dm::DMStag, rows, diag = 1; x = nothing, b = nothing)

[`zero_rows_local!`](@ref) with the rows given as stencils of `dm`, which `J` was
created from. Collective: a rank that owns none of the rows passes an empty vector.

# External Links
$(doc_external("DMStag/DMStagStencilToIndexLocal"))
"""
function zero_rows_local!(
    J::AbstractPetscMat{PetscLib},
    dm::DMStag{PetscLib, N},
    rows::AbstractVector{LibPETSc.DMStagStencil},
    diag = 1;
    x = nothing,
    b = nothing,
) where {PetscLib, N}
    T = inttype(PetscLib)
    rows_0b = LibPETSc.DMStagStencilToIndexLocal(
        getlib(PetscLib), dm, T(N), T(length(rows)), stencil_vector(rows),
    )
    return zero_rows_local!(J, rows_0b, diag; x, b)
end

"""
    LibPETSc.IS(dm::DMStag, loc => dof, ...)
    LibPETSc.IS(dm::DMStag, pairs::AbstractVector{<:Pair})

The index set, in the global numbering, of every point of `dm` at the given
locations and components: each pair is a `LibPETSc.DMStagStencilLocation` and a
0-based component, as [`stencil`](@ref) takes them. The caller owns the result. It
suits [`set_fieldsplit_is!`](@ref):

```julia
flow = LibPETSc.IS(dm, face_location(dm, 1) => 0, face_location(dm, 2) => 0,
                   element_location(dm) => 0)
set_fieldsplit_is!(pc(ksp), "flow", flow)
```

# External Links
$(doc_external("DMStag/DMStagCreateISFromStencils"))
"""
LibPETSc.IS(dm::DMStag, first::Pair, rest::Pair...) = LibPETSc.IS(dm, [first, rest...])

function LibPETSc.IS(dm::DMStag{PetscLib}, pairs::AbstractVector{<:Pair}) where {PetscLib}
    T = inttype(PetscLib)
    stencils = LibPETSc.DMStagStencil[
        LibPETSc.DMStagStencil(loc, zero(T), zero(T), zero(T), T(dof)) for (loc, dof) in pairs
    ]
    return LibPETSc.DMStagCreateISFromStencils(getlib(PetscLib), dm, T(length(stencils)), stencils)
end

# The C calls take a `Vector`; any other vector is copied into one
stencil_vector(v::Vector{LibPETSc.DMStagStencil}) = v
stencil_vector(v::AbstractVector{LibPETSc.DMStagStencil}) = collect(v)
value_vector(::Type{PetscLib}, v::Vector) where {PetscLib} =
    eltype(v) === scalartype(PetscLib) ? v : Vector{scalartype(PetscLib)}(v)
value_vector(::Type{PetscLib}, v::AbstractVector) where {PetscLib} =
    Vector{scalartype(PetscLib)}(v)

# ============================================================================
#   Coordinates
# ============================================================================

"""
    with_product_coordinates(f, dm::DMStag)

Call `f(x)`, `f(x, y)` or `f(x, y, z)` with the local coordinate arrays of `dm`,
one per axis, and return what `f` returns. It needs product coordinates, which
[`set_uniform_coordinates!`](@ref) sets. The arrays are PETSc's own and are read
only: they are handed back when `f` returns or throws, and writing to them, or
using them after `f`, is undefined.

Each array is indexed `[i, slot]`: `i` is the 1-based element index, ghost
elements included, as [`ghost_corners`](@ref) gives it, and `slot` is 1 for the
coordinate of the element's lower face and 2 for its centre.

```julia
with_product_coordinates(dm) do x, y
    x[i, 1], x[i, 2]   # x of the lower face and of the centre of element i
end
```

# External Links
$(doc_external("DMStag/DMStagGetProductCoordinateArraysRead"))
$(doc_external("DMStag/DMStagRestoreProductCoordinateArraysRead"))
"""
function with_product_coordinates(f, dm::DMStag{PetscLib, N}) where {PetscLib, N}
    lib = getlib(PetscLib)
    x, y, z = LibPETSc.DMStagGetProductCoordinateArraysRead(lib, dm)
    try
        N == 1 && return f(x)
        N == 2 && return f(x, y)
        return f(x, y, z)
    finally
        LibPETSc.DMStagRestoreProductCoordinateArraysRead(lib, dm, x, y, z)
    end
end
