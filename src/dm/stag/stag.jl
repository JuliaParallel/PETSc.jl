# Helper to convert points_per_proc tuples to PetscInt
to_petscint_tuple(t::Tuple, PetscInt) = map(arr -> PetscInt.(arr), t)

"""
    da = DMStag(
        petsclib::PetscLib
        comm::MPI.Comm,
        boundary_type::NTuple{D, DMBoundaryType},
        global_dim::NTuple{D, Integer},
        dof_per_node::NTuple{1 + D, Integer};
        stencil_width::Integer = 1,
        stencil_type = :box,
        points_per_proc::Tuple,
        processors::Tuple,
        set_from_options::Bool = true,
        setup::Bool = true,
        prefix::AbstractString = "",
        options...
    )

Creates a `D`-dimensional distributed staggered array.

The Tuple `dof_per_node` specifies how many degrees of freedom are at all the
staggerings in the order:
 - 1D: `(vertex, element)`
 - 2D: `(vertex, edge, element)`
 - 3D: `(vertex, edge, face, element)`

If keyword argument `points_per_proc[k] isa Vector{petsclib.PetscInt}` then this
specifies the points per processor in dimension `k`.

If keyword argument `processors[k] isa Integer` then this specifies the number of
processors used in dimension `k`; ignored when `D == 1`.

`stencil_width` is the number of ghost points on each side. `stencil_type` is
`:box`, `:star` or `:none` (or a `DMStagStencilType`).

With `set_from_options = true`, [`set_from_options!`](@ref) is applied with the
remaining keywords `options...` as the options database, under `prefix`. With
`setup = true`, [`setup!`](@ref) is called.

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
    dof_per_node::NTuple{N1, Integer};
    stencil_width::Integer = 1,
    stencil_type = :box,
    points_per_proc::Union{Tuple, Nothing} = nothing,
    processors = nothing,
    set_from_options::Bool = true,
    setup::Bool = true,
    prefix::AbstractString = "",
    options...,
) where {PetscLib, N, N1}
    check_initialized(getlib(PetscLib))
    1 <= N <= 3 || throw(ArgumentError("a DMStag has 1, 2 or 3 dimensions, got $N"))
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
    return finish_dm!(da, comm, options; prefix, set_from_options, setup)
end


"""
    DMStag(dm::DMStag, dof_per_node; set_from_options = true, setup = true, prefix = "", options...)

A `DMStag` compatible with `dm` (same dimension, communicator and layout) but
with the degrees of freedom per stratum given by `dof_per_node`. The keywords
work as in the other `DMStag` constructor.

# External Links
$(doc_external("DMSTAG/DMStagCreateCompatibleDMStag"))
"""
function DMStag(
    dm::DMStag{PetscLib, N},
    dof_per_node::Union{NTuple{2, Integer}, NTuple{3, Integer}, NTuple{4, Integer}};
    set_from_options::Bool = true,
    setup::Bool = true,
    prefix::AbstractString = "",
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
    return finish_dm!(dmnew, comm(dm), options; prefix, set_from_options, setup)
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

function Base.size(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    M, Ny, P = LibPETSc.DMStagGetGlobalSizes(PetscLib, dm)
    gsize = (Int(M), Int(Ny), Int(P))
    return ntuple(i -> gsize[i], Val(N))
end
