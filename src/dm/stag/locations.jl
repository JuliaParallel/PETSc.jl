# ============================================================================
#   Locations and stencils
# ============================================================================

"""
    slot::Int = dof_slot(dm::DMStag, loc::LibPETSc.DMStagStencilLocation, dof::Integer)

Returns the location `slot` for a degree of freedom `dof` at a given stencil location `loc` in the DMStag `dm`.
Note that the returned `slot` is 1-based for Julia compatibility.
`dof` is PETSc's 0-based component number at that location, as in [`stencil`](@ref).
"""
function dof_slot(dm::DMStag{PetscLib}, loc::LibPETSc.DMStagStencilLocation, dof::Integer) where {PetscLib}
    slot = LibPETSc.DMStagGetLocationSlot(getlib(PetscLib), dm, loc, PetscLib.PetscInt(dof))
    return slot+1
end


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

"""
    on_lower_side(loc::LibPETSc.DMStagStencilLocation, axis::Integer)

Whether `loc` sits on the lower side of its element along `axis`: `true` for
`DMSTAG_LEFT`, `DMSTAG_DOWN_LEFT` and every other location with `LEFT` along axis
1, `DOWN` along axis 2 or `BACK` along axis 3; `false` for `DMSTAG_ELEMENT` and for
the upper and centre sides. Points on the lower side of `axis` have one more index
along it than there are elements, since the upper boundary is stored as the lower
side of element `n + 1`.

It takes no DM, so it cannot check `axis` against the dimension: any `axis` in
`1:3` is accepted, and a 2D location such as `DMSTAG_LEFT` is not on the lower side
along axis 3.

```julia
PETSc.on_lower_side(face_location(dm, 1), 1)   # true
PETSc.on_lower_side(face_location(dm, 1), 2)   # false
```

See also [`face_location`](@ref), [`stencil`](@ref).

# External Links
$(doc_external("DMStag/DMStagStencilLocation"))
"""
function on_lower_side(loc::LibPETSc.DMStagStencilLocation, axis::Integer)
    1 <= axis <= 3 || throw(ArgumentError("axis $axis is not an axis of a DMStag"))
    loc == LibPETSc.DMSTAG_NULL_LOCATION &&
        throw(ArgumentError("DMSTAG_NULL_LOCATION is on no side"))
    # the locations count up in base 3, axis 1 fastest: 0 lower, 1 centre, 2 upper
    return div(Int(loc) - 1, 3^(axis - 1)) % 3 == 0
end

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
