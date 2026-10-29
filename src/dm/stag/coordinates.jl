# ============================================================================
#   Coordinates
# ============================================================================

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
