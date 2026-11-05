# Matrix indexing - set single value
function Base.setindex!(m::AbstractPetscMat{PetscLib}, val, i::Integer, j::Integer) where {PetscLib}
    checkbounds(m, i, j)
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), PetscInt[i - 1], PetscInt(1), PetscInt[j - 1],
                          PetscScalar[val], LibPETSc.INSERT_VALUES)
    return m
end

# Matrix indexing - set multiple values with vectors of indices
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractMatrix, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    checkbounds(m, rows, cols)
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    petsc_rows = indices_0b(PetscInt, rows)
    petsc_cols = indices_0b(PetscInt, cols)

    size(vals) == (length(petsc_rows), length(petsc_cols)) || throw(
        DimensionMismatch(
            "block is $(size(vals)) but the index ranges are " *
            "$(length(petsc_rows))x$(length(petsc_cols))",
        ),
    )

    # MatSetValues reads its value array row by row, so the block is transposed
    # before it is flattened: Julia lays a matrix out column by column.
    petsc_vals = vec(permutedims(PetscScalar.(vals)))

    nrows = PetscInt(length(petsc_rows))
    ncols = PetscInt(length(petsc_cols))
    LibPETSc.MatSetValues(PetscLib, m, nrows, petsc_rows, ncols, petsc_cols, petsc_vals, LibPETSc.INSERT_VALUES)
    return m
end

# Matrix indexing - set row with vector of values
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractVector, i::Integer, cols::AbstractVector{<:Integer}) where {PetscLib}
    checkbounds(m, i, cols)
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    petsc_cols = indices_0b(PetscInt, cols)
    length(vals) == length(petsc_cols) || throw(
        DimensionMismatch("tried to assign $(length(vals)) values to $(length(petsc_cols)) columns"),
    )
    petsc_vals = c_vector(PetscScalar, vals)
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), PetscInt[i - 1], PetscInt(length(petsc_cols)),
                          petsc_cols, petsc_vals, LibPETSc.INSERT_VALUES)
    return m
end

# Matrix indexing - set column with vector of values
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractVector, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    checkbounds(m, rows, j)
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    petsc_rows = indices_0b(PetscInt, rows)
    length(vals) == length(petsc_rows) || throw(
        DimensionMismatch("tried to assign $(length(vals)) values to $(length(petsc_rows)) rows"),
    )
    petsc_vals = c_vector(PetscScalar, vals)
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(length(petsc_rows)), petsc_rows, PetscInt(1),
                          PetscInt[j - 1], petsc_vals, LibPETSc.INSERT_VALUES)
    return m
end

# Matrix indexing - get single value
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, j::Integer) where {PetscLib}
    checkbounds(m, i, j)
    PetscInt = inttype(PetscLib)
    values = Vector{scalartype(PetscLib)}(undef, 1)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), PetscInt[i - 1], PetscInt(1), PetscInt[j - 1], values)
    return values[1]
end

# Matrix indexing - get block of values
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    checkbounds(m, rows, cols)
    PetscInt = inttype(PetscLib)
    petsc_rows = indices_0b(PetscInt, rows)
    petsc_cols = indices_0b(PetscInt, cols)
    nrows = length(petsc_rows)
    ncols = length(petsc_cols)

    # PETSc returns the block row by row; Julia stores a matrix column by column.
    values = Vector{scalartype(PetscLib)}(undef, nrows * ncols)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(nrows), petsc_rows, PetscInt(ncols), petsc_cols, values)
    return index_shaped(permutedims(reshape(values, ncols, nrows)), rows, cols)
end

# Matrix indexing - get row values at specified columns
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, cols::AbstractVector{<:Integer}) where {PetscLib}
    checkbounds(m, i, cols)
    PetscInt = inttype(PetscLib)
    petsc_cols = indices_0b(PetscInt, cols)
    values = Vector{scalartype(PetscLib)}(undef, length(petsc_cols))
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), PetscInt[i - 1], PetscInt(length(petsc_cols)), petsc_cols, values)
    return index_shaped(values, cols)
end

# Matrix indexing - get column values at specified rows
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    checkbounds(m, rows, j)
    PetscInt = inttype(PetscLib)
    petsc_rows = indices_0b(PetscInt, rows)
    values = Vector{scalartype(PetscLib)}(undef, length(petsc_rows))
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(length(petsc_rows)), petsc_rows, PetscInt(1), PetscInt[j - 1], values)
    return index_shaped(values, rows)
end

# Matrix indexing - get entire row
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, ::Colon) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, i, 1:ncols)
end

# Matrix indexing - get entire column  
function Base.getindex(m::AbstractPetscMat{PetscLib}, ::Colon, j::Integer) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, 1:nrows, j)
end

# Matrix indexing - get all values (use with caution for large matrices!)
function Base.getindex(m::AbstractPetscMat{PetscLib}, ::Colon, ::Colon) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, 1:nrows, 1:ncols)
end


function Base.setindex!(
    M::AbstractPetscMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
     ms_i = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(i[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(i[2] - 1),
        PetscInt(i[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(i[4] - 1),
    )
    ms_j = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(j[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(j[2] - 1),
        PetscInt(j[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(j[4] - 1),
    )
    set_values!(M, [ms_i], [ms_j], [PetscScalar(val)]; mode = INSERT_VALUES)
    return M
end

"""
    add_index!(M::AbstractPetscMat, val, i::CartesianIndex{N}, j::CartesianIndex{N})

Add `val` to the entry of `M` at grid row `i` and grid column `j`, for a matrix
created from a DM such as `PetscMat(da)`. The indices are 1-based global grid
indices `(x, y, z, component)`, of which the first `N` are given and the rest
are taken as 1, so a 2D DMDA with several dofs per point needs
`CartesianIndex(x, y, 1, component)`. `M[i, j] = val` with the same indices
inserts instead of adding; PETSc needs an [`assemble!`](@ref) between inserting
and adding, and before `M` is read.

```jldoctest
julia> petsclib = PETSc.petsclibs[1];

julia> PETSc.initialize(petsclib)

julia> da = PETSc.DMDA(petsclib, MPI.COMM_SELF, (PETSc.DM_BOUNDARY_NONE,), (3,), 1);

julia> A = PETSc.PetscMat(da);

julia> PETSc.add_index!(A, 1.0, CartesianIndex(2), CartesianIndex(2));

julia> PETSc.add_index!(A, 2.5, CartesianIndex(2), CartesianIndex(2));

julia> PETSc.assemble!(A);

julia> A[2, 2]
3.5

julia> PETSc.destroy!(A); PETSc.destroy!(da);
```

# External Links
$(doc_external("Mat/MatSetValuesStencil"))
"""
function add_index!(
    M::AbstractPetscMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    ms_i = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(i[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(i[2] - 1),
        PetscInt(i[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(i[4] - 1),
    )
    ms_j = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(j[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(j[2] - 1),
        PetscInt(j[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(j[4] - 1),
    )
    set_values!(M, [ms_i], [ms_j], [PetscScalar(val)]; mode = ADD_VALUES)
    return M
end
