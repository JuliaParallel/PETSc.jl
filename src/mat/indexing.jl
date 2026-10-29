# Matrix indexing - set single value
function Base.setindex!(m::AbstractPetscMat{PetscLib}, val, i::Integer, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc (Julia uses 1-based)
    row = PetscInt(i - 1)
    col = PetscInt(j - 1)
    value = PetscScalar(val)
    
    # Use MatSetValues for single entry
    # MatSetValues(mat, m, idxm, n, idxn, y, INSERT_VALUES)
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), [row], PetscInt(1), [col], [value], LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - set multiple values with vectors of indices  
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractMatrix, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc
    petsc_rows = PetscInt[r - 1 for r in rows]
    petsc_cols = PetscInt[c - 1 for c in cols]

    size(vals) == (length(rows), length(cols)) || throw(
        DimensionMismatch(
            "block is $(size(vals)) but the index ranges are " *
            "$(length(rows))x$(length(cols))",
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
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    row = PetscInt(i - 1)
    petsc_cols = PetscInt[c - 1 for c in cols]
    petsc_vals = PetscScalar.(vals)
    
    ncols = PetscInt.(length(petsc_cols))
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), [row], ncols, petsc_cols, petsc_vals, LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - set column with vector of values  
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractVector, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    petsc_rows = PetscInt[r - 1 for r in rows]
    col = PetscInt(j - 1)
    petsc_vals = PetscScalar.(vals)
    
    nrows = PetscInt.(length(petsc_rows))
    LibPETSc.MatSetValues(PetscLib, m, nrows, petsc_rows, PetscInt(1), [col], petsc_vals, LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - get single value
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc (Julia uses 1-based)
    row = PetscInt(i - 1)
    col = PetscInt(j - 1)
    
    # Use MatGetValues for single entry
    values = Vector{PetscScalar}(undef, 1)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), [row], PetscInt(1), [col], values)
    
    return values[1]
end

# Matrix indexing - get block of values
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc
    petsc_rows = PetscInt[r - 1 for r in rows]
    petsc_cols = PetscInt[c - 1 for c in cols]
    
    # Use MatGetValues for block of entries
    nrows = PetscInt.(length(petsc_rows))
    ncols = PetscInt.(length(petsc_cols))

    # PETSc returns values in row-major order
    values = Vector{PetscScalar}(undef, nrows * ncols)
    LibPETSc.MatGetValues(PetscLib, m, nrows, petsc_rows, ncols, petsc_cols, values)
    
    # Reshape to Julia matrix (column-major)
    result = Matrix{PetscScalar}(undef, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            # Convert from PETSc row-major to Julia column-major
            petsc_idx = (i-1) * ncols + (j-1) + 1
            result[i, j] = values[petsc_idx]
        end
    end
    
    return result
end

# Matrix indexing - get row values at specified columns
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    row = PetscInt(i - 1)
    petsc_cols = PetscInt[c - 1 for c in cols]

    ncols = PetscInt.(length(petsc_cols))
    values = Vector{PetscScalar}(undef, ncols)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), [row], ncols, petsc_cols, values)
    
    return values
end

# Matrix indexing - get column values at specified rows
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    petsc_rows = PetscInt[r - 1 for r in rows]
    col = PetscInt(j - 1)
    
    nrows = length(petsc_rows)
    values = Vector{PetscScalar}(undef, nrows)
    LibPETSc.MatGetValues(PetscLib, m, nrows, petsc_rows, PetscInt(1), [col], values)
    return values
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
