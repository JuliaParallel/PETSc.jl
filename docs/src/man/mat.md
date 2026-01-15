# Mat

PETSc matrices (`Mat`) provide sparse and dense matrix storage with efficient parallel operations. They are essential for discretizing PDEs and setting up linear/nonlinear systems.

## Overview

PETSc matrices support:
- **Sparse formats**: AIJ (CSR), BAIJ (block CSR), and more
- **Dense format**: For small matrices or dense operations
- **Parallel distribution**: Row-based distribution across MPI processes
- **Matrix-free operations**: Via MatShell for custom operators

## Creating Matrices

### Sparse Matrices (AIJ/CSR Format)

```julia
# Create sparse matrix with estimated non-zeros per row
A = MatSeqAIJ(petsclib, num_rows, num_cols, nnz_per_row)

# From Julia SparseMatrixCSC
using SparseArrays
S = sprand(100, 100, 0.1)
A = MatCreateSeqAIJ(petsclib, MPI.COMM_SELF, S)

# With varying non-zeros per row
nnz = PetscInt[5, 3, 4, ...]  # One value per row
A = MatSeqAIJ(petsclib, num_rows, num_cols, nnz)
```

### Dense Matrices

```julia
# Wrap a Julia matrix (no copy)
julia_mat = rand(10, 10)
A = MatSeqDense(petsclib, julia_mat)
```

### From DM Objects

```julia
# Create matrix with sparsity pattern from DM
A = DMCreateMatrix(dm)
```

### Matrix Shell (Matrix-Free)

```julia
# Create a shell matrix with custom mult operation
A = MatShell(petsclib, m, n, mult_function, context)
```

## Setting Values

```julia
# Set individual element (0-based internally, 1-based in Julia)
A[i, j] = value

# Use setvalues! for efficient batch insertion
setvalues!(A, rows, cols, values, INSERT_VALUES)

# For stencil-based assembly
setvalues!(A, stencil_row, stencil_col, value, INSERT_VALUES)
```

## Assembly

Matrices must be assembled after setting values:

```julia
# Set all values first
A[1, 1] = 2.0
A[1, 2] = -1.0
# ...

# Then assemble
assemble!(A)
```

## Common Operations

```julia
size(A)              # Get (rows, cols)
ownershiprange(A)    # Get rows owned by this process
setup!(A)            # Complete matrix setup
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["mat.jl"]
```
