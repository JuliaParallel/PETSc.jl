# Mat

PETSc matrices (`Mat`) provide sparse and dense matrix storage with efficient parallel operations. They are essential for discretizing PDEs and setting up linear/nonlinear systems.

## Overview

PETSc matrices support:
- **Sparse formats**: AIJ (CSR), BAIJ (block CSR), and more
- **Dense format**: For small matrices or dense operations
- **Parallel distribution**: Row-based distribution across MPI processes
- **Matrix-free operations**: Via MatShell for custom operators

`PetscMat` is the one constructor: it replaces v0.4's `MatSeqAIJ`,
`MatSeqDense`, `MatCreateSeqAIJ`, `MatSeqAIJWithArrays` and `MatAIJ`
([naming conventions](naming.md), §6). The old spellings still work in v0.5 and
warn once.

## Creating Matrices

### Sparse Matrices (AIJ/CSR Format)

```julia
# Create sparse matrix with estimated non-zeros per row
A = PetscMat(petsclib, num_rows, num_cols, nnz_per_row)

# From Julia SparseMatrixCSC
using SparseArrays
S = sprand(100, 100, 0.1)
A = PetscMat(petsclib, MPI.COMM_SELF, S)

# With varying non-zeros per row
nnz = petsclib.PetscInt[5, 3, 4]  # One value per row
A = PetscMat(petsclib, num_rows, num_cols, nnz)
```

### Dense Matrices

```julia
# Wrap a Julia matrix (no copy)
julia_mat = rand(10, 10)
A = PetscMat(petsclib, julia_mat)
```

### From DM Objects

```julia
# Create matrix with sparsity pattern from DM
A = LibPETSc.DMCreateMatrix(petsclib, dm)
```

### Matrix Shell (Matrix-Free)

```julia
# Create a shell matrix with custom mult operation, y = mult_function(y, x)
A = PETSc.MatShell(petsclib, mult_function, MPI.COMM_SELF, local_rows, local_cols)
```

## Setting Values

```julia
# Set individual element. `setindex!` is 1-based, and converts for you
A[i, j] = value

# `set_values!` is the bulk route, and it keeps PETSc's 0-based indices: the
# `_0b` in the parameter names says so (naming conventions, §12.1)
PETSc.set_values!(A, rows_0b, cols_0b, values, PETSc.INSERT_VALUES)

# The same call with `MatStencil` rows and columns, for stencil-based assembly
PETSc.set_values!(A, stencils_0b, stencils_0b, values, PETSc.INSERT_VALUES)
```

## Assembly

Matrices must be assembled after setting values:

```julia
# Set all values first
A[1, 1] = 2.0
A[1, 2] = -1.0
# ...

# Then assemble
PETSc.assemble!(A)
```

## Common Operations

```julia
size(A)                    # Get (rows, cols)
PETSc.ownership_range(A)   # The 1-based rows owned by this process
PETSc.setup!(A)            # Complete matrix setup
PETSc.destroy!(A)          # Release it
```

`ownership_range(A)` is 1-based and takes no second argument. v0.4's
`ownershiprange(A, false)` still works and warns; it is a `MethodError` in v0.6
([naming conventions](naming.md), §12.1).

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["mat.jl"]
```
