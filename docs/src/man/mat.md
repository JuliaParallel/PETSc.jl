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
# A shell matrix whose multiply is mult_function(y, x), filling y
A = PETSc.MatShell(mult_function, petsclib, MPI.COMM_SELF, local_rows, local_cols)

# With the transpose, the diagonal (for pc_type = jacobi) and z = y + Ax as well
A = PETSc.MatShell(mult_function, petsclib, MPI.COMM_SELF, local_rows, local_cols;
                   mult_transpose = (y, x) -> ..., diagonal = d -> ..., mult_add = (z, x, y) -> ...)
```

With `mult_transpose`, `transpose(A) * x` and `A' * x` work on the shell.

## Setting Values

```julia
# Set individual element. `setindex!` is 1-based, and converts for you
A[i, j] = value

# `set_values!` is the bulk route, and it keeps PETSc's 0-based indices: the
# `_0b` in the parameter names says so (naming conventions, §12.1)
PETSc.set_values!(A, rows_0b, cols_0b, values; mode = PETSc.INSERT_VALUES)

# The same call with `MatStencil` rows and columns, for stencil-based assembly
PETSc.set_values!(A, stencils_0b, stencils_0b, values; mode = PETSc.INSERT_VALUES)
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

A `PetscMat` (and a `MatShell`) is an `AbstractMatrix` of the library's scalar type. `A[i, j]` is one C call and reaches only the rows this rank owns; PETSc raises an error for any other row.

`A == B` compares the entries as Julia matrices do: an entry stored in only one of the two counts as zero, so an explicitly stored zero makes no difference, and `NaN != NaN`. It is collective. [`PETSc.equal(A, B)`](@ref PETSc.equal) is PETSc's own `MatEqual`, which is faster but also compares which entries are stored. `A ≈ B` compares Frobenius norms, as for Julia matrices. `collect(A)` and `Matrix(A)` give every rank the whole matrix as a dense `Matrix`, so every rank must call them and each needs memory for all entries. `sum(A)` uses `MatGetRowSum`. `hash(A)` throws, so key a `Dict` of matrices with an `IdDict`. Other generic reductions (`maximum`, `sum(f, A)`) still read the matrix entry by entry: correct on one rank but slow, and on several ranks they stop at the first row another rank owns.

An assembled `PetscMat` supports `copy(A)`, `copyto!(B, A)`, `A + B`, `A - B`, `-A` and `α * A`. Each returns a new matrix and is collective; `A ± B` has the union of the two nonzero patterns. On more than one rank the result has no finalizer, so call `destroy!` on it on every rank.

`ownership_range(A)` is 1-based and takes no second argument. v0.4's
`ownershiprange(A, false)` still works and warns; it is a `MethodError` in v0.6
([naming conventions](naming.md), §12.1).

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["mat/mat.jl", "mat/indexing.jl", "mat/assembly.jl", "mat/linalg.jl", "mat/shell.jl"]
```
