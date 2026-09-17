# Vec

PETSc vectors (`Vec`) are the fundamental building blocks for storing solution data, right-hand sides, and other distributed arrays. PETSc.jl provides a Julia-friendly interface that makes `Vec` objects behave like native Julia arrays.

## Overview

PETSc vectors support:
- **Distributed parallel storage**: Split across MPI processes
- **Sequential storage**: For serial computations
- **Ghost points**: For communication in stencil operations
- **Julia array interface**: Use familiar indexing and broadcasting syntax

## Creating Vectors

### Sequential Vectors

```julia
# Create a sequential vector of length n
v = PetscVec(petsclib, n)

# Wrap an existing Julia array (no copy)
julia_array = zeros(100)
v = PetscVec(petsclib, julia_array)

# `destroy!` releases it; a vector on MPI.COMM_SELF also gets a finalizer
PETSc.destroy!(v)
```

### From DM Objects

```julia
# Create global and local vectors from a DM
gvec = PETSc.global_vec(dm)
lvec = PETSc.local_vec(dm)
```

`PetscVec` replaces v0.4's `VecSeq` and `as_petsc_vec`: construction goes
through the type ([naming conventions](naming.md), §6). The old spellings still
work in v0.5 and warn once.

## Julia Array Interface

PETSc vectors implement the Julia array interface:

```julia
v[1] = 1.0           # Set single element
v[1:10] .= 2.0       # Set range
x = v[5]             # Get element
length(v)            # Get length
size(v)              # Get size tuple
```

## Assembly

After setting values, vectors must be assembled:

```julia
v[1] = 1.0
v[2] = 2.0
PETSc.assemble!(v)  # Finalize vector assembly
```

## Ghost Point Updates

For vectors with ghost points (from DMDA/DMStag):

```julia
# Update ghost values from neighboring processes
PETSc.ghost_update!(vec, PETSc.INSERT_VALUES, PETSc.SCATTER_FORWARD)

# Or use begin/end for non-blocking:
PETSc.ghost_update_begin!(vec, PETSc.INSERT_VALUES, PETSc.SCATTER_FORWARD)
# ... do other work ...
PETSc.ghost_update_end!(vec, PETSc.INSERT_VALUES, PETSc.SCATTER_FORWARD)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["vec.jl"]
```
