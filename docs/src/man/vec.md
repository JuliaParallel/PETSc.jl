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

## Linear Algebra

The vector operations of `LinearAlgebra` and the arithmetic operators call PETSc directly. All of them are collective, so every rank must make the call:

```julia
using LinearAlgebra

w = copy(v)            # VecDuplicate + VecCopy
dot(x, y)              # xᴴy, the global result on every rank
axpy!(α, x, y)         # y = αx + y
axpby!(α, x, β, y)     # y = αx + βy
rmul!(v, α)            # v = αv
z = x + y              # a new vector; so are x - y, -x, α * x and x / α
```

`copy` and the operators return a new `PetscVec` with the layout of their first vector argument. On more than one rank it has no finalizer, so call `destroy!` on it on every rank. Broadcasting is different for now: `y .= α .* x .+ y` updates `y` in place and works on any number of ranks, while an out-of-place broadcast such as `2 .* x` returns a Julia `Vector` and needs an undistributed vector.

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
PETSc.ghost_update!(vec; mode = PETSc.INSERT_VALUES, scatter = PETSc.SCATTER_FORWARD)

# Or use begin/end for non-blocking:
PETSc.ghost_update_begin!(vec; mode = PETSc.INSERT_VALUES, scatter = PETSc.SCATTER_FORWARD)
# ... do other work ...
PETSc.ghost_update_end!(vec; mode = PETSc.INSERT_VALUES, scatter = PETSc.SCATTER_FORWARD)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["vec/vec.jl", "vec/array_interface.jl", "vec/local_arrays.jl", "vec/ghost.jl", "vec/linalg.jl"]
```
