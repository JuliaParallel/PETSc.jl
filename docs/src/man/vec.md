# Vec

PETSc vectors (`Vec`) are the fundamental building blocks for storing solution data, right-hand sides, and other distributed arrays. PETSc.jl provides a Julia-friendly interface that makes `Vec` objects behave like native Julia arrays.

## Overview

PETSc vectors support:
- **Distributed parallel storage**: Split across MPI processes
- **Sequential storage**: For serial computations
- **Ghost points**: For communication in stencil operations
- **Julia array interface**: a `PetscVec` is an `AbstractVector`, with bulk indexing, broadcasting and reductions

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

A `PetscVec` is an `AbstractVector` of the library's scalar type, so `eltype`, `size`, `length`, `axes` and the generic `AbstractVector` methods work as for any vector:

```julia
v[1:3] = [1.0, 2.0, 3.0]   # set a range (then assemble!)
v[2]                       # one entry
a = v[1:3]                 # a range, as a Vector
length(v)                  # global length, an Int
size(v)                    # (length(v),)
collect(v)                 # every entry, as a Vector
```

`v[i]` and `v[i] = x` are one C call each and reach only the entries **owned by this rank**; PETSc raises an error for any other. So on one rank a vector behaves like a `Vector`, and on several ranks code that touches only its own entries works. A loop over every entry is correct but slow, since each step is a C call: use [`with_local_array!`](@ref) for the entries each rank owns, which hands you a plain Julia array.

```julia
PETSc.with_local_array!(v) do a
    a[1] = 1.0
end
```

Reductions are collective and give every rank the global answer. `sum`, `maximum`, `minimum`, `dot`, `norm` and `isapprox` call PETSc directly. `any`, `all`, `count` and `mapreduce` (so `sum(f, v)`, `prod`, `extrema` and `reduce`) run the function on the entries each rank owns and combine the results. `collect(v)`, `Vector(v)` and `copyto!(::Vector, v)` give every rank the whole vector, so every rank must call them, and each needs memory for all `length(v)` entries. `v == w` compares the entries, also collectively. `hash(v)` throws, so key a `Dict` of vectors with an `IdDict`.

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

`copy` and the operators return a new `PetscVec` with the layout of their first vector argument. On more than one rank it has no finalizer, so call `destroy!` on it on every rank. Broadcasting follows the same rule: `y .= α .* x .+ y` updates `y` in place, and an out-of-place broadcast such as `2 .* x` returns a new `PetscVec` with the layout of the first `PetscVec` in the expression. Both run on the entries each rank owns and work on any number of ranks. A broadcast whose entries are not numbers (`x .> 0` gives `Bool`) throws, since a `PetscVec` holds only the library's scalar type: broadcast over `x[:]` for those.

## Assembly

After setting values, vectors must be assembled:

```julia
v[1:2] = [1.0, 2.0]
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
