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
v = VecSeq(petsclib, n)

# Wrap an existing Julia array (no copy)
julia_array = zeros(100)
v = VecSeq(petsclib, julia_array)
```

### From DM Objects

```julia
# Create global and local vectors from a DM
global_vec = DMGlobalVec(dm)
local_vec = DMLocalVec(dm)
```

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
assemble!(v)  # Finalize vector assembly
```

## Ghost Point Updates

For vectors with ghost points (from DMDA/DMStag):

```julia
# Update ghost values from neighboring processes
ghostupdate!(vec, INSERT_VALUES, SCATTER_FORWARD)

# Or use begin/end for non-blocking:
ghostupdatebegin!(vec, INSERT_VALUES, SCATTER_FORWARD)
# ... do other work ...
ghostupdateend!(vec, INSERT_VALUES, SCATTER_FORWARD)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["vec.jl"]
```
