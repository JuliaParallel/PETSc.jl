# Utilities (PetscRandom, PetscLayout, Subcommunicators, ...) - Low-level Interface

Small helper classes used throughout PETSc:

- **PetscRandom**: parallel random number generators
- **PetscLayout**: the parallel row distribution of vectors and matrices
- **PetscSubcomm** / **PetscShmComm**: splitting a communicator, shared-memory communicators
- **PetscBag**: option-settable parameter structs
- **PetscToken**: string tokenizer
- **PetscSegBuffer**, **PetscHeap**, **PetscIntStack**: growable buffers, heaps and stacks
- **PetscKDTree**, **PetscGridHash**: spatial search structures used by `DMPlex`/`DMSwarm`
- **PetscOmpCtrl**, **PetscMatlabEngine**, **PetscBench**: OpenMP control, MATLAB engine, benchmarks

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Random numbers in [0, 10)
rnd = LibPETSc.PetscRandomCreate(petsclib, MPI.COMM_SELF)
LibPETSc.PetscRandomSetType(petsclib, rnd, LibPETSc.PETSCRANDER48)
LibPETSc.PetscRandomSetInterval(petsclib, rnd, 0.0, 10.0)
x = LibPETSc.PetscRandomGetValue(petsclib, rnd)          # PetscScalar
LibPETSc.PetscRandomDestroy(petsclib, rnd)

# Parallel layout of 100 rows
layout = LibPETSc.PetscLayoutCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.PetscLayoutSetSize(petsclib, layout, 100)
LibPETSc.PetscLayoutSetUp(petsclib, layout)
rstart, rend = LibPETSc.PetscLayoutGetRange(petsclib, layout)   # rows owned by this rank
LibPETSc.PetscLayoutDestroy(petsclib, layout)

PETSc.finalize(petsclib)
```

## Function Reference

### PetscRandom

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscRandom_wrappers.jl"]
Order   = [:function]
```

### PetscLayout

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscLayout_wrappers.jl"]
Order   = [:function]
```

### PetscSubcomm and PetscShmComm

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Petsccomm_wrappers.jl"]
Order   = [:function]
```

### PetscBag

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscBag_wrappers.jl"]
Order   = [:function]
```

### PetscToken

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscToken_wrappers.jl"]
Order   = [:function]
```

### PetscSegBuffer

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscSegBuffer_wrappers.jl"]
Order   = [:function]
```

### PetscHeap

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscHeap_wrappers.jl"]
Order   = [:function]
```

### PetscIntStack

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscIntStack_wrappers.jl"]
Order   = [:function]
```

### PetscKDTree

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscKDTree_wrappers.jl"]
Order   = [:function]
```

### PetscGridHash

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscGridHash_wrappers.jl"]
Order   = [:function]
```

### PetscOmpCtrl

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscOmpCtrl_wrappers.jl"]
Order   = [:function]
```

### PetscMatlabEngine

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscMatlabEngine_wrappers.jl"]
Order   = [:function]
```

### PetscBench

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscBench_wrappers.jl"]
Order   = [:function]
```
