# PetscObject, Logging and Devices - Low-level Interface

- **PetscObject**: operations shared by every PETSc object (names, options prefixes,
  reference counting, composed objects, `PetscObjectGetComm`). Any Julia handle
  (`PetscVec`, `PetscMat`, `KSP`, ...) is accepted where a `PetscObject` is expected.
- **PetscLog**: performance logging with stages and events (`-log_view`).
- **PetscDevice**: GPU device and stream management (`PetscDeviceContext`).
- **PetscContainer**, **PetscFunctionList**, **PetscDLLibrary**: attaching user data to
  objects, registering implementations and loading shared libraries.

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

v = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
LibPETSc.PetscObjectSetName(petsclib, v, "rhs")
name = LibPETSc.PetscObjectGetName(petsclib, v)           # "rhs"
cls  = LibPETSc.PetscObjectGetClassName(petsclib, v)      # "Vec"
comm = LibPETSc.PetscObjectGetComm(petsclib, v)           # an MPI.Comm

# Logging: register a stage and an event, then view the summary at the end
LibPETSc.PetscLogDefaultBegin(petsclib)
stage = LibPETSc.PetscLogStageRegister(petsclib, "Assembly")
classid = LibPETSc.PetscClassIdRegister(petsclib, "MyApp")
event = LibPETSc.PetscLogEventRegister(petsclib, "MyAssembly", classid)
LibPETSc.PetscLogStagePush(petsclib, stage)
# ... work ...
LibPETSc.PetscLogStagePop(petsclib)
LibPETSc.PetscLogView(petsclib, LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib))

LibPETSc.VecDestroy(petsclib, v)
PETSc.finalize(petsclib)
```

`PetscLogEventBegin`/`PetscLogEventEnd` are C macros and have no wrapper; time your own code
with stages, or pass `-log_view` on the command line to get PETSc's own event timings.

## Function Reference

### PetscObject

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscObject_wrappers.jl"]
Order   = [:function]
```

### PetscLog

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscLog_wrappers.jl"]
Order   = [:function]
```

### PetscDevice

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDevice_wrappers.jl"]
Order   = [:function]
```

### PetscContainer

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscContainer_wrappers.jl"]
Order   = [:function]
```

### PetscFunctionList

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscFunctionList_wrappers.jl"]
Order   = [:function]
```

### PetscDLLibrary

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDLLibrary_wrappers.jl"]
Order   = [:function]
```
