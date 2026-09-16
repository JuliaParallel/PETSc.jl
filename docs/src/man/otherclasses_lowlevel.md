# PF, PetscPartitioner, PetscRegressor, PetscDA, Characteristic - Low-level Interface

- **PF**: mathematical functions applied pointwise to arrays and vectors (`PFApply`)
- **PetscPartitioner**: graph partitioners used by `DMPlexDistribute` (simple, ParMETIS, PTScotch, ...)
- **PetscRegressor**: regression models (linear, ridge, lasso) fitted with Tao
- **PetscDA**: data assimilation (ensemble methods), new in PETSc 3.25
- **Characteristic**: the method of characteristics for semi-Lagrangian advection on a `DMDA`

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# A pointwise function R^1 -> R^1: the identity
pf = LibPETSc.PFCreate(petsclib, MPI.COMM_SELF, 1, 1)
LibPETSc.PFSetType(petsclib, pf, LibPETSc.PFIDENTITY, C_NULL)
x = [1.0, 2.0, 3.0]
y = LibPETSc.PFApply(petsclib, pf, 3, x)     # y == x
LibPETSc.PFDestroy(petsclib, pf)

# The partitioner used when distributing a DMPlex
part = LibPETSc.PetscPartitionerCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.PetscPartitionerSetType(petsclib, part, LibPETSc.PETSCPARTITIONERSIMPLE)
LibPETSc.PetscPartitionerDestroy(petsclib, part)

PETSc.finalize(petsclib)
```

## Function Reference

### PF

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PF_wrappers.jl"]
Order   = [:function]
```

### PetscPartitioner

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscPartitioner_wrappers.jl"]
Order   = [:function]
```

### PetscRegressor

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscRegressor_wrappers.jl"]
Order   = [:function]
```

### PetscDA

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDA_wrappers.jl"]
Order   = [:function]
```

### Characteristic

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Characteristic_wrappers.jl"]
Order   = [:function]
```
