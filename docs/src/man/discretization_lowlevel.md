# Discretization (PetscFE, PetscDS, PetscSpace, PetscDualSpace, PetscFV) - Low-level Interface

These classes describe finite-element and finite-volume discretizations on a `DMPlex`:

- **PetscSpace**: the approximation space (polynomials of a given degree, tensor/sum spaces)
- **PetscDualSpace**: the dual space that defines the degrees of freedom (Lagrange nodes, ...)
- **PetscFE**: a finite element = space + dual space + quadrature; attached to a DM field
- **PetscDS**: the discrete system: residual/Jacobian point functions, boundary conditions and
  constants for all fields of a DM (see [DMPlex](dmplex.md) for the Julia callback helpers)
- **PetscFV** and **PetscLimiter**: finite-volume discretizations and slope limiters
- **PetscConvEst**: convergence-rate estimation by mesh refinement

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

dim, ncomp, simplex = 2, 1, LibPETSc.PETSC_TRUE
# P2 Lagrange element with a quadrature order chosen by PETSc (qorder = -1 -> PETSC_DETERMINE)
fe = LibPETSc.PetscFECreateLagrange(petsclib, MPI.COMM_SELF, dim, ncomp, simplex, 2, -1)
LibPETSc.PetscObjectSetName(petsclib, fe, "velocity")

nb = LibPETSc.PetscFEGetDimension(petsclib, fe)          # basis functions per cell
sp = LibPETSc.PetscFEGetBasisSpace(petsclib, fe)         # owned by the FE
deg = LibPETSc.PetscSpaceGetDegree(petsclib, sp)         # (degree, maxdegree)

LibPETSc.PetscFEDestroy(petsclib, fe)
PETSc.finalize(petsclib)
```

Elements are attached to a DM with `DMSetField` / `DMCreateDS`; the high-level helpers
`PETSc.fe_create_default`, `PETSc.set_field!` and `PETSc.create_ds!` on the
[DMPlex](dmplex.md) page wrap this workflow.

## Function Reference

### PetscFE

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscFE_wrappers.jl"]
Order   = [:function]
```

### PetscDS

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDS_wrappers.jl"]
Order   = [:function]
```

### PetscSpace

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscSpace_wrappers.jl"]
Order   = [:function]
```

### PetscDualSpace

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDualSpace_wrappers.jl"]
Order   = [:function]
```

### PetscFV

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscFV_wrappers.jl"]
Order   = [:function]
```

### PetscLimiter

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscLimiter_wrappers.jl"]
Order   = [:function]
```

### PetscConvEst

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscConvEst_wrappers.jl"]
Order   = [:function]
```
