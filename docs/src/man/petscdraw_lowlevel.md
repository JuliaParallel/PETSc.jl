# PetscDraw (Graphics) - Low-level Interface

`PetscDraw` is PETSc's simple graphics layer (X11 windows, image files, or a null device) with
line graphs (`PetscDrawLG`), scatter plots (`PetscDrawSP`), histograms (`PetscDrawHG`), bar
charts (`PetscDrawBar`) and axes. `PetscViewerDrawOpen` on the
[PetscViewer](petscviewer_lowlevel.md) page uses it to plot vectors and matrix structure.

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# "" selects the default display; the type comes from -draw_type (x, image, null, ...)
draw = LibPETSc.PetscDrawCreate(petsclib, MPI.COMM_SELF, "", "Line graph",
                                Cint(0), Cint(0), Cint(400), Cint(300))
LibPETSc.PetscDrawSetType(petsclib, draw, LibPETSc.PETSC_DRAW_NULL)   # headless
LibPETSc.PetscDrawSetFromOptions(petsclib, draw)

lg = LibPETSc.PetscDrawLGCreate(petsclib, draw, 1)      # one curve
for i in 0:9
    LibPETSc.PetscDrawLGAddPoint(petsclib, lg, [Float64(i)], [sin(i)])
end
LibPETSc.PetscDrawLGDraw(petsclib, lg)
LibPETSc.PetscDrawFlush(petsclib, draw)

LibPETSc.PetscDrawLGDestroy(petsclib, lg)
LibPETSc.PetscDrawDestroy(petsclib, draw)
PETSc.finalize(petsclib)
```

## Function Reference

### PetscDraw

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscDraw_wrappers.jl"]
Order   = [:function]
```
