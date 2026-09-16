# PetscViewer - Low-level Interface

The PetscViewer component provides flexible I/O capabilities for visualizing and saving PETSc objects, including vectors, matrices, and other data structures. Viewers support multiple output formats for analysis, debugging, and post-processing.

## Overview

PETSc viewers enable:
- **Text output**: ASCII formatted data for debugging
- **Binary I/O**: Efficient storage and checkpointing
- **Visualization**: Integration with visualization tools (VTK, HDF5, MATLAB)
- **Monitoring**: Runtime inspection of solver progress
- **Logging**: Recording solver statistics and performance data

Available viewer types:
- **PETSCVIEWERASCII**: Human-readable text output
- **PETSCVIEWERBINARY**: Platform-independent binary format
- **PETSCVIEWERVTK**: VTK format for ParaView, VisIt
- **PETSCVIEWERHDF5**: HDF5 hierarchical data format
- **PETSCVIEWERDRAW**: X-window graphics (2D plots, contours)
- **PETSCVIEWERSOCKET**: Network streaming to MATLAB, Python
- **PETSCVIEWERMATLAB**: MATLAB-compatible output

## Basic Usage

```julia
using PETSc, MPI

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a viewer for ASCII output to stdout
viewer = LibPETSc.PetscViewerCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.PetscViewerSetType(petsclib, viewer, LibPETSc.PETSCVIEWERASCII)
LibPETSc.PetscViewerFileSetMode(petsclib, viewer, LibPETSc.FILE_MODE_WRITE)

# View a vector
# LibPETSc.VecView(petsclib, vec, viewer)

# View a matrix  
# LibPETSc.MatView(petsclib, mat, viewer)

# Cleanup
LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Convenience Functions

For commonly used viewers, PETSc.jl provides convenience functions:

```julia
using PETSc, MPI

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Get stdout viewer (single process)
viewer_stdout_self = LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)

# Get stdout viewer (all processes)
viewer_stdout_world = LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)

# Get stderr viewer (single process)
viewer_stderr_self = LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)

# Get stderr viewer (all processes)
viewer_stderr_world = LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()

# Use them to view objects
# LibPETSc.VecView(petsclib, vec, viewer_stdout_self)
# LibPETSc.MatView(petsclib, mat, viewer_stderr_world)
```

## Output to Files

### ASCII File Output

```julia
# Create ASCII file viewer
viewer = LibPETSc.PetscViewerASCIIOpen(petsclib, LibPETSc.PETSC_COMM_SELF, "output.txt")

# Set format (optional)
LibPETSc.PetscViewerPushFormat(petsclib, viewer, LibPETSc.PETSC_VIEWER_ASCII_MATLAB)

# View object
# LibPETSc.MatView(petsclib, mat, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

### Binary File Output

```julia
# Create binary viewer for checkpointing
viewer = LibPETSc.PetscViewerBinaryOpen(petsclib, MPI.COMM_WORLD, "checkpoint.dat",
                                        LibPETSc.FILE_MODE_WRITE)

# Save vector
# LibPETSc.VecView(petsclib, vec, viewer)

# Save matrix
# LibPETSc.MatView(petsclib, mat, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

### Loading from Binary Files

```julia
# Open for reading
viewer = LibPETSc.PetscViewerBinaryOpen(petsclib, MPI.COMM_WORLD, "checkpoint.dat",
                                        LibPETSc.FILE_MODE_READ)

# Load vector
vec = LibPETSc.VecCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.VecLoad(petsclib, vec, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Visualization Formats

### VTK Output

```julia
# Create VTK viewer for ParaView/VisIt
viewer = LibPETSc.PetscViewerVTKOpen(petsclib, MPI.COMM_WORLD, "solution.vtu",
                                     LibPETSc.FILE_MODE_WRITE)

# View DM-based solution
# LibPETSc.DMView(petsclib, dm, viewer)
# LibPETSc.VecView(petsclib, solution, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

### HDF5 Output

```julia
# Create HDF5 viewer for hierarchical data
viewer = LibPETSc.PetscViewerHDF5Open(petsclib, MPI.COMM_WORLD, "data.h5",
                                      LibPETSc.FILE_MODE_WRITE)

# Organize data in groups
LibPETSc.PetscViewerHDF5PushGroup(petsclib, viewer, "/timestep_001")
# LibPETSc.VecView(petsclib, vec, viewer)
LibPETSc.PetscViewerHDF5PopGroup(petsclib, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Standard Viewers

PETSc provides predefined viewers:

```julia
# Standard output
LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)

# Standard error
LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)
LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)

# Example: view to stdout
# LibPETSc.VecView(petsclib, vec, LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib))
```

## Format Options

Control output detail with `PetscViewerPushFormat`:

- **PETSC_VIEWER_DEFAULT**: Standard format
- **PETSC_VIEWER_ASCII_MATLAB**: MATLAB-compatible format
- **PETSC_VIEWER_ASCII_DENSE**: Dense matrix format
- **PETSC_VIEWER_ASCII_INFO**: Summary information only
- **PETSC_VIEWER_ASCII_INFO_DETAIL**: Detailed information

## Draw Viewer (Graphics)

For interactive 2D visualization:

```julia
# Create draw viewer (X-window)
# display "" selects the default display
viewer = LibPETSc.PetscViewerDrawOpen(petsclib, LibPETSc.PETSC_COMM_SELF, "", "Plot",
                                      Cint(0), Cint(0), Cint(600), Cint(600))

# View vector as bar chart
# LibPETSc.VecView(petsclib, vec, viewer)

# View matrix structure
# LibPETSc.MatView(petsclib, mat, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Socket Viewer (MATLAB/Python)

Stream data to external tools:

```julia
# Create socket viewer
viewer = LibPETSc.PetscViewerSocketOpen(petsclib, MPI.COMM_WORLD, "localhost", Cint(5000))

# Send data
# LibPETSc.VecView(petsclib, vec, viewer)

LibPETSc.PetscViewerDestroy(petsclib, viewer)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Monitoring Convergence

Viewers are used with KSP/SNES monitors:

```julia
# Monitor KSP residuals with a Julia callback (see the KSP page for the callback signature)
# LibPETSc.KSPMonitorSet(petsclib, ksp, monitor_cfunction, C_NULL, C_NULL)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscViewer_wrappers.jl"]
Order   = [:function]
```

### Standard viewers (hand-written)

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/extra_wrappers.jl"]
Order   = [:function]
Filter  = t -> startswith(string(nameof(t)), "PETSC_VIEWER")
```
