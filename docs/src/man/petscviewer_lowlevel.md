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
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Create a viewer for ASCII output to stdout
viewer = LibPETSc.PetscViewerCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.PetscViewerSetType(petsclib, viewer, "ascii")  # String convenience wrapper
LibPETSc.PetscViewerFileSetMode(petsclib, viewer, LibPETSc.FILE_MODE_WRITE)

# View a vector
# LibPETSc.VecView(petsclib, vec, viewer)

# View a matrix  
# LibPETSc.MatView(petsclib, mat, viewer)

# Cleanup - wrap in Ref since PetscViewerDestroy expects Ptr{PetscViewer}
viewer_ref = Ref(viewer)
LibPETSc.PetscViewerDestroy(petsclib, viewer_ref)
```

## Convenience Functions

For commonly used viewers, PETSc.jl provides convenience functions:

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Get stdout viewer (single process)
viewer_stdout_self = LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)

# Get stdout viewer (all processes)
viewer_stdout_world = LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)

# Get stderr viewer (single process)
viewer_stderr_self = LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)

# Get stderr viewer (all processes)
viewer_stderr_world = LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)

# Use them to view objects
# LibPETSc.VecView(petsclib, vec, viewer_stdout_self)
# LibPETSc.MatView(petsclib, mat, viewer_stderr_world)
```

## Output to Files

### ASCII File Output

```julia
# Create ASCII file viewer
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerASCIIOpen(petsclib, LibPETSc.PETSC_COMM_SELF, "output.txt", viewer)

# Set format (optional)
LibPETSc.PetscViewerPushFormat(petsclib, viewer[], LibPETSc.PETSC_VIEWER_ASCII_MATLAB)

# View object
# LibPETSc.MatView(petsclib, mat, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

### Binary File Output

```julia
# Create binary viewer for checkpointing
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerBinaryOpen(petsclib, MPI.COMM_WORLD, "checkpoint.dat", 
                               LibPETSc.FILE_MODE_WRITE, viewer)

# Save vector
# LibPETSc.VecView(petsclib, vec, viewer[])

# Save matrix
# LibPETSc.MatView(petsclib, mat, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

### Loading from Binary Files

```julia
# Open for reading
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerBinaryOpen(petsclib, MPI.COMM_WORLD, "checkpoint.dat",
                               LibPETSc.FILE_MODE_READ, viewer)

# Load vector
vec = LibPETSc.VecCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.VecLoad(petsclib, vec, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

## Visualization Formats

### VTK Output

```julia
# Create VTK viewer for ParaView/VisIt
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerVTKOpen(petsclib, MPI.COMM_WORLD, "solution.vtu",
                            LibPETSc.FILE_MODE_WRITE, viewer)

# View DM-based solution
# LibPETSc.DMView(petsclib, dm, viewer[])
# LibPETSc.VecView(petsclib, solution, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

### HDF5 Output

```julia
# Create HDF5 viewer for hierarchical data
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerHDF5Open(petsclib, MPI.COMM_WORLD, "data.h5",
                             LibPETSc.FILE_MODE_WRITE, viewer)

# Organize data in groups
LibPETSc.PetscViewerHDF5PushGroup(petsclib, viewer[], "/timestep_001")
# LibPETSc.VecView(petsclib, vec, viewer[])
LibPETSc.PetscViewerHDF5PopGroup(petsclib, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
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
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerDrawOpen(petsclib, LibPETSc.PETSC_COMM_SELF, C_NULL, "Plot", 
                              0, 0, 600, 600, viewer)

# View vector as bar chart
# LibPETSc.VecView(petsclib, vec, viewer[])

# View matrix structure
# LibPETSc.MatView(petsclib, mat, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

## Socket Viewer (MATLAB/Python)

Stream data to external tools:

```julia
# Create socket viewer
viewer = Ref{LibPETSc.PetscViewer}()
LibPETSc.PetscViewerSocketOpen(petsclib, MPI.COMM_WORLD, "localhost", 5000, viewer)

# Send data
# LibPETSc.VecView(petsclib, vec, viewer[])

LibPETSc.PetscViewerDestroy(petsclib, viewer)
```

## Monitoring Convergence

Viewers are used with KSP/SNES monitors:

```julia
# Monitor KSP residuals (automatic viewer to stdout)
# LibPETSc.KSPMonitorSet(petsclib, ksp, LibPETSc.KSPMonitorDefault, 
#                        LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib), C_NULL)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscViewer_wrappers.jl"]
Order   = [:function]
```
