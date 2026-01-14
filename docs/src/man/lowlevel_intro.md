# Low-Level Interface (LibPETSc)

PETSc.jl provides two ways to interact with the PETSc library:

1. **High-level interface**: Julia-friendly wrappers that handle memory management, use 1-based indexing, and provide convenient syntax (e.g., `A[1,2] = 3.0`)
2. **Low-level interface** (`LibPETSc`): Direct access to nearly the entire PETSc C API (~13,000+ functions)

This guide focuses on the low-level interface.

## When to Use the Low-Level Interface

Use the low-level `LibPETSc` interface when:

- You need PETSc functionality not yet wrapped in the high-level interface
- You're porting existing PETSc C/Fortran code to Julia
- You need fine-grained control over PETSc objects
- You want to access advanced or experimental PETSc features
- Performance-critical code requires eliminating Julia wrapper overhead

For most common tasks (creating vectors/matrices, solving linear systems), the high-level interface is recommended.

## Basic Usage Pattern

The low-level interface follows the PETSc C API closely. All functions require a `petsclib` parameter as the first argument:

```julia
using PETSc

# Get a library instance (usually you'll use the first one)
petsclib = PETSc.petsclibs[1]

# Initialize PETSc
PETSc.initialize(petsclib)

# Create and use PETSc objects
vec = LibPETSc.VecCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.VecSetSizes(petsclib, vec, 10, 10)
LibPETSc.VecSetType(petsclib, vec, "seq")  # Set vector type
LibPETSc.VecSetFromOptions(petsclib, vec)

# Work with the vector
LibPETSc.VecSet(petsclib, vec, 1.0)
println("Vector size: ", LibPETSc.VecGetSize(petsclib, vec))

# Clean up
LibPETSc.VecDestroy(petsclib, vec)
PETSc.finalize(petsclib)
```

## Key Concepts

### 1. The `petsclib` Parameter

Every low-level function requires a `petsclib` instance. This specifies which PETSc library variant to use (e.g., Float64 vs Float32, Int32 vs Int64):

```julia
# Default libraries available (uses prebuilt PETSc binaries)
petsclib = PETSc.petsclibs[1]  # Usually Float64, Int64

# Access the scalar and integer types
PetscScalar = petsclib.PetscScalar  # The floating-point type (e.g., Float64)
PetscInt = petsclib.PetscInt        # The integer type (e.g., Int64)
PetscReal = petsclib.PetscReal      # The real type (real part of PetscScalar)
```

#### Using a Custom PETSc Build

You can link to your own custom build of PETSc instead of using the prebuilt binaries. This is useful when you need:
- Specific configure options not available in the default build
- Custom optimizations for your hardware
- Debug builds for troubleshooting
- Integration with specific external packages

To use a custom PETSc installation:

```julia
using PETSc

# Create a custom library instance pointing to your PETSc installation
petsclib = PETSc.SetPetscLib("/path/to/your/libpetsc.so"; 
                             PetscScalar=Float64, 
                             PetscInt=Int64)

# Initialize and use as normal
PETSc.initialize(petsclib)

# ... your code using petsclib ...

PETSc.finalize(petsclib)
```

**Important notes for custom builds:**
- The dynamic library path should point to `libpetsc.so` (Linux), `libpetsc.dylib` (macOS), or `libpetsc.dll` (Windows)
- The `PetscScalar` and `PetscInt` types must match how your PETSc was configured
- Your custom PETSc must be compatible with the MPI version used by `MPI.jl`
- You can check available precompiled libraries with `[PETSc.petsclibs...]`

### 2. Zero-Based Indexing

**Important**: The low-level interface uses 0-based indexing (C convention), while Julia uses 1-based indexing:

```julia
# Low-level (0-based)
indices = PetscInt[0, 1, 2]  # First three elements
LibPETSc.VecSetValues(petsclib, vec, 3, indices, values, INSERT_VALUES)

# High-level (1-based)
vec[1] = value  # First element
```

### 3. Error Handling

Low-level functions return `PetscErrorCode`. Use the `@chk` macro to check for errors:

```julia
using PETSc.LibPETSc: @chk

err = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
@chk err  # Throws an error if PETSc returned non-zero
```

### 4. String Convenience Wrappers

Many PETSc `SetType` functions accept C string pointers. For convenience, PETSc.jl provides Julia `String` overloads:

```julia
# String convenience wrapper (recommended)
LibPETSc.MatSetType(petsclib, mat, "seqaij")
LibPETSc.VecSetType(petsclib, vec, "seq")
LibPETSc.KSPSetType(petsclib, ksp, "gmres")
LibPETSc.SNESSetType(petsclib, snes, "newtonls")
LibPETSc.PCSetType(petsclib, pc[], "ilu")
LibPETSc.TSSetType(petsclib, ts, "bdf")
LibPETSc.TaoSetType(petsclib, tao, "lmvm")
LibPETSc.DMSetType(petsclib, dm, "da")
LibPETSc.PetscViewerSetType(petsclib, viewer, "ascii")

# Equivalent low-level C pointer syntax (not recommended unless necessary)
ptr = Base.unsafe_convert(Ptr{Int8}, pointer(Vector{UInt8}("seqaij\0")))
LibPETSc.MatSetType(petsclib, mat, ptr)
```

The string wrappers handle the C string conversion internally, making the code cleaner and more Julia-friendly.

Most wrapper functions already include error checking, but when calling C functions directly, use `@chk`.

### 4. Memory Management

PETSc objects created with `Create` functions must be destroyed with corresponding `Destroy` functions:

```julia
# Create
mat = LibPETSc.MatCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.MatSetSizes(petsclib, mat, m_local, n_local, m_global, n_global)

# ... use mat ...

# Destroy when done
LibPETSc.MatDestroy(petsclib, mat)
```

For serial (single-process) objects, the high-level interface handles this automatically via finalizers.

### 5. Assembly

After setting values in vectors or matrices, you must call assembly functions:

```julia
# Set values
LibPETSc.VecSetValues(petsclib, vec, n, indices, values, INSERT_VALUES)

# Assemble
LibPETSc.VecAssemblyBegin(petsclib, vec)
LibPETSc.VecAssemblyEnd(petsclib, vec)
```

## Common Patterns

### Creating and Filling a Vector

```julia
petsclib = PETSc.petsclibs[1]
PetscInt = petsclib.PetscInt
PetscScalar = petsclib.PetscScalar

# Create vector
vec = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
LibPETSc.VecSetSizes(petsclib, vec, 10, 10)
LibPETSc.VecSetType(petsclib, vec, "seq")  # Set vector type
LibPETSc.VecSetFromOptions(petsclib, vec)

# Set values (0-based indices!)
indices = PetscInt[0, 1, 2, 3, 4]
values = PetscScalar[1.0, 2.0, 3.0, 4.0, 5.0]
LibPETSc.VecSetValues(petsclib, vec, 5, indices, values, INSERT_VALUES)

# Assemble
LibPETSc.VecAssemblyBegin(petsclib, vec)
LibPETSc.VecAssemblyEnd(petsclib, vec)

# View (print to stdout)
LibPETSc.VecView(petsclib, vec, C_NULL)

# Clean up
LibPETSc.VecDestroy(petsclib, vec)
```

### Creating a Sparse Matrix

```julia
# Create matrix
mat = LibPETSc.MatCreate(petsclib, MPI.COMM_SELF)
LibPETSc.MatSetSizes(petsclib, mat, 5, 5, 5, 5)
LibPETSc.MatSetType(petsclib, mat, "seqaij")  # String convenience wrapper
LibPETSc.MatSetUp(petsclib, mat)

# Set values (0-based indexing!)
row = PetscInt[0]
cols = PetscInt[0, 1]
vals = PetscScalar[2.0, -1.0]
LibPETSc.MatSetValues(petsclib, mat, 1, row, 2, cols, vals, INSERT_VALUES)

# Assemble
LibPETSc.MatAssemblyBegin(petsclib, mat, LibPETSc.MAT_FINAL_ASSEMBLY)
LibPETSc.MatAssemblyEnd(petsclib, mat, LibPETSc.MAT_FINAL_ASSEMBLY)

# View
LibPETSc.MatView(petsclib, mat, C_NULL)

# Clean up
LibPETSc.MatDestroy(petsclib, mat)
```

### Solving a Linear System

```julia
# Create KSP solver
ksp = LibPETSc.KSPCreate(petsclib, MPI.COMM_SELF)
LibPETSc.KSPSetOperators(petsclib, ksp, mat, mat)
LibPETSc.KSPSetFromOptions(petsclib, ksp)

# Solve Ax = b
LibPETSc.KSPSolve(petsclib, ksp, b, x)

# Get convergence info
reason = Ref{PetscInt}()
LibPETSc.KSPGetConvergedReason(petsclib, ksp, reason)

iterations = Ref{PetscInt}()
LibPETSc.KSPGetIterationNumber(petsclib, ksp, iterations)

println("Converged in $(iterations[]) iterations, reason: $(reason[])")

# Clean up
LibPETSc.KSPDestroy(petsclib, ksp)
```

## Mixing High-Level and Low-Level Interfaces

You can mix both interfaces. High-level objects provide `.ptr` field to access the underlying C pointer:

```julia
# Create with high-level interface
vec_high = VecSeq(petsclib, 10)

# Use with low-level interface
LibPETSc.VecSet(petsclib, vec_high.ptr, 5.0)

# Or use the object directly (if it's a compatible type)
LibPETSc.VecView(petsclib, vec_high, C_NULL)
```

## Finding Functions

The low-level interface provides wrappers for most PETSc functions. To find a function:

1. **Check the PETSc manual**: https://petsc.org/release/docs/
2. **Use Julia's help system**: Type `?LibPETSc.FunctionName`
3. **Browse the documentation**: See the reference pages for each PETSc class (Vec, Mat, KSP, etc.)
4. **Autocomplete**: In the REPL, type `LibPETSc.Vec` and press Tab to see all Vec functions

## Type Conventions

Low-level functions use these type patterns:

```julia
# PETSc objects (opaque pointers)
CVec = Ptr{Cvoid}        # Vector
CMat = Ptr{Cvoid}        # Matrix
CDM = Ptr{Cvoid}         # DM (domain management)
CKSP = Ptr{Cvoid}        # KSP (linear solver)
CSNES = Ptr{Cvoid}       # SNES (nonlinear solver)

# PETSc data types (depend on petsclib)
PetscInt                 # Integer type (Int32 or Int64)
PetscScalar              # Scalar type (Float64, Float32, ComplexF64, etc.)
PetscReal                # Real type (real part of PetscScalar)
```

## Reference Pages

Detailed documentation for low-level functions by category:

- [Vec (Vectors)](@ref vec_lowlevel.md) - ~293 functions for vector operations
- [Mat (Matrices)](@ref mat_lowlevel.md) - ~756 functions for matrix operations
- [KSP (Linear Solvers)](@ref ksp_lowlevel.md) - ~256 functions for iterative linear solvers
- [SNES (Nonlinear Solvers)](@ref snes_lowlevel.md) - ~333 functions for nonlinear solvers

## Getting Help

If you encounter issues:

1. Check the [PETSc documentation](https://petsc.org/release/docs/)
2. Review the [examples](https://github.com/JuliaParallel/PETSc.jl/tree/main/examples)
3. Ask questions on [Julia Discourse](https://discourse.julialang.org/) with the `petsc` tag
4. Open an issue on [GitHub](https://github.com/JuliaParallel/PETSc.jl/issues)
