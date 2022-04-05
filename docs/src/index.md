# PETSc.jl

 [PETSc.jl](https://github.com/JuliaParallel/PETSc.jl) is a Julia wrapper for the Portable, Extensible Toolkit for Scientific Computation [PETSc](https://petsc.org/) package, which allows solving ordinary and partial differential equations in parallel on laptops or massively parallel high-performance systems.

 The use of Julia greatly simplifies the code that developers have to write, while allowing to employ Julia features such as automatic differentiation. The Julia wrapper also comes with a pre-built library, which greatly simplifies the process of getting your first code working in parallel, on different operating systems. In many cases, the Julia code is significantly shorter than its C counterpart.

 This wrapper mimics the PETSc-functionality as closely as possible, but remains work in progress (meaning that not everything has been translated yet). See the official [user guide](https://petsc.org/release/overview/) if you want to learn more about PETSc in general. For Julia-specific examples, have a look at our [examples](https://github.com/JuliaParallel/PETSc.jl/tree/main/examples) or [tests](https://github.com/JuliaParallel/PETSc.jl/tree/main/test).

This package includes a low-level, automatically-generated wrapper layer, upon which a higher-level interface is built.


## The High-Level Interface

The high level interface is designed to be familiar and convenient for Julia users, but exposes only a small portion of the functionality
of the underlying PETSc library.  This interface should be considered *unstable*, as its implementation involves design decisions which are likely to be revisited - note the low version number
of this package if considering relying on these features.

For example, with this interface, PETSc's [KSP](https://petsc.org/release/docs/manual/ksp) linear solvers (including Krylov methods) can be used in a way similar to solvers from other Julia packages. See the example in [Getting started](@ref) and the API in [KSP](@ref).

## The Low-Level Interface

The low-level interface covers more of the PETSc API, but may be awkward to work with and likely requires
previous experience with PETSc to use effectively. It is automatically generated with [Clang.jl](https://github.com/JuliaInterop/Clang.jl).

The high-level interface described in [KSP](@ref) creates a [KSP](https://petsc.org/release/docs/manual/ksp) linear solver object via the low-level interface to [`KSPSolve()`](https://petsc.org/release/docs/manualpages/KSP/KSPCreate.html), with the use of constructs which require knowledge of PETSc's nature as a [C library](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/). Expert users are of course free to directly use the low level interface, as in this simple example which directly calls [`PetscGetVersionNumber()`](https://petsc.org/release/docs/manualpages/PetscGetVersionNumber.html).

```@example
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

petsclib = PETSc.petsclibs[1]
PetscInt = petsclib.PetscInt
PetscErrorCode = Cint

PETSc.initialize(petsclib)
major = Ref{PetscInt}(0)
minor = Ref{PetscInt}(0)
subminor = Ref{PetscInt}(0)
release = Ref{PetscInt}(0)
error_code = ccall(
      (:PetscGetVersionNumber, petsclib.petsc_library),
      PetscErrorCode,
      (Ptr{PetscInt}, Ptr{PetscInt}, Ptr{PetscInt}, Ptr{PetscInt}),
      major, minor, subminor, release
      )

version = (major[], minor[], subminor[])
println("PETSc $(version[1]).$(version[2]).$(version[3])")
PETSc.finalize(petsclib)
```
