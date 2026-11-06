---
title: 'PETSc.jl: Julia bindings for the Portable, Extensible Toolkit for Scientific Computation'
tags:
  - Julia
  - PETSc
  - partial differential equations
  - parallel computing
  - MPI
  - GPU
  - nonlinear solvers
authors:
  - name: Boris J. P. Kaus
    orcid: 0000-0002-0247-8660
    corresponding: true
    affiliation: 1
  - name: Filippo Ferrazzini
    affiliation: 1
  - name: Viral B. Shah
    orcid: 0000-0001-9602-4012
    affiliation: 2
  - name: Jeremy E. Kozdon
    orcid: 0000-0002-2493-4292
    affiliation: 3
  - name: Valentin Churavy
    orcid: 0000-0002-9033-165X
    affiliation: 1
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
  - name: Erik Schnetter
    orcid: 0000-0002-4518-9017
    affiliation: 4
  - name: Simon Byrne
    orcid: 0000-0001-8048-6810
    affiliation: 5
affiliations:
  - name: Johannes Gutenberg University Mainz, Germany
    index: 1
  - name: Julia Computing, USA
    index: 2
  - name: NextSilicon, USA
    index: 3
  - name: Perimeter Institute, Canada
    index: 4
  - name: NVIDIA, USA
    index: 5
date: 6 November 2026
bibliography: paper.bib
---

# Summary

PETSc, the Portable, Extensible Toolkit for Scientific Computation [@petsc-user-ref; @petsc-efficient], is one of the most widely used libraries for the parallel solution of partial differential equations. It provides distributed vectors and sparse matrices, Krylov solvers and preconditioners, Newton-type nonlinear solvers, time integrators, optimization, and structured and unstructured mesh management, all on top of MPI and with growing GPU support. PETSc is written in C, and application codes have traditionally been written in C, C++ or Fortran, or in Python through petsc4py [@petsc4py].

PETSc.jl makes PETSc available from Julia [@julia]. It ships pre-built PETSc binaries for Linux, macOS and Windows, so that `pkg> add PETSc` gives a working parallel solver stack without compiling anything, and it can equally be pointed at a PETSc built on a supercomputer. Two layers are provided. A low-level layer, `LibPETSc`, is generated from PETSc's own sources and covers about 6300 functions of the PETSc 3.25 C API with Julia-side type information and docstrings. A high-level layer exposes the most used objects, `Vec`, `Mat`, `KSP`, `SNES`, `TS`, `DMDA`, `DMStag` and `DMPlex`, with Julia semantics: a PETSc vector is an `AbstractVector`, a PETSc matrix an `AbstractMatrix`, entries are set with `A[i, j] = v`, a linear system is solved with `ksp \ b`, and residual and Jacobian routines are ordinary Julia functions or closures. Julia's automatic differentiation [@forwarddiff] and GPU packages [@cudajl; @kernelabstractions] can be used inside those callbacks.

# Statement of need

Julia has a rich ecosystem for differential equations [@diffeq; @nonlinearsolve], but its native solvers were not designed for distributed-memory parallelism at the scale of thousands of MPI ranks, and reimplementing the preconditioners, multigrid hierarchies and mesh infrastructure that PETSc has accumulated over three decades is not realistic. At the same time, PETSc users who write their physics in C or Fortran lose the productivity of a high-level language: residual functions are verbose, automatic differentiation of Jacobians requires external tools, and prototyping a new discretisation is slow.

PETSc.jl bridges the two. A researcher writes the problem-specific code, residuals, Jacobians, boundary conditions and material models, in Julia, and delegates the parallel linear algebra, the Newton and Krylov iterations, the time stepping and the mesh distribution to PETSc. The package was developed for geodynamic modelling, where Stokes and thermo-mechanical problems on staggered and unstructured grids are solved on hundreds to tens of thousands of cores, but nothing in it is specific to that field. The examples directory contains Julia translations of PETSc's own tutorial examples, including the DMPlex finite-element Stokes examples (`ex62`, `ex69`), a staggered-grid Stokes solver, and time-dependent problems.

Compared with petsc4py, PETSc.jl gives the same reach into the C library without a separate build step, and its callbacks compile to native code, so the residual evaluations that PETSc calls at every quadrature point or grid node carry no interpreter overhead. Compared with wrapping PETSc through Julia's `ccall` by hand, which the previous generation of PETSc.jl did for a few hundred functions, the generated layer keeps the whole API available and tracks new PETSc releases with a regeneration script.

# Design

**Two layers.** `LibPETSc` is generated once per PETSc release by a modified version of the Python tooling PETSc uses to produce its own Fortran bindings. Every generated wrapper takes the PETSc library as its first argument, checks the returned error code and throws a `PetscError` with the C function name and PETSc's message, returns output arguments as a tuple or `NamedTuple` instead of taking `Ref`s, and rejects a destroyed or null handle with an `ArgumentError` instead of a segmentation fault. The high-level layer is hand-written and follows a documented naming convention that maps PETSc's `VecAXPY` style onto Julia's `axpy!` style.

**Several libraries in one process.** PETSc is compiled for one scalar type and one integer width. PETSc.jl ships four builds, real and complex in single and double precision, loads them side by side, and carries the library in the type of every object, so dispatch selects the right C function at compile time and a `Float32` and a `ComplexF64` solver can coexist.

**Object lifetime.** PETSc objects are reference-counted C handles and have to be destroyed on every MPI rank in the same order. PETSc.jl distinguishes owned from borrowed handles, destroys objects created on a single rank from a finalizer that is drained on the next PETSc call (so no PETSc function runs on a garbage-collector thread), and requires an explicit `destroy!` on several ranks, where collective destruction cannot be left to the garbage collector. Using a handle after `destroy!` throws.

**Callbacks.** Residual, Jacobian, monitor and convergence-test functions are Julia closures. The package installs one C-callable trampoline per callback type and passes the closure as the PETSc context, so a user function receives PETSc vectors as Julia arrays, can capture any state, and may throw: an exception raised inside a callback is caught at the C boundary and rethrown by the Julia call that started the solve. The pointwise functions that PETSc's finite-element machinery calls at every quadrature point are written with macros that produce allocation-free closures over PETSc's memory.

**Options.** Every PETSc option can be given as a keyword, `KSP(A; ksp_type = "cg", pc_type = "gamg")`, as a command-line string, or through PETSc's options database, with the usual PETSc precedence, so a solver can be reconfigured from a batch script without changing code.

**Memory and GPUs.** `with_local_array!` checks a vector's local storage out as a Julia array for the duration of a block and returns it to PETSc afterwards. Package extensions for CUDA.jl and AMDGPU.jl make the same call return a `CuArray` or `ROCArray` when the vector lives on the device, so a KernelAbstractions kernel can fill a PETSc vector in place and PETSc's CUDA and HIP vector and matrix types, and the solvers that run on them, are used unchanged. Host code keeps working on device vectors; PETSc moves the data.

**HPC.** The pre-built binaries link MPI through MPItrampoline [@mpitrampoline], so the same binaries run against the vendor MPI of a cluster after MPI.jl [@mpijl] has been pointed at it, and no PETSc build is required on the machine. A locally built PETSc is selected with one call and stored as a preference.

# Example

A 1D Laplacian solved with conjugate gradients and a Jacobi preconditioner:

```julia
using PETSc
petsclib = PETSc.petsclibs[1]            # Float64, Int64
PETSc.initialize(petsclib)

n = 10
A = PETSc.PetscMat(petsclib, n, n, 3)    # at most 3 nonzeros per row
for i in 1:n
    A[i, i] = 2
    i > 1 && (A[i, i - 1] = -1)
    i < n && (A[i, i + 1] = -1)
end
PETSc.assemble!(A)

ksp = PETSc.KSP(A; ksp_type = "cg", pc_type = "jacobi", ksp_rtol = 1e-10)
x = ksp \ ones(n)

PETSc.destroy!(ksp); PETSc.destroy!(A)
PETSc.finalize(petsclib)
```

The same program runs under `mpiexec -n 4` without changes: the matrix is then distributed, and each rank sets the rows it owns.

# Performance

The Julia layer adds a fixed cost per call that does not grow with problem size; this is enforced by a test that checks the allocations of every high-level call at two sizes. \autoref{fig:scaling} shows a weak-scaling study of PETSc's `ex45` tutorial, a 3D Laplacian solved with conjugate gradients and geometric multigrid, on the LUMI-C supercomputer, from 64 to 32768 MPI ranks with a grid of $4097^3$ points at the largest size. The Julia version using the pre-built binaries, the Julia version linked against a locally built PETSc, and the original C program produce identical iteration counts and errors. The time spent inside `KSPSolve` is within a few seconds of the C version up to 4096 ranks; the remaining gap in total time is Julia-side setup, in particular mesh construction, assembly and just-in-time compilation, which is paid once per run.

![Weak scaling of the 3D Laplacian example (`ex45`) on LUMI-C, comparing PETSc.jl with pre-built binaries, PETSc.jl with a locally built PETSc, and the native C program.\label{fig:scaling}](weak_scaling_ex45_LUMI.png)

# Acknowledgements

PETSc.jl was started by Jared Crean and Viral Shah and rewritten by Jeremy Kozdon, Valentin Churavy and Simon Byrne. The current version was developed by Boris Kaus and Filippo Ferrazzini. Funding came from the European Research Council Consolidator Grant MAGMA (grant agreement 771143), the GPU4GEO and dGPU4GEO projects of the Platform for Advanced Scientific Computing (PASC), a Gutenberg Research Council Fellowship of the Johannes Gutenberg University Mainz, the DFG project MOSAIC, and the project DEGREE of the German Federal Ministry of Research, Technology and Space (BMFTR). The scaling runs used the LUMI supercomputer, owned by the EuroHPC Joint Undertaking and hosted by CSC (Finland).

# References
