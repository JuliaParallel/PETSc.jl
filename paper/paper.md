---
title: 'PETSc.jl: A Julia interface for the Portable, Extensible Toolkit for Scientific Computation'
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
    orcid: 0009-0006-8207-2498
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

Many problems in science and engineering can be described by partial differential equations that have to be solved for millions to billions of unknowns on parallel computers. PETSc, the Portable, Extensible Toolkit for Scientific Computation [@petsc-user-ref; @petsc-efficient], is one of the most widely used libraries for this. It provides distributed vectors and sparse matrices, linear and nonlinear solvers, preconditioners, time integrators and mesh management, built on MPI and with growing GPU support. PETSc is written in C and is used from C, C++, Fortran, or Python through petsc4py [@petsc4py].

PETSc.jl makes PETSc available from the Julia language [@julia]. It builds on [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl), pre-built PETSc binaries that we maintain for all major operating systems, so that `pkg> add PETSc` gives a working parallel solver stack without compiling anything. Users write the problem-specific parts of a simulation, such as residuals, Jacobians and boundary conditions, as ordinary Julia functions, and PETSc does the parallel linear algebra, the solves and the mesh handling. PETSc vectors and matrices behave like Julia arrays, and Julia's automatic differentiation [@forwarddiff] and GPU packages [@cudajl; @kernelabstractions] can be used inside user code. The same program runs on a laptop and, using `MPI`, on thousands of cores.

# Statement of need

Julia has a rich ecosystem for differential equations [@diffeq; @nonlinearsolve], but its native solvers were not designed for distributed-memory parallelism at the scale of thousands of MPI ranks, and reimplementing the preconditioners, multigrid hierarchies and mesh infrastructure that PETSc has accumulated over three decades is challenging. At the same time, PETSc users who write their physics in C or Fortran lose the productivity of a high-level language: residual functions are verbose, automatic differentiation of Jacobians requires external tools, and prototyping a new discretisation is slow. Likewise, rewriting compute kernels for GPU's is challenging particularly given the wide range of programming languages for that.

PETSc.jl lets a researcher write the physics in Julia and delegate the parallel linear algebra, the Newton and Krylov iterations, the time stepping and the parallel mesh distribution to PETSc. We initially focussed on geodynamic applications, where codes such as LaMEM [@lamem] use PETSc directly to solve (nonlinear) thermo-mechanical problems on staggered and unstructured grids on hundreds to tens of thousands of cores, but nothing in PETSc.jl is specific to that field. Examples include Julia translations of PETSc's own tutorials, among them the DMStag, DMDA structured (staggered) grid finite difference examples, or DMPlex finite-element examples with multigrid preconditioners.

# State of the field

petsc4py [@petsc4py] is the established high-level Python interface to PETSc. Compared with it, PETSc.jl comes with precompiled binaries, and its callbacks compile to native code, so the residual evaluations that PETSc calls at every grid node or quadrature point carry no interpreter overhead and run natively on GPUs.

Within Julia, GridapPETSc, the PETSc backend of the Gridap finite-element framework [@gridap; @gridapdistributed], as well as PetscWrap and PetscCall, wrap the parts of PETSc they need, and all of them obtain the library from PETSc_jll. Earlier versions of PETSc.jl likewise wrapped a few hundred functions by hand. Hand-written wrappers fall behind each PETSc release and leave most of the API out of reach. Rather than start another package, we rewrote PETSc.jl around bindings generated for the complete API, and we maintain PETSc_jll, so that improvements to the binaries benefit every Julia package that uses PETSc.

# Software design

**Low and high-level layers.** The low-level layer, `LibPETSc`, is automatically generated for each PETSc release by a modified version of the Python tooling PETSc uses for its own Fortran bindings, and covers about 6300 functions of the PETSc 3.25 C API. The generated wrappers check error codes and throw a `PetscError` carrying PETSc's message, return output arguments as tuples instead of taking pointers, and reject a destroyed or null handle instead of crashing. A hand-written, more Julia-like, high-level layer covers the most used objects, `Vec`, `Mat`, `KSP`, `SNES`, `TS`, `DMDA`, `DMStag` and `DMPlex`, with Julia semantics: a PETSc vector is an `AbstractVector`, entries are set with `A[i, j] = v`, and a linear system is solved with `ksp \ b`. Generation keeps the whole API reachable and turns a new PETSc release into a script run; hand-written code is reserved for the part of the API where idiomatic Julia matters most.

**Pre-built binaries.** Configuring and compiling PETSc together with its external packages is often the hardest step in using it. We therefore maintain [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl), which is cross-compiled with Julia's BinaryBuilder [@binarybuilder] for Linux, macOS and FreeBSD on x86-64 and ARM processors and for Windows on x86-64, together with a CUDA-enabled build for Linux. It comes in eight variants, real or complex scalars, single or double precision, and 32- or 64-bit indices, which are linked against various combinations of the sparse direct solvers MUMPS, SuperLU_DIST and SuiteSparse, the hypre algebraic multigrid preconditioners, HDF5 input and output, and the mesh generators Triangle and TetGen. The binaries are not tied to Julia: each release is published as one tarball per platform with PETSc's headers and libraries, and the libraries it links against are distributed in the same way, so C, C++ or Fortran codes can use them without building PETSc.

**Several libraries in one process.** PETSc is compiled for one scalar type and one integer width. PETSc.jl loads four variants side by side, real and complex in single and double precision, with the index width chosen by a preference, and carries the library in the type of every object. Dispatch thus selects the right C function at compile time, and a `Float32` and a `ComplexF64` solver can coexist.

**Object lifetime.** PETSc objects shared across MPI ranks must be destroyed collectively, in the same order on every rank, which Julia's garbage collector cannot guarantee. PETSc.jl therefore frees single-rank objects from a finalizer, deferred to the next PETSc.jl call that creates or destroys an object, or to `finalize`, so that no PETSc function runs on a garbage-collector thread, and requires an explicit `destroy!` for objects on multiple ranks.

**Callbacks and options.** Residual, Jacobian and monitor functions are Julia closures. For each callback type, PETSc.jl registers one fixed C function with PETSc and passes the closure as PETSc's user context; when PETSc calls the C function, it retrieves the closure and calls it. A user function receives PETSc vectors as Julia arrays, can capture any state, and may throw errors. In this way, it is possible to use Julia automatic differentiation tools within a callback. Every PETSc option can also be given as a keyword, `KSP(A; ksp_type = "cg", pc_type = "gamg")`, through the `PETSC_OPTIONS` environment variable, or on the command line when the script passes its arguments on, as in `PETSc.initialize(petsclib; options = ARGS)`, so solvers can be reconfigured from a batch script as in C.

**GPUs.** `with_local_array!` passes a vector's local storage as a Julia array to a `do` block and hands it back to PETSc afterwards. With the CUDA.jl extension, the array is a `CuArray` when the vector lives on the device; the AMDGPU.jl extension does the same with a `ROCArray`, for a locally built PETSc with HIP support, as PETSc_jll has no HIP build yet. In both cases, a [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl) kernel [@kernelabstractions] can fill a PETSc vector in place, and PETSc's GPU solvers are used unchanged. This thus forms a simple way to write GPU kernels that are highly portable.

**HPC.** The PETSc_jll binaries can be deployed directly on a supercomputer, without compiling PETSc there. PETSc_jll is built for the standard MPI application binary interface introduced in MPI 5.0 [@mpi50], as well as for MPICH, Open MPI and MPItrampoline [@mpitrampoline], and loads the build matching the MPI that MPI.jl [@mpijl] is configured for. Pointing MPI.jl at the vendor MPI of a cluster, through the standard ABI or through MPItrampoline, is therefore enough to run the same binaries there. A locally built PETSc is selected with one call and stored as a preference.

# Research impact statement

\autoref{fig:scaling} shows a weak-scaling study of PETSc's `ex45` tutorial, a 3D Laplacian solved with conjugate gradients and geometric multigrid, on the LUMI-C supercomputer, from 64 to 32768 MPI ranks with $4097^3$ grid points at the largest size. The Julia version with the pre-built binaries, the Julia version linked against a locally built PETSc, and the original C program give identical iteration counts and errors, and the time spent inside `KSPSolve` is within a few seconds of C up to 4096 ranks. The remaining difference in total time is Julia-side setup and just-in-time compilation, paid once per run. A test checks that the cost the Julia layer adds to each call does not grow with problem size.

![Weak scaling of the 3D Laplacian example (`ex45`) on LUMI-C, comparing PETSc.jl with pre-built binaries, PETSc.jl with a locally built PETSc, and the native C program. [NOTE: To be replaced with newer results on Alps, as this was done with an earlier version of the package]\label{fig:scaling}](weak_scaling_ex45_LUMI.png)

PETSc.jl is used as a solver backend elsewhere in the Julia ecosystem: LinearSolve.jl and NonlinearSolve.jl [@nonlinearsolve] offer PETSc solvers through package extensions, and SafePETSc builds on it. PETSc_jll is a dependency of the geodynamics code LaMEM [@lamem], of GridapPETSc and GridapDistributed's PETSc wrappers, and of PetscWrap and PetscCall, which therefore all run on the binaries described here.

# AI usage disclosure
A large part of this package, including the automatically generated wrappers, was written before we started using AI tools. Generative AI (Claude by Anthropic, through Claude Code) has been used since January 2026 in developing PETSc.jl and PETSc_jll build recipes, by refactoring the automatic wrappers, for code, tests and documentation, and to draft this paper. All AI-assisted changes were reviewed by the authors and merged through pull requests that run the test suite, including the examples, on Linux, macOS and Windows, and against a PETSc built from source. The authors checked and edited the text of this paper and take responsibility for its content.

# Acknowledgements

PETSc.jl was started by Jared Crean and Viral Shah and rewritten by Jeremy Kozdon, Boris Kaus, Valentin Churavy and Simon Byrne with most recent contributions by Boris Kaus and Filippo Ferrazzini. Funding came from the European Research Council Consolidator Grant MAGMA (grant agreement 771143), the GPU4GEO and dGPU4GEO projects of the Platform for Advanced Scientific Computing (PASC), a Gutenberg Research Council Fellowship of the Johannes Gutenberg University Mainz, the DFG project MOSAIC, and the project DEGREE of the German Federal Ministry of Research, Technology and Space (BMFTR). The scaling runs used the LUMI supercomputer, owned by the EuroHPC Joint Undertaking and hosted by CSC (Finland).

# References
