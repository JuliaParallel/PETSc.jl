# Design notes



- PETSc can only be built for a single `PetscScalar` type. A workaround is to build multiple PETSc libraries for all supported scalar types (`Float32`, `Float64`, `Complex{Float64}`)
  - Appears that `Complex{Float32}` not supported (https://github.com/JuliaParallel/PETSc.jl/blob/old/deps/build_petscs.jl#L128).
    * TODO: check if still the case
  - Need JLL support for this (https://github.com/JuliaPackaging/Yggdrasil/issues/1527)
  - Define macro to `@eval` all functions which use `ccall` for different scalar types.
  - All PETSc types and methods which involve `ccall` need a `PetscScalar` parameter.
  * TODO: GPU support: need separate ones for CUDA as well?


- Lazily initialize each library if an object of that parameter is constructed.
  - Also initialize MPI if not already initialized
    * TODO: check what thread level is required?
  - Add `atexit` hook to finalize PETSc (this should be okay with MPI `atexit`, due to LIFO)
  - Disable the PETSc signal handler


- A Julia object matching each PETSc object (`Vec`, `Mat`, `KSP`, etc.). 
  - These will typically have a `ptr` as the first field, which lets us use the `cconvert`/`unsafe_convert` trick to pass pointer by value/reference.
  - Most (all?) objects will have a `comm` field, for the MPI communicator
  - Objects which wrap Julia objects will also need a reference to those objects to prevent GC.


- For convenience, attach finalizers to call `destroy` for single-process ("sequential") objects (`VecSeq`, `MatSeqXXX`, or any others where `comm = MPI.COMM_SELF`). 
  - We can't attach finalizers for distributed objects (i.e. `VecMPI`), as `destroy` needs to be called collectively on all MPI ranks.
  - Appears to be safe for users to call `destroy` manually if finalizer already defined 
    * TODO: check this with PETSc devs
  - Unclear how to handle objects that are contained within others, e.g. `PC` from `KSPGetPC`, `KSP` from `SNESGetKSP`, etc.
    * There appears to be some sort of reference counting, unclear if this is valid.
      `PetscObjectReference` / `PetscObjectDereference`

- For PETSc objects which are equivalent to Julia objects (e.g. `VecSeq` : `Vector{PetscScalar}`, `MatSeqDense` : `Matrix{PetscScalar}`), use `XXXCreateSeqWithArray` methods so that they can share same memory.
    * TODO: check PETSc guarantees on accessing the Julia objects directly.
  - For other objects (`MatSeqAIJ`), for now we let PETSc manage memory (may want to re-evaluate this later)
  - Define conversion routines to wrap with `Seq` objects where possible.
  - Define convenience versions of functions which take/return Julia `Vector`s, e.g. `y = KSP(M) \ x` where `y` and `x` are `Vector`s.


- For specifying object options, there are 2 possible approaches:
    (a) use `PetscOptions` objects (key-value pairs) to capture keyword args, which can be pushed and popped to the global options, then use `XXXSetFromOption` methods, e.g. `KSP(mat, ksp_atol=1e-8)`
    (b) use C setter functions (e.g. `KSPSetTolerances`)
  - for now, we go with (a).
    - it's easier
    - not all options are available via C setters, e.g. `mg_coarse_XXX`/`mg_levels_XXX` options
  - ideally we would create a more "object-oriented" interface: e.g. each preconditioner would be a different Julia type, but this doesn't yet seem possible.


- For cases where PETSc needs to call Julia functions (`MatShell`, `SNES`), PETSc provides a mechanism to pass through a context pointer. We can use this to pass through a pointer to the object itself via `pointer_from_objref`.
  * Can we pass `NULL` to vec/matrix args? What does that do?
  - What should the callback interface look like?
  - How to handle errors from within callbacks?


- TODO: Error handling:
  - https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscPushErrorHandler.htmls