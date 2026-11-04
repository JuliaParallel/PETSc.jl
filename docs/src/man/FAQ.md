# Frequently Asked Questions


## 1. Can I use my own PETSc library?
Yes. You do need to compile PETSc as a dynamic (shared) library, and `MPI.jl` must be configured to be compatible with the MPI used to compile PETSc. If you run this on an HPC system and don't know the configuration options, compile one of the PETSc examples and run it with `-log_view`; the output lists all configuration options used on that machine.

Please note that the version of PETSc should be compatible with the version used for the wrappers.

Point PETSc.jl at it with `PETSc.set_library!`; see [Using a custom PETSc build](@ref) for the steps, how to revert, and a one-off session without changing the stored preference.

## 2. Help, my code crashes?
That is very possible. If you provide a *short* minimum working example (MWE), feel free to open an issue on the github repo, so we can check it out.

## 3. What about the garbage collector in Julia?
Using the GC in combination with MPI code is a tricky business. The users are therefore responsible to free PETSc objects in the code, as shown in the various examples. We have some help for that:
1. `PETSc.with_destroy` destroys a group of objects when a block ends, also when it throws, in reverse order, so every rank does it the same way:
   ```julia
   PETSc.with_destroy() do list
       dm  = PETSc.track!(list, PETSc.DMStag(petsclib, comm, ...))
       ksp = PETSc.track!(list, PETSc.KSP(dm))
       ...
   end
   ```
   An object that should outlive the block is taken out with `PETSc.untrack!(list, obj)`.
2. The script `scripts/audit.jl` in the PETSc.jl repository: `julia scripts/audit.jl path/to/your/file.jl` scans the file and tries to guess whether objects are destroyed.
3. Initialize with `log_view = true`: at `finalize`, PETSc lists the objects created and destroyed per class, so a leak shows as a count mismatch.

Objects on a one-rank communicator also have a finalizer, so the garbage collector frees them if you forget. The finalizer does not call PETSc: Julia runs it on whichever thread triggered the collection, and PETSc is not thread safe, so a destroy there could overlap a PETSc call of another thread. It queues the object instead, and the next PETSc.jl call that creates or destroys an object, or `PETSc.finalize`, destroys what has piled up, on that thread. A program that creates many serial objects and then stays away from PETSc for a long time keeps them until its next such call.

## 4. Is it compatible with GPUs?
The precompiled `PETSc_jll` binaries are currently not compatible with GPUs. You can, however, use your own PETSc compilation and things should be fine. 

## 5. I really like this package!
Thanks a lot, we appreciate if you give us a star on github!

## 6. How do I cite this in publications?
We are planning to submit a JOSS paper once the current release is sufficiently stable. Meanwhile, you can cite the [zenodo release](https://doi.org/10.5281/zenodo.18274809).