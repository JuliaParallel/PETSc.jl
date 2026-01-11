# Frequently Asked Questions


##### Can I use my own PETSc library?
Yes, see the function `SetPetscLib`. You do need to compile this library as a dynamic library, and `MPI.jl` must be configured to be compatible with the MPI library used to compile PETSc. If you run this on a HPC system, and you don't know exactly which options were used, you can compile one of the PETSc examples and run it with the `-log_view` command line option. At the end of the simulation, it will give you the configuration options used on that machine.  

Please note that the version of PETSc should be compatible 

##### Help, my code crashes?
That is very possible. If you provide a *short* minimum working example (MWE), feel free to open an issue on the github repo, so we can check it out.

##### What about the garbage compiler in julia?
Using the GC in combination with MPI code is a tricky business. The users are therefore responsible to free PETSc objects in the code, as shown in the varous examples. We have some help for that:
1. The julia function `PETSc.audit_petsc_file("path/to/your/file.jl")` which scans your julia file and tries to guess whether objects are destroyed.
2. You can initialize the PETSc library with `log_view=true`. At the end of the code, it will give an  

##### Is it compatible with GPU's?
The precompiled `PETSc_jll` binaries are currently not compatible with GPU's. You can, however, use your own PETSc compilation and things should be fine. 
