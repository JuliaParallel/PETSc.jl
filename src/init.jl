
@for_libpetsc begin
    function initialized(::Type{$PetscScalar})
        r_flag = Ref{PetscBool}()
        @chk ccall((:PetscInitialized, $libpetsc), PetscErrorCode, (Ptr{PetscBool},), r_flag)
        return r_flag[] == PETSC_TRUE
    end
    function initialize(::Type{$PetscScalar})
        if !initialized($PetscScalar)
            MPI.Initialized() || MPI.Init()
            @chk ccall((:PetscInitializeNoArguments, $libpetsc), PetscErrorCode, ())

            # disable signal handler
            @chk ccall((:PetscPopSignalHandler, $libpetsc), PetscErrorCode, ())

            atexit(() -> finalize($PetscScalar))
        end
        return nothing
    end
    function finalized(::Type{$PetscScalar})
        r_flag = Ref{PetscBool}()
        @chk ccall((:PetscFinalized, $libpetsc), PetscErrorCode, (Ptr{PetscBool},), r_flag)
        return r_flag[] == PETSC_TRUE
    end
    function finalize(::Type{$PetscScalar} )
        if !finalized($PetscScalar)
            @chk ccall((:PetscFinalize, $libpetsc), PetscErrorCode, ())
        end
        return nothing
    end 
end

"""
    PETSc.initialize([T])

Initialize PETSc for scalar type `T`, if not already initialized. If no `T` is
provided, then it will be initialized for all supported scalar types.

Additionally:
 - This will initialize MPI if it has not already been initialized.
 - It will disable the PETSc signal handler (via
   [`PetscPopSignalHandler`](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscPopSignalHandler.html))
 - Add an [`atexit`](https://docs.julialang.org/en/v1/base/base/#Base.atexit)
   hook to call [`PETSc.finalize`](@ref).
"""
function initialize()
    map(initialize, scalar_types)
    return nothing
end


"""
    PETSc.finalize([T])

Finalize PETSc for scalar type `T`. If no `T` is provided, then it will be finalized for all supported scalar types.

It is generally not necessary to call this function directly, as it is added as an [`atexit`](https://docs.julialang.org/en/v1/base/base/#Base.atexit) hook in [`PETSc.initialize`](@ref).
"""
function finalize()
    map(finalize, scalar_types)
    return nothing
end
