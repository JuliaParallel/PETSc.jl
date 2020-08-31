
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

function initialize()
    map(initialize, scalar_types)
    return nothing
end
function finalize()
    map(finalize, scalar_types)
    return nothing
end
