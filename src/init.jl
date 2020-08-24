function initialize()
    @chk ccall((:PetscInitializeNoArguments, libpetsc), PetscErrorCode, ())
end

function finalize()
    @chk ccall((:PetscFinalize, libpetsc), PetscErrorCode, ())
end
