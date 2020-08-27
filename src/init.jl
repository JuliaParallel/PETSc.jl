
@for_libpetsc begin
    function initialize(::Type{$PetscScalar})
        @chk ccall((:PetscInitializeNoArguments, $libpetsc), PetscErrorCode, ())
    end
    function finalize(::Type{$PetscScalar})
        @chk ccall((:PetscFinalize, $libpetsc), PetscErrorCode, ())
    end 
end


function initialize()
    for T in scalar_types
        initialize(T)
    end
end
function finalize()
    for T in scalar_types
        finalize(T)
    end
end