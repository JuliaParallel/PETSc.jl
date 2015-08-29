# PETSc error codes and error handling

typealias PetscErrorCode Cint # from petscsys.h

# wrap this around all PETSc ccalls, in order to catch exceptions
function chk(errnum)
    if errnum != 0
        throw(PetscError(errnum))
    end
    nothing
end

immutable PetscError <: Exception
    errnum::PetscErrorCode
    msg::String
    PetscError(n::Integer, m="") = new(convert(PetscErrorCode, n), 
                                       PetscErrorMessage(n))
end

# return error message string corresponding to errnum
function PetscErrorMessage(errnum)
    msg = Ptr{Uint8}[ C_NULL ]
    ccall((:PetscErrorMessage,petsc), PetscErrorCode,
              (Cint, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}),
              errnum, msg, C_NULL)
    msg[1] == C_NULL ? "(unknown)" : bytestring(msg[1])
end

# we want Petsc to return errors to us, rather than using its own
# error handlers, so that we can catch error codes and throw exceptions
ccall((:PetscPushErrorHandler,petsc), PetscErrorCode, (Ptr{Void},Ptr{Void}),
      cglobal((:PetscIgnoreErrorHandler,petsc)), C_NULL)
