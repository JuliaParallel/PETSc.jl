# PETSc error codes and error handling
export chk, PetscError, PetscErrorMessage
#typealias PetscErrorCode Cint # from petscsys.h

# wrap this around all PETSc ccalls, in order to catch exceptions
function chk(errnum)
  if errnum != 0
    throw(PetscError(errnum))
  end
  errnum
end

immutable PetscError <: Exception
  errnum::PetscErrorCode
  msg::AbstractString

  function PetscError(n::Integer, msg::AbstractString="")
    if isempty(msg)
      new(PetscErrorCode(n), PetscErrorMessage(n))
    else
      new(PetscErrorMessage(n), msg)
    end
  end
end

# return error message string corresponding to errnum
function PetscErrorMessage(errnum)
  msg_ptr = Ref{Ptr{UInt8}}(C_NULL)
  # use the first petsc library for all cases
  println(STDERR, "errnum = ", errnum)
  C.PetscErrorMessage(C.petsc_type[1], errnum, msg_ptr, Ref{Ptr{UInt8}}(C_NULL))
  # error num strings are stored in a constant table so
  # the pointer does not have to be free'd here
  return bytestring(msg_ptr[])
end
