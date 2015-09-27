# PETSc error codes and error handling

#typealias PetscErrorCode Cint # from petscsys.h

# wrap this around all PETSc ccalls, in order to catch exceptions
function chk(errnum)
    if errnum != 0
         println(STDERR, "Error number $errnum has occured")
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
#    println("entered PetscErrorMessage")
    arr = Array(Ptr{Uint8}, 1)
    arr_ptr = pointer(arr)
    # use the first petsc library for all cases
    C.PetscErrorMessage(C.petsc_type[1], errnum,  arr_ptr, Ref{Ptr{Uint8}}(0))
#    println("retrieved error message from petsc") 
    str = bytestring(arr[1])
#    println("error message = ", str)
#    msg[1] ==  ? "(unknown)" : bytestring(msg[1])
end

# we want Petsc to return errors to us, rather than using its own
# error handlers, so that we can catch error codes and throw exceptions
# need to do this for all Petsc versions
for i=1:C.numlibs
  libname = C.petsc_libs[i]
  val = @eval(cglobal((:PetscIgnoreErrorHandler, C.$libname)))
  C.PetscPushErrorHandler(C.petsc_type[i], val, C_NULL)
end

