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
    msg::AbstractString
    PetscError(n::Integer, m="") = new(convert(PetscErrorCode, n),
                                       PetscErrorMessage(n))
end

# return error message string corresponding to errnum
function PetscErrorMessage(errnum)
#    println("entered PetscErrorMessage")
    arr = Array(Ptr{UInt8}, 1)
    arr_ptr = pointer(arr)
    # use the first petsc library for all cases
    C.PetscErrorMessage(C.petsc_type[1], errnum,  arr_ptr, Ref{Ptr{UInt8}}(0))
#    println("retrieved error message from petsc")
    str = bytestring(arr[1])
#    println("error message = ", str)
#    msg[1] ==  ? "(unknown)" : bytestring(msg[1])
end


