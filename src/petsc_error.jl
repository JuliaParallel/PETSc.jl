# have Petsc report errors to Julia

function error_handler(comm::comm_type, line::Int32, func::Ptr{Uint8}, file::Ptr{Uint8}, n::PetscErrorCode, p::PetscErrorType, mess::Ptr{Uint8}, ctx::Ptr{Void})
# receives the error call from Petsc

func_string = bytestring(func)
file_string = bytestring(file)

if p == PETSC_ERROR_INITIAL
  p_string = "Initial Petsc Error"
elseif p == PETSC_ERROR_REPEAT
  p_string = "Repeat Petsc Error"
elseif  p == PETSC_ERROR_IN_CXX
  p_string == "CXX Petsc Error"
else
  pstring = "Unknown PetscErrorType"
end

mess_string = bytestring(mess)
print_with_color(:red, STDERR, string("\n### ERROR: PETSc Internal Error ###\n"))
print_with_color(:red, STDERR, string("Error in function ", func_string, ", file ", file_string, "\n"))
print_with_color(:red, STDERR, string("Error number: ", n, ", Error type: ", p_string, "\n"))
print_with_color(:red, STDERR, string("Error message: \n", mess_string, "\n"))
print_with_color(:red, STDERR, string("### Error Message Finished ###\n"))

bt = backtrace()
s = sprint(io->Base.show_backtrace(io, bt))
print_with_color(:red, STDERR, string("backtrace: ", s))
print(STDERR, "\n\n")


return PetscErrorCode(0)

end

# tell Petsc abut the Error handler
cfunc = cfunction(error_handler, PetscErrorCode, (comm_type, Int32, Ptr{Uint8}, Ptr{Uint8}, PetscErrorCode, PetscErrorType, Ptr{Uint8}, Ptr{Void}) )
ctx = C_NULL
ierr = ccall( (:PetscPushErrorHandler, petsc), PetscErrorCode, (Ptr{Void}, Ptr{Void}), cfunc, ctx)


function chkerrq(i::PetscErrorCode)
  if i != 0
    error("Petsc Error Code: $i")
  end

  return nothing
end


