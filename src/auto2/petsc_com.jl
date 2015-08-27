# common functions (initialization, destruction etsc)

function __init__()
  
  PetscInitialize([])
end

function PetscInitialize(args)
  PetscInitialize(args,"","")
end

function PetscInitialize(args,filename,help)
  # argument list starts with program name
  args = ["julia";args];
  println("typeof(args) = ", typeof(args))
  #
  #   If the user forgot to PetscFinalize() we do it for them, before restarting PETSc
  #

#=
 init = PetscInitialized()
 if (init != 0)
#    gc() # call garbage collection to force all PETSc objects be destroy that are queued up for destruction
#    err = ccall( (:PetscFinalize,  libpetsclocation),Int32,());if (err != 0) return err; end
    err = PetscFinalize()

    if err != 0
      return err
    end
  end
=#
  # convert arguments to Cstring
  arr = Array(ASCIIString, length(args))
  for i = 1:length(args)
    arr[i] = Base.cconvert(Cstring, args[i])
    # = cstring(args[i])
  end
#  ptrs = _jl_pre_exec(arr)
#  err = ccall(Libdl.dlsym(libpetsc, :PetscInitializeNoPointers),Int32,(Int32,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{Uint8}), length(ptrs), ptrs,cstring(filename),cstring(help));

  println("typeof(arr) = ", typeof(arr))
  println("typeof(filename) = ", typeof(bytestring(filename)))
  println("typeof(help) = ", typeof(bytestring(help)))

  err = C.PetscInitializeNoPointers(Float64, Int32(length(arr)), arr, bytestring(filename), bytestring(help))

  if err != 0
    println("non zero return status: ", err)
  end
#  err = ccall( (:PetscInitializeNoPointers,  libpetsclocation),Int32,(Int32,Ptr{Ptr{Uint8}},Cstring,Cstring), length(arr), arr,filename,help);
  return err
end


function petsc_sizeof(t::C.PetscDataType)
    s = Array(Csize_t, 1)
    chk(PetscDataTypeGetSize(C.petsc_type[1], t, s))
    s[1]
end


