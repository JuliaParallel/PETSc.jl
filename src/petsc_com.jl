# common functions (initialization, destruction etsc)

function __init__()

  # initialize MPI first so Petsc won't finalize it
  if !MPI.Initialized()
    MPI.Init()
  end
    
  for i in [Float64, Float32, Complex128]
    PetscInitialize(i)
    println("initialized Petsc library ", i)
  end

  # we want Petsc to return errors to us, rather than using its own
  # error handlers, so that we can catch error codes and throw exceptions
  # need to do this for all Petsc versions
  for i=1:3
    libname = C.petsc_libs[i]
    val = @eval(cglobal((:PetscIgnoreErrorHandler, C.$libname)))
    C.PetscPushErrorHandler(C.petsc_type[i], val, C_NULL)
  end

  # register atexit finalizers for the Petsc libraries
  atexit(() -> C.PetscFinalize(Float64))
  atexit(() -> C.PetscFinalize(Float32))
  atexit(() -> C.PetscFinalize(Complex128))
end


function PetscInitialize(lib::DataType)  
  PetscInitialize(lib, [])
end


function PetscInitialize(lib::DataType, args)
  PetscInitialize(lib, args,"","")
end

function PetscInitialize(lib::DataType, args,filename,help)
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
  arr = Array(UTF8String, length(args))  # convert to C string
  arr2 = Array(Ptr{Uint8}, length(args))  # get array of pointers
  for i = 1:length(args)
    arr[i] = Base.cconvert(Cstring, args[i])
    arr2[i] = pointer(arr[i])
    # = cstring(args[i])
  end
#  ptrs = _jl_pre_exec(arr)
#  err = ccall(Libdl.dlsym(libpetsc, :PetscInitializeNoPointers),Int32,(Int32,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{Uint8}), length(ptrs), ptrs,cstring(filename),cstring(help));

  println("typeof(arr) = ", typeof(arr))
  println("typeof(filename) = ", typeof(bytestring(filename)))
  println("typeof(help) = ", typeof(bytestring(help)))

  val = length(arr)
  val_ptr = Int32[val]
  arr2_tmp = typeof(arr2)[arr2]
  arr2_ptr = pointer(arr2)
  arr2_ptr_tmp = typeof(arr2_ptr)[arr2_ptr]
  arr3_ptr = pointer(arr2_ptr_tmp)

  err = C.PetscInitialize(lib, val_ptr, arr3_ptr, bytestring(filename), bytestring(help))

  if err != 0
    println("non zero return status: ", err)
  end
#=
  # convert returned argc, argv to Julia representation
  val = val_ptr[1]
  println("number of returned arguments = ", val)
  arr2_ret = pointer_to_array(argv_ptr, val)  # own = true?
  arr_ret = Array(Cstring, val)
  for i=1:val
    arr_ret[i] = bytestring(arr2_ret[i])
  end
  
  println("returned argv = ", arr_ret)
=#
#  err = ccall( (:PetscInitializeNoPointers,  libpetsclocation),Int32,(Int32,Ptr{Ptr{Uint8}},Cstring,Cstring), length(arr), arr,filename,help);
  return err
end


function petsc_sizeof(t::C.PetscDataType)
    s = Array(Csize_t, 1)
    chk(PetscDataTypeGetSize(C.petsc_type[1], t, s))
    s[1]
end


