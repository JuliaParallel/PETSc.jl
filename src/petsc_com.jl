# common functions (initialization, destruction etsc)

function __init__()

  # initialize MPI first so Petsc won't finalize it
  if !MPI.Initialized()
    MPI.Init()
  end

  for i=1:3
    if have_petsc[i]
      PetscInitialize(petsc_type[i])
    end
  end

  # we want Petsc to return errors to us, rather than using its own
  # error handlers, so that we can catch error codes and throw exceptions
  # need to do this for all Petsc versions

  for i=1:3
    if have_petsc[i]
      libname = C.petsc_libs[i]
      val = @eval(cglobal((:PetscIgnoreErrorHandler, C.$libname)))
      C.PetscPushErrorHandler(C.petsc_type[i], val, C_NULL)
    end
  end

  # register atexit finalizers for the Petsc libraries
  atexit() do
    for i=1:3
      if have_petsc[i]
        C.PetscFinalize(petsc_type[i])
      end
    end
  end
end


PetscInitialize(lib::DataType) = PetscInitialize(lib, [""], "", "")
PetscInitialize(lib::DataType, args) = PetscInitialize(lib, args, "", "")

function PetscInitialize(lib::DataType, args, filename, help)
  # argument list starts with program name
  args = filter(a->!isempty(a), args)
  args = map(string, vcat("julia", args))

  nargs = Int32[length(args)]

  # args are passed in as an Ptr{Ptr{Ptr{UInt8}}}
  args_ptr = map(pointer, args)
  args_ptr_ptr = [pointer(args_ptr)]
  args_ptr_ptr_ptr = pointer(args_ptr_ptr)

  err = C.PetscInitialize(lib, nargs, args_ptr_ptr_ptr, String(filename), String(help))
  if err != 0
    error("PETSc initialization error, code: $err")
  end
  return err
end

function petsc_sizeof(T::C.PetscDataType)
  sz = Array(Csize_t, 1)
  chk(C.PetscDataTypeGetSize(C.petsc_type[1], T, sz))
  sz[1]
end

function PetscFinalized(T::Type)
  ret = Array(PetscBool, 1)
  C.PetscFinalized(T, ret)
  ret[1] != 0
end
