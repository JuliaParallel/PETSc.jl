#
#   Allows the PETSc dynamic library to be used from Julia. (http::/julialang.org)
#     PETSc must be configured with --with-shared-libraries
#     You can build with or without MPI, but cannot run on more than one process
#
#  Load the PETSc dynamic library
#

module PETSc


using MPI
include("petsc_constants.jl")


# export all names
export PetscInitialize, getPETSC_COMM_SELF, PetscView, PetscIS, PetscDestroy, PetscISSetType, PetscISGetSizePetscISGetIndices, PetscDataTypeFromString, PetscDataTypeGetSize

export KSP, KSPSetOperators, KSPSetFromOptions, KSPSolve,  KSPSetUp


# -------------------------------------
function echodemo(filename)
  f = open(filename)
  h = readall(f)
  close(f)
  pos = 0
  while (pos <= length(h))
    (ex,pos)=parse(h,pos)
    str = string(ex.args)
    if _jl_have_color
      print("\033[1m\033[30m")
    end
    #  the next line doesn't work right for multiple line commands like for loops
    println(str[2:strlen(str)-1])
    if _jl_have_color
      print(_jl_answer_color())
      println(" ") # force a color change, otherwise answer color is not used
    end
    e = eval(ex)
    println(" ")
  end
end

# -------------------------------------


#=
PETSC_INSERT_VALUES = 1;
PETSC_ADD_VALUES    = 2;
PETSC_COPY_VALUES   = 0;

PETSC_NORM_1         = 0;
PETSC_NORM_2         = 1;
PETSC_NORM_FROBENIUS = 2;
PETSC_NORM_INFINITY  = 3;
PETSC_NORM_MAX       = PETSC_NORM_INFINITY;
=#

# -------------------------------------
#    These are the Julia interface methods for all the visible PETSc functions. Julia datatypes for PETSc objects simply contain the C pointer to the
#    underlying PETSc object.
#
# -------------------------------------
function PetscInitialize()
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
  init = 0;
  err = ccall( (:PetscInitialized,  libpetsclocation),Int32,(Ptr{Int32},),&init);
  if (init != 0)
    gc() # call garbage collection to force all PETSc objects be destroy that are queued up for destruction
    err = ccall( (:PetscFinalize,  libpetsclocation),Int32,());if (err != 0) return err; end
  end
  arr = Array(ByteString, length(args))
  for i = 1:length(args)
    arr[i] = Base.cconvert(Cstring, args[i])
    # = cstring(args[i])
  end
#  ptrs = _jl_pre_exec(arr)
#  err = ccall(Libdl.dlsym(libpetsc, :PetscInitializeNoPointers),Int32,(Int32,Ptr{Ptr{Uint8}},Ptr{Uint8},Ptr{Uint8}), length(ptrs), ptrs,cstring(filename),cstring(help));

  err = ccall( (:PetscInitializeNoPointers,  libpetsclocation),Int32,(Int32,Ptr{Ptr{Uint8}},Cstring,Cstring), length(arr), arr,filename,help);
  return err
end

function PetscFinalize()
  gc() # call garbage collection to force all PETSc objects be destroy that are queued up for destruction
  return ccall( (:PetscFinalize,  libpetsclocation),Int32,());
end

function getPETSC_COMM_SELF()
  comm = Array(Int64, 1)
  err = ccall( (:PetscGetPETSC_COMM_SELF,  libpetsclocation),Int32,(Ptr{Int64},),comm);
#  return comm[1]
   return MPI.COMM_WORLD.val
end

function PetscDataTypeFromString(name::AbstractString)
    ptype = Array(Cint, 1)
    found = Array(PetscBool, 1)
    ccall((:PetscDataTypeFromString,petsc),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}), name, ptype, found)

    return ptype[1], convert(Bool, found[1])
end


function PetscDataTypeGetSize(dtype::PetscDataType)
    datasize = Array(Csize_t, 1)
    ccall((:PetscDataTypeGetSize,petsc),PetscErrorCode,(PetscDataType,Ptr{Csize_t}), dtype, datasize)

    return datasize[1]
end






# -------------------------------------
#
abstract PetscObject

function PetscView(obj::PetscObject)
  PetscView(obj,0)
end

type PetscIS <: PetscObject
  pobj::Int64
  function PetscIS(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    is = Array(Int64,1)
    err = ccall( (:ISCreate,  libpetsclocation),Int32,(Int64,Ptr{Int64}),comm.val,is);if (err != 0) return err;end
    is = new(is[1])
    finalizer(is,PetscDestroy)
    # does not seem to be called immediately when is is no longer visible, is it called later during garbage collection?
    return is
  end
end

  function PetscDestroy(is::PetscIS)
    if (is.pobj != 0) then
      err = ccall( ( :ISDestroy,  libpetsclocation),Int32,(Ptr{Int64},), &is.pobj);
    end
    is.pobj = 0
    println("ISDestroy called")
    return 0
  end

  function PetscIS(indices::Array{Int64})
    is = PetscIS()
    err = ccall( ( :ISSetType,  libpetsclocation),Int32,(Int64,Cstring), is.pobj,"general");
    err = ccall( (:ISGeneralSetIndices,  libpetsclocation),Int32,(Int64,Int32,Ptr{Int32},Int32),is.pobj,length(indices),convert(Array{Int32},indices),PETSC_COPY_VALUES)
    return is
  end

  function PetscISSetType(vec::PetscIS,name)
    err = ccall( (:ISSetType,  libpetsclocation),Int32,(Int64,Cstring), vec.pobj,name);
  end

  function PetscView(obj::PetscIS,viewer)
   err = ccall( ( :ISView,  libpetsclocation),Int32,(Int64,Int64),obj.pobj,0);
  end

  function PetscISGetSize(obj::PetscIS)
    n = Array(Int32,1)
    err = ccall( ( :ISGetSize,  libpetsclocation),Int32,(Int64,Ptr{Int32}), obj.pobj,n);
    return n[1]
  end

  function PetscISGetIndices(obj::PetscIS)
    len = PetscISGetSize(obj)
    indices = Array(Int32,len);
    err = ccall( (:ISGetIndicesCopy,  libpetsclocation),Int32,(Int64,Ptr{Int32}),obj.pobj,indices);
    indices = indices + 1
    return indices
  end

# -------------------------------------
#

include("petsc_vec.jl")

# -------------------------------------
include("petsc_mat.jl")

type KSP <: PetscObject
  pobj::Ptr{Void}

  function KSP(comm::MPI_Comm)
      ptr = Array(Ptr{Void}, 1)
      ierr = ccall((:KSPCreate,petsc),PetscErrorCode,(Int64,Ptr{Void}),comm.val, ptr)
      @assert(ierr == 0)
      obj = new(ptr[1])
#      finalizer(obj, PetscDestroy)
      return obj
  end
end


function PetscDestroy(ksp::KSP)
    err = ccall((:KSPDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),&ksp.pobj)
    println("KSPDestroy called")
#    sleep(5)
end


function KSPSetOperators(ksp::KSP,Amat::PetscMat,Pmat::PetscMat)
   err = ccall((:KSPSetOperators,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, Amat.pobj, Pmat.pobj)
end


function KSPSetFromOptions(ksp::KSP)
    ccall((:KSPSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},),ksp.pobj)
end

function KSPSolve(ksp::KSP, b::PetscVec, x::PetscVec)
    err = ccall((:KSPSolve,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, b.pobj, x.pobj)
end

function KSPSetUp(ksp::KSP)
    err = ccall((:KSPSetUp,petsc),PetscErrorCode,(Ptr{Void},), ksp.pobj)
end




end  # end module


