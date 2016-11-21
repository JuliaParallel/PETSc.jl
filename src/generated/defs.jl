using MPI

export Scalar
export comm_type
export MPI_Comm
export PetscErrorCode
export PetscBool
export PetscInt # defined in depfile
export have_petsc, petsc_libs, petsc_type
const depfile = joinpath(dirname(@__FILE__), "..", "..", "deps", "deps.jl")
isfile(depfile) || error("PETSc not properly installed. Please run Pkg.build(\"PETSc\")")
include(depfile)

const petsc_libs = [:petscRealDouble, :petscRealSingle, :petscComplexDouble]
const petsc_type = [Float64, Float32, Complex128]

typealias Scalar Union{Float32, Float64, Complex128}

const MPI_COMM_SELF = MPI.COMM_SELF
typealias MPI_Comm Union{MPI.CComm, MPI.Comm}
typealias comm_type MPI.CComm

# some auxiliary functions used by ccall wrappers
# get an array of pointers to UInt8s that is the same shape as
# the Symbol array
# does *not* allocate the pointers
function Symbol_get_before(sym_arr)
  ptr_arr = Array(Ptr{UInt8}, length(sym_arr))
#  println("ptr_arr = ", ptr_arr)
#  for i=1:length(sym_arr)
#    println("ptr_arr[$i] = ", ptr_arr[i])
#  end

  return pointer(ptr_arr), ptr_arr
end

# turn array of strings (UInt8 *) and  puts them into Symbol array
function Symbol_get_after(ptr, sym_arr)
  ptr_arr = pointer_to_array(ptr, length(sym_arr))

  for i=1:length(sym_arr)
    sym_arr[i] = string(ptr_arr[i])
  end

end

function Symbol_set_before(sym_arr)
  str_arr = similar(sym_arr, UTF8String)

  for i=1:length(str_arr)
    str_arr[i] = string(sym_arr[i])
  end

end

# TODO: auto-generate these
function PetscObjectComm(::Type{Float64},arg1::Ptr)
   ccall((:PetscObjectComm,petscRealDouble),comm_type,(Ptr{Void},),arg1)
end
function PetscObjectComm(::Type{Float32},arg1::Ptr)
   ccall((:PetscObjectComm,petscRealSingle),comm_type,(Ptr{Void},),arg1)
end
function PetscObjectComm(::Type{Complex128},arg1::Ptr)
   ccall((:PetscObjectComm,petscComplexDouble),comm_type,(Ptr{Void},),arg1)
end
