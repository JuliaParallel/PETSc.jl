using MPI

export Scalar
export comm_type
export MPI_Comm
export PetscErrorCode
export PetscBool
export PetscInt # defined in depfile

const depfile = joinpath(dirname(@__FILE__), "..", "..", "deps", "deps.jl")
isfile(depfile) || error("PETSc not properly installed. Please run Pkg.build(\"PETSc\")")
include(depfile)

const petsc_libs = [:petscRealDouble, :petscRealSingle, :petscComplexDouble]
const petsc_type = [Float64, Float32, Complex128]

typealias Scalar Union{Float32, Float64, Complex128}

const MPI_COMM_SELF = MPI.COMM_SELF
typealias MPI_Comm MPI.Comm
typealias comm_type typeof(MPI.COMM_WORLD.val)

# some auxiliary functions used by ccall wrappers
function symbol_get_before(sym_arr)
  ptr_arr = Array(Ptr{UInt8}, length(sym_arr))
  println("ptr_arr = ", ptr_arr)
  for i=1:length(sym_arr)
    println("ptr_arr[$i] = ", ptr_arr[i])
  end

  return pointer(ptr_arr), ptr_arr
end

function symbol_get_after(ptr, sym_arr)
  ptr_arr = pointer_to_array(ptr, length(sym_arr))

  for i=1:length(sym_arr)
    println("ptr_arr[$i] = ", ptr_arr[i])
    sym_arr[i] = bytestring(ptr_arr[i])
  end

end

function symbol_set_before(sym_arr)
  str_arr = similar(sym_arr, UTF8String)

  for i=1:length(str_arr)
    str_arr[i] = bytestring(sym_arr[i])
  end

end
