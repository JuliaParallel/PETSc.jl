using MPI

export Scalar
export comm_type
export MPI_Comm
export PetscErrorCode
export PetscBool
export PetscInt

include("lib_locations.jl")
#=
#global const PETSC_DIR = ENV["PETSC_DIR"]
#global const PETSC_ARCH = ENV["PETSC_ARCH"]

pkgdir = Pkg.dir("PETSc")
petsc_version = "3.6.0"
petsc_arch = "arch-linux2-c-debug"
#global const libpetsclocation = string(PETSC_DIR, "/", PETSC_ARCH, "/lib/", "libpetsc")
#global const petsc1 = libpetsclocation
global const petsc1 = joinpath(pkgdir, "deps/petsc-$petsc_version/$petsc_arch", "lib/libpetsc")
println("petsc1 = ", petsc1)
global const petsc2 = petsc1
global const petsc3 = petsc1
global const numlibs = 1  # number of libraries actually present
#global const libpetsc = Libdl.dlopen(libpetsclocation)
=#

numlibs = 1


global const petsc_libs = [:petscRealDouble, :petscRealSingle, :petscComplexDouble]
global const petsc_type = [Float64, Float32, Complex128]

typealias Scalar Union{Float32, Float64, Complex128}

MPI_COMM_SELF = MPI.COMM_SELF
typealias MPI_Comm MPI.Comm
typealias comm_type typeof(MPI.COMM_WORLD.val)
typealias PetscInt Int64

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



