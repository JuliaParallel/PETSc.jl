using MPI

export Scalar
export comm_type
export MPI_Comm
export PetscErrorCode
export PetscBool

global const PETSC_DIR = ENV["PETSC_DIR"]
global const PETSC_ARCH = ENV["PETSC_ARCH"]

global const libpetsclocation = string(PETSC_DIR, "/", PETSC_ARCH, "/lib/", "libpetsc")
global const petsc1 = libpetsclocation # for compatability with auto generated wrappers
global const petsc2 = libpetsclocation
global const petsc3 = libpetsclocation
global const numlibs = 1  # number of libraries actually present
#global const libpetsc = Libdl.dlopen(libpetsclocation)

global const petsc_libs = [:petsc1, :petsc2, :petsc3]
global const petsc_type = [Float64, Float32, Complex128]

typealias Scalar Union(Float32, Float64, Complex128)

MPI_COMM_SELF = MPI.COMM_SELF
typealias MPI_Comm MPI.Comm
typealias comm_type typeof(MPI.COMM_WORLD.val)
