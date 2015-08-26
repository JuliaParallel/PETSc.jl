using MPI


global const PETSC_DIR = ENV["PETSC_DIR"]
global const PETSC_ARCH = ENV["PETSC_ARCH"]

global const libpetsclocation = string(PETSC_DIR, "/", PETSC_ARCH, "/lib/", "libpetsc")
global const petsc1 = libpetsclocation # for compatability with auto generated wrappers
#global const libpetsc = Libdl.dlopen(libpetsclocation)




MPI_COMM_SELF = MPI.COMM_SELF
typealias MPI_Comm MPI.Comm
typealias comm_type typeof(MPI.COMM_WORLD.val)
