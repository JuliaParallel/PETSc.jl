module PETSc

const petsc = haskey(ENV, "PETSC_DIR") ? joinpath(ENV["PETSC_DIR"], get(ENV, "PETSC_ARCH", ""), "lib", "libpetsc") : "libpetsc"

include("error.jl")
include("types.jl")

# Initialize PETSc using command-line ARGS
function PetscInitialize(ARGS...)
    argc = Cint[length(ARGS) + 1]
    sARGS = unshift!(UTF8String[ bytestring(string(a)) for a in ARGS ],"julia")
    args = Ptr{Uint8}[ convert(Ptr{Uint8},s) for s in sARGS ]
    pargs = Ptr{Ptr{Uint8}}[ convert(Ptr{Ptr{Uint8}}, args) ]
    chk(ccall((:PetscInitialize, petsc), PetscErrorCode,
              (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              argc, pargs, C_NULL, C_NULL))
end
PetscInitialize() # TODO: support command-line arguments?

PetscFinalize() = chk(ccall((:PetscFinalize, petsc), PetscErrorCode, ()))

const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,petsc), Ptr{Void}))

include("vec.jl")
include("mat.jl")

comm(o::Union(Vec,Mat)) = ccall((:PetscObjectComm,petsc), Ptr{Void}, (Ptr{Void},), o)

end
