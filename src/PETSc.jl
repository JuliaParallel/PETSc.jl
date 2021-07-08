module PETSc

using MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

using Libdl

include("const.jl")
include("startup.jl")
include("lib.jl")
include("init.jl")
include("ref.jl")
include("viewer.jl")
include("options.jl")
include("vec.jl")
include("mat.jl")
include("matshell.jl")
include("ksp.jl")
include("pc.jl")
include("snes.jl")

end
