# __precompile__(false)

module PETSc

using MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

include("LibPETSc.jl")
using .LibPETSc
export LibPETSc

using Libdl

function _petsc_link(fname)
"""
[`$fname`](https://petsc.org/release/docs/manualpages/$fname.html)
"""
end

function _doc_external(fname)
"""
- PETSc Manual: $(_petsc_link(fname))
"""
end

include("init.jl")
include("viewer.jl")
include("options.jl")
include("vec.jl")       # not yet autowrapped!
include("mat.jl")
include("matshell.jl")

include("sys.jl")

##include("startup.jl")  # can be removed (later)
##include("lib.jl")      # can be removed (later)
##include("ref.jl")      # can be removed (later)

#=
include("dm.jl")
include("dmda.jl")
include("dmstag.jl")
include("ksp.jl")
include("pc.jl")
include("snes.jl")
=#


end
