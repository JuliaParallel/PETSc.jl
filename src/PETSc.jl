# __precompile__(false)

module PETSc

using MPI, LinearAlgebra, SparseArrays, OffsetArrays

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

#=
include("init.jl")
include("utils.jl")
include("viewer.jl")
include("options.jl")
include("vec.jl")           # not yet wrapped!
include("mat.jl")           # not yet wrapped!
include("matshell.jl")      # not yet wrapped!
include("dm.jl")            # partly wrapped, no tests yet
include("dmda.jl")          # not yet wrapped!
include("dmstag.jl")        # mostly wrapped and tested
include("pc.jl")            # to be fixed/wrapped
include("ksp.jl")           # part is wrapped
include("snes.jl")          # not yet wrapped!
include("sys.jl")

##include("startup.jl")  # can be removed (later)
##include("lib.jl")      # can be removed (later)
##include("ref.jl")      # can be removed (later)

=#

end
