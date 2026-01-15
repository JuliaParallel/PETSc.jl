# __precompile__(false)

module PETSc

using MPI, LinearAlgebra, SparseArrays, OffsetArrays

MPI.Initialized() || MPI.Init()

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

include("LibPETSc.jl")
using .LibPETSc
export LibPETSc
export audit_petsc_file
export set_petsclib

using Libdl

include("init.jl")
include("vec.jl")       
include("mat.jl")          
include("options.jl")
include("ksp.jl")
include("snes.jl")          
include("dm.jl")          
include("sys.jl")
include("dmda.jl")          
include("dmstag.jl")

# String convenience wrappers for SetType functions
include("string_wrappers.jl")       
include("string_wrappers_extra.jl")


include("audit.jl")



#=
include("utils.jl")
include("viewer.jl")

include("matshell.jl")      # not yet wrapped!
include("dm.jl")            # partly wrapped, no tests yet
include("dmda.jl")          # not yet wrapped!
include("pc.jl")            # to be fixed/wrapped
include("ksp.jl")           # part is wrapped
include("sys.jl")

##include("startup.jl")  # can be removed (later)
##include("lib.jl")      # can be removed (later)
##include("ref.jl")      # can be removed (later)

=#

end
