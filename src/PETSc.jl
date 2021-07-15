module PETSc

using MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

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
include("dm.jl")
include("dmda.jl")
include("ksp.jl")
include("pc.jl")
include("snes.jl")
include("sys.jl")

end
