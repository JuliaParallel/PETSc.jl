module PETSc
using MPI

include("LibPETSc.jl")
using .LibPETSc

using MPI, LinearAlgebra, SparseArrays

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
# include("ref.jl")
include("vec.jl")
include("mat.jl")
include("matshell.jl")
include("dm.jl")
include("dmda.jl")
# include("dmstag.jl")
include("ksp.jl")
# include("pc.jl")
# include("snes.jl")
include("sys.jl")

end
