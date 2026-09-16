# __precompile__(false)

module PETSc

using MPI, LinearAlgebra, SparseArrays, OffsetArrays, Preferences

MPI.Initialized() || MPI.Init()

function petsc_link(fname)
"""
[`$fname`](https://petsc.org/release/docs/manualpages/$fname.html)
"""
end

function doc_external(fname)
"""
- PETSc Manual: $(petsc_link(fname))
"""
end

# `_doc_external` is interpolated into several thousand docstrings in
# src/autowrapped/, which the generator in wrapping/ emits verbatim. Keeping the
# old spelling as an alias leaves those files untouched by the rename.
const _doc_external = doc_external

include("LibPETSc.jl")
using .LibPETSc

# The export list is the v0.4 one (§13 replaces it in a later step), minus
# `HostBackend`, which was exported but never defined: host memory is `nothing`,
# not a backend type. The exported names that were renamed stay exported through
# their shims in src/deprecations.jl.
export LibPETSc
export audit_petsc_file
export set_petsclib
export set_library!, unset_library!, library_info
export AbstractPetscMemBackend, AbstractPETScMemBackend
export determine_memtype
export get_petsc_arrays, restore_petsc_arrays
export dmda_star_fd_coloring

using Libdl


include("init.jl")
include("vec.jl")       
include("mat.jl")          
include("options.jl")
include("ts.jl")
include("ksp.jl")
include("snes.jl")          
include("dm.jl")          
include("sys.jl")
include("dmda.jl")
include("dmstag.jl")
include("dmplex.jl")

include("audit_names.jl")   # generated from scripts/renames.jl
include("audit.jl")

include("deprecations.jl")  # generated from scripts/renames.jl
include("public_names.jl")  # generated from scripts/renames.jl



#=
include("utils.jl")
include("viewer.jl")

include("matshell.jl")      # not yet wrapped!
include("dm.jl")            # partly wrapped, no tests yet
include("dmda.jl")          # not yet wrapped!
include("pc.jl")            # to be fixed/wrapped
include("ksp.jl")           # part is wrapped
include("sys.jl")

##include("lib.jl")      # can be removed (later)
##include("ref.jl")      # can be removed (later)

=#

end
