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

"""
    doc_borrowed()

One sentence for the docstring of a reader that hands back a PETSc object.

The handle is borrowed: it belongs to the object it was asked of, carries no
finalizer, and `destroy!` on it is a no-op (see `owns`). `scripts/api_surface.jl
--sweeps` greps for a call to this helper, so a reader returning a PETSc object
without one is reported.
"""
function doc_borrowed()
"""
!!! note "Borrowed handle"
    The returned object is owned by the object it was asked of, not by the
    caller: it carries no finalizer and must not be destroyed. `destroy!` on it
    is a no-op ([`owns`](@ref)).
"""
end

"""
    doc_callback(returns = "")

The paragraph every callback setter's docstring carries: how long the closure
lives, what its return value means and where an exception comes out. 
`returns` describes a return value that carries a meaning; by default it is ignored.
"""
function doc_callback(returns::AbstractString = "")
    ret = isempty(returns) ?
        "its return value is ignored, except that a nonzero `Integer` still fails the call and warns until v0.6" :
        returns
"""
!!! note "Callback"
    The closure is kept with the PETSc object, not with the wrapper passed
    here, so a borrowed handle is fine to pass. It succeeds by returning and
    fails by throwing; $ret. An exception it throws comes out of the `solve!`,
    `step!` or `setup!` that ran it.
"""
end

# `_doc_external` is interpolated into several thousand docstrings in
# src/autowrapped/, which the generator in wrapping/ emits verbatim. Keeping the
# old spelling as an alias leaves those files untouched by the rename.
const _doc_external = doc_external

include("LibPETSc.jl")
using .LibPETSc

# Only types, construction entry points and `petsclibs` are exported (§13); every
# verb and accessor stays qualified. The v0.4 list exported twelve functions and
# no types at all; those names are still reachable qualified (`PETSc.set_library!`)
# and, where they were renamed, through the shims in src/deprecations.jl. The rest
# of the API is marked with `public` in src/public_names.jl (§13.1).
export LibPETSc
export DMDA, DMStag, DMPlex
export PetscVec, PetscMat, PetscOptions
export KSP, SNES, TS
export petsclibs

using Libdl


include("init.jl")
include("callbacks.jl")
include("vec.jl")       
include("mat.jl")          
include("options.jl")
include("ts.jl")
include("ksp.jl")
include("pc.jl")
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
