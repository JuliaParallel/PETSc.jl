# these are julia-specific PETSc functions, which makes PETSc a bit easier to use
# a (nearly) full wrapper is in wrapping/

import .LibPETSc: AbstractPetscVec, PetscVec, CVec

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscVec{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc Vec (null pointer)")
        return
    else
        # Get size using PETSc function (commented out until VecGetSize is available)
        si = LibPETSc.VecGetSize(PetscLib,v)
        ty = LibPETSc.VecGetType(PetscLib,v)
        print(io, "PETSc $ty Vec; length=$si")
    end
    return nothing
end

# =============================================================================
# Multiple dispatch to make PetscVec behave like Julia Vector
# =============================================================================

# Array interface - size and length
Base.size(v::PetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetSize(PetscLib,v)
Base.length(v::PetscVec{PetscLib}) where {PetscLib} = prod(size(v))
Base.lastindex(v::PetscVec{PetscLib}) where {PetscLib} = length(v)
type(m::PetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetType(PetscLib,m)

function Base.getindex(v::PetscVec{PetscLib}, i::Integer) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]))
    return val[1]
end

# Range indexing
function Base.getindex(v::PetscVec{PetscLib}, r::AbstractRange) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1))
    return val
end

function Base.setindex!(v::PetscVec{PetscLib}, val, i::Integer) where {PetscLib} 
     PetscInt = inttype(PetscLib)
     LibPETSc.VecSetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]), [val], PETSc.INSERT_VALUES)
     return nothing
end

function Base.setindex!(v::PetscVec{PetscLib}, vals, r::AbstractRange) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    LibPETSc.VecSetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1), vals, PETSc.INSERT_VALUES)
    return nothing
end

Base.fill!(v::PetscVec{PetscLib}, val) where {PetscLib} = LibPETSc.VecSet(PetscLib,v, val)

# Iterator interface
Base.iterate(v::PetscVec{PetscLib})  where {PetscLib} = iterate(v, 1)
function Base.iterate(v::PetscVec{PetscLib}, state)  where {PetscLib}
    if state > length(v)
        return nothing
    end
    return (v[state], state + 1)
end

"""
    assemble!(A::PetscVec) 

Assembles a PETSc vector after setting values.
"""
function assemble!(A::PetscVec{PetscLib}) where {PetscLib}
    LibPETSc.VecAssemblyBegin(PetscLib, A)
    LibPETSc.VecAssemblyEnd(PetscLib, A)
end



destroy(m::PetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecDestroy(PetscLib,m)
