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

"""
    VecPtr(petsclib, v::CVec, own)

Container type for a PETSc Vec that is just a raw pointer.

If the `own` then the finalizer is set on the vector; calling `destroy` when
`!own` is a no-op.
"""
mutable struct VecPtr{PetscLib} <:
               AbstractPetscVec{PetscLib}
    ptr::CVec
    own::Bool
end
function VecPtr(
    petsclib::PetscLib,
    ptr::CVec,
    own,
) where {PetscLib <: PetscLibType}
    v = VecPtr{PetscLib}(ptr, own)
    #comm = getcomm(v)

    #comm = MPI.Comm()
    comm = LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)

    if own && MPI.Comm_size(comm) == 1
        finalizer(destroy, v)
    end
    return v
end
VecPtr(::Type{PetscLib}, x...) where {PetscLib <: PetscLibType} = VecPtr(getlib(PetscLib), x...)



# =============================================================================
# Multiple dispatch to make PetscVec behave like Julia Vector
# =============================================================================

# Treat PETSc vectors as 1-D array-like objects for broadcasting/indexing
Base.ndims(::Type{<:AbstractPetscVec}) = 1
Base.IndexStyle(::Type{<:AbstractPetscVec}) = IndexLinear()
Base.axes(v::AbstractPetscVec) = (Base.OneTo(length(v)),)
Base.BroadcastStyle(::Type{<:AbstractPetscVec}) = Broadcast.DefaultArrayStyle{1}()

# Array interface - size and length
Base.size(v::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetSize(PetscLib,v)
Base.length(v::AbstractPetscVec{PetscLib}) where {PetscLib} = prod(size(v))
Base.lastindex(v::AbstractPetscVec{PetscLib}) where {PetscLib} = length(v)
Base.similar(v::AbstractPetscVec{PetscLib}) where {PetscLib} =  LibPETSc.VecDuplicate(getlib(PetscLib), v)
type(m::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetType(PetscLib,m)

function Base.getindex(v::AbstractPetscVec{PetscLib}, i::Integer) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]))
    return val[1]
end

# Range indexing
function Base.getindex(v::AbstractPetscVec{PetscLib}, r::AbstractRange) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1))
    return val
end

# Get all values
function Base.getindex(v::AbstractPetscVec{PetscLib}, ::Colon) where {PetscLib}
    n  = length(v)
    return getindex(v, 1:n)
end

Base.isapprox(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}; kwargs...) where {PetscLib} = all(v[:] .≈ w[:]; kwargs...)

function Base.setindex!(v::AbstractPetscVec{PetscLib}, val, i::Integer) where {PetscLib} 
     PetscInt = inttype(PetscLib)
     PetscScalar = PETSc.scalartype(PetscLib)
     LibPETSc.VecSetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]), [PetscScalar(val)], PETSc.INSERT_VALUES)
     return nothing
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, vals, r::AbstractRange) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    LibPETSc.VecSetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1), vals, PETSc.INSERT_VALUES)
    return nothing
end

Base.fill!(v::AbstractPetscVec{PetscLib}, val) where {PetscLib} = LibPETSc.VecSet(PetscLib,v, PetscLib.PetscScalar(val))

# Broadcasting assignment support (dest is not an AbstractArray)
function Base.copyto!(dest::AbstractPetscVec{PetscLib}, bc::Base.Broadcast.Broadcasted) where {PetscLib}
    # Evaluate the broadcasted RHS
    rhs = Base.materialize(bc)
    if rhs isa Number
        # Fast path for scalar RHS
        Base.fill!(dest, rhs)
        return dest
    end
    # Array-like RHS: shape must match
    axes(dest) == axes(rhs) || throw(DimensionMismatch("broadcast axes $(axes(rhs)) do not match destination axes $(axes(dest))"))
    @inbounds for i in eachindex(rhs)
        dest[i] = rhs[i]
    end
    return dest
end

Base.materialize!(dest::AbstractPetscVec, bc::Base.Broadcast.Broadcasted) = (Base.copyto!(dest, bc); dest)

# Iterator interface
Base.iterate(v::AbstractPetscVec{PetscLib})  where {PetscLib} = iterate(v, 1)
function Base.iterate(v::AbstractPetscVec{PetscLib}, state)  where {PetscLib}
    if state > length(v)
        return nothing
    end
    return (v[state], state + 1)
end

"""
    assemble!(A::PetscVec) 

Assembles a PETSc vector after setting values.
"""
function assemble!(A::AbstractPetscVec{PetscLib}) where {PetscLib}
    LibPETSc.VecAssemblyBegin(PetscLib, A)
    LibPETSc.VecAssemblyEnd(PetscLib, A)
end



destroy(m::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecDestroy(PetscLib,m)


"""
    unsafe_localarray(vec::AbstractVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc `vec`

Use `read=false` if the array is write-only; `write=false` if read-only.

!!! note
    `Base.finalize` should be called on the `Array` before the data can be used.

# External Links
$(_doc_external("Vec/VecGetArray"))
$(_doc_external("Vec/VecGetArrayWrite"))
$(_doc_external("Vec/VecGetArrayRead"))
$(_doc_external("Vec/VecRestoreArray"))
$(_doc_external("Vec/VecRestoreArrayWrite"))
$(_doc_external("Vec/VecRestoreArrayRead"))
"""
function unsafe_localarray(
    vec::AbstractPetscVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}
    if write && read
        v = LibPETSc.VecGetArray(PetscLib, vec)
    elseif write
        v = LibPETSc.VecGetArrayWrite(PetscLib, vec)
    elseif read
        v = LibPETSc.VecGetArrayRead(PetscLib, vec)
    else
        error("either read or write should be true")
    end
    
    if write && read
        finalizer(v) do v
            LibPETSc.VecRestoreArray(PetscLib, vec, v)
            return nothing
        end
    elseif write
        finalizer(v) do v
            LibPETSc.VecRestoreArrayWrite(PetscLib, vec, v)         
            return nothing
        end
    elseif read
        finalizer(v) do v
            LibPETSc.VecRestoreArrayRead(PetscLib, vec,  v)
            return nothing
        end
    end
    return v
end


"""
    withlocalarray!(
        f!,
        vecs::NTuple{N, AbstractVec};
        read::Union{Bool, NTuple{N, Bool}} = true,
        write::Union{Bool, NTuple{N, Bool}} = true,
    )

Convert `x` to an `Array{PetscScalar}` using [`unsafe_localarray`](@ref) and
apply the function `f!`.

Use `read=false` if the array is write-only; `write=false` if read-only.

# Examples
```julia-repl
julia> withlocalarray!(x; write=true) do x
   @. x .*= 2
end

julia> withlocalarray!(
           x,
           y;
           read = (false, true),
           write = (true, false)
       ) do x, y
   @. x .= 2 .+ y
end

!!! note
    `Base.finalize` is automatically called on the array.
"""
function withlocalarray!(
    f!,
    vecs::NTuple{N, AbstractPetscVec};
    read::Union{Bool, NTuple{N, Bool}} = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
) where {N}
    read isa NTuple{N, Bool} || (read = ntuple(_ -> read, N))
    write isa NTuple{N, Bool} || (write = ntuple(_ -> write, N))

    arrays = map(vecs, read, write) do v, r, w
        unsafe_localarray(v; read = r, write = w)
    end
    val = f!(arrays...)
    map(arrays) do array
        Base.finalize(array)
    end
    return val
end
withlocalarray!(f!, vecs...; kwargs...) = withlocalarray!(f!, vecs; kwargs...)
