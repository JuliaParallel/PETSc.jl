# AbstractVec
#   - VecSeq: wrap
#   - VecMPI (TODO)
#   - VecGhost (TODO)
# for the MPI variants we won't be able to attach finalizers, as destroy needs to be called collectively.

const CVec = Ptr{Cvoid}

abstract type AbstractVec{T} <: AbstractVector{T} end
scalartype(::AbstractVec{T}) where {T} = T

# allows us to pass XXVec objects directly into CVec ccall signatures
Base.cconvert(::Type{CVec}, obj::AbstractVec) = obj.ptr
# allows us to pass XXVec objects directly into Ptr{CVec} ccall signatures
Base.unsafe_convert(::Type{Ptr{CVec}}, obj::AbstractVec) =
    convert(Ptr{CVec}, pointer_from_objref(obj))


"""
    VecSeq(v::Vector)

A standard, sequentially-stored serial PETSc vector, wrapping the Julia vector `v`.

This reuses the array `v` as storage, and so `v` should not be `resize!`-ed or otherwise have its length modified while the PETSc object exists.

This should only be need to be called for more advanced uses, for most simple usecases, users should be able to pass `Vector`s directly and have the wrapping performed automatically

# External Links
$(_doc_external("Vec/VecCreateSeqWithArray"))
"""
mutable struct VecSeq{T} <: AbstractVec{T}
    ptr::CVec
    array::Vector{T}
end

"""
    Vec(v::CVec)

Container for an abstract PETSc vector

# External Links
$(_doc_external("Vec/Vec"))
"""
mutable struct Vec{T} <: AbstractVec{T}
    ptr::CVec
end

Base.eltype(::Type{V}) where {V<:AbstractVec{T}} where T = T
Base.eltype(v::AbstractVec{T}) where {T} = T
Base.size(v::AbstractVec) = (length(v),)
Base.parent(v::AbstractVec) = v.array

# this allows setting V[1:2] = 3:4 on a PetscVec (more convenient)
function Base.setindex!(v::AbstractVec, val, I)
    v.array[I]=val
end
Base.getindex(v::AbstractVec, I) = v.array[I]


@for_libpetsc begin
    function VecSeq(comm::MPI.Comm, X::Vector{$PetscScalar}; blocksize=1)
        @assert initialized($petsclib)
        v = VecSeq(C_NULL, X)
        @chk ccall((:VecCreateSeqWithArray, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
                comm, blocksize, length(X), X, v)
        finalizer(destroy, v)
        return v
    end
    function destroy(v::AbstractVec{$PetscScalar})
        finalized($petsclib) ||
        @chk ccall((:VecDestroy, $libpetsc), PetscErrorCode, (Ptr{CVec},), v)
        return nothing
    end
    function Base.length(v::AbstractVec{$PetscScalar})
        r_sz = Ref{$PetscInt}()
        @chk ccall((:VecGetSize, $libpetsc), PetscErrorCode,
          (CVec, Ptr{$PetscInt}), v, r_sz)
        return r_sz[]
    end
    function LinearAlgebra.norm(v::AbstractVec{$PetscScalar}, normtype::NormType=NORM_2)
        r_val = Ref{$PetscReal}()
        @chk ccall((:VecNorm, $libpetsc), PetscErrorCode,
                   (CVec, NormType, Ptr{$PetscReal}),
                   v, normtype,r_val)
        return r_val[]
    end

    function assemblybegin(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyBegin, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end
    function assemblyend(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyEnd, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end
    function assemble(V::AbstractVec{$PetscScalar})
        assemblybegin(V)
        assemblyend(V)
    end

    function ownershiprange(vec::AbstractVec{$PetscScalar})
        r_lo = Ref{$PetscInt}()
        r_hi = Ref{$PetscInt}()
        @chk ccall((:VecGetOwnershipRange, $libpetsc), PetscErrorCode,
          (CVec, Ptr{$PetscInt}, Ptr{$PetscInt}), vec, r_lo, r_hi)
        r_lo[]:(r_hi[]-$PetscInt(1))
    end

    function view(vec::AbstractVec{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, getcomm(vec)))
        @chk ccall((:VecView, $libpetsc), PetscErrorCode,
                    (CVec, CPetscViewer),
                vec.ptr, viewer);
        return nothing
    end

    #    This function just overwrites the same function multiple times
    #    because there is no way to dispatch to the right version based on the
    #    function signature (which operates on Ptr{Void}). Disable for now
    #    since it is not used anywhere.
    #=
    function localsize(cv::CVec)
        r_sz = Ref{$PetscInt}()
        @chk ccall((:VecGetLocalSize, $libpetsc), PetscErrorCode,
            (CVec, Ptr{$PetscInt}), cv, r_sz)
        return r_sz[]
    end
    =#

    function setvalues!(vec::AbstractVec{$PetscScalar},idxs,vals, insertmode::InsertMode)
        idxs = Vector(Int32,idxs);
        vals = Vector(vals);

        @chk ccall((:VecSetValues, $libpetsc), PetscErrorCode,
                 (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar},InsertMode), vec, length(idxs), idxs, vals, insertmode)
        return nothing

    end

    function unsafe_localarray(::Type{$PetscScalar}, cv::CVec; read::Bool=true, write::Bool=true)
        r_pv = Ref{Ptr{$PetscScalar}}()
        if write
            if read
                @chk ccall((:VecGetArray, $libpetsc), PetscErrorCode,
                    (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
            else
                @chk ccall((:VecGetArrayWrite, $libpetsc), PetscErrorCode,
                    (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
            end
        else
            @chk ccall((:VecGetArrayRead, $libpetsc), PetscErrorCode,
                (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
        end
        r_sz = Ref{$PetscInt}()

        @chk ccall((:VecGetLocalSize, $libpetsc), PetscErrorCode,
            (CVec, Ptr{$PetscInt}), cv, r_sz)
        v = unsafe_wrap(Array, r_pv[], r_sz[]; own = false)

        if write
            if read
                finalizer(v) do v
                    @chk ccall((:VecRestoreArray, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                    return nothing
                end
            else
                finalizer(v) do v
                    @chk ccall((:VecRestoreArrayWrite, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                    return nothing
                end
            end
        else
            finalizer(v) do v
                @chk ccall((:VecRestoreArrayRead, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                return nothing
            end
        end
        return v
    end

    function Base.fill!(v::AbstractVec{$PetscScalar}, x)
        @chk ccall((:VecSet, $libpetsc),
                   PetscErrorCode,
                   (CVec, $PetscScalar),
                   v, $PetscScalar(x))
        return v
    end

    function Base.setindex!(
        v::AbstractVec{$PetscScalar},
        val,
        i::Integer,
    )
        @chk ccall(
            (:VecSetValues, $libpetsc),
            PetscErrorCode,
            (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
            v,
            1,
            Ref{$PetscInt}(i - 1),
            Ref{$PetscScalar}(val),
            INSERT_VALUES,
        )

        return val
    end

    function Base.getindex(
        v::AbstractVec{$PetscScalar},
        i::Integer,
    )
        vals = [$PetscScalar(0)]
        @chk ccall(
            (:VecGetValues, $libpetsc),
            PetscErrorCode,
            (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
            v,
            1,
            Ref{$PetscInt}(i - 1),
            vals,
        )

        return vals[1]
    end
end

"""
    unsafe_localarray(PetscScalar, ptr:CVec; read=true, write=true)
    unsafe_localarray(ptr:AbstractVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc data.

Use `read=false` if the array is write-only; `write=false` if read-only.

!!! note
    `Base.finalize` should be called on the `Array` before the data can be used.
"""
unsafe_localarray

unsafe_localarray(v::AbstractVec{T}; kwargs...) where {T} =
    unsafe_localarray(T, v.ptr; kwargs...)

"""
    map_unsafe_localarray!(f!, x::AbstractVec{T}; read=true, write=true)

Convert `x` to an `Array{T}` and apply the function `f!`.

Use `read=false` if the array is write-only; `write=false` if read-only.

# Examples
```julia-repl
julia> map_unsafe_localarray(x; write=true) do x
   @. x .*= 2
end

!!! note
    `Base.finalize` is automatically called on the array.
"""
function map_unsafe_localarray!(f!, v::AbstractVec{T}; kwargs...) where {T}
    array = unsafe_localarray(T, v.ptr; kwargs...)
    f!(array)
    Base.finalize(array)
end



function Base.show(io::IO, ::MIME"text/plain", vec::AbstractVec)
    _show(io, vec)
end

VecSeq(X::Vector{T}; kwargs...) where {T} = VecSeq(MPI.COMM_SELF, X; kwargs...)
AbstractVec(X::AbstractVector) = VecSeq(X)


"""
    ownershiprange(vec::AbstractVec)

The range of indices owned by this processor, assuming that the vectors are laid out with the first n1 elements on the first processor, next n2 elements on the second, etc. For certain parallel layouts this range may not be well defined.

Note: unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(_doc_external("Vec/VecGetOwnershipRange"))
"""
ownershiprange

"""
    setvalues!(
        vector::AbstractVec{PetscScalar},
        indices::Vector{PetscInt},
        vals::Vector{PetscScalar},
        mode::InsertMode;
        num_vals = length(ind)
    )

Insert a set of values into the `vector`. Equivalent to one of the following
depending on the `mode`
```julia
vector[indices[1:num_vals]] .= vals[1:num_vals]
vector[indices[1:num_vals]] .+= vals[1:num_vals]
```

!!! warning
    `indices` should use 0-based indexing!

# External Links
$(_doc_external("Vec/VecSetValues"))
"""
function setvalues!(::AbstractVec) end

@for_libpetsc function setvalues!(
    vec::AbstractVec{$PetscScalar},
    inds::Vector{$PetscInt},
    vals::Vector{$PetscScalar},
    mode::InsertMode;
    num_vals = length(inds)
)
    @chk ccall(
         (:VecSetValues, $libpetsc),
         PetscErrorCode,
         (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
         vec,
         num_vals,
         inds,
         vals,
         mode,
    )
    return vals
end

"""
    getvalues!(
        vector::AbstractVec{PetscScalar},
        indices::Vector{PetscInt},
        vals::Vector{PetscScalar};
        num_vals = length(inds)
    )

Get a set of values from the `vector`. Equivalent to one of the following
```julia
vals[1:num_vals] .= vector[indices[1:num_vals]]
```

!!! warning
    `indices` should use 0-based indexing!

# External Links
$(_doc_external("Vec/VecGetValues"))
"""
function getvalues!(::AbstractVec) end

@for_libpetsc function getvalues!(
    vec::AbstractVec{$PetscScalar},
    inds::Vector{$PetscInt},
    vals::Vector{$PetscScalar};
    num_vals = length(inds)
)
    @chk ccall(
         (:VecGetValues, $libpetsc),
         PetscErrorCode,
         (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
         vec,
         num_vals,
         inds,
         vals,
    )
    return vals
end
