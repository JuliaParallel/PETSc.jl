# this creates a Julia wrapper around a PETSc array
# that holds a pointer to the underlying PETSc data
# and a Julia array that references the data, so we 
# can use Julia array operations on it

"""
    PetscArray{T,N} <: AbstractArray{T,N}
Behaves like a normal Julia array but holds a pointer to PETSc data, which can be destroyed
"""
struct PetscArray{T,N} <: AbstractArray{T,N}
    data::Array{T,N}  # julia array that holds the data
    ptr::Vector  # ptr to the underlying PETSc data (as a vector so we can put it to zero)
end

# ============================================================================
# Array interface - make PetscArray behave like a normal Julia array
# ============================================================================

# Size and shape
Base.size(A::PetscArray) = size(A.data)
Base.axes(A::PetscArray) = axes(A.data)
Base.length(A::PetscArray) = length(A.data)
Base.eltype(::Type{PetscArray{T,N}}) where {T,N} = T
Base.ndims(::Type{PetscArray{T,N}}) where {T,N} = N

# Indexing - linear and cartesian
Base.getindex(A::PetscArray, i::Int) = getindex(A.data, i)
Base.getindex(A::PetscArray, I::Vararg{Int,N}) where {N} = getindex(A.data, I...)
Base.getindex(A::PetscArray, I...) = getindex(A.data, I...)

Base.setindex!(A::PetscArray, v, i::Int) = setindex!(A.data, v, i)
Base.setindex!(A::PetscArray, v, I::Vararg{Int,N}) where {N} = setindex!(A.data, v, I...)
Base.setindex!(A::PetscArray, v, I...) = setindex!(A.data, v, I...)

# Iteration
Base.iterate(A::PetscArray) = iterate(A.data)
Base.iterate(A::PetscArray, state) = iterate(A.data, state)

# Conversion and similarity
Base.similar(A::PetscArray, ::Type{T}, dims::Dims) where {T} = similar(A.data, T, dims)
Base.similar(A::PetscArray, ::Type{T}) where {T} = similar(A.data, T)
Base.similar(A::PetscArray) = similar(A.data)

Base.copy(A::PetscArray) = PetscArray(copy(A.data), A.ptr)
Base.copyto!(dest::PetscArray, src::PetscArray) = (copyto!(dest.data, src.data); dest)
Base.copyto!(dest::Array, src::PetscArray) = copyto!(dest, src.data)
Base.copyto!(dest::PetscArray, src::Array) = (copyto!(dest.data, src); dest)

# Broadcasting
Base.BroadcastStyle(::Type{<:PetscArray}) = Broadcast.ArrayStyle{Array}()
Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Array}}, ::Type{ElType}) where {ElType} = similar(Array{ElType}, axes(bc))

# Show methods
Base.show(io::IO, ::MIME"text/plain", A::PetscArray) = show(io, MIME("text/plain"), A.data)
Base.show(io::IO, A::PetscArray) = show(io, A.data)

# Mathematical operations
Base.:(+)(A::PetscArray, B::PetscArray) = A.data + B.data
Base.:(+)(A::PetscArray, B::Array) = A.data + B
Base.:(+)(A::Array, B::PetscArray) = A + B.data

Base.:(-)(A::PetscArray, B::PetscArray) = A.data - B.data
Base.:(-)(A::PetscArray, B::Array) = A.data - B
Base.:(-)(A::Array, B::PetscArray) = A - B.data
Base.:(-)(A::PetscArray) = -A.data

Base.:(*)(A::PetscArray, B::PetscArray) = A.data * B.data
Base.:(*)(A::PetscArray, B::Array) = A.data * B
Base.:(*)(A::Array, B::PetscArray) = A * B.data
Base.:(*)(A::PetscArray, x::Number) = A.data * x
Base.:(*)(x::Number, A::PetscArray) = x * A.data

Base.:(/)(A::PetscArray, x::Number) = A.data / x
Base.:(\)(A::PetscArray, B::PetscArray) = A.data \ B.data
Base.:(\)(A::PetscArray, B::Array) = A.data \ B
Base.:(\)(A::Array, B::PetscArray) = A \ B.data

# Comparison
Base.:(==)(A::PetscArray, B::PetscArray) = A.data == B.data
Base.:(==)(A::PetscArray, B::Array) = A.data == B
Base.:(==)(A::Array, B::PetscArray) = A == B.data

Base.isapprox(A::PetscArray, B::PetscArray; kwargs...) = isapprox(A.data, B.data; kwargs...)
Base.isapprox(A::PetscArray, B::Array; kwargs...) = isapprox(A.data, B; kwargs...)
Base.isapprox(A::Array, B::PetscArray; kwargs...) = isapprox(A, B.data; kwargs...)

# Array operations
Base.fill!(A::PetscArray, x) = (fill!(A.data, x); A)
Base.sum(A::PetscArray; dims=:) = sum(A.data; dims=dims)
Base.sum(f, A::PetscArray) = sum(f, A.data)
Base.prod(A::PetscArray; dims=:) = prod(A.data; dims=dims)
Base.minimum(A::PetscArray; dims=:) = minimum(A.data; dims=dims)
Base.maximum(A::PetscArray; dims=:) = maximum(A.data; dims=dims)
Base.extrema(A::PetscArray; dims=:) = extrema(A.data; dims=dims)

Base.transpose(A::PetscArray) = transpose(A.data)
Base.adjoint(A::PetscArray) = adjoint(A.data)
Base.permutedims(A::PetscArray, perm) = permutedims(A.data, perm)

# Reshape and views
Base.reshape(A::PetscArray, dims::Dims) = reshape(A.data, dims)
Base.reshape(A::PetscArray, dims::Int...) = reshape(A.data, dims...)
Base.vec(A::PetscArray) = vec(A.data)
Base.view(A::PetscArray, I...) = view(A.data, I...)

# Type conversions
Base.convert(::Type{Array{T,N}}, A::PetscArray{T,N}) where {T,N} = A.data
Base.convert(::Type{Array}, A::PetscArray) = A.data
Base.Array(A::PetscArray) = A.data

# First, last, etc.
Base.first(A::PetscArray) = first(A.data)
Base.last(A::PetscArray) = last(A.data)
Base.firstindex(A::PetscArray) = firstindex(A.data)
Base.lastindex(A::PetscArray) = lastindex(A.data)
Base.eachindex(A::PetscArray) = eachindex(A.data)

# Linear algebra support
import LinearAlgebra
LinearAlgebra.norm(A::PetscArray, p::Real=2) = LinearAlgebra.norm(A.data, p)
LinearAlgebra.dot(A::PetscArray, B::PetscArray) = LinearAlgebra.dot(A.data, B.data)
LinearAlgebra.dot(A::PetscArray, B::Array) = LinearAlgebra.dot(A.data, B)
LinearAlgebra.dot(A::Array, B::PetscArray) = LinearAlgebra.dot(A, B.data)

