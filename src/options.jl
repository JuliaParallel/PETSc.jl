# interface to the PETSc options database(s), by providing an OPTIONS[T]
# dictionary-like object that is analogous to the Julia Base.ENV object
# for environment variables.

type Options{T<:Scalar} <: Associative{String,String}; end
const OPTIONS = [T => Options{T}() for T in C.petsc_type]
export OPTIONS, withoptions

typealias SymOrStr Union{AbstractString,Symbol}

function Base.setindex!{T}(::Options{T}, v, k::SymOrStr)
  chk(C.PetscOptionsSetValue(T, string('-',k), string(v)))
  return v
end

function Base.setindex!{T}(::Options{T}, v::Void, k::SymOrStr)
  chk(C.PetscOptionsClearValue(T, string('-',k)))
  return v
end

# PETSc complains if you don't use an option, remove this?
# allow OPTIONS[k]=v to set options for all PETSc scalar types simultaneously
function Base.setindex!(o::typeof(OPTIONS), v, k::SymOrStr)
  for opts in values(o)
    opts[k] = v
  end
  return v
end

const _optionstr = Array(UInt8, 1024)
function Base.get{T}(::Options{T}, k::SymOrStr, def)
  b = Ref{PetscBool}()
  chk(C.PetscOptionsGetString(T, Cstring(Ptr{UInt8}(C_NULL)), string('-',k),
                              pointer(_optionstr), Csize_t(length(_optionstr)),
                              b))
  return b[] != 0 ? unsafe_string(pointer(_optionstr)) : def
end

function Base.haskey{T}(::Options{T}, k::SymOrStr)
  b = Ref{PetscBool}()
  chk(C.PetscOptionsHasName(T, Cstring(Ptr{UInt8}(C_NULL)), string('-',k), b))
  return b[] != 0
end

Base.similar(::Options) = Dict{ByteString,ByteString}()

Base.pop!(o::Options, k::SymOrStr) = (v = o[k]; o[k] = nothing; v)
Base.pop!(o::Options, k::SymOrStr, def) = haskey(o,k) ? pop!(o,k) : def
Base.delete!(o::Options, k::SymOrStr) = (o[k] = nothing; o)
Base.delete!(o::Options, k::SymOrStr, def) = haskey(o,k) ? delete!(o,k) : def

# need to override show: default show function doesn't work because
# there seems to be no way to iterate over the PETSc options database (grr).
Base.show{T}(io::IO, ::Options{T}) = print(io, "PETSc{$T} options database")

# temporarily set some options, call f, and then unset them; like withenv
function withoptions{T<:Scalar}(f::Function, ::Type{T}, keyvals)
  old = Dict{SymOrStr,Any}()
  o = OPTIONS[T]
  for (key,val) in keyvals
    old[key] = get(o,key,nothing)
    val !== nothing ? (o[key]=val) : delete!(o, key)
  end
  try f()
  finally
    for (key,val) in old
      val !== nothing ? (o[key]=val) : delete!(o, key)
    end
  end
end
withoptions{T<:Scalar,K<:SymOrStr}(f::Function, ::Type{T}, keyvals::Pair{K}...) = withoptions(f, T, keyvals)
withoptions{T<:Scalar}(f::Function, ::Type{T}) = f() # handle empty keyvals case; see julia#10853
