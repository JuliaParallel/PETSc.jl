
const CPetscOptions = Ptr{Cvoid}

#TODO: should it be <: AbstractDict{String,String}?
abstract type AbstractOptions{T}
end

"""
    GlobalOptions{T}()

The PETSc global options database.
"""
struct GlobalOptions{T} <: AbstractOptions{T}
end
Base.cconvert(::Type{CPetscOptions}, obj::GlobalOptions) = C_NULL

mutable struct Options{T} <: AbstractOptions{T}
    ptr::CPetscOptions
end
Base.cconvert(::Type{CPetscOptions}, obj::Options) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPetscOptions}}, obj::Options) =
    convert(Ptr{CPetscOptions}, pointer_from_objref(obj))

scalartype(::Options{T}) where {T} = T

@for_libpetsc begin
    function Options{$PetscScalar}()
        @assert initialized($petsclib)
        opts = Options{$PetscScalar}(C_NULL)
        @chk ccall((:PetscOptionsCreate, $libpetsc), PetscErrorCode, (Ptr{CPetscOptions},), opts)
        finalizer(destroy, opts)
        return opts
    end
    function destroy(opts::Options{$PetscScalar})
        finalized($petsclib) ||
        @chk ccall((:PetscOptionsDestroy, $libpetsc), PetscErrorCode, (Ptr{CPetscOptions},), opts)
        return nothing
    end

    function Base.push!(::GlobalOptions{$PetscScalar}, opts::Options{$PetscScalar})
        @chk ccall((:PetscOptionsPush, $libpetsc), PetscErrorCode, (CPetscOptions,), opts)
        return nothing
    end
    function Base.pop!(::GlobalOptions{$PetscScalar})
        @chk ccall((:PetscOptionsPop, $libpetsc), PetscErrorCode, ())
        return nothing
    end
    function Base.setindex!(opts::AbstractOptions{$PetscScalar}, val, key)
        @chk ccall((:PetscOptionsSetValue, $libpetsc), PetscErrorCode,
            (CPetscOptions, Cstring, Cstring),
            opts, string('-',key), (val === true || isnothing(val)) ? C_NULL : string(val))
    end

    function view(opts::AbstractOptions{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib))
        @chk ccall((:PetscOptionsView, $libpetsc), PetscErrorCode,
                  (CPetscOptions, CPetscViewer),
                  opts, viewer);
        return nothing
    end
end

"""
    Options{T}(kw => arg, ...)

Create a new PETSc options database.
"""
function Options{T}(ps::Pair...) where {T}
    opts = Options{T}()
    for (k,v) in ps
        if  v!=false
            opts[k] = v
        end
    end
    return opts
end

Base.show(io::IO, opts::AbstractOptions) = _show(io, opts)

"""
    with(f, opts::Options)

Call `f()` with the [`Options`](@ref) `opts` set temporarily (in addition to any global options).
"""
function with(f, opts::Options{T}) where {T}
  global_opts = GlobalOptions{T}()
  push!(global_opts, opts)
  try
    f()
  finally
    pop!(global_opts)
  end
end
