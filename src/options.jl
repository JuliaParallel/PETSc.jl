
const CPetscOptions = Ptr{Cvoid}

mutable struct Options{T}
    ptr::CPetscOptions
end
Base.cconvert(::Type{CPetscOptions}, obj::Options) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPetscOptions}}, obj::Options) =
    convert(Ptr{CPetscOptions}, pointer_from_objref(obj))

struct GlobalOptions{T}
end

@for_libpetsc begin
    function Options{$PetscScalar}()
        opts = Options{$PetscScalar}(C_NULL)
        @chk ccall((:PetscOptionsCreate, $libpetsc), PetscErrorCode, (Ptr{CPetscOptions},), opts)
        return opts
    end
    function Base.setindex!(opts::Options{$PetscScalar}, val, key)
        @chk ccall((:PetscOptionsSetValue, $libpetsc), PetscErrorCode,
            (CPetscOptions, Cstring, Cstring), 
            opts, string('-',key), val == true ? C_NULL : string(val))
    end

    function destroy(opts::Options{$PetscScalar})
        @chk ccall((:PetscOptionsDestroy, $libpetsc), PetscErrorCode, (Ptr{CPetscOptions},), opts)
    end

    function Base.push!(::GlobalOptions{$PetscScalar}, opts::Options{$PetscScalar})
        @chk ccall((:PetscOptionsPush, $libpetsc), PetscErrorCode, (CPetscOptions,), opts)
        return nothing
    end
    function Base.pop!(::GlobalOptions{$PetscScalar})
        @chk ccall((:PetscOptionsPop, $libpetsc), PetscErrorCode, ())
        return nothing
    end
end
function Options{T}(ps::Pair...) where {T}
    opts = Options{T}()
    for (k,v) in ps
        opts[k] = v
    end
    return opts
end

function with(f, opts::Options{T}) where {T}
  global_opts = GlobalOptions{T}()
  push!(global_opts, opts)
  try
    f()
  finally
    pop!(global_opts)
  end
end