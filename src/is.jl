# index sets and vector scatters

###########################################################################
export IS # index sets
# Note: we expose a 1-base Julian index interface, but internally
# PETSc's indices are 0-based.
#TODO: support block versions
type IS{T}
  p::C.IS{T}
  function IS(p::C.IS{T})
    o = new(p)
    finalizer(o, PetscDestroy)
    return o
  end
end

comm{T}(a::IS{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function PetscDestroy{T}(o::IS{T})
  PetscFinalized(T) || C.ISDestroy(Ref(o.p))
end

# internal constructor, takes array of zero-based indices:
function IS_{T<:Scalar}(::Type{T}, idx::Array{PetscInt}; comm::MPI.Comm=MPI.COMM_SELF)
  is_c = Ref{C.IS{T}}()
  chk(C.ISCreateGeneral(comm, length(idx), idx, C.PETSC_COPY_VALUES, is_c))
  return IS{T}(is_c[])
end

IS{I<:Integer, T<:Scalar}(::Type{T}, idx::AbstractArray{I}; comm::MPI.Comm=MPI.COMM_SELF) =
  IS_(T, PetscInt[i-1 for i in idx]; comm=comm)

function IS{I<:Integer, T<:Scalar}(::Type{T}, idx::Range{I}; comm::MPI.Comm=MPI.COMM_SELF)
  is_c = Ref{C.IS{T}}()
  chk(C.ISCreateStride(comm, length(idx), start(idx)-1, step(idx), is_c))
  return IS{T}(is_c[])
end

#function ISBlock{I<:Integer, T<:Scalar}(::Type{T}, idx::

function Base.copy{T}(i::IS{T})
  is_c = Ref{C.IS{T}}()
  chk(C.ISDuplicate(i.p, is_c))
  return IS{T}(is_c[])
end

function Base.length(i::IS)
  len = Ref{PetscInt}()
  chk(C.ISGetSize(i.p, len))
  return Int(len[])
end

function lengthlocal(i::IS)
  len = Ref{PetscInt}()
  chk(C.ISGetLocalSize(i.p, len))
  return Int(len[])
end

import Base.==
function =={T}(i::IS{T}, j::IS{T})
  b = Ref{PetscBool}()
  chk(C.ISEqual(i.p, j.p, b))
  return b[] != 0
end

function Base.sort!(i::IS)
  chk(C.ISSort(i.p))
  return i
end
Base.sort(i::IS) = sort!(copy(i))

function Base.issorted(i::IS)
  b = Ref{PetscBool}()
  chk(C.ISSorted(i.p, b))
  return b[] != 0
end

function Base.union{T}(i::IS{T}, j::IS{T})
  is_c = Ref{C.IS{T}}()
  chk(C.ISExpand(i.p, j.p, is_c))
  return IS{T}(is_c[])
end

function Base.setdiff{T}(i::IS{T}, j::IS{T})
  is_c = Ref{C.IS{T}}()
  chk(C.ISDifference(i.p, j.p, is_c))
  return IS{T}(is_c[])
end

function Base.extrema(i::IS)
  min = Ref{PetscInt}()
  max = Ref{PetscInt}()
  chk(C.ISGetMinMax(i.p, min, max))
  return (Int(min[])+1, Int(max[])+1)
end
Base.minimum(i::IS) = extrema(i)[1]
Base.maximum(i::IS) = extrema(i)[2]

function Base.convert{T<:Integer}(::Type{Vector{T}}, idx::IS)
  pref = Ref{Ptr{PetscInt}}()
  chk(C.ISGetIndices(idx.p, pref))
  inds = Int[i+1 for i in pointer_to_array(pref[], lengthlocal(idx))]
  chk(C.ISRestoreIndices(idx.p, pref))
  return inds
end
Base.Set(i::IS) = Set(Vector{Int}(i))

export set_blocksize, get_blocksize

function set_blocksize(is::IS, bs::Integer)
  chk(C.ISSetBlockSize(is.p, bs))
end

function get_blocksize(is::IS)
  bs = Ref{PetscInt}()
  chk(C.ISGetBlockSize(is.p, bs))
  return Int(bs[])
end

function petscview{T}(is::IS{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.ISView(is.p, viewer))
end


###############################################################################
# we expose a 1 based API, but internally ISLoalToGlobalMappings are zero based

export ISLocalToGlobalMapping

type ISLocalToGlobalMapping{T}
  p::C.ISLocalToGlobalMapping{T}
  function ISLocalToGlobalMapping(p::C.ISLocalToGlobalMapping{T})
    o = new(p)
    finalizer(o, PetscDestroy)
    return o
  end
end


# zero based, not exported
function _ISLocalToGlobalMapping{T}(::Type{T}, indices::AbstractArray{PetscInt}, bs=1; comm=MPI_COMM_WORLD, copymode=C.PETSC_COPY_VALUES)

  isltog = Ref{C.ISLocalToGlobalMapping}()
  chk(C.ISLocalToGlobalMappingCreate(comm, bs, length(indices), indices, copymode, isltog))

  return ISLocalToGlobalMapping(isltog[])
end

# one based, exported
#TODO: add a data argument to ISLocalToGlobalMapping, to store intermediate
# array for copymode = don't copy
function ISLocalToGlobalMapping{T, I <: Integer}(::Type{T}, indices::AbstractArray{I}, bs=1; comm=MPI_COMM_WORLD, copymode=C.PETSC_COPY_VALUES)

  indices_0 = PetscInt[ i-1 for i in indices]
  return _ISLocalToGlobalMapping(t, indices, bs=bs, comm=comm, copymode=copymode)

end

function ISLocalToGlobalMapping{T}(is::IS{T})

  isltog = Ref{C.ISLocalToGlobalMapping{T}}()
  chk(C.ISLocalToGlobalMappingCreateIS(is.p, isltog))
  return ISLocalToGlobalMapping{T}(isltog[])
end


comm{T}(a::ISLocalToGlobalMapping{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function PetscDestroy{T}(o::ISLocalToGlobalMapping{T})
  PetscFinalized(T) || C.ISLocalToGlobalMappingDestroy(Ref(o.p))
end

function petscview{T}(is::ISLocalToGlobalMapping{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.ISLocalToGlobalMappingView(is.p, viewer))
end


