# index sets and vector scatters

###########################################################################
export IS # index sets
# Note: we expose a 1-base Julian index interface, but internally
# PETSc's indices are 0-based.

type IS{T}
  p::C.IS{T}
  function IS(p::C.IS{T})
    o = new(p)
    finalizer(o, ISDestroy)
    return o
  end
end

function ISDestroy{T}(o::IS{T})
  PetscFinalized(T) || C.ISDestroy(Ref(o.p))
end

# internal constructor, takes array of zero-based indices:
function IS_(T::DataType, idx::Array{PetscInt}; comm=MPI.COMM_SELF)
  is_c = Ref{C.IS{T}}()
  chk(C.ISCreateGeneral(comm, length(idx), idx, C.PETSC_COPY_VALUES, is_c))
  return IS{T}(is_c[])
end

IS{I<:Integer}(T::DataType, idx::AbstractArray{I}; comm=MPI.COMM_SELF) =
  IS_(T, PetscInt[i-1 for i in idx]; comm=comm)

function IS{I<:Integer}(T::DataType, idx::Range{I}; comm=MPI.COMM_SELF)
  is_c = Ref{C.IS{T}}()
  chk(C.ISCreateStride(comm, length(idx), start(idx)-1, step(idx), is_c))
  return IS{T}(is_c[])
end

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

###########################################################################
export VecScatter

# describes a scatter operation context (input/output index sets etc.)
type VecScatter{T}
  p::C.VecScatter{T}
  function VecScatter(p::C.VecScatter{T})
    o = new(p)
    finalizer(o, VecScatterDestroy)
    return o
  end
end

function VecScatterDestroy{T}(o::IS{T})
  PetscFinalized(T) || C.VecScatterDestroy(Ref(o.p))
end

function VecScatter{T}(x::Vec{T}, ix::IS{T}, y::Vec{T}, iy::IS{T})
  scatter_c = Ref{C.VecScatter{T}}()
  chk(C.VecScatterCreate(x.p, ix.p, y.p, iy.p, scatter_c))
  return VecScatter{T}(scatter_c[])
end

function Base.copy{T}(i::VecScatter{T})
  vs_c = Ref{C.VecScatter{T}}()
  chk(C.VecScatterCopy(i.p, vs_c))
  return VecScatter{T}(vs_c[])
end

###########################################################################
export scatter!

function scatter!{T}(scatter::VecScatter{T}, x::Vec{T}, y::Vec{T}; imode=C.INSERT_VALUES, smode=C.SCATTER_FORWARD)
  chk(C.VecScatterBegin(scatter.p, x.p, y.p, imode, smode))
  yield() # do async computations while messages are in transit
  chk(C.VecScatterEnd(scatter.p, x.p, y.p, imode, smode))
  return y
end

function scatter!{T,I1,I2}(x::Vec{T}, ix::AbstractVector{I1},
                           y::Vec{T}, iy::AbstractVector{I2};
                          imode=C.INSERT_VALUES, smode=C.SCATTER_FORWARD)
  scatter = VecScatter(x, IS(T, ix, comm=x.comm),
                       y, IS(T, iy, comm=y.comm))
  scatter!(scatter, x, y; imode=imode, smode=smode)
end
