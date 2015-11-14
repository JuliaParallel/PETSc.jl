# index sets and vector scatters
export IS, VecScatter, scatter!

type IS{T}
  p::C.IS{T}

  # use default inner constructor
end

function IS(T::DataType, idx::AbstractArray{PetscInt}; comm=MPI.COMM_SELF)
  is_c = Ref{C.IS{T}}()
#  println("is_c = ", is_c)
#  println("typeof(is_c) = ", typeof(is_c))
  chk(C.ISCreateGeneral(comm, length(idx), idx, C.PETSC_COPY_VALUES, is_c))
  return IS{T}(is_c[])
end

type VecScatter{T}
  p::C.VecScatter{T}

  # use default innter constructor
end

function VecScatter{T}(x::Vec{T}, ix::IS{T}, y::Vec{T}, iy::IS{T})
  scatter_c = Ref{C.VecScatter{T}}()
  chk(C.VecScatterCreate(x.p, ix.p, y.p, iy.p, scatter_c))
  return VecScatter{T}(scatter_c[])
end

function scatter!{T}(scatter::VecScatter{T}, x::Vec{T}, y::Vec{T}; imode=C.INSERT_VALUES, smode=C.SCATTER_FORWARD)

  chk(C.VecScatterBegin(scatter.p, x.p, y.p, imode, smode))
  chk(C.VecScatterEnd(scatter.p, x.p, y.p, imode, smode))
end
