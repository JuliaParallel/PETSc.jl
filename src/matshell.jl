"""
    MatShell{T}(obj, m, n)

Create a `m×n` PETSc shell matrix object wrapping `obj`.

If `obj` is a `Function`, then the multiply action `obj(y,x)`; otherwise it calls `mul!(y, obj, x)`.
This can be changed by `PETSc._mul!`.

"""
mutable struct MatShell{T,A} <: AbstractMat{T}
  ptr::CMat
  comm::MPI.Comm
  obj::A
end


struct MatOp{T,Op}
end


function _mul!(y,mat::MatShell{T,F},x) where {T, F<:Function}
  mat.obj(y, x)
end

function _mul!(y,mat::MatShell{T},x) where {T}
  LinearAlgebra.mul!(y, mat.obj, x)
end

MatShell{T}(obj, m, n) where {T} = MatShell{T}(obj, MPI.COMM_SELF, m, n, m, n)


@for_libpetsc begin
  function MatShell{$PetscScalar}(obj::A, comm::MPI.Comm, m, n, M, N) where {A}
    mat = MatShell{$PetscScalar,A}(C_NULL, comm, obj)
    # we use the MatShell object itsel
    ctx = pointer_from_objref(mat)
    @chk ccall((:MatCreateShell, $libpetsc), PetscErrorCode,
      (MPI.MPI_Comm,$PetscInt,$PetscInt,$PetscInt,$PetscInt,Ptr{Cvoid},Ptr{CMat}),
      comm, m, n, M, N, ctx, mat)

    mulptr = @cfunction(MatOp{$PetscScalar, MATOP_MULT}(), $PetscInt, (CMat, CVec, CVec))
    @chk ccall((:MatShellSetOperation, $libpetsc), PetscErrorCode, (CMat, MatOperation, Ptr{Cvoid}), mat, MATOP_MULT, mulptr)
    return mat
  end

  function (::MatOp{$PetscScalar, MATOP_MULT})(M::CMat,cx::CVec,cy::CVec)::$PetscInt
    #try
      r_ctx = Ref{Ptr{Cvoid}}()
      @chk ccall((:MatShellGetContext, $libpetsc), PetscErrorCode, (CMat, Ptr{Ptr{Cvoid}}), M, r_ctx)
      ptr = r_ctx[]
      mat = unsafe_pointer_to_objref(ptr)

      x = unsafe_localarray($PetscScalar, cx)
      y = unsafe_localarray($PetscScalar, cy)

      _mul!(y,mat,x)

      Base.finalize(y)
      Base.finalize(x)
      return $PetscInt(0)
    # catch e
    #   if e isa PetscErrorCode
    #     return e.code
    #   else
    #     return 
    #   end
    # end
  end
    
end
