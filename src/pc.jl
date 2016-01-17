export PC, petscview

# preconditioner context
type PC{T}
  p::C.PC{T}
  function PC(p::C.PC{T})
    o = new(p)
    finalizer(o, PetscDestroy)
    return o
  end
end

comm{T}(a::PC{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function PetscDestroy{T}(o::PC{T})
  PetscFinalized(T) || C.PCDestroy(Ref(o.p))
end

function petscview{T}(o::PC{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.PCView(o.p, viewer))
end

"""
    PC(A::Mat, PA=A, kws...)

Create a preconditioner (PC) context, given the matrix `A` of the
linear system to be solved, and optionally a different
matrix `PA` from which the preconditioner is constructed.

The remaining keywords specify zero or more additional options:
* `pc_type="a"`: use preconditioning algorithm `a`
* `pc_use_amat=true`: use Amat (instead of Pmat) to define preconditioner in nested inner solves
* ... additional options that depend on the preconditioner type ...
"""
function PC{T}(::Type{T}; comm::MPI.Comm=MPI.COMM_SELF, kws...)
  pc_c = Ref{C.PC{T}}()
  chk(C.PCCreate(comm, pc_c))
  withoptions(T, kws) do
    chk(C.PCSetFromOptions(pc_c[]))
  end
  return PC{T}(pc_c[])
end

function PC{T}(A::Mat{T}, PA::Mat{T}=A; kws...)
  pc_c = Ref{C.PC{T}}()
  chk(C.PCCreate(comm(A), pc_c))
  pc = pc_c[]
  chk(C.PCSetOperators(pc, A.p, PA.p))
  withoptions(T, kws) do
    chk(C.PCSetFromOptions(pc))
  end
  return PC{T}(pc)
end
