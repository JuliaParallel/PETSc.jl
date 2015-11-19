# provide some most commonly used options, leave rest as low level
# common options: Orthogonilization type, KSP type, PC type
# use PC context created as part of KSP

export KSP

type KSP{T, MType}
  pksp::C.KSP{T}
  ppc::C.PC{T}
  own_pc::Bool  # is the PC owned by the ksp context
  A::Mat{T, MType}
end

comm{T}(a::KSP{T}) = comm(A)

function KSP{T, MType}(A::Mat{T, MType}, pc_mat::Mat{T, MType}=A; kws...)
  ksp_arr = Array(C.KSP{T}, 1)
  pc_arr = Array(C.PC{T}, 1)

  chk(C.KSPCreate(comm(A), ksp_arr))
  ksp = ksp_arr[1]
  C.KSPGetPC(ksp, pc_arr)
  pc = pc_arr[1]

  rank = MPI.Comm_rank(comm(A))
  chk(C.KSPSetOperators(ksp, A.p, pc_mat.p))
  withoptions(T, kws) do
    chk(C.KSPSetFromOptions(ksp))
  end

  # todo set tolerances from kws...

  # todo: finalizer

  return KSP{T, MType}(ksp, pc, true, A)
end


function KSPDestroy{T}(ksp::KSP{T})
  if !PetscFinalized(T)
    ksp.own_pc && C.PCDestroy(Ref(ksp.ppc))
    C.KSPDestroy(Ref(ksp.pksp))
  end
end

# can use options databse instead
function settolerances{T}(ksp::KSP{T}; rtol=1e-8, abstol=1e-12, dtol=1e5, maxits=size(ksp.A, 1))

  C.KSPSetTolerances(ksp, rtol, abstol, dtol, maxits)
end

# x = A \ b
function Base.A_ldiv_B!{T}(ksp::KSP{T}, b::Vec{T}, x::Vec{T})
  # if solving multiple rhs with the same matrix A,
  # the preconditioner is re-used automatically
  # if A changes, the preconditioner is recomputed

  # assemble the matrix
  AssemblyBegin(ksp.A, C.MAT_FINAL_ASSEMBLY)
  AssemblyEnd(ksp.A, C.MAT_FINAL_ASSEMBLY)

  # assemble the vector
  AssemblyBegin(b)
  AssemblyEnd(b)

  AssemblyBegin(x)
  AssemblyEnd(x)

  chk(C.KSPSolve(ksp.pksp, b.p, x.p))

  reason_arr = Array(Cint, 1)
  chk(C.KSPGetConvergedReason(ksp.pksp, reason_arr))
  if reason_arr[1] < 0
    println(STDERR, "Warning: KSP Solve did not converge")
  end

  return x
end

import Base: \
(\){T}(ksp::KSP{T}, b::Vec{T}) = A_ldiv_B!(ksp, b, similar(b, size(ksp.A, 2)))
