# provide some most commonly used options, leave rest as low level
# common options: Orthogonilization type, KSP type, PC type
# use PC context created as part of KSP

type KSP{T, MType}
  pksp::C.KSP{T}
  ppc::C.PC{T}
  own_pc::Bool  # is the PC owned by the ksp context
  A::Mat{T, MType}
end

function KSP{T, MType}(A::Mat{T, MType}, pc_mat::Mat{T, MType}=lhs_mat)
  ksp_arr = Array(C.KSP{T}, 1)
  pc_arr = Array(C.PC{T}, 1)

  C.KSPCreate(lhs_mat.comm, ksp_arr)
  ksp = ksp_arr[1]
  C.KSPGetPC(ksp, pc_arr)
  pc = pc_arr[1]

  KSPSetOperators(ksp, A.p, pc_mat.p)

  return KSP{T, MType}(ksp, pc, true, A)
end


function KSPDestroy(ksp::KSP)

  tmp = Array(PetscBool, 1)
  C.PetscFinalized(eltype(mat), tmp)
   
  if tmp[1] == 0  # if petsc has not been finalized yet
    if !ksp.own_pc
      C.PCDestroy([ksp.ppc])
    end

    C.KSPDestroy([ksp.pksp])
  end

   # if Petsc has been finalized, let the OS deallocate the memory
end



function settolerances{T}(ksp::KSP{T}; rtol=1e-8, abstol=1e-12, dtol=1e5, maxits=size(ksp.A, 1))

  C.KSPSetTolerances(ksp, rtol, abstol, dtol, maxits)
end

function solve(ksp::KSP{T}, b::Vec{T}, x::Vec{T})
# perform the solve
# users should specify all the options they want to use
# before calling this function
# if solving multiple rhs with the same matrix A,
# the preconditioner is resued automatically
# if A changes, the preconditioner is recomputed

  KSPSetFromOption(ksp.pksp)
#  KSPSetUp(ksp)   # this is called by KSPSolve if needed
                   # decreases logging accurace of setup operations
  KSPSolve(ksp.pksp, b.p, x.p)

  reason_arr = Array(PetscInt, 1)
  C.KSPGetConvergedReason(ksp, reason) 
  reason = reason_arr[1]



  if reason < 0
    println(STDERR, "Warning: KSP Solve did not converge")
  end

end

function solve(ksp::KSP{T}, b::Vec{T})
  x = similar(b)
  solve(ksp, b, x)
end




  
  
