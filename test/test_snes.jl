using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

# This shows some examples of how SNES can be used. 
# there are 2 ways to use it:
#    1) you will have julia vectors (but PETSc matrices) within your user-routine
#       That is fine for sequential runs.
#    2) you will receive pointers to PETSc vectors within your user routines.
#       That is important for parallel simulations, where global residual
#       and solution vectors are passed to the user-routines, but we work with
#       local vectors


# Case 1) Julia vectors

x = ones(2)
# this is an example of a case where julia vectors are passed to the routines
function fn_julia!(fx::Vector, x::Vector, args...)
  fx[1] = x[1]^2 + x[1]*x[2] - 3
  fx[2] = x[1]*x[2] + x[2]^2 - 6
end

J = zeros(2,2)
PJ = PETSc.MatSeqDense(J)
function update_jac_julia!(x::Vector, args...)
    J[1,1] = 2x[1] + x[2]
    J[1,2] = x[1]
    J[2,1] = x[2]
    J[2,2] = x[1] + 2x[2] 
end

S_julia = PETSc.SNES{Float64}(MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none", ksp_monitor=false, ksp_converged_reason=true)
PETSc.setfunction!(S_julia, fn_julia!, PETSc.VecSeq(zeros(2)))
PETSc.setjacobian!(S_julia, update_jac_julia!, PJ, PJ)
x_julia = PETSc.VecSeq([2.0, 3.0]);
b_julia = PETSc.VecSeq([0.0, 0.0]);
PETSc.solve!(x_julia, S_julia, b_julia)
@test x_julia.array ≈ [1.0,2.0] rtol=1e-4



# Case 2) PETSc vectors
# this is an example of a case where pointers to PETSc vectors are passed to the routines
function fn!(cfx, cx, user_ctx)
  # We could do Global->Local here on cfx/cx, provided a pointer to the local
  #  vector is available in user_ctx
  x  = PETSc.unsafe_localarray(Float64, cx;  write=false)   # read array
  fx = PETSc.unsafe_localarray(Float64, cfx; read=false)    # write array
 
  fx[1] = x[1]^2 + x[1]*x[2] - 3
  fx[2] = x[1]*x[2] + x[2]^2 - 6
  
  Base.finalize(fx)
  Base.finalize(x)
end


function update_jac!(cx, J1::PETSc.AbstractMat{Float64}, args...)
    x  = PETSc.unsafe_localarray(Float64, cx;  write=false)
    @show J1, typeof(J1)
    J1[1,1] = 2x[1] + x[2]
    J1[1,2] = x[1]
    J1[2,1] = x[2]
    J1[2,2] = x[1] + 2x[2] 

    Base.finalize(x)
    PETSc.assemble(J1)          
end

# structure with which we can pass data to the user-routines above
mutable struct Data
    vec
    julia
end

julia_vec = 0;  # we want pointers to local vectors 
S = PETSc.SNES{Float64}(MPI.COMM_SELF, julia_vec;  
                ksp_rtol=1e-4, 
                pc_type="none", 
                ksp_monitor=true, 
                ksp_converged_reason=true)


data        = Data([100;2], 1)
S.user_ctx  = data;      # we can pack anything we need in this struct


PETSc.setfunction!(S, fn!, PETSc.VecSeq(zeros(2)))
PJ = PETSc.MatSeqDense(zeros(2,2))
PETSc.setjacobian!(S, update_jac!, PJ, PJ)
x  = PETSc.VecSeq([2.0, 3.0]);
b  = PETSc.VecSeq([0.0, 0.0]);
PETSc.solve!(x, S, b)
@test x.array ≈ [1.0,2.0] rtol=1e-4


# cleanup
#finalize(x)
#finalize(J);
#PETSc.destroy(x0);
#PETSc.destroy(x);
#PETSc.destroy(b);
#PETSc.destroy(PJ);
#PETSc.destroy(S);




#include("test_dmstag.jl")



#PETSc.finalize()