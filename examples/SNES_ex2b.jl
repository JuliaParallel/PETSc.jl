# This implements src/snes/examples/tutorials/ex2.c from PETSc using the PETSc.jl package, using SNES
#
# This is the same as SNES_ex2.jl, except that we show how automatic differentiation can be used to
# compute the jacobian. That implies that the user does not have to provide a hand-coded jacobian.
# Note that the way we compute the jacobian here is not very efficient as we do not use the sparsity structure of the 
# matrix; see the other examples for faster implementions. 
#
# Newton method to solve u'' + u^{2} = f, sequentially.

using PETSc, MPI, LinearAlgebra, SparseArrays, UnicodePlots, ForwardDiff
using Test

if ~MPI.Initialized()
    MPI.Init()
end

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
comm = MPI.COMM_SELF


"""
    Computes initial guess 
"""
function FormInitialGuess!(x)
    for i=1:length(x)
        x[i] = 0.50;
    end
end

"""
    Computes the residual f, given solution vector x
"""
function FormResidual!(cf,cx, args...)
    if typeof(cx) <: Ptr{Nothing}
        # When this routine is called from PETSc, cx is a pointer to a global vector
        # That's why we have to transfer it first to 
        x   =   PETSc.unsafe_localarray(PETSc.scalartype(petsclib),cx, write=false)
    else
        x   = cx;
    end
    if typeof(cf) <: Ptr{Nothing}
        f   =   PETSc.unsafe_localarray(PETSc.scalartype(petsclib),cf, write=true)
    else
        f   = cf;
    end
    n       =   length(x);
    xp      =   LinRange(0.0,1.0, n);
    F       =   6.0.*xp .+ (xp .+1.e-12).^6.0;      # define source term function
    
    dx      =   1.0/(n-1.0);
    f[1]    =   x[1] - 0.0;
    for i=2:n-1
        f[i] = (x[i-1] - 2.0*x[i] + x[i+1])/dx^2 + x[i]*x[i] - F[i]
    end
    f[n]    =   x[n] - 1.0;
    Base.finalize(x)
    Base.finalize(f)

end

"""
    Wrapper which makes it easier to compute the jacobian using automatic differntiation
"""
function  ForwardDiff_res(x)

    f   = zero(x)               # vector of zeros, of same type as x
    FormResidual!(f,x);

    return f;
end

"""
    Computes the jacobian, given solution vector x
"""
function FormJacobian!(cx, args...)

    if typeof(cx) <: Ptr{Nothing}
        x   =   PETSc.unsafe_localarray(PETSc.scalartype(petsclib),cx, write=false)
    else
        x   =   cx;
    end
    J        =  args[1];        # preconditioner = args[2], in case we want it to be different from J

    # Use AD to compute jacobian; by transferring x into sparse, the output will be sparse
    J_julia  =  ForwardDiff.jacobian(ForwardDiff_res,sparse(x));

    J       .=  J_julia;   # copy to petsc format
end


# ==========================================
# Main code 


# Compute initial solution
n   =   101;
x   =   zeros(n);

FormInitialGuess!(x);

# Compute initial jacobian using a julia structure to obtain the nonzero structure
# Note that we can also obtain this structure in a different manner
Jstruct  = zeros(n,n);
FormJacobian!(x, Jstruct);                              # jacobian in julia form
Jsp      =   sparse(Float64.(abs.(Jstruct) .> 0))       # sparse julia, with 1.0 in nonzero spots
PJ       =   PETSc.MatSeqAIJWithArrays(petsclib, comm, Jsp);  # transfer to PETSc format

# Setup SNES
x_s = PETSc.VecSeq(petsclib, comm, x);                  # solution vector
res = PETSc.VecSeq(petsclib, comm, zeros(size(x)));     # residual vector

S = PETSc.SNES(petsclib,MPI.COMM_SELF; 
        snes_rtol=1e-12, 
        snes_monitor=nothing,
        snes_converged_reason=nothing);
PETSc.setfunction!(S, FormResidual!, res)
PETSc.setjacobian!(S, FormJacobian!, PJ, PJ)

# solve
PETSc.solve!(x_s, S);

# Extract & plot solution
x_sol = x_s[:];                  # convert solution to julia format
FormResidual!(res[:],x_sol)      # just for checking, compute residual

@show norm(res[:])

PETSc.finalize(petsclib)

# plot solution in REPL
lineplot(LinRange(0,1,n),x_sol,xlabel="width",ylabel="solution")
