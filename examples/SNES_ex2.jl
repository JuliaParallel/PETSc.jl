# This implements src/snes/examples/tutorials/ex2.c from PETSc using the PETSc.jl package, using SNES
#
# This solves the equations sequentially
# 
# Newton method to solve u'' + u^{2} = f, sequentially.

using PETSc, MPI, LinearAlgebra, SparseArrays, UnicodePlots

if ~MPI.Initialized()
    MPI.Init()
end

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
comm = MPI.COMM_WORLD
"""    
    FormInitialGuess!(x)
Computes initial guess 
"""
function FormInitialGuess!(x)
    for i in eachindex(x)
        x[i] = 0.50;
    end
    return nothing
end

""" 
    F = SetInitialArrays(n)
Computes rhs forcing function 
""" 
function SetInitialArrays(n)
    h =  1.0/(n-1.0)
    F = zeros(n);
    xp = 0.0;
    for i=1:n 
        v    = 6.0*xp + (xp+1.e-12)^6.0; 
        F[i] = v;
        xp   = xp+h;
    end

    return F
end

"""
    FormResidual!(cf,cx, args...)
Computes the residual `f`, given solution vector `x`
"""
function FormResidual!(f,snes, x)
    n       =   length(x);
    xp      =   LinRange(0.0,1.0, n);
    F       =   6.0.*xp .+ (xp .+1.e-12).^6.0;      # define source term function
    
    dx      =   1.0/(n-1.0);
    f[1]    =   x[1] - 0.0;
    for i=2:n-1
        f[i] = (x[i-1] - 2.0*x[i] + x[i+1])/dx^2 + x[i]*x[i] - F[i]
    end
    f[n]    =   x[n] - 1.0;

    return 0
end

"""
    Computes the jacobian, given solution vector x
"""
function FormJacobian!(J, snes, x)
    n = length(x)
    dx  =   1.0/(n-1.0);

    # interior points (hand-coded jacobian)
    for i=2:n-1
        J[i,i-1] = 1.0/dx^2;
        J[i,i  ] = -2.0/dx^2 + 2.0*x[i];
        J[i,i+1] = 1.0/dx^2;
    end

    # boundary points
    J[1,1] = 1.0;
    J[n,n] = 1.0;
    if !isnothing(snes)
        PETSc.assemble!(J)
    end
    return 0
end


# ==========================================
# Main code 


# Compute initial solution
n   =   21;
F   =   SetInitialArrays(n);
x   =   zeros(n);

FormInitialGuess!(x);

# Compute initial jacobian using a julia structure to obtain the nonzero structure
# Note that we can also obtain this structure in a different manner
Jstruct  = zeros(n,n);
FormJacobian!(Jstruct, nothing, x);                                              # jacobian in julia form
Jsp      =   sparse(Float64.(abs.(Jstruct) .> 0))                       # sparse julia, with 1.0 in nonzero spots
PJ       =   PETSc.MatSeqAIJWithArrays(petsclib, comm, Jsp);  # transfer to PETSc format

# Setup snes
x_s = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, length(x), x)    # solution vector
res = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, length(F), F)    # residual vector
b   = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, length(F), F)    # residual vector

S = PETSc.SNES(petsclib,comm; 
        snes_rtol=1e-12, 
        snes_monitor=true,
        snes_converged_reason=false);

# Set functions for residual and jacobian computations
PETSc.setfunction!(S, FormResidual!, res)
PETSc.setjacobian!(S, FormJacobian!, PJ)

# solve
PETSc.solve!(x_s, S);

# Extract & plot solution
x_sol = x_s[:];                  # convert solution to julia format
FormResidual!(res,S, x_s)

@show norm(res[:])

# cleanup
PETSc.destroy(x_s)
PETSc.destroy(res)
PETSc.destroy(b)
PETSc.destroy(PJ)
PETSc.destroy(S)
PETSc.finalize(petsclib)

# plot solution in REPL
lineplot(LinRange(0,1,n),x_sol,xlabel="width",ylabel="solution")

