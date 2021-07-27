# This implements src/snes/examples/tutorials/ex2.c from PETSc using the PETSc.jl package, using SNES
#
# Note that yhe PETSc.jl package does currently not support MPI-parallel cases (even when PETSC_jll does support it)
#
# Newton method to solve u'' + u^{2} = f, sequentially.

using PETSc, MPI, LinearAlgebra, SparseArrays, UnicodePlots

if ~MPI.Initialized()
    MPI.Init()
end

PETSc.initialize()

```
    Computes initial guess 
```
function FormInitialGuess!(x)
    for i=1:length(x)
        x[i] = 0.50;
    end
end

```
    Computes rhs forcing function 
```
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

```
    Computes the residual f, given solution vector x
```
function FormResidual!(cf,cx, args...)
    if typeof(cx) <: Vector{Float64}
        x = cx;
    else
        x   =   PETSc.unsafe_localarray(Float64,cx)
    end
    if typeof(cf) <: Vector{Float64}
        f = cf;
    else
        f   =   PETSc.unsafe_localarray(Float64,cf)
    end
    n   =   length(x);
    print("Length of x is ", n, "\n");
    print("f is ", f,"\n");
    xp  =   LinRange(0.0,1.0, n);
    F   .=  6.0.*xp .+ (xp .+1.e-12).^6.0;      # define source term function
    
    dx      =   1.0/(n-1.0);
    f[1]    = x[1] - 0.0;
    for i=2:n-1
        f[i] = (x[i-1] - 2.0*x[i] + x[i+1])/dx^2 + x[i]*x[i] - F[i]
    end
    f[n]    = x[n] - 1.0;
    Base.finalize(x)
    Base.finalize(f)

end

```
    Computes the jacobian, given solution vector x
```
function FormJacobian!(x, args...)

    #x = PETSc.unsafe_localarray(Float64,cx)
    J   =   args[1];        # preconditioner = args[2], in case we want it to be different from J
    n   =   length(x);
    dx  =   1.0/(n-1.0);
    
    # interior points
    for i=2:n-1
        J[i,i-1] = 1.0/dx^2;
        J[i,i  ] = -2.0/dx^2 + 2.0*x[i];
        J[i,i+1] = 1.0/dx^2;
    end

    # boundary points
    J[1,1] = 1.0;
    J[n,n] = 1.0;
  
    if typeof(J) <: PETSc.AbstractMat
        PETSc.assemble(J);  # finalize assembly
    end
    Base.finalize(x)
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
FormJacobian!(x, Jstruct);                              # jacobian in julia form
Jsp      =   sparse(Float64.(abs.(Jstruct) .> 0))       # sparse julia, with 1.0 in nonzero spots
PJ       =   PETSc.MatSeqAIJ(Jsp);                      # transfer to 

# Setup snes
x_s = PETSc.VecSeq(x);                  # solution vector
res = PETSc.VecSeq(zeros(size(x)));     # residual vector

S = PETSc.SNES{Float64}(MPI.COMM_SELF; 
        snes_rtol=1e-12, 
        snes_monitor=nothing,
        snes_converged_reason=nothing);
PETSc.setfunction!(S, FormResidual!, res)
PETSc.setjacobian!(S, FormJacobian!, PJ, PJ)

# solve
PETSc.solve!(x_s, S);

# Extract & plot solution
x_sol = x_s.array;                  # convert solution to julia format
FormResidual!(res.array,x_sol)      # just for checking, compute residual
@show norm(res.array)


lineplot(LinRange(0,1,n),x_sol,xlabel="width",ylabel="solution")

#PETSc.finalize()
