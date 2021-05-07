# This implements src/snes/examples/tutorials/ex2.c from PETSc using the PETSc.jl package, using SNES
#
# This is the same as SNES_ex2b.j, except that we show how automatic differentiation can be used to
# compute the jacobian. 
#
# Newton method to solve u'' + u^{2} = f, sequentially.

using PETSc, MPI, LinearAlgebra, SparseArrays, Plots, ForwardDiff

if ~MPI.Initialized()
    MPI.Init()
end

PETSc.initialize()

```
    Computes initial guess 
```
function FormInitialGuess!(dm,x)

    x_Local = PETSc.DMCreateLocalVector(dm);
    x_array = PETSc.DMStagVecGetArray(dm,x_Local);
    x_array[1:2,1:2] .= 1;
    PETSc.DMLocalToGlobal(dm, x_Local, PETSc.INSERT_VALUES, x)
end

```
    Computes the residual f, given solution vector x
```
function FormResidual!(f,x)
    n       =   length(x);
    xp      =   LinRange(0.0,1.0, n);
    F       =   6.0.*xp .+ (xp .+1.e-12).^6.0;      # define source term function
    
    dx      =   1.0/(n-1.0);
    f[1]    =   x[1] - 0.0;
    for i=2:n-1
        f[i] = (x[i-1] - 2.0*x[i] + x[i+1])/dx^2 + x[i]*x[i] - F[i]
    end
    f[n]    =   x[n] - 1.0;

end

```
    Wrapper which makes it easier to compute the jacobian using automatic differntiation
```
function  ForwardDiff_routine(x)

    f   = zero(x)               # vector of zeros, of same type as e
    FormResidual!(f,x);

    return f;
end

```
    This copies a julia sparse matrix to PETSc MatSeqAIJ format
```
function Mat_JuliaToPETSc!(J::PETSc.MatSeqAIJ, J_julia::SparseMatrixCSC)

    for i = 1:size(J_julia,1)
        col = J_julia[i,:];
        row = ones(Int32,length(col.nzind))*i;
        for j=1:length(col.nzind)
            J[i, col.nzind[j]] = col.nzval[j];
        end
    end
    PETSc.assemble(J);  # finalize assembly

end

```
    Computes the jacobian, given solution vector x
```
function FormJacobian!(x, args...)

    J        =   args[1];        # preconditioner = args[2], in case we want it to be different from J

    # Use AD to compute jacobian; by transferring x into sparse, the output will be sparse
    J_julia  =   ForwardDiff.jacobian(ForwardDiff_routine,sparse(x));

    if typeof(J) <: PETSc.AbstractMat
        Mat_JuliaToPETSc!(J, J_julia);  # transfer julia sparse matrix 2 petsc
    else
        J .= J_julia;
    end
end


# ==========================================
# Main code 


# Compute initial solution
nx   =   5;
x0   =   0;
xend =   1;

# create dmstag for solution and setup
dm  = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,nx,1,1,PETSc.DMSTAG_STENCIL_BOX,1);
# creat uniform coordinates
PETSc.DMStagSetUniformCoordinates(dm, x0, xend);
#determine boundary type
bnd =   PETSc.DMStagGetBoundaryTypes(dm);

a = 1.0; b = 2.0; c = 1.0;
if bnd == PETSc.DM_BOUNDARY_PERIODIC
    b = a;
    c = 0.0;
end


x       = PETSc.DMCreateGlobalVector(dm);
x_Local = PETSc.DMCreateLocalVector(dm);
x_array = PETSc.DMStagVecGetArray(dm,x_Local);
#need coordinate stuff
start,n,nExtra = PETSc.DMStagGetCorners(dm);
iu = PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_LEFT, 0);
ip = PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_ELEMENT, 0);

dmForcing = PETSc.DMStagCreateCompatibleDMStag(dm,1,0)
f         = PETSc.DMCreateGlobalVector(dmForcing);
fLocal    = PETSc.DMCreateLocalVector(dmForcing);
f        .= c;
fLocal   .= c;

A   = PETSc.DMCreateMatrix(dm);
rhs = PETSc.DMCreateGlobalVector(dm);
#FormInitialGuess!(dm,x);

# Compute initial jacobian using a julia structure to obtain the nonzero structure
# Note that we can also obtain this structure in a different manner
#Jstruct  = zeros(n,n);
#FormJacobian!(x, Jstruct);                              # jacobian in julia form
#Jsp      =   sparse(Float64.(abs.(Jstruct) .> 0))       # sparse julia, with 1.0 in nonzero spots
#PJ       =   PETSc.MatSeqAIJ(Jsp);                      # transfer to PETSc (initialize matrix with correct nonzero pattern)

# Setup SNES
#x_s = PETSc.VecSeq(x);                  # solution vector
#res = PETSc.VecSeq(zeros(size(x)));     # residual vector

#S = PETSc.SNES{Float64}(MPI.COMM_SELF; 
#        snes_rtol=1e-12, 
#        snes_monitor=true, 
#        snes_converged_reason=true);
#PETSc.setfunction!(S, FormResidual!, res)
#PETSc.setjacobian!(S, FormJacobian!, PJ, PJ)

# solve
#PETSc.solve!(x_s, S);

# Extract & plot solution
#x_sol = x_s.array;                  # convert solution to julia format
#FormResidual!(res.array,x_sol)      # just for checking, compute residual
#@show norm(res.array)

#plot(LinRange(0,1,n),x_sol,xlabel="width",ylabel="solution")

#PETSc.finalize()
