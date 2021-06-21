# This shows how to solve the 2D incompressible Stokes equations using SNES solvers,
# using a staggered grid discretization and a Velocity-Pressure formulation
#   
#   Governing equations:
#
#              dVx/dx + dVz/dz = 0                              |   Mass balance
#   -dP/dx + dTxx/dx + dTxz/dz = 0                              |   Horizontal force balance
#   -dP/dz + dTxz/dx + dTzz/dz = rho*g                          |   Vertical force balance
#   
#   with:
#       Exx = dVx/dx, Ezz=dVz/dz, Exy= 0.5*(dVx/dz + dVz/dx)    |   Strain rate definition
#
#      | Txx |   | 2*eta   0    0  |  | Exx |
#      | Tzz | = |   0   2*eta  0  |  | Ezz |                   |   Linear viscous (isotropic) rheology
#      | Txz |   |   0     0   eta |  | Exz |        
#
# This exemple also uses the Automatic differenciation package ForwardDiff

using PETSc, MPI
using Plots
using ForwardDiff, SparseArrays, SparseDiffTools

PETSc.initialize()



function FormRes!(f_g, x_g, user_ctx)
    
    # Copy global to local vectors
    PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    PETSc.DMGlobalToLocal(user_ctx.dm, f_g, PETSc.INSERT_VALUES,  user_ctx.f_l) 

    # Retrieve arrays from the local vectors
    ArrayLocal_x     =   PETSc.DMStagVecGetArrayRead(user_ctx.dm,   user_ctx.x_l);      # array with all local x-data (solution array)
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm,       user_ctx.f_l);      # array with all local residual

    # Compute local residual 
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx)

    # Finalize local arrays
    Base.finalize(ArrayLocal_x)
    Base.finalize(ArrayLocal_f)

    # Copy local into global residual vector
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.f_l, PETSc.INSERT_VALUES, f_g) 

end

function  ForwardDiff_res(x, user_ctx)
    f   = zero(x)               # vector of zeros, of same type as x (local vector)

    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all local x-data
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, f);        # array with all local residual
    
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

    # As the residual vector f is linked with ArrayLocal_f, we don't need to
    # pass ArrayLocal_f back to f

    return f;
end

function ComputeSparsityPatternJacobian_automatic(x_l, user_ctx)
    # This computes the sparsity pattern and coloring of the jacobian automatically
    # This will works for any equation but is slow @ high resolutions 
    
    f_Residual  =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional arguments into the routine
    J_julia     =   ForwardDiff.jacobian(f_Residual,x_l*0 .+ 1);# Compute a "blank" jacobian  

    # employ sparse structure to compute jacobian - to be moved inside routines
    jac         =   sparse(J_julia);
    colors      =   matrix_colors(jac)          # the number of nonzeros per row

    return jac, colors
end

function  f(out, x, user_ctx)

    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);          # array with all local x-data (solition array)
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, out);        # array with all local residual
    
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

    return nothing
end

function FormJacobian!(cx_g, J, P, user_ctx)
    # This requires several steps:
    #
    #   1) Extract local vector from global solution (x) vector
    #   2) Compute local jacobian from the residual routine (note that
    #       this routine requires julia vectors as input)
    #   3) Store the matrix inside the global Petsc matrix J

    # Extract the local vector
    PETSc.DMGlobalToLocal(user_ctx.dm, cx_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    x               =   PETSc.unsafe_localarray(Float64, user_ctx.x_l.ptr;  write=false, read=true);

    # Check the sparsity pattern

    if isnothing(user_ctx.jac)
        error("You first have to define the sparsity pattern of the jacobian")
    else
        jac     =   user_ctx.jac;
        colors  =   user_ctx.colors;
    end
    out         =   similar(x);

    f_Res           =   ((out,x)->f(out, x, user_ctx));  # pass additional arguments into the routine

    if isnothing(user_ctx.jac_cache)
        # Allocate data required for jacobian computations
        user_ctx.jac_cache = ForwardColorJacCache(f_Res,x; colorvec=colors, sparsity = jac);
    end

    # Compute Jacobian using automatic differenciation
    forwarddiff_color_jacobian!(jac, f_Res, x, user_ctx.jac_cache);


    # Store global Jacobian inside Petsc Matrix
    ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);
    J              .=   jac[ind,ind];

    # Store Julia matrix and coloring
    user_ctx.jac    =   jac;
    user_ctx.colors =   colors;


   return jac, ind
end

# Define a struct that holds data we need in the local SNES routines below   
mutable struct Data
    # DMs and vectors
    dm
    dmCoeff
    coeff_l
    x_l
    f_l
    # physical parameters
    eta1
    eta2
    rho1
    rho2
    gz
    kappa
    # dimensions
    dx
    dz
    xlim
    zlim
    # jacobian and sparsity pattern
    jac
    jac_cache
    colors
end
user_ctx = Data(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 


function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
    # Compute the local residual. The vectors include ghost points 

    # Compute shear stresses
    Txx,Tzz,Txz = ComputeStresses!(user_ctx, ArrayLocal_x);

    # Extracting arrays from the residual and solution vectors
    P       = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_ELEMENT, 0);
    Vx      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_LEFT, 0);
    Vz      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_DOWN, 0);

    f_p     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_ELEMENT, 0);
    f_x     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_LEFT, 0);
    f_z     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_DOWN, 0);

    # Extracting physical parameters
    RhoC    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_DOWN_LEFT, 1);

    # Gettins sizes
    dx = user_ctx.dx;
    dz = user_ctx.dz;

    # Getting global indices for the local vectors
    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    # Horizontal force balance f(x)

    # Boundary conditions (Dirichlet with no normal velocity at boundary)
    f_x[ix[1]    ,:] .= Vx[ix[1]    ,:] .- 0.0;
    f_x[ix[end]+1,:] .= Vx[ix[end]+1,:] .+ 0.0;

    # -dP/dx + dTxx/dx + dTxz/dz = 0
    f_x[ix[2]:ix[end],iz] .= .-(P[ix[2]:ix[end],iz] .- P[ix[1]:ix[end-1],iz])./dx .+
                               (Txx[2:end,:]        .- Txx[1:end-1,:])       ./dx .+
                               (Txz[2:end-1,2:end]  .- Txz[2:end-1,1:end-1]) ./dz;

    # Vertical force balance f(z)

    # Boundary conditions (Dirichlet with no normal velocity at boundary)
    f_z[:,iz[1]    ] .= Vz[:,iz[1]    ] .- 0.0;
    f_z[:,iz[end]+1] .= Vz[:,iz[end]+1] .+ 0.0;

    # -dP/dz + dTxz/dx + dTzz/dz = rho*g
    f_z[ix,iz[2]:iz[end]] .= .-(P[ix,iz[2]:iz[end]] .- P[ix,iz[1]:iz[end-1]])./dz .+
                               (Tzz[:,2:end]        .- Tzz[:,1:end-1])       ./dz .+
                               (Txz[ 2:end,2:end-1] .- Txz[1:end-1,2:end-1]) ./dx .+
                               (RhoC[ix.+1,iz[2]:iz[end]] .+ RhoC[ix,iz[2]:iz[end]  ]) .* 0.5 .* user_ctx.gz;

    # Mass balance f(p)
    # dVx/dx + dVz/dz = 0
    kappa = user_ctx.kappa; # penalty term to help solver
    f_p[ix,iz] .= (Vx[ix.+1,iz].-Vx[ix,iz])./dx .+ (Vz[ix,iz.+1].-Vz[ix,iz])./dz .+ 1/kappa*P[ix,iz];

    # Cleanup
    Base.finalize(Vx);    Base.finalize(Vz);       
    Base.finalize(P);     Base.finalize(RhoC);
    Base.finalize(f_x);   Base.finalize(f_z);  
    Base.finalize(Txx);   Base.finalize(Tzz);       
    Base.finalize(Txz);                                    
end

function ComputeStresses!(user_ctx, ArrayLocal_x)
    
    # Compute strain rates
    Exx,Ezz,Exz = ComputeStrainRates!(user_ctx, ArrayLocal_x);

    # Getting Eta at the center and corner of the cells
    EtaE    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_ELEMENT, 0);
    EtaC    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_DOWN_LEFT, 0);
    
    # Getting indices for center nodes (not ghost)
    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points
    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    # Compute shear stresses
    Txx    = 2 .* EtaE[ix,iz] .* Exx;
    Tzz    = 2 .* EtaE[ix,iz] .* Ezz;
    Txz    = 2 .* EtaC[sx[1]:sn[1]+1,sx[2]:sn[2]+1] .* Exz;

    return Txx,Tzz,Txz

    Base.finalize(Exx);   Base.finalize(Ezz);       
    Base.finalize(Exz);  
    
end

function ComputeStrainRates!(user_ctx, ArrayLocal_x)

    # Getting velocity vectors
    Vx      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_LEFT, 0);
    Vz      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_DOWN, 0);
    
    # Gettind dimensions
    dx = user_ctx.dx;
    dz = user_ctx.dz;

    # Getting indices for center nodes (not ghost)
    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points
    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    # Set ghost points (free slip)
    Vx[:,iz[1]-1]   .=  Vx[:,iz[1]];                # ghost points on left and right size (here: free slip)
    Vx[:,iz[end]+1] .=  Vx[:,iz[end]];
    Vz[ix[1]-1,:]   .=  Vz[ix[1],:];
    Vz[ix[end]+1,:] .=  Vz[ix[end],:];

    # Compute deviatoric Strain rates
    DivV  = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .+ (Vz[ix,iz.+1].-Vz[ix,iz])./dz;
    Exx   = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .- 1/3 .* DivV;
    Ezz   = (Vz[ix,iz.+1].-Vz[ix,iz])./dz .- 1/3 .* DivV;
    Exz   = 0.5.*((Vx[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vx[sx[1]:sn[1].+1,sx[2].-1:sn[2]])./dz .+
                  (Vz[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vz[sx[1].-1:sn[1],sx[2]:sn[2].+1])./dx);
    
    return Exx,Ezz,Exz

end

function PopulateCoefficientData!(ctx)

    # Create coefficient local vector and array (to store eta and rho)
    user_ctx.coeff_l    =   PETSc.DMCreateLocalVector(user_ctx.dmCoeff);
    coeff_array         =   PETSc.DMStagVecGetArray(user_ctx.dmCoeff,user_ctx.coeff_l);

    # Create Coordinates dm, local vector and array
    dm_coord  = PETSc.DMGetCoordinateDM(user_ctx.dmCoeff);
    vec_coord = PETSc.DMGetCoordinatesLocal(user_ctx.dmCoeff);
    coord     = PETSc.DMStagVecGetArray(dm_coord, vec_coord);
    start,n,nExtra = PETSc.DMStagGetCorners(user_ctx.dmCoeff);

    # Get the correct entries for each of our variables in local element-wise storage
    iec = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_DOWN_LEFT, 0);  # location eta corner
    irc = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_DOWN_LEFT, 1);  # location rho corner
    iee = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_ELEMENT, 0);    # location eta element
    ixc = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_DOWN_LEFT, 0);          # location coord corner
    ixe = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_ELEMENT, 0);            # location coord element

    # Fill coefficient vectors with suitable value (determined by x-coordinate in function GetPhase())

    # Element nodes (eta)

    Coeff          = coeff_array[:,:,iee+1];
    X_coord        = coord[:,:,ixe+1];
    index1         = findall(x -> GetPhase(user_ctx,x,1),X_coord);
    Coeff[index1] .= user_ctx.eta1;

    index2         = findall(x -> GetPhase(user_ctx,x,2),X_coord);
    Coeff[index2] .= user_ctx.eta2;

    coeff_array[:,:,iee+1] .= Coeff;

    # Corner nodes (rho and eta)

    CoeffE         = coeff_array[:,:,iec+1];
    CoeffR         = coeff_array[:,:,irc+1];
    X_coord        = coord[:,:,ixc+1];
    index1         = findall(x -> GetPhase(user_ctx,x,1),X_coord);
    CoeffE[index1].= user_ctx.eta1;
    CoeffR[index1].= user_ctx.rho1;

    index2         = findall(x -> GetPhase(user_ctx,x,2),X_coord);
    CoeffE[index2].= user_ctx.eta2;
    CoeffR[index2].= user_ctx.rho2;

    coeff_array[:,:,iec+1] .= CoeffE;
    coeff_array[:,:,irc+1] .= CoeffR;


end

function GetPhase(ctx,x,n)
    # Divide domain in phases according to x-coordinate (similar to SolCx benchmark)
    if x < (ctx.xlim[2]-ctx.xlim[1])/2
        if n == 1 return true else return false end
    else
        if n == 1 return false else return true end
    end
end

function SetVecX!(user_ctx,x_g)
    # Set first guess for vector x to 0
    user_ctx.x_l .= 0; 
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g);
end


# Main Solver
nx, nz           =   32, 32;                  # number of nodes is x and z direction
user_ctx.xlim    =   [0,2];                   # x and z dimensions
user_ctx.zlim    =   [0,1];
xlim             =   user_ctx.xlim;
zlim             =   user_ctx.zlim;
user_ctx.dx      =   (xlim[2]-xlim[1])/nx;   # x-resolution
user_ctx.dz      =   (zlim[2]-zlim[1])/nz;   # z-resolution
user_ctx.eta1    =   1;                      # viscosity phase 1
user_ctx.eta2    =   2;                      # viscosity phase 2
user_ctx.rho1    =   3;                      # density phase 1
user_ctx.rho2    =   4;                      # density phase 2
user_ctx.gz      =   -1;                     # gravity magnitude
user_ctx.kappa   =   1e2;                    # incompressible penaly term

# Create Solution and coefficient DMs
user_ctx.dm           =   PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_GHOSTED,PETSc.DM_BOUNDARY_GHOSTED,nx,nz,1,1,0,1,1,PETSc.DMSTAG_STENCIL_BOX,1);  # V edge, P Element
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dm, xlim[1], xlim[2], zlim[1], zlim[2]);            # set coordinates
user_ctx.dmCoeff      =   PETSc.DMStagCreateCompatibleDMStag(user_ctx.dm,2,0,1,0);   # rho and eta on VERTEX, eta on ELEMENT
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dmCoeff, xlim[1], xlim[2], zlim[1], zlim[2]);

# Populate phases
PopulateCoefficientData!(user_ctx);

# Create solution and residual vectors
x_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm);
f_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm);
user_ctx.x_l    =   PETSc.DMCreateLocalVector(user_ctx.dm);
user_ctx.f_l    =   PETSc.DMCreateLocalVector(user_ctx.dm);

# Create Petsc matrix for Jacobian
J               =   PETSc.DMCreateMatrix(user_ctx.dm);

# Set first guess values for solution vector
SetVecX!(user_ctx,x_g);

# Compute jacobian sparcity pattern
user_ctx.jac, user_ctx.colors = ComputeSparsityPatternJacobian_automatic(user_ctx.x_l.array, user_ctx);

#Setting up SNES
S = PETSc.SNES{Float64}(MPI.COMM_SELF, 0; 
        snes_rtol=1e-12, 
        snes_monitor=true, 
        snes_max_it = 500,
        snes_monitor_true_residual=true, 
        snes_converged_reason=true);
S.user_ctx  =       user_ctx;


PETSc.setfunction!(S, FormRes!, f_g)
PETSc.setjacobian!(S, FormJacobian!, J, J)

#Solve
PETSc.solve!(x_g, S);

# Copy solution to local vector
PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l);

# Extract solution
Vx  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_LEFT,      0); 
Vz  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_DOWN,      0); 
P   =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_ELEMENT,   0); 

# Getting indices for center nodes (not ghost)
sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm);
ix      =   sx[1]:sn[1];           
iz      =   sx[2]:sn[2];

# Getting rid of ghost points
Vx  =   Vx[sx[1]:sn[1]+1,iz];
Vz  =   Vz[ix,sx[2]:sn[2]+1];
P   =    P[ix,iz];

# Compute central average velocities
Vx_cen = (Vx[2:end,:]  + Vx[1:end-1,:])/2;
Vz_cen = (Vz[:,2:end]  + Vz[:,1:end-1])/2;

#get coordinates
dm_coord  = PETSc.DMGetCoordinateDM(user_ctx.dmCoeff);
vec_coord = PETSc.DMGetCoordinatesLocal(user_ctx.dmCoeff);

XCoord_e = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_ELEMENT         ,      0); # location coord corner
ZCoord_e = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_ELEMENT         ,      1); # location coord corner
XCoord_c = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_DOWN_LEFT,      0);   # location coord element
ZCoord_c = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_DOWN_LEFT,      1);   # location coord element

XCoord_e = XCoord_e[ix,iz];
ZCoord_e = ZCoord_e[ix,iz];
XCoord_c = XCoord_c[sx[1]:sn[1]+1,sx[2]:sn[2]+1];
ZCoord_c = ZCoord_c[sx[1]:sn[1]+1,sx[2]:sn[2]+1];

xe_1D    = XCoord_e[:,1];
ze_1D    = ZCoord_e[1,:];
xc_1D    = XCoord_c[:,1];
zc_1D    = ZCoord_c[1,:];

# Plot
p1 = heatmap(xe_1D,ze_1D, P', xlabel="Width", ylabel="Depth", title="Pressure",aspect_ratio = 1)
p2 = heatmap(xe_1D,zc_1D, Vz', xlabel="Width", ylabel="Depth", title="Vz",aspect_ratio = 1)
p3 = heatmap(xc_1D,ze_1D, Vx', xlabel="Width", ylabel="Depth", title="Vx",aspect_ratio = 1)
quiver!(XCoord_e[:],ZCoord_e[:],quiver=(Vx_cen[:]*0.02,Vz_cen[:]*0.02), color=:white)
plot(p1, p2, p3, layout = (3, 1), legend = false)