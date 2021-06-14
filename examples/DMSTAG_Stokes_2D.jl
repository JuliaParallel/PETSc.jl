# This is an example of a 1D viscoelastic porosity wave as described in 
# Vasyliev et al. Geophysical Research Letters (25), 17. p. 3239-3242
# 
# It simulates how a pulse of magma migrates upwards in the Earth, which 
# can be described by a set of coupled nonlinear PDE's
#
# This example only requires the specification of a residual routine; automatic
# differentiation is used to generate the jacobian.  

using PETSc, MPI
using Plots
using ForwardDiff, SparseArrays, SparseDiffTools

PETSc.initialize()



function FormRes!(f_g, x_g, user_ctx)

    # Note that in PETSc, cx and cfx are pointers to global vectors. 
    
    # Copy global to local vectors
    PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    PETSc.DMGlobalToLocal(user_ctx.dm, f_g, PETSc.INSERT_VALUES,  user_ctx.f_l) 

    # Retrieve arrays from the local vectors
    ArrayLocal_x     =   PETSc.DMStagVecGetArrayRead(user_ctx.dm,   user_ctx.x_l);      # array with all local x-data
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm,       user_ctx.f_l);      # array with all local residual
    
    #print("going through FormRes \n")

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

    #@show typeof(x) typeof(f)
    
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

    # As the residual vector f is linked with ArrayLocal_f, we don't need to
    # pass ArrayLocal_f back to f

    return f;
end

function ComputeSparsityPatternJacobian_automatic(x_l, user_ctx)
    # This computes the sparsity pattern and coloring of the jacobian automatically
    # This will works for any equation but is slow @ high resolutions 
    
    f_Residual  =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional arguments into the routine
    J_julia     =   ForwardDiff.jacobian(f_Residual,x_l*0 .+ 1);  

    # employ sparse structure to compute jacobian - to be moved inside routines
    jac         =   sparse(J_julia);
    colors      =   matrix_colors(jac)          # the number of nonzeros per row

    return jac, colors
end

function  f(out, x, user_ctx)

    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all local x-data
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, out);        # array with all local residual
    
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

    return nothing
end

#function  f(out, x, user_ctx)
#
#    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all local x-data
#    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, out);        # array with all local residual
#    
#    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);
#
#    return nothing
#end

function FormJacobian!(cx_g, J, P, user_ctx)
    # This requires several steps:
    #
    #   1) Extract local vector from global solution (x) vector
    #   2) Compute local jacobian from the residual routine (note that
    #       this routine requires julia vectors as input)

    # Extract the local vector
    PETSc.DMGlobalToLocal(user_ctx.dm, cx_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    x               =   PETSc.unsafe_localarray(Float64, user_ctx.x_l.ptr;  write=false, read=true);

    #@show typeof(x) typeof(user_ctx.x_l)

    if isnothing(user_ctx.jac)
        error("You first have to define the sparsity pattern of the jacobian")
    else
        jac     =   user_ctx.jac;
        colors  =   user_ctx.colors;
    end
    out         =   similar(x);

    f_Res           =   ((out,x)->f(out, x, user_ctx));

    if isnothing(user_ctx.jac_cache)
        # Allocate data required for jacobian computations
        user_ctx.jac_cache = ForwardColorJacCache(f_Res,x; colorvec=colors, sparsity = jac);
    end

    forwarddiff_color_jacobian!(jac, f_Res, x, user_ctx.jac_cache);

    
    #n                =   size(P,1)
    #J               .=   sparse(J_julia[1:n,1:n]);
    ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);    # extr

    # @show size(J) size(J_julia[ind,ind])
    J              .=   jac[ind,ind];

    user_ctx.jac    =   jac;
    user_ctx.colors =   colors;
    #J               =   J_julia[ind,ind];


   return jac, ind
end

# Define a struct that holds data we need in the local SNES routines below   
mutable struct Data
    dm
    dmCoeff
    dmEps
    dmTau
    coeff_l
    Eps_l
    Tau_l
    x_l
    f_l
    eta1
    eta2
    rho1
    rho2
    gz
    dt
    dx
    dz
    xlim
    zlim
    jac
    jac_cache
    colors
end
user_ctx = Data(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 


function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
    # Compute the local residual. The vectors include ghost points 

    #print("going through ComputeLocalResidual \n")

    Txx,Tzz,Txz = ComputeStresses!(user_ctx, ArrayLocal_x);

    P       = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_ELEMENT, 0);
    Vx      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_LEFT, 0);
    Vz      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_DOWN, 0);

    f_p     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_ELEMENT, 0);
    f_x     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_LEFT, 0);
    f_z     = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,     PETSc.DMSTAG_DOWN, 0);

    RhoC    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_DOWN_LEFT, 1);

    #Txx,Tzz,Txz = ComputeStresses!(user_ctx, ArrayLocal_x);

    #print("Exiting computeStresses \n")

    dx = user_ctx.dx;
    dz = user_ctx.dz;

    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    # Force balance f(x)

    f_x[ix[1]    ,:] .= Vx[ix[1]    ,:] .- 0.0;    #Dirichlet
    f_x[ix[end]+1,:] .= Vx[ix[end]+1,:] .+ 0.0;

    f_x[ix[2]:ix[end],iz] .= .-(P[ix[2]:ix[end],iz] .- P[ix[1]:ix[end-1],iz])./dx .+
                               (Txx[2:end,:]        .- Txx[1:end-1,:])       ./dx .+
                               (Txz[2:end-1,2:end]  .- Txz[2:end-1,1:end-1]) ./dz;

    # Force balance f(z)

    f_z[:,iz[1]    ] .= Vz[:,iz[1]    ] .- 0.0;    #Dirichlet
    f_z[:,iz[end]+1] .= Vz[:,iz[end]+1] .+ 0.0;

    f_z[ix,iz[2]:iz[end]] .= .-(P[ix,iz[2]:iz[end]] .- P[ix,iz[1]:iz[end-1]])./dz .+
                               (Tzz[:,2:end]        .- Tzz[:,1:end-1])       ./dz .+
                               (Txz[ 2:end,2:end-1] .- Txz[1:end-1,2:end-1]) ./dx .+
                               (RhoC[ix.+1,iz[2]:iz[end]] .+ RhoC[ix,iz[2]:iz[end]  ]) .* 0.5 .* user_ctx.gz;

    # Mass balance f(p)

    f_p[ix,iz] .= (Vx[ix.+1,iz].-Vx[ix,iz])./dx .+ (Vz[ix,iz.+1].-Vz[ix,iz])./dz;

    #println(f_x,"\n\n",f_z,"\n\n",f_p,"\n\n")

    #print("f_x = ",f_x," \n f_z = ",f_z," \n f_p = ",f_p,"\n ArrayLocal_f = ",ArrayLocal_f)

    # Cleanup
    #Base.finalize(Vx);    Base.finalize(Vz);       
    #Base.finalize(P);     Base.finalize(RhoC);
    #Base.finalize(Txx);   Base.finalize(Tzz);       
    #Base.finalize(Txz);   Base.finalize(RhoC);                                         
end

function ComputeStresses!(user_ctx, ArrayLocal_x)

    #print("going through ComputeStresses \n")
    
    Exx,Ezz,Exz = ComputeStrainRates!(user_ctx, ArrayLocal_x);

    #print("Exiting ComputeStrainrates \n")

    EtaE    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_ELEMENT, 0);
    EtaC    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_DOWN_LEFT, 0);
    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    Txx    = 2 .* EtaE[ix,iz] .* Exx;
    Tzz    = 2 .* EtaE[ix,iz] .* Ezz;
    Txz    = 2 .* EtaC[sx[1]:sn[1]+1,sx[2]:sn[2]+1] .* Exz;

    #println(Txx,"\n\n",Tzz,"\n\n",Txz,"\n\n")

    #println(EtaE,"\n\n",EtaC,"\n\n")

    return Txx,Tzz,Txz

    Base.finalize(Exx);   Base.finalize(Ezz);       
    Base.finalize(Exz);  
    
end

function ComputeStrainRates!(user_ctx, ArrayLocal_x)

    #print("going through ComputeStrainRates \n")

    dx = user_ctx.dx;
    dz = user_ctx.dz;

    Vx      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_LEFT, 0);
    Vz      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_DOWN, 0);
    sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

    ix      =   sx[1]:sn[1];           
    iz      =   sx[2]:sn[2];

    Vx[:,iz[1]-1]   .=  Vx[:,iz[1]];                # ghost points on left and right size (here: free slip)
    Vx[:,iz[end]+1] .=  Vx[:,iz[end]];
    Vz[ix[1]-1,:]   .=  Vz[ix[1],:];
    Vz[ix[end]+1,:] .=  Vz[ix[end],:];
    #@show typeof(Exx) typeof(Vx) typeof(ArrayLocal_x)
    DivV                              = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .+ (Vz[ix,iz.+1].-Vz[ix,iz])./dz;
    
    #Exx[ix,iz]                       .= Ezz[ix,iz]; #(Vx[ix.+1,iz].-Vx[ix,iz])./dx .- 1/3 .* DivV;
    #Ezz[ix,iz]                       .= 1; #(Vz[ix,iz.+1].-Vz[ix,iz])./dz .- 1/3 .* DivV;
    #Exz[sx[1]:sn[1]+1,sx[2]:sn[2]+1] .= 0.5.*((Vx[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vx[sx[1]:sn[1].+1,sx[2].-1:sn[2]])./dz .+
    #                                          (Vz[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vz[sx[1].-1:sn[1],sx[2]:sn[2].+1])./dx);

    Exx = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .- 1/3 .* DivV;
    Ezz = (Vz[ix,iz.+1].-Vz[ix,iz])./dz .- 1/3 .* DivV;
    Exz = 0.5.*((Vx[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vx[sx[1]:sn[1].+1,sx[2].-1:sn[2]])./dz .+
                (Vz[sx[1]:sn[1].+1,sx[2]:sn[2].+1].-Vz[sx[1].-1:sn[1],sx[2]:sn[2].+1])./dx);

    #print("Exx = ",Exx," \n Ezz = ",Ezz,"Exz = ",Exz)
    #println(Exx,"\n\n",Ezz,"\n\n",Exz)
    
    return Exx,Ezz,Exz

end


function SetInitialPerturbations(user_ctx, x_g)
    # Computes the initial perturbations as in the paper

    # Retrieve coordinates from DMStag
    DMcoord     =   PETSc.DMGetCoordinateDM(user_ctx.dm)
    vec_coord   =   PETSc.DMGetCoordinatesLocal(user_ctx.dm);
    Coord       =   PETSc.DMStagVecGetArray(DMcoord, vec_coord);
    Z_cen       =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,Coord, PETSc.DMSTAG_ELEMENT,    0); # center (has 1 extra)
    Z           =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,Coord, PETSc.DMSTAG_LEFT,       0)
    user_ctx.dz =   Z[2]-Z[1];

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); 
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_LEFT,    0); 
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_ELEMENT, 0); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_LEFT,    0); 
    
    Phi0 =1; dPhi1=8; dPhi2=1; z1=10; z2=40; lambda=1
    dPe1        =   dPhi1/user_ctx.De; 
    dPe2        =   dPhi2/user_ctx.De;
    
    Phi        .=   Phi0 .+ dPhi1.*exp.( -(Z .- z1).^2.0) +  dPhi2.*exp.( -(Z .- z2).^2.0);
    Pe         .=          -dPe1 .*exp.( -(Z .- z1).^2.0) -   dPe2.*exp.( -(Z .- z2).^2.0);

    Phi_old    .=   Phi0 .+ dPhi1.*exp.( -(Z .- z1).^2.0) +  dPhi2.*exp.( -(Z .- z2).^2.0);
    Pe_old     .=          -dPe1 .*exp.( -(Z .- z1).^2.0) -   dPe2.*exp.( -(Z .- z2).^2.0);

    # Copy local into global residual vector
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g); 

    # send back coordinates
    return Z, Z_cen
end

function PopulateCoefficientData!(ctx)

    user_ctx.coeff_l    =   PETSc.DMCreateLocalVector(user_ctx.dmCoeff);
    coeff_array         =   PETSc.DMStagVecGetArray(user_ctx.dmCoeff,user_ctx.coeff_l);

    dm_coord  = PETSc.DMGetCoordinateDM(user_ctx.dmCoeff);
    vec_coord = PETSc.DMGetCoordinatesLocal(user_ctx.dmCoeff);
    coord     = PETSc.DMStagVecGetArray(dm_coord, vec_coord);
    start,n,nExtra = PETSc.DMStagGetCorners(user_ctx.dmCoeff);

    # Get the correct entries for each of our variables in local element-wise storage
    iec = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_DOWN_LEFT, 0);       # location eta corner
    irc = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_DOWN_LEFT, 1);       # location rho corner
    iee = PETSc.DMStagGetLocationSlot(user_ctx.dmCoeff, PETSc.DMSTAG_ELEMENT, 0);         # location eta element
    ixc = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_DOWN_LEFT, 0); # location coord corner
    ixe = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_ELEMENT, 0);   # location coord element

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
    if x < (ctx.xlim[2]-ctx.xlim[1])/2
        if n == 1 return true else return false end
    else
        if n == 1 return false else return true end
    end
end

function SetVecX!(user_ctx,x_g)
    user_ctx.x_l .= 0; 
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g);
end


# Main Solver
nx               =   32;
nz               =   32;
user_ctx.xlim    =   [0,1];
user_ctx.zlim    =   [0,1];
xlim             =   user_ctx.xlim;
zlim             =   user_ctx.zlim;
user_ctx.dx      =   (xlim[2]-xlim[1])/nx;
user_ctx.dz      =   (zlim[2]-zlim[1])/nz;
user_ctx.eta1    =   1;          # viscosity matrix
user_ctx.eta2    =   2;          # viscosity anomaly
user_ctx.rho1    =   3;          # viscosity matrix
user_ctx.rho2    =   4;          # viscosity anomaly
user_ctx.gz      =   -1;
#user_ctx.dt     =   2e-2;       # Note that the timestep has to be tuned a bit depending on the   
user_ctx.dm           =   PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_GHOSTED,PETSc.DM_BOUNDARY_GHOSTED,nx,nz,1,1,0,1,1,PETSc.DMSTAG_STENCIL_BOX,1);  # V edge, P Element
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dm, xlim[1], xlim[2], zlim[1], zlim[2]);            # set coordinates
user_ctx.dmCoeff      =   PETSc.DMStagCreateCompatibleDMStag(user_ctx.dm,2,0,1,0);   # rho and eta on VERTEX, eta on ELEMENT
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dmCoeff, xlim[1], xlim[2], zlim[1], zlim[2]);
user_ctx.dmEps =   PETSc.DMStagCreateCompatibleDMStag(user_ctx.dm,1,0,2,0);   # Exz on VERTEX, Exx and Ezz on ELEMENT
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dmEps, xlim[1], xlim[2], zlim[1], zlim[2]);
user_ctx.dmTau =   PETSc.DMStagCreateCompatibleDMStag(user_ctx.dm,1,0,2,0);   # Txz on VERTEX, Txx and Tzz on ELEMENT
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dmTau, xlim[1], xlim[2], zlim[1], zlim[2]);

PopulateCoefficientData!(user_ctx);

x_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm);
f_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm);
user_ctx.x_l    =   PETSc.DMCreateLocalVector(user_ctx.dm);
user_ctx.Eps_l  =   PETSc.DMCreateLocalVector(user_ctx.dmEps);
user_ctx.Tau_l  =   PETSc.DMCreateLocalVector(user_ctx.dmTau);
#user_ctx.xold_l =   PETSc.DMCreateLocalVector(user_ctx.dm)
#user_ctx.xold_g =   PETSc.DMCreateGlobalVector(user_ctx.dm)
user_ctx.f_l    =   PETSc.DMCreateLocalVector(user_ctx.dm);
J               =   PETSc.DMCreateMatrix(user_ctx.dm);                  # Jacobian from DMStag


# initial non-zero structure of jacobian
#Z, Z_cen        =   SetInitialPerturbations(user_ctx, x_g)
SetVecX!(user_ctx,x_g);

user_ctx.jac, user_ctx.colors = ComputeSparsityPatternJacobian_automatic(user_ctx.x_l.array, user_ctx);

#x0              =   PETSc.VecSeq(rand(size(x_g,1)));
#FormRes!(f_g, x0, user_ctx);
#J_julia,ind     =   FormJacobian!(x0, J, J, user_ctx);

FormRes!(f_g, x_g, user_ctx);
J_julia,ind     =   FormJacobian!(x_g, J, J, user_ctx);

#@show J

sol = J\f_g     # julia vector
x_g .= sol;     # copy to PETSc vecror

# Copy solution to local vector
PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 

# Extract solution
Vx  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_LEFT,      0); 
Vz  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_DOWN,      0); 
P   =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_ELEMENT,   0); 

sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

ix      =   sx[1]:sn[1];           
iz      =   sx[2]:sn[2];

Vx  =   Vx[sx[1]:sn[1]+1,iz];
Vz  =   Vz[ix,sx[2]:sn[2]+1];
P   =    P[ix,iz];

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
heatmap(xe_1D,ze_1D, P', xlabel="Width", ylabel="Depth", title="Pressure")
heatmap(xe_1D,zc_1D, Vz', xlabel="Width", ylabel="Depth", title="Vz")
heatmap(xc_1D,ze_1D, Vx', xlabel="Width", ylabel="Depth", title="Vx")


#S = PETSc.SNES{Float64}(MPI.COMM_SELF, 0; 
#        snes_rtol=1e-12, 
#        snes_monitor=true, 
#        snes_max_it = 500,
#        snes_monitor_true_residual=true, 
#        snes_converged_reason=true);
#S.user_ctx  =       user_ctx;
#

#SetInitialPerturbations(user_ctx, x_g)

#PETSc.setfunction!(S, FormRes!, f_g)
#PETSc.setjacobian!(S, FormJacobian!, J, J)

#FormRes!(f_g, x_g, user_ctx);   
#FormJacobian!(x_g, J, J, user_ctx);

# Preparation of visualisation
#ENV["GKSwstype"]="nul"; 
#if isdir("viz_out")==true
#    rm("viz_out", recursive=true)
#end
#mkdir("viz_out") 
#loadpath = "./viz_out/"; anim = Animation(loadpath,String[])


#time = 0.0;
#it   = 1;
#while time<25.0
#    global time, Z, Z_cen, it
#
#    # Solve one (nonlinear) timestep
#sol = J\f_g     # julia vector
#
#    # Update old local values
#    user_ctx.xold_g  =  x_g;
#    PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
#    PETSc.DMGlobalToLocal(user_ctx.dm, user_ctx.xold_g,  PETSc.INSERT_VALUES,  user_ctx.xold_l) 
#
#    # Update time
#    time += user_ctx.dt;
#    it   += 1;
#
#    if mod(it,20)==0  # Visualisation
#        # Extract values and plot 
#        Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); 
#        Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_LEFT,    0); 
#
#        p1 = plot(Phi[1:end-1], Z_cen[1:end-1],  ylabel="Z", xlabel="ϕ",  xlims=( 0.0, 2.0), label=:none, title="De=$(user_ctx.De)"); 
#        p2 = plot(Pe,           Z    ,                       xlabel="Pe", xlims=(-0.1, 0.1), label=:none, title="$(round(time;sigdigits=3))"); 
#        
#        plot(p1, p2, layout=(1,2)); frame(anim)
#
#        Base.finalize(Pe);  Base.finalize(Phi)
#    end
#
#    println("Timestep $it, time=$time")
#    
#end
#
#gif(anim, "Example_1D.gif", fps = 15)   # create a gif animation
#