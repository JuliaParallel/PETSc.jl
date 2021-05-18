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
using ForwardDiff, SparseArrays

PETSc.initialize()



function FormRes!(cfx_g, cx_g, user_ctx)

    # Note that in PETSc, cx and cfx are pointers to global vectors. 
    
    # Copy global to local vectors
    PETSc.DMGlobalToLocal(user_ctx.dm, cx_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    PETSc.DMGlobalToLocal(user_ctx.dm, cfx_g, PETSc.INSERT_VALUES,  user_ctx.f_l) 

    # Retrieve arrays from the local vectors
    ArrayLocal_x     =   PETSc.DMStagVecGetArrayRead(user_ctx.dm,   user_ctx.x_l);      # array with all local x-data
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm,       user_ctx.f_l);      # array with all local residual
    
    # Compute local residual 
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx)

    # Finalize local arrays
    Base.finalize(ArrayLocal_x)
    Base.finalize(ArrayLocal_f)

    # Copy local into global residual vector
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.f_l, PETSc.INSERT_VALUES, cfx_g) 

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

function FormJacobian!(cx_g, J, P, user_ctx)
    # This requires several steps:
    #
    #   1) Extract local vector from global solution (x) vector
    #   2) Compute local jacobian from the residual routine (note that
    #       this routine requires julia vectors as input)

    # Extract the local vector
    PETSc.DMGlobalToLocal(user_ctx.dm, cx_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    x               =   PETSc.unsafe_localarray(Float64, user_ctx.x_l.ptr;  write=false, read=true)

    f_Residual      =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional arguments into the routine

    J_julia         =   ForwardDiff.jacobian(f_Residual,x);  

    # @show J_julia, size(J_julia)
    n                =   size(P,1)
    J               .=   sparse(J_julia[1:n,1:n]);       

   return J_julia
end

# Define a struct that holds data we need in the local SNES routines below   
mutable struct Data
    dm
    dmCoeff
    coeff_g
    eta1
    eta2
    rho1
    rho2
    g
    dt
    dz
    xlim
    zlim
end
user_ctx = Data(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 


function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
    # Compute the local residual. The vectors include ghost points 
    n           =   2.0
    m           =   3.0
    dt          =   user_ctx.dt;
    dz          =   user_ctx.dz;
    De          =   user_ctx.De;

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   0); 
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_ELEMENT,   0); 
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_LEFT,      0); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_LEFT,      0); 
    
    res_Phi     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   0); 
    res_Pe      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_LEFT,      0); 
    
    # compute the FD stencil
    nPhi        =   length(Phi)-1;    # array length
    
    # Porosity residual @ center points
    i               =   2:nPhi-1
    res_Phi[1]      =   Phi[1]   -1.0;                                  # left BC
    res_Phi[nPhi]   =   Phi[nPhi]-1.0;                                  # right BC
    res_Phi[i]      =   (Phi[i] - Phi_old[i])/dt +   ((Phi[i .+ 0].^n) .* ( (Pe[i .+ 1] - Pe[i     ])/dz .+ 1.0)  
                                                    - (Phi[i .- 1].^n) .* ( (Pe[i     ] - Pe[i .- 1])/dz .+ 1.0))/dz
   
    # Pressure update @ nodal points
    nP              =    length(Pe);    
    i               =    2:nP-1
    res_Pe[1]       =    Pe[1]  - 0.;                                   # left BC
    res_Pe[nP]      =    Pe[nP] - 0.;                                   # right BC
    res_Pe[i]       =    De.*(Pe[i]-Pe_old[i])/dt -   ((Phi[i .+ 0].^n) .* ( (Pe[i .+ 1] - Pe[i     ])/dz .+ 1.0)  
                                                     - (Phi[i .- 1].^n) .* ( (Pe[i     ] - Pe[i .- 1])/dz .+ 1.0))/dz +    (Phi[i].^m)   .* Pe[i];

    # Cleanup
    Base.finalize(Phi);    Base.finalize(Phi_old);       
    Base.finalize(Pe);     Base.finalize(Pe_old);       
                                           
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
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g) 

    # send back coordinates
    return Z, Z_cen
end

function PopulateCoefficientData!(ctx)

    user_ctx.coeff_g    =   PETSc.DMCreateGlobalVector(user_ctx.dmCoeff);
    coeff_l             =   PETSc.DMCreateLocalVector(user_ctx.dmCoeff);
    coeff_array         =   PETSc.DMStagVecGetArray(user_ctx.dmCoeff,coeff_l);

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

    PETSc.DMLocalToGlobal(user_ctx.dmCoeff,coeff_l, PETSc.INSERT_VALUES, user_ctx.coeff_g);


end

function GetPhase(ctx,x,n)
    if x < (ctx.xlim[2]-ctx.xlim[1])/2
        if n == 1 return true else return false end
    else
        if n == 1 return false else return true end
    end
end


# Main Solver
nx               =   5;
nz               =   5;
user_ctx.xlim    =   [0,1];
user_ctx.zlim    =   [0,1];
xlim             =   user_ctx.xlim;
zlim             =   user_ctx.zlim;
dx               =   (xlim[2]-xlim[1])/nx;
dz               =   (zlim[2]-zlim[1])/nz;
user_ctx.eta1    =   1;          # viscosity matrix
user_ctx.eta2    =   2;          # viscosity anomaly
user_ctx.rho1    =   3;          # viscosity matrix
user_ctx.rho2    =   4;          # viscosity anomaly
user_ctx.g       =   -1;
#user_ctx.dt     =   2e-2;       # Note that the timestep has to be tuned a bit depending on the   
user_ctx.dm      =   PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,nx,nz,1,1,0,1,1,PETSc.DMSTAG_STENCIL_BOX,1);  # V edge, P Element
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dm, xlim[1], xlim[2], zlim[1], zlim[2]);            # set coordinates
user_ctx.dmCoeff =   PETSc.DMStagCreateCompatibleDMStag(user_ctx.dm,2,0,1,0);   # rho and eta on VERTEX, eta on ELEMENT
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dmCoeff, xlim[1], xlim[2], zlim[1], zlim[2]);  

PopulateCoefficientData!(user_ctx);

#x_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)
#f_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)
#user_ctx.x_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
#user_ctx.xold_l =   PETSc.DMCreateLocalVector(user_ctx.dm)
#user_ctx.xold_g =   PETSc.DMCreateGlobalVector(user_ctx.dm)
#user_ctx.f_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
#J               =   PETSc.DMCreateMatrix(user_ctx.dm);                  # Jacobian from DMStag


# initial non-zero structure of jacobian
#Z, Z_cen        =   SetInitialPerturbations(user_ctx, x_g)

#x0              =   PETSc.VecSeq(rand(size(x_g,1)));
#J_julia         =   FormJacobian!(x0.ptr, J, J, user_ctx)


#S = PETSc.SNES{Float64}(MPI.COMM_SELF, 0; 
#        snes_rtol=1e-12, 
#        snes_monitor=true, 
#        snes_max_it = 500,
#        snes_monitor_true_residual=true, 
#        snes_converged_reason=true);
#S.user_ctx  =       user_ctx;


#SetInitialPerturbations(user_ctx, x_g)

#PETSc.setfunction!(S, FormRes!, f_g)
#PETSc.setjacobian!(S, FormJacobian!, J, J)

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
#    PETSc.solve!(x_g, S);
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