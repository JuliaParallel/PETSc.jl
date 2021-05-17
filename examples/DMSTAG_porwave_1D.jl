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
using SparseArrays, SparseDiffTools, ForwardDiff

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

function  f(out, x, user_ctx)

    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all local x-data
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

    # Extract the local vector
    PETSc.DMGlobalToLocal(user_ctx.dm, cx_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    x               =   PETSc.unsafe_localarray(Float64, user_ctx.x_l.ptr;  write=false, read=true)


    if isnothing(user_ctx.jac)
        # Compute sparsity pattern of jacobian. This is relatvely slow, but only has to be done once.
        # Theoretically, more efficient tools for this exists (jacobian_sparsity in the SparsityDetection.jl package),
        # but they don't seem to work with the PETSc approach we use. Therefore we employ      
        f_Residual  =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional arguments into the routine
        J_julia     =   ForwardDiff.jacobian(f_Residual,x);  

        # employ sparse structure to compute jacobian - to be moved inside routines
        jac         =   sparse(J_julia);
        colors      =   matrix_colors(jac)          # the number of nonzeros per row

    else
        jac     =   user_ctx.jac;
        colors  =   user_ctx.colors;

    end
    out         =   similar(x);
        
   f_Res           =   ((out,x)->f(out, x, user_ctx));        # pass additional arguments into the routine
   forwarddiff_color_jacobian!(jac, f_Res, x, colorvec = colors)

    ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);    # extr
    J              .=   jac[ind,ind];    

    user_ctx.jac    =   jac;
    user_ctx.colors =   colors;
    
   return jac, ind
end

# Define a struct that holds data we need in the local SNES routines below   
mutable struct Data
    dm
    x_l
    xold_l
    xold_g
    f_l
    dt
    dz
    De
    jac
    colors
end
user_ctx = Data(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 


function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
    # Compute the local residual. The vectors include ghost points 
    n           =   3.0
    m           =   2.0
    dt          =   user_ctx.dt;
    dz          =   user_ctx.dz;
    De          =   user_ctx.De;

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   0); 
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_ELEMENT,   0); 
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   1); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_ELEMENT,   1); 
    
    res_Phi     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   0); 
    res_Pe      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   1); 
    
    # compute the FD stencil
    sx, sn          =     PETSc.DMStagGetCentralNodes(dm);          # indices of (center/element) points, not including ghost values.
    
    # Porosity residual @ center points
    iz                  =   sx[1]:sn[1];            # Phi is on center points
    i                   =   iz[2:end-1];
    res_Phi[iz[1]]      =   Phi[iz[1]]   - 1.0;                                  # left BC
    res_Phi[iz[end]]    =   Phi[iz[end]] - 1.0;                                  # right BC
    res_Phi[i]          =   (Phi[i] - Phi_old[i])/dt + De.*(Pe[i]-Pe_old[i])/dt + (Phi[i].^m)   .* Pe[i]
   
    # Pressure update @ nodal points
    iz                  =   sx[1]:sn[1];          # Pe is on center points as well (dof=2)
    i                   =   iz[2:end-1];
    res_Pe[iz[1]]       =   Pe[iz[1]]   - 0.;                                   # left BC
    res_Pe[iz[end]]     =   Pe[iz[end]] - 0.;                                   # right BC
    res_Pe[i]           =   De.*(Pe[i]-Pe_old[i])/dt -   (   ((0.5*(Phi[i .+ 1] + Phi[i .+ 0])).^n) .* ( (Pe[i .+ 1] - Pe[i     ])/dz .+ 1.0)  
                                                          -  ((0.5*(Phi[i .- 1] + Phi[i .+ 0])).^n) .* ( (Pe[i     ] - Pe[i .- 1])/dz .+ 1.0))/dz +    (Phi[i].^m)   .* Pe[i];

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
    user_ctx.dz =   Z_cen[2]-Z_cen[1];

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); 
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 1); 
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_ELEMENT, 0); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_ELEMENT, 1); 
    
    Phi0 =1.0; dPhi1=8.0; dPhi2=1.0; z1=0.0; z2=40.0; lambda=1.0
    dPe1        =   dPhi1/user_ctx.De; 
    dPe2        =   dPhi2/user_ctx.De;
    
    Phi        .=   Phi0 .+ dPhi1.*exp.( -((Z_cen .- z1).^2.0)/lambda^2) +  dPhi2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);
    Pe         .=          -dPe1 .*exp.( -((Z_cen .- z1).^2.0)/lambda^2) -   dPe2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);

    Phi_old    .=   Phi0 .+ dPhi1.*exp.( -((Z_cen .- z1).^2.0)/lambda^2) +  dPhi2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);
    Pe_old     .=          -dPe1 .*exp.( -((Z_cen .- z1).^2.0)/lambda^2) -   dPe2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);

    # Copy local into global residual vector
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g) 

    # send back coordinates
    return Z_cen
end


# Main Solver
nx              =   2001;
L               =   150;
user_ctx.De     =   1e-2;        # Deborah number
user_ctx.dt     =   1e-3;        # Note that the timestep has to be tuned a bit depending on De in order to obtain convergence     
user_ctx.dm     =   PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,nx,0,2);  # both Phi and Pe are on center points 
PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dm, -20, L)            # set coordinates
x_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)


f_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)
user_ctx.x_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
user_ctx.xold_l =   PETSc.DMCreateLocalVector(user_ctx.dm)
user_ctx.xold_g =   PETSc.DMCreateGlobalVector(user_ctx.dm)
user_ctx.f_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
J               =   PETSc.DMCreateMatrix(user_ctx.dm);                  # Jacobian from DMStag


# initial non-zero structure of jacobian
Z_cen           =   SetInitialPerturbations(user_ctx, x_g)

x0              =   PETSc.VecSeq(rand(size(x_g,1)));
J_julia,ind     =   FormJacobian!(x0.ptr, J, J, user_ctx)


S = PETSc.SNES{Float64}(MPI.COMM_SELF, 0; 
        snes_rtol=1e-12, 
        snes_monitor=true, 
        snes_max_it = 500,
        snes_monitor_true_residual=true, 
        snes_converged_reason=true);
S.user_ctx  =       user_ctx;


SetInitialPerturbations(user_ctx, x_g)

PETSc.setfunction!(S, FormRes!, f_g)
PETSc.setjacobian!(S, FormJacobian!, J, J)

# Preparation of visualisation
ENV["GKSwstype"]="nul"; 
if isdir("viz_out")==true
    rm("viz_out", recursive=true)
end
mkdir("viz_out") 
loadpath = "./viz_out/"; anim = Animation(loadpath,String[])


time = 0.0;
it   = 1;
while time<25.0
    global time, Z, Z_cen, it

    # Solve one (nonlinear) timestep
    PETSc.solve!(x_g, S);

    # Update old local values
    user_ctx.xold_g  =  x_g;
    PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    PETSc.DMGlobalToLocal(user_ctx.dm, user_ctx.xold_g,  PETSc.INSERT_VALUES,  user_ctx.xold_l) 

    # Update time
    time += user_ctx.dt;
    it   += 1;

    if mod(it,200)==0  # Visualisation
        # Extract values and plot 
        Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); 
        Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 1); 

        p1 = plot(Phi[1:end-1], Z_cen[1:end-1],  ylabel="Z", xlabel="ϕ",  xlims=( 0.0, 6.0), label=:none, title="De=$(user_ctx.De)"); 
        p2 = plot(Pe[1:end-1],  Z_cen[1:end-1]    ,                       xlabel="Pe", xlims=(-1.25, 1.25), label=:none, title="$(round(time;sigdigits=3))"); 
       # p1 = plot(Phi[1:end-1], Z_cen[1:end-1],  ylabel="Z", xlabel="ϕ",   label=:none, title="De=$(user_ctx.De)"); 
       # p2 = plot(Pe[1:end-1],  Z_cen[1:end-1]    ,                       xlabel="Pe",  label=:none, title="$(round(time;sigdigits=3))"); 
        
        plot(p1, p2, layout=(1,2)); frame(anim)

        Base.finalize(Pe);  Base.finalize(Phi)
    end

    println("Timestep $it, time=$time")
    
end

gif(anim, "Example_1D.gif", fps = 15)   # create a gif animation
