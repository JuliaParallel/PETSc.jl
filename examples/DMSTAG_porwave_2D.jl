# EXCLUDE FROM TESTING
# NOTE: This is temporarily not working until we merge the DMSTAG routines with the new Clang branch
#
#
# This is an example of a 2D viscoelastic porosity wave as described in 
# Vasyliev et al. Geophysical Research Letters (25), 17. p. 3239-3242
# 
#
# This example only requires the specification of a residual routine; automatic
# differentiation is used to generate the jacobian. See the equivalent 1D version
# as well. 
# This example also shows how to employ ghost points to set lateral boundary conditions  

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
        error("You first have to define the sparsity pattern of the jacobian")
    else
        jac     =   user_ctx.jac;
        colors  =   user_ctx.colors;
    end
    out         =   similar(x);
        
    f_Res           =   ((out,x)->f(out, x, user_ctx));        # pass additional arguments into the routine

    if isnothing(user_ctx.jac_cache)
        # Allocate data required for jacobian computations
        user_ctx.jac_cache = ForwardColorJacCache(f_Res,x; colorvec=colors, sparsity = jac);
    end

    @time forwarddiff_color_jacobian!(jac, f_Res, x, user_ctx.jac_cache);

    ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);    # extr
    J              .=   jac[ind,ind];    


    user_ctx.jac    =   jac;
    user_ctx.colors =   colors;
    
   return jac, ind
end


# Define a struct that holds data we need in the residual SNES routines below   
mutable struct Data_PorWav2D
    dm
    x_l
    xold_l
    xold_g
    f_l
    dt
    dx
    dz
    De
    jac
    colors
    jac_cache
end
user_ctx = Data_PorWav2D(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 

function ComputeSparsityPatternJacobian(x_l, user_ctx)
    # This computes the sparsity pattern of our jacobian by hand. 
    # That is signficantly faster than the automatic method, yet you will have to analyze your residual routine for it.
    #
    # As you will see, there are however any similaries between the two routines, so it is usually not all that difficult
    # to define this

    AddInd!(i, array) =  append!(i,vec(array));   # small helper routine for readability of the lines below

    n               =   length(x_l);
    ind_x           =   Vector(1:n);
    ind_f           =   Vector(1:n);
    ArrayLocal_x    =   PETSc.DMStagVecGetArray(user_ctx.dm, ind_x);        # array with all local x-data
    ArrayLocal_f    =   PETSc.DMStagVecGetArray(user_ctx.dm, ind_f);        # array with all local residual
    
    Phi             =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   0); 
    Pe              =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   1); 
    res_Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   0); 
    res_Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   1); 
    
    # compute the FD stencil
    sx, sn              =   PETSc.DMStagGetCentralNodes(user_ctx.dm);          # indices of (center/element) points, not including ghost values.    

    # Porosity residual @ center points
    ix                  =   sx[1]:sn[1];           
    iz                  =   sx[2]:sn[2];      

    #res_Phi[ix,iz]      =   (Phi[ix,iz] - Phi_old[ix,iz])/dt + De.*(Pe[ix,iz]-Pe_old[ix,iz])/dt + (Phi[ix,iz].^m) .* Pe[ix,iz]
    i,j                 =   Vector{Int64}(), Vector{Int64}()
    AddInd!(i,res_Phi[ix,iz]), AddInd!(j,Phi[ix,iz]);   #   (Phi[ix,iz])/dt term
    AddInd!(i,res_Phi[ix,iz]), AddInd!(j,Pe[ix,iz]);    #   new non-zeros because of (Phi[ix,iz].^m) .* Pe[ix,iz] term
   
    # Pressure update @ nodal points
    ix                  =   sx[1]:sn[1];                # lateral BC's are set using ghost points
    iz                  =   sx[2]+1:sn[2]-1;            # constant Pe on top and bottom, so this is only for center points
    Pe[ix[1]-1,:]       =   Pe[ix[1],:];                # ghost points on left and right size (here: flux free)
    Pe[ix[end]+1,:]     =   Pe[ix[end],:];

    # BC's
    AddInd!(i, res_Pe[:,iz[1]-1]),      AddInd!(j,Pe[:,iz[1]-1]);       # bottom BC
    AddInd!(i, res_Pe[:,iz[end]+1]),    AddInd!(j,Pe[:,iz[end]+1]);     # top BC
    
    
#    res_Pe[ix,iz]       =   De*( Pe[ix,iz] - Pe_old[ix,iz])/dt +
#    (Phi[ix,iz].^m).* Pe[ix,iz]     -
#    ((   ((0.5.*(Phi[ix,iz .+ 1] + Phi[ix    ,iz ])).^n) .* ( (Pe[ix,iz .+ 1] - Pe[ix, iz     ])/dz .+ 1.0)  -
#         ((0.5.*(Phi[ix,iz .- 1] + Phi[ix    ,iz ])).^n) .* ( (Pe[ix,iz     ] - Pe[ix, iz .- 1])/dz .+ 1.0))/dz)  
# -  ((   ((0.5.*(Phi[ix,iz     ] + Phi[ix.+ 1,iz ])).^n) .* ( (Pe[ix.+ 1,iz ] - Pe[ix, iz     ])/dx       )  -
#         ((0.5.*(Phi[ix,iz     ] + Phi[ix.- 1,iz ])).^n) .* ( (Pe[ix    ,iz ] - Pe[ix.- 1,iz  ])/dx       ))/dx)
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j,Phi[ix,iz]);          # new part of (Phi[ix,iz].^m).* Pe[ix,iz] 
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j,Phi[ix,iz .+ 1 ]);    # Phi[ix,iz .+ 1]
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j,Phi[ix,iz .- 1 ]);    # 
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j,Phi[ix .- 1, iz]);    # 
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j,Phi[ix .+ 1, iz]);    # 
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j, Pe[ix,iz]      );    # De*( Pe[ix,iz] - Pe_old[ix,iz])/dt
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j, Pe[ix .+ 1, iz]);    # Pe[ix.+ 1,iz ]
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j, Pe[ix .- 1, iz]);    # Pe[ix.- 1,iz ]
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j, Pe[ix, iz .- 1]);    # Pe[ix ,iz - 1]
    AddInd!(i, res_Pe[ix,iz]),  AddInd!(j, Pe[ix, iz .+ 1]);    # Pe[ix ,iz + 1]
 
    jac                 =   sparse(i,j,ones(Float64,length(i[:])),n,n); # create the sparse jacobian
    colors              =   matrix_colors(jac)                          # the number of nonzeros per row
    return jac, colors
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

function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
    # Compute the local residual. The vectors include ghost points 
    n           =   3.0
    m           =   2.0
    dt          =   user_ctx.dt;
    dz          =   user_ctx.dz;
    dx          =   user_ctx.dx;
    De          =   user_ctx.De;

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   0); 
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_ELEMENT,   0); 
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x   ,   PETSc.DMSTAG_ELEMENT,   1); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,user_ctx.xold_l,   PETSc.DMSTAG_ELEMENT,   1); 
    
    res_Phi     =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   0); 
    res_Pe      =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f,      PETSc.DMSTAG_ELEMENT,   1); 
    
    # compute the FD stencil
    sx, sn              =   PETSc.DMStagGetCentralNodes(dm);          # indices of (center/element) points, not including ghost values.
    
    # Porosity residual @ center points
    ix                  =   sx[1]:sn[1];           
    iz                  =   sx[2]:sn[2];      
    res_Phi[ix,iz]      =   (Phi[ix,iz] - Phi_old[ix,iz])/dt + De.*(Pe[ix,iz]-Pe_old[ix,iz])/dt + (Phi[ix,iz].^m) .* Pe[ix,iz]
    
    # Pressure update @ nodal points
    ix                  =   sx[1]:sn[1];                # lateral BC's are set using ghost points
    iz                  =   sx[2]+1:sn[2]-1;            # constant Pe on top and bottom, so this is only for center points
    Pe[ix[1]-1,:]       =   Pe[ix[1],:];                # ghost points on left and right size (here: flux free)
    Pe[ix[end]+1,:]     =   Pe[ix[end],:];
    res_Pe[:,iz[1]-1]   =   Pe[:,iz[1]-1]   .- 0.;      # bottom BC
    res_Pe[:,iz[end]+1] =   Pe[:,iz[end]+1] .- 0.;      # top BC

  
    res_Pe[ix,iz]       =   De*( Pe[ix,iz] - Pe_old[ix,iz])/dt +
                                (Phi[ix,iz].^m).* Pe[ix,iz]     -
                                ((   ((0.5.*(Phi[ix,iz .+ 1] + Phi[ix    ,iz ])).^n) .* ( (Pe[ix,iz .+ 1] - Pe[ix, iz     ])/dz .+ 1.0)  -
                                     ((0.5.*(Phi[ix,iz .- 1] + Phi[ix    ,iz ])).^n) .* ( (Pe[ix,iz     ] - Pe[ix, iz .- 1])/dz .+ 1.0))/dz)  -
                                ((   ((0.5.*(Phi[ix,iz     ] + Phi[ix.+ 1,iz ])).^n) .* ( (Pe[ix.+ 1,iz ] - Pe[ix, iz     ])/dx       )  -
                                     ((0.5.*(Phi[ix,iz     ] + Phi[ix.- 1,iz ])).^n) .* ( (Pe[ix    ,iz ] - Pe[ix.- 1,iz  ])/dx       ))/dx)
                               
        
    # remark to the lines above: these are long, multi line, expressions. For code readability, I chopped them up over several lines. 
    # In that case you MUST put the +/- sign at the end of every line. If not, some lines will not be evaluated but the code will still run, which is tricky to spot


    # Cleanup
    Base.finalize(Phi);  Base.finalize(Phi_old);       
    Base.finalize(Pe);   Base.finalize(Pe_old);       

end


function SetInitialPerturbations(user_ctx, x_g)
    # Computes the initial perturbations as in the paper

    # Retrieve coordinates from DMStag
    DMcoord     =   PETSc.DMGetCoordinateDM(user_ctx.dm)
    vec_coord   =   PETSc.DMGetCoordinatesLocal(user_ctx.dm);
    Coord       =   PETSc.DMStagVecGetArray(DMcoord, vec_coord);
    X_cen       =   Coord[:,:,1]; 
    Z_cen       =   Coord[:,:,2]; 
    
    user_ctx.dz =   Z_cen[2,2]-Z_cen[2,1];
    user_ctx.dx =   X_cen[1,2]-X_cen[2,2];

    Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); # center
    Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 1); # center
    
    Phi_old     =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_ELEMENT, 0); 
    Pe_old      =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.xold_l,  PETSc.DMSTAG_ELEMENT, 1); 
    
    
    Phi0 =1.0; dPhi1=8.0; dPhi2=1.0; z1=0.0; z2=40.0; x1=0.0; x2=0.0; lambda=1.0
    dPe1        =   dPhi1/user_ctx.De; 
    dPe2        =   dPhi2/user_ctx.De;
    
    Phi        .=   Phi0 .+ dPhi1.*exp.( -((Z_cen .- z1).^2.0)/lambda^2) +  dPhi2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);
    Pe         .=          -dPe1 .*exp.( -((Z_cen .- z1).^2.0)/lambda^2) -   dPe2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);

    Phi_old    .=   Phi0 .+ dPhi1.*exp.( -((Z_cen .- z1).^2.0)/lambda^2) +  dPhi2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);
    Pe_old     .=          -dPe1 .*exp.( -((Z_cen .- z1).^2.0)/lambda^2) -   dPe2.*exp.( -((Z_cen .- z2).^2.0)/lambda^2);

    # Copy local into global residual vector
    PETSc.DMLocalToGlobal(user_ctx.dm,user_ctx.x_l, PETSc.INSERT_VALUES, x_g) 

    # send back coordinates (mainly for plotting)
    return X_cen, Z_cen
end


# Main solver
nx, nz          =   100, 500;
W, L            =   100,150;

user_ctx.De     =   1e-1;        # Deborah number
user_ctx.dt     =   1e-5;       # Note that the timestep has to be tuned a bit depending on the     

# Both Pe and Phi are @ defined centers
user_ctx.dm =   PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_GHOSTED,PETSc.DM_BOUNDARY_NONE,nx,nz,1,1,0,0,2,PETSc.DMSTAG_STENCIL_BOX,1)


PETSc.DMStagSetUniformCoordinatesExplicit(user_ctx.dm, -W/2, W/2, -20, L) # set coordinates
x_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)


f_g             =   PETSc.DMCreateGlobalVector(user_ctx.dm)
user_ctx.x_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
user_ctx.xold_l =   PETSc.DMCreateLocalVector(user_ctx.dm)
user_ctx.xold_g =   PETSc.DMCreateGlobalVector(user_ctx.dm)
user_ctx.f_l    =   PETSc.DMCreateLocalVector(user_ctx.dm)
J               =   PETSc.DMCreateMatrix(user_ctx.dm);                  # Jacobian from DMStag


# initial non-zero structure of jacobian
X_cen, Z_cen    =   SetInitialPerturbations(user_ctx, x_g)

# plotting stuff initial setup
Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT,   0); 
Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT,   1); 
Phi_pl      =   Phi[2:end-1,2:end-1];
Pe_pl       =   Pe[2:end-1,2:end-1]
xc_1D,zc_1D =   X_cen[2:end-1,2], Z_cen[2,2:end-1];

#heatmap(xc_1D,zc_1D,Phi_pl', xlabel="Width", ylabel="Depth", title="Phi")
heatmap(xc_1D,zc_1D,Pe_pl', xlabel="Width", ylabel="Depth", title="Pe")


# Compute sparsity & coloring
# do this in case you solve new equations:
#@time user_ctx.jac, user_ctx.colors = ComputeSparsityPatternJacobian_automatic(user_ctx.x_l.array, user_ctx);   

# much faster, but requires reformulation in case you solve new equations:
@time user_ctx.jac, user_ctx.colors = ComputeSparsityPatternJacobian(user_ctx.x_l.array, user_ctx);    


S = PETSc.SNES{Float64}(MPI.COMM_SELF; 
        snes_rtol=1e-12, 
        snes_monitor=true, 
        snes_max_it = 500,
        snes_monitor_true_residual=true, 
        snes_converged_reason=true);
S.user_ctx  =       user_ctx;


if true


#SetInitialPerturbations(user_ctx, x_g)

PETSc.setfunction!(S, FormRes!, f_g)
PETSc.setjacobian!(S, FormJacobian!, J, J)



# Preparation of visualisation
ENV["GKSwstype"]="nul"; 
if isdir("viz2D_out")==true
    rm("viz2D_out", recursive=true)
end
mkdir("viz2D_out") 
loadpath = "./viz2D_out/"; anim = Animation(loadpath,String[])


t    = 0.0;
it   = 1;
while t < 25.0
    global t, Z, Z_cen, it

    # Solve one (nonlinear) timestep
    PETSc.solve!(x_g, S);

    # Update old local values
    user_ctx.xold_g  =  x_g;
    PETSc.DMGlobalToLocal(user_ctx.dm, x_g,  PETSc.INSERT_VALUES,  user_ctx.x_l) 
    PETSc.DMGlobalToLocal(user_ctx.dm, user_ctx.xold_g,  PETSc.INSERT_VALUES,  user_ctx.xold_l) 

    # Update time
    t    += user_ctx.dt;
    it   += 1;

    if mod(it,20)==0  # Visualisation
        # Extract values and plot 
        Phi         =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 0); 
        Pe          =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l,     PETSc.DMSTAG_ELEMENT, 1); 

        #p1 = plot(Phi[2,2:end-1], zc_1D,  ylabel="Z", xlabel="ϕ",  xlims=( 0.0, 2.0), label=:none, title="De=$(user_ctx.De)"); 
        #p2 = plot(Pe[2,2:end-1],  zc_1D ,                       xlabel="Pe", xlims=(-0.1, 0.1), label=:none, title="$(round(time;sigdigits=3))"); 
        p1 = plot(Phi[2,2:end-1], zc_1D,  ylabel="Z", xlabel="ϕ", label=:none, title="De=$(user_ctx.De)"); 
        p2 = plot(Pe[2,2:end-1],  zc_1D ,                       xlabel="Pe",  label=:none, title="$(round(t;sigdigits=3))"); 
        
        #p1 = heatmap(xc_1D,zc_1D,Phi_pl', xlabel="Width", ylabel="Depth", title="Phi")
        #p2 = heatmap(xc_1D,zc_1D,Pe_pl', xlabel="Width", ylabel="Depth", title="Pe")

        plot(p1, p2, layout=(1,2)); frame(anim)

        Base.finalize(Pe);  Base.finalize(Phi)
    end

    println("Timestep $it, time=$t")
    
end

gif(anim, "Example_2D.gif", fps = 15)   # create a gif animation
end