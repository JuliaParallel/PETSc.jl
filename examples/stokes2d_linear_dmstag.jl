# EXCLUDE FROM TESTING
# This shows how to solve the 2D incompressible Stokes equations using SNES solvers,
# using a staggered grid discretization and a Velocity-Pressure formulation
#   
#   Governing equations:
#              dVx/dx + dVz/dz = αP                             |   Mass balance (with numerical compressibility)
#   -dP/dx + dTxx/dx + dTxz/dz = 0                              |   Horizontal force balance
#   -dP/dz + dTxz/dx + dTzz/dz = rho*g                          |   Vertical force balance
#   
#   with:
#       Exx = dVx/dx, Ezz = dVz/dz, Exz = 0.5*(dVx/dz + dVz/dx) |   Strain rate definition
#
#      | Txx |   | 2*eta   0     0   |  | Exx |
#      | Tzz | = |   0   2*eta   0   |  | Ezz |                 |   Linear viscous (isotropic) rheology
#      | Txz |   |   0     0   2*eta |  | Exz |
#
# This example also uses the Automatic differentiation package ForwardDiff

using PETSc, MPI
using ForwardDiff, SparseArrays, LinearAlgebra, CairoMakie

petsclib = first(PETSc.petsclibs);
PETSc.initialized(petsclib) || PETSc.initialize(petsclib)


"""
    rVx,rVz,rP = local_residuals(Vx_l::Array{_T,2}, Vz_l::Array{_T,2}, P_l::Array{_T,2}, params::NamedTuple) 

Computes the local residual for a small patch of Vx,Vz,P
"""
function local_residuals(Vx_l::Matrix, Vz_l::Matrix, P_l::Matrix, params::NamedTuple) 
    # strainrates
    Exx_l = diff(Vx_l,dims=1) ./ params.Δx;            # dVx/dx
    Ezz_l = diff(Vz_l,dims=2) ./ params.Δz;            # dVz/dz
    Exz_l = 0.5*(diff(Vz_l[:,2:end],dims=1) ./ params.Δx + diff(Vx_l[2:end,:],dims=2) ./ params.Δz);             # dVz/dz

    # 2nd invariant of strainrate @ center (interpolate xz ->xx)
    # NOTE: I think this is using the correct Exx/Ezz points, but 
    # better to double check
    #Eii = sqrt(Exx_l[2,2]^2 + Ezz_l[2,2]^2 + 2*sum(Exz_l.^2)/4);
            
    # once we have that we can compute Tii
    #Tii = 2*params.ηl * Eii
    #η   = Tii / (2*Eii + 1e-16)   # avoid division by zero

    # stresses, assuming linear viscosity
    Txx_l = 2*params.η_center[1:2,:] .* Exx_l;        # Txx = 2*eta*Exx
    Tzz_l = 2*params.η_center[:,1:2] .* Ezz_l;        # Tzz = 2*eta*Ezz
    Txz_l = 2*params.η_vertex .* Exz_l;        # Txz = 2*eta*Exz

    # compute derivatives
    dPdx_l   = diff(  P_l[:,2:end  ],dims=1) ./ params.Δx;
    dTxxdx_l = diff(Txx_l[:,2:end-1],dims=1) ./ params.Δx;
    dTxzdz_l = diff(Txz_l[1:end-1,:],dims=2) ./ params.Δz;
            
    # x-momentum
    rVx     = -dPdx_l[1] + dTxxdx_l[1] + dTxzdz_l[1];

    # compute derivatives
    dPdz_l   = diff(P_l[2:end,:],dims=2) ./ params.Δz;
    dTzzdz_l = diff(Tzz_l[2:end-1,:],dims=2) ./ params.Δz;
    dTxzdx_l = diff(Txz_l[:,1:end-1],dims=1) ./ params.Δx;
            
    # z-momentum
    rVz = -dPdz_l[1] + dTzzdz_l[1] + dTxzdx_l[1] - params.ρ;
            
    # mass conservation
    rP = Exx_l[2,2] + Ezz_l[2,2] + (1/params.κ) .* P_l[2,2];
    #rP = Exx_l[2,2] + Ezz_l[2,2] ;
    
    return rVx, rVz, rP
end

"""
    res = local_residual_vec(x_in::Vector, params::NamedTuple) 

Residual function that calls `local_residuals` and can be used with ForwardDiff.
It also sets ghost point values if needed. 
"""
function local_residual_vec(x_in::Vector, params, I, N)
    nVx  = (3,3)
    nVz  = (3,3)
    nP   = (2,2)
    Vx_l = reshape(x_in[1:9],   nVx)
    Vz_l = reshape(x_in[10:18], nVz)
    P_l  = reshape(x_in[19:22], nP )

    # we need to deal with BC's here; there may be better ways to do that (e.g, adjusting )
    if I[2]==1 
        if BC.bottom == :free_slip
            Vx_l[:,1] = Vx_l[:,2]   # free slip
        end
    end
    if I[2]==N[2]
        if BC.top == :free_slip
            Vx_l[:,3] = Vx_l[:,2]   # free slip
        end
    end
    if I[1]==1 
        if BC.left == :free_slip
            Vz_l[1,:] = Vz_l[:,2]   # free slip
        end
    end
    if I[1]==N[1] 
        if BC.right == :free_slip
            Vz_l[3,:] = Vz_l[2,:]   # free slip
        end
    end
                
    rVx,rVz,rP = local_residuals(Vx_l, Vz_l, P_l, params)

    return [rVx,rVz,rP]
end

function FormRes!(r_g, snes, x_g, user_ctx)
    dm = PETSc.getDM(snes)
    
    # Copy global to local vectors
    LibPETSc.VecSet(petsclib, user_ctx.r_l, 0.0) # set residual to zero before accumulating contributions
    LibPETSc.VecSet(petsclib, r_g, 0.0) # set residual to zero before accumulating contributions
    
    LibPETSc.VecSet(petsclib, user_ctx.x_l, 0.0) # set solution to zero before accumulating contributions
    PETSc.dm_global_to_local!(x_g, user_ctx.x_l, dm)

    # ghost point values
    set_ghostpoint_values!(user_ctx.x_l, dm, petsclib, user_ctx.BC )
    
    # get coordinates
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    # get the corners 
    corners       = PETSc.getcorners_dmstag(dm)

    LibPETSc.VecSet(petsclib,user_ctx.r_l, 0.0) # set residual to zero before accumulating contributions
    Xlocal = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.x_l)
    Rlocal = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.r_l)
    
    P  = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

    rP  = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    rVx = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    rVz = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

    Δx, Δz = user_ctx.dx, user_ctx.dz;

    # Get coefficient array for density and eta (stored in user_ctx.coeff_l)
    coeff_array = LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dmCoeff, user_ctx.coeff_l)
    rho_vz     = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_LEFT, 0)])
    eta_center = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_ELEMENT,   0)])
    eta_vertex = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_DOWN_LEFT, 0)])

#=
    # ----
    # just checking: do this by hand
    Exx = diff(Vx[1:end,1:end],dims=1) ./ Δx;             # dVx/dx
    Ezz = diff(Vz[1:end,1:end],dims=2) ./ Δz;             # dVz/dz
    Exz = 0.5 .* (diff(Vx[2:end,:],dims=2) ./ Δz .+         # 0.5*(dVx/dz + dVz/dx)
                  diff(Vz[:,2:end],dims=1) ./ Δx);

    Txx = 2 .* user_ctx.eta1 .* Exx;    # Txx = 2*eta*Exx
    Tzz = 2 .* user_ctx.eta1 .* Ezz;    # Tzz = 2*eta*Ezz

    # We actually have to do some padding around this to take BC's properly into account
    Txz = 2 .* user_ctx.eta1 .* Exz;    # Txz = 2*eta*Exz

                ρ = rho_vz[ix, iy];
    dP_dx = diff(P[2:end-1,2:end-1],dims=1) ./ Δx;
    dTxx_dx = diff(Txx[2:end,2:end-1],dims=1) ./ Δx;
    dTxz_dz = diff(Txz[2:end-1,1:end],dims=2) ./ Δz;
    rM1 = -dP_dx + dTxx_dx + dTxz_dz        

    # 2nd momentum balance
    dP_dz   = diff(P[2:end-1,2:end-1],dims=2) ./ Δz;
    dTzz_dz = diff(Tzz[2:end-1,2:end],dims=2) ./ Δz;
    dTxz_dx = diff(Txz[1:end,2:end-1],dims=1) ./ Δx;
    
    ρg = zeros(size(dP_dz))
    for i=1:size(ρg,1), j=1:size(ρg,2)
        x,z = X_coord[i+1,1], Z_coord[j,2];         # coordinate of Vz points
        if GetPhase(user_ctx,x,1)
            ρg[i,j] = user_ctx.rho1;
        else
            ρg[i,j] = user_ctx.rho2;  
        end
    end


    rM2 = -dP_dz + dTzz_dz + dTxz_dx + ρg    # to be added:  gravity term

    # conservation of mass
    rMass = Exx[2:end,2:end-1] + Ezz[2:end-1,2:end] + (1/user_ctx.kappa) .* P[2:end-1,2:end-1];
    =#

    # ----

    #=
    # Staggered grid layout (3×4 cells): 

    Vx   P   Vx   P   Vx   P   Vx   P  
         
    o - Vz - o - Vz - o - Vz - o   Vz 
    |        |        |        |       
    Vx   P   Vx   P   Vx   P   Vx   P  
    |        |        |        |       
    o - Vz - o - Vz - o - Vz - o   Vz 
    |        |        |        |       
    Vx   P   Vx   P   Vx   P   Vx   P 
    |        |        |        |      
    o - Vz - o - Vz - o - Vz - o   Vz 
    |        |        |        |       
    Vx   P   Vx   P   Vx   P   Vx   P  
    |        |        |        |       
    o - Vz - o - Vz - o - Vz - o   Vz 
    |        |        |        |       
    Vx   P   Vx   P   Vx   P   Vx   P 
    |        |        |        |       
    o - Vz - o - Vz - o - Vz - o   Vz 
    
    Note: (1) the grid points outside the mesh are added by default to make the arrays
              regular in size and can be used to set ghost point boundary conditions.
          (2) Adding ghost points adds layers of "P" cells around it, so adding one ghostpoints
              would make the P grid be 4x3 "real" points, but "6x5" including additional points (at left and bottom).
    =#

    # collect the residuals
    for ix=corners.lower[1] : corners.upper[1]
        for iy=corners.lower[2] : corners.upper[2]
            Vx_l = [Vx[ix-1, iy-1] Vx[ix-1, iy] Vx[ix-1, iy+1];
                    Vx[ix  , iy-1] Vx[ix  , iy] Vx[ix  , iy+1];
                    Vx[ix+1, iy-1] Vx[ix+1, iy] Vx[ix+1, iy+1]];
            Vz_l = [Vz[ix-1, iy-1] Vz[ix-1, iy] Vz[ix-1, iy+1];
                    Vz[ix  , iy-1] Vz[ix  , iy] Vz[ix  , iy+1];
                    Vz[ix+1, iy-1] Vz[ix+1, iy] Vz[ix+1, iy+1]];
            P_l  = [P[ix-1, iy-1]  P[ix-1, iy];
                    P[ix  , iy-1]  P[ix,   iy] ];
            
            η_center = eta_center[ix-1 : ix+1,  iy-1 : iy+1 ]
            η_vertex = eta_vertex[ix   : ix+1,  iy   : iy+1 ]
            #η_vertex = eta_vertex[ix-1 : ix  ,  iy-1 : iy   ]
            params = (Δx=Δx, Δz=Δz, κ=user_ctx.kappa, ρ=rho_vz[ix, iy], η_center=η_center, η_vertex=η_vertex); 

            rVx_l,rVz_l,rP_l = local_residuals(Vx_l, Vz_l, P_l, params)

            # momentum residuals
            if iy>corners.lower[2] && iy<corners.upper[2]
                rVz[ix, iy] = rVz_l;
            else 
                rVz[ix, iy] = 0.0; # set BC
            end
            if ix>corners.lower[1] && ix<corners.upper[1]
                rVx[ix, iy] = rVx_l;
            else 
                rVx[ix, iy] = 0.0; # set BC
            end
            # mass conservation
            rP[ix, iy] = rP_l
        end
    end

    # restore arrays
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, user_ctx.r_l, Rlocal)
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, user_ctx.x_l, Xlocal)

    # restore coefficient array
    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dmCoeff, user_ctx.coeff_l, coeff_array)

    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dm, X_coord,Z_coord,nothing)

    # Copy local into global residual vector
    PETSc.dm_local_to_global!(user_ctx.r_l, r_g, dm)
    
    return 0
end


function FormJacobian!(J, snes, x_g, user_ctx)
    dm = PETSc.getDM(snes)

    # Extract the local vector
    PETSc.dm_global_to_local!(x_g, user_ctx.x_l, dm, PETSc.INSERT_VALUES)

    # get coordinates
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    # get the corners 
    corners = PETSc.getcorners_dmstag(dm)
    Xlocal  = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.x_l)
    
    P  = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

    Δx, Δz = user_ctx.dx, user_ctx.dz;

    # Get coefficient array for density and eta (stored in user_ctx.coeff_l)
    coeff_array = LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dmCoeff, user_ctx.coeff_l)
    rho_vz    = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_LEFT, 0)])
    iee = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_ELEMENT,   0)
    iec = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_DOWN_LEFT, 0)
    eta_center = @view(coeff_array[:,:,iee])
    eta_vertex = @view(coeff_array[:,:,iec])

    # collect the residuals 
    for ix=corners.lower[1] : corners.upper[1] 
        for iy=corners.lower[2] : corners.upper[2]
            Vx_l = [Vx[ix-1, iy-1] Vx[ix-1, iy] Vx[ix-1, iy+1];
                    Vx[ix  , iy-1] Vx[ix  , iy] Vx[ix  , iy+1];
                    Vx[ix+1, iy-1] Vx[ix+1, iy] Vx[ix+1, iy+1]];
            Vz_l = [Vz[ix-1, iy-1] Vz[ix-1, iy] Vz[ix-1, iy+1];
                    Vz[ix  , iy-1] Vz[ix  , iy] Vz[ix  , iy+1];
                    Vz[ix+1, iy-1] Vz[ix+1, iy] Vz[ix+1, iy+1]];
            P_l  = [P[ix-1, iy-1]  P[ix-1, iy];
                    P[ix  , iy-1]  P[ix,   iy] ];


            # location of all points in stencil format (following the ordering above)
            iix = ix - 1       
            iiy = iy - 1
            i_vec = [# Vx
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix-1,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix  ,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix+1,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix-1,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix  ,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix+1,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix-1,iiy+1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix  ,iiy+1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, iix+1,iiy+1,0,0),
                       
                    # Vz
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix-1,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix  ,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix+1,iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix-1,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix  ,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix+1,iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix-1,iiy+1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix  ,iiy+1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN, iix+1,iiy+1,0,0),

                    # P
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, iix-1, iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, iix  , iiy-1,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, iix-1, iiy  ,0,0),
                    LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, iix  , iiy  ,0,0)
                ]  

            # local parameters needed in the local residual routine
            η_center = eta_center[ix-1 : ix+1,  iy-1 : iy+1 ]
            η_vertex = eta_vertex[ix   : ix+1,  iy   : iy+1 ]
            #η_vertex = eta_vertex[ix-1 : ix  ,  iy-1 : iy   ]
            params = (Δx=Δx, Δz=Δz, κ=user_ctx.kappa, ρ=rho_vz[ix, iy], η_center=η_center, η_vertex=η_vertex, BC=user_ctx.BC); 

            # use AD to compute the coefficients of the local coefficients (note that in this case they are trivial)
            r_local = (x) -> local_residual_vec(x, params, (ix, iy), (corners.upper[1], corners.upper[2]))

            x_in = vcat(Vx_l[:], Vz_l[:], P_l[:])
            x_in = rand(Float64,size(x_in))
            Jder = ForwardDiff.jacobian(r_local, x_in)
            
            iP  = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT,iix,iiy,0,0)
            iVx = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,   iix,iiy,0,0)
            iVz = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,   iix,iiy,0,0)
            
            iP_coeff= PETSc.LibPETSc.DMStagStencil[]
            iVx_coeff = PETSc.LibPETSc.DMStagStencil[]
            iVz_coeff = PETSc.LibPETSc.DMStagStencil[]
            iP_val  = Float64[]
            iVx_val = Float64[]
            iVz_val = Float64[]
            for i=1:length(i_vec)
                #if Jder[3,i] != 0.0
                    push!(iP_coeff, i_vec[i])
                    push!(iP_val,   Jder[3,i])
                #end
                #if Jder[1,i] != 0.0
                    push!(iVx_coeff, i_vec[i])
                    push!(iVx_val,   Jder[1,i])
                #end
                #if Jder[2,i] != 0.0
                    push!(iVz_coeff, i_vec[i])
                    push!(iVz_val,   Jder[2,i])
                #end
            end

            if ix == 1 
                # Dirichlet BC on left boundary for Vx
                iVx_coeff = [iVx]
                iVx_val = Float64[1.0]
            end
            if iy == 1 
                # Dirichlet BC on left boundary for Vz
                iVz_coeff = [iVz]
                iVz_val = Float64[1.0]
            end

            # mass conservation   
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iP], length(iP_coeff), iP_coeff, iP_val, PETSc.INSERT_VALUES)
            
            # x-momentum
            if ix <=corners.upper[1]  
                LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], length(iVx_coeff), iVx_coeff, iVx_val, PETSc.INSERT_VALUES)
            end
            if iy <=corners.upper[2]  
                # z-momentum
                LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], length(iVz_coeff), iVz_coeff, iVz_val, PETSc.INSERT_VALUES)
            end
            
        end
    end

    # dirichlet BC's on right & top boundaries for Vx,Vz 
    for iy=corners.lower[2] : corners.upper[2] 
        ix = corners.upper[1] + 1
        iix = ix - 1
        iiy = iy - 1
        iVx = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], 1, [iVx], [1.0], PETSc.INSERT_VALUES)

        ix = corners.lower[1] 
        iVx = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], 1, [iVx], [1.0], PETSc.INSERT_VALUES)
    end

    for ix=corners.lower[1] : corners.upper[1] 
        iy = corners.upper[2] + 1
        iix = ix - 1
        iiy = iy - 1
        iVz = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], 1, [iVz], [1.0], PETSc.INSERT_VALUES)

        iy = corners.lower[2] 
        iVz = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], 1, [iVz], [1.0], PETSc.INSERT_VALUES)
    end

    PETSc.assemble!(J);

    # restore coefficient array
    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dmCoeff, user_ctx.coeff_l, coeff_array)

    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dm, X_coord,Z_coord,nothing)
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, user_ctx.x_l, Xlocal)

    # Store Julia matrix and coloring
    user_ctx.jac    =   J;
    #user_ctx.colors =   colors;

    return 0
end

# Define a struct that holds data we need in the local SNES routines below   
mutable struct Data_Stokes2D
    # DMs and vectors
    dm
    dmCoeff
    coeff_l
    x_l # local solution vector
    r_l # local residual vector
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
    # boundary conditions
    BC
    # jacobian and sparsity pattern
    jac
    jac_cache
    colors
end
user_ctx = Data_Stokes2D(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 

function PopulateCoefficientData!(user_ctx)
    # NOTE: I don't think that ghost nodes are propely set this way; to be checked (also why the arrays are larger than needed)

    ghost_corners = PETSc.getcorners_dmstag(user_ctx.dmCoeff)
    corners = PETSc.getcorners_dmstag(user_ctx.dmCoeff)
  
    user_ctx.coeff_l    =   PETSc.DMLocalVec(user_ctx.dmCoeff); # DO WE NEED TO DO Global2Local before this in parallel?
    coeff_array         =   LibPETSc.DMStagVecGetArray(petsclib,user_ctx.dmCoeff,user_ctx.coeff_l);

    # Get 1D coordinate arrays of the DM that contain the coordinates of the vertex and center points
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dmCoeff)
    @show X_coord
    @show Z_coord, typeof(Z_coord)
    @show X_coord[corners.lower[1], 1]
    
    # Get the correct entries for each of our variables in local element-wise storage
    iec = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_DOWN_LEFT, 0);  # location eta corner
    iee = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_ELEMENT,   0);  # location eta element
    irc = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_LEFT,      0);  # location rho @ vz points
    
    # Helper views
    η_center = view(coeff_array,:,:,iee);
    η_vertex = view(coeff_array,:,:,iec);
    rho_vz   = view(coeff_array,:,:,irc);

    #nx,nz = size(user_ctx.dmCoeff)[1:2]

    # set properties at centers
    for i=corners.lower[1]:corners.upper[1], j=corners.lower[2]:corners.upper[2]
        x,z = X_coord[i,2], Z_coord[j,2];  # coordinate of center points
        if GetPhase(user_ctx,x,z) == 1
            η_center[i,j] = user_ctx.eta1;
        else
            η_center[i,j] = user_ctx.eta2;
        end
    end
    #η_center[:,0] = η_center[:,1]
    #η_center[0,:] = η_center[1,:]
    #η_center[:,nz+1] = η_center[:,nz]
    #η_center[nx+1,:] = η_center[nx,:]

    # Vertexes
    for i=corners.lower[1]:corners.upper[1]+1, j=corners.lower[2]:corners.upper[2]+1
        x,z = X_coord[i,1], Z_coord[j,1];  # coordinate of vertex points
        if GetPhase(user_ctx,x,z) == 1
            η_vertex[i,j] = user_ctx.eta1;
        else
            η_vertex[i,j] = user_ctx.eta2;
        end
    end
    #η_vertex[:,0] = η_vertex[:,1]
    #η_vertex[0,:] = η_vertex[1,:]
    
    # Vz points
    for i=corners.lower[1]:corners.upper[1], j=corners.lower[2]:corners.upper[2]+1
        x,z = X_coord[i,2], Z_coord[j,1];  # coordinate of vertex points
        if GetPhase(user_ctx,x,z) == 1
            rho_vz[i,j] = user_ctx.rho1;
        else
            rho_vz[i,j] = user_ctx.rho2;  
        end
    end

    # restore arrays
    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dmCoeff, X_coord,Z_coord,nothing)
    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dmCoeff,user_ctx.coeff_l, coeff_array)
    # Note: Do not call Base.finalize on arrays obtained from PETSc; they are managed automatically

    return nothing
end

function GetPhase(ctx,x,z)

    # Returns phase 1 on left side, phase 2 on right side
    phase = 1
    if x < (ctx.xlim[2]-ctx.xlim[1])/2
    #    phase = 2
    end

    if (x-0.5)^2 + (z-0.5)^2 < 0.1^2
        phase = 2
    end


    return phase
end

function set_initial_solution!(x_g, user_ctx, petsclib)
    
    # Sets 
    LibPETSc.VecSet(petsclib, x_g, 0.0) # set global solution vector to zero (including ghost points)
    
    corners       = PETSc.getcorners_dmstag(user_ctx.dm)
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    LibPETSc.VecSet(petsclib,user_ctx.x_l, 0.0)
    X_write = LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dm, user_ctx.x_l)
    
    # add views (makes the code below more readable)
    P  = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_DOWN,   0)]);

    # Loop over P points (center)
    for ix=corners.lower[1] : corners.upper[1] 
        for iy=corners.lower[2] : corners.upper[2] 
            x,z = X_coord[ix,2], Z_coord[iy,2];         # coordinate of center points
            P[ix, iy] = π^2 * cos(π*x)*cos(π*z);  # Set DOF at element center
        end
    end

    # Vx points
    for ix=corners.lower[1] : corners.upper[1] 
        for iy=corners.lower[2] : corners.upper[2] 
            x,z = X_coord[ix,1], Z_coord[iy,2];         # coordinate of center points
            Vx[ix, iy] = π*sin(π*x)*cos(π*z);     # Set DOF at Vx points
        end
    end

    # Vz points (not including ghost points)
    for ix=corners.lower[1] : corners.upper[1] 
        for iy=corners.lower[2] : corners.upper[2] 
            x,z = X_coord[ix,2], Z_coord[iy,1];         # coordinate of center points
            
            # Set DOF at Vz points
            Vz[ix, iy] = -π*cos(π*x)*sin(π*z);
        end
    end

    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dm, user_ctx.x_l, X_write)

    #Base.finalize(X_write)  # Not needed; PETSc manages the array
    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dm, X_coord,Z_coord,nothing)

    # set ghost point values
    set_ghostpoint_values!(user_ctx.x_l, user_ctx.dm, petsclib, user_ctx.BC )

    PETSc.dm_local_to_global!(user_ctx.x_l,x_g, user_ctx.dm)

    return nothing
end

# sets ghostpoint values for a DMStag vector for free slip BCs left
function set_ghostpoint_values!(x_l, dm, petsclib, BC )
    corners       = PETSc.getcorners_dmstag(dm)
    ghost_corners = PETSc.getghostcorners_dmstag(dm)
    
    X_write = LibPETSc.DMStagVecGetArray(petsclib, dm, x_l)
    
    # add views (makes the code below more readable)
    #P  = @view(X_write[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(X_write[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(X_write[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

    # set ghost points for Vx on bottom & top boundary (free slip)
    # NOTE: this may require modifications in parallel as we only want to set ghostpoints on the physical boundary
    for ix=corners.lower[1] : corners.upper[1] 
        if BC.bottom == :free_slip
        # (τxz=0, implies that dVx/dz=0 at boundary)
            Vx[ix, ghost_corners.lower[2]] = Vx[ix, corners.lower[2]]
        end
        if BC.top == :free_slip
            Vx[ix, ghost_corners.upper[2]] = Vx[ix, corners.upper[2]]
        end
    end

    # set ghost points for Vz on left & right boundary (free slip)
    for iy=corners.lower[2] : corners.upper[2] 
        # (τxz=0, implies that dVz/dx=0 at boundary)
        if BC.left == :free_slip
            Vz[ghost_corners.lower[1], iy] = Vz[corners.lower[1], iy]
        end
        if BC.right == :free_slip
            Vz[ghost_corners.upper[1], iy] = Vz[corners.upper[1], iy]
        end
    end
    
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, x_l, X_write)
    
    return nothing
end

"""
    Vx,Vy,P,X,Xc, ρ_vz = extract_solution_julia(x_g, user_ctx)
Extracts the solution and returns them as julia arrays (without ghost points). 
"""
function extract_solution_julia(x_g::LibPETSc.PetscVec{PetscsLib}, user_ctx) where PetscsLib
    petsclib = PETSc.getlib(PetscsLib)
    PETSc.dm_global_to_local!(x_g, user_ctx.x_l, user_ctx.dm)
    
    corners     = PETSc.getcorners_dmstag(user_ctx.dm)
    X_write     = LibPETSc.DMStagVecGetArray(petsclib   , user_ctx.dm, user_ctx.x_l)
    coeff_array = LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dmCoeff, user_ctx.coeff_l)

    # add views (makes the code below more readable)
    P    = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx   = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz   = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_DOWN,   0)]);
    ρ    = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_LEFT, 0)])
    ηc   = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_ELEMENT, 0)])
    ηv   = @view(coeff_array[:,:,PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_DOWN_LEFT, 0)])
    
    nx,ny = corners.size[1:2]

    # extract solution not including ghost points
    Vx_sol = copy(Vx[1:nx+1, 1:ny  ])
    Vz_sol = copy(Vz[1:nx  , 1:ny+1])
    P_sol  = copy( P[1:nx  , 1:ny  ])
    ρ_vz   = copy( ρ[1:nx  , 1:ny+1])
    η_c    = copy( ηc[1:nx  , 1:ny  ])
    η_v    = copy( ηv[1:nx+1, 1:ny+1])

    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dm, user_ctx.x_l, X_write)

    X_coord,Y_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    x, xc = X_coord[1:nx+1,1], X_coord[1:nx,2];
    y, yc = Y_coord[1:nx+1,1], Y_coord[1:nx,2];
    
    return Vx_sol, Vz_sol, P_sol, (x,y), (xc,yc), (; η_c, η_v, ρ_vz)
end 


# Read command-line options 
cli_opts = PETSc.parse_options(ARGS)           # NamedTuple
cli = Dict{Symbol,Any}()
for (k, v) in pairs(cli_opts)
    cli[k] = v
end

# Re-pack merged options into a NamedTuple for keyword splatting
opts = (; [(k => cli[k]) for k in keys(cli)]...)


# Main Solver
nx, nz           =   64, 64;                  # number of nodes is x and z direction
#nx, nz           =   10, 10;                  # number of nodes is x and z direction

user_ctx.xlim    =   [0.0,1.0];                   # x and z dimensions
user_ctx.zlim    =   [0.0,1.0];
xlim             =   user_ctx.xlim;
zlim             =   user_ctx.zlim;
user_ctx.dx      =   (xlim[2]-xlim[1])/nx;   # x-resolution
user_ctx.dz      =   (zlim[2]-zlim[1])/nz;   # z-resolution
user_ctx.eta1    =   1;                      # viscosity phase 1
user_ctx.eta2    =   1e3;                      # viscosity phase 2
user_ctx.rho1    =   1;                      # density phase 1
user_ctx.rho2    =   2;                      # density phase 2
user_ctx.gz      =   -1;                     # gravity magnitude
user_ctx.kappa   =   1e8;                    # incompressible penalty term
BC               =   (left=:free_slip, right=:free_slip, top = :free_slip, bottom = :free_slip);  # boundary conditions
user_ctx.BC      =   BC;                     # boundary conditions

# Create Solution and coefficient DMs
comm        = MPI.COMM_WORLD
dofVertex   = 0
dofEdge     = 1
dofCenter   = 1
user_ctx.dm = PETSc.DMStag(petsclib, comm,
                      (PETSc.DM_BOUNDARY_GHOSTED, PETSc.DM_BOUNDARY_GHOSTED),
                      #(PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
                      (nx, nz),
                      (dofVertex, dofEdge, dofCenter),
                      1,
                      PETSc.DMSTAG_STENCIL_BOX;
                      opts...)

# set coordinates
PETSc.setuniformcoordinates_stag!( user_ctx.dm, (xlim[1],zlim[1]), (xlim[2],zlim[2]))

# create coefficient DM (for rho and eta)
user_ctx.dmCoeff      =   LibPETSc.DMStagCreateCompatibleDMStag(petsclib, user_ctx.dm,1,1,1,0);   # rho and eta on VERTEX, eta on ELEMENT
PETSc.setuniformcoordinates_stag!(user_ctx.dmCoeff, (xlim[1],zlim[1]), (xlim[2],zlim[2]))

# Populate phases
PopulateCoefficientData!(user_ctx);

#=
# Create solution and residual vectors
x_g             =   PETSc.DMGlobalVec(user_ctx.dm);
r_g             =   PETSc.DMGlobalVec(user_ctx.dm);
LibPETSc.VecSet(petsclib,x_g, 0.0)
LibPETSc.VecSet(petsclib,r_g, 0.0)
user_ctx.x_l    =   PETSc.DMLocalVec(user_ctx.dm);
user_ctx.r_l    =   PETSc.DMLocalVec(user_ctx.dm);
LibPETSc.VecSet(petsclib,user_ctx.x_l, 0.0)
LibPETSc.VecSet(petsclib,user_ctx.r_l, 0.0)

# Create PETSc matrix for Jacobian
J               =   LibPETSc.DMCreateMatrix(petsclib,user_ctx.dm);

# Setting up SNES
snes = PETSc.SNES(petsclib,comm; 
        snes_rtol=1e-12, 
        snes_monitor=true, 
        snes_max_it = 500,
        #snes_type = "ksponly",
        snes_monitor_true_residual=true, 
        snes_view=true,
        snes_linesearch_monitor=false,
        #snes_linesearch_view=true,
        #snes_linesearch_type="basic",
        ksp_type = "preonly",
        pc_type = "lu",
        snes_maxit=10,
        snes_converged_reason=true, help=false,
        opts...);

snes.user_ctx  =       user_ctx;       # crashes

PETSc.setDM!(snes, user_ctx.dm)

# Set first guess values for solution vector
set_initial_solution!(x_g, user_ctx, petsclib)

PETSc.setfunction!(snes, FormRes!, r_g)
PETSc.setjacobian!(snes, FormJacobian!, J, J)

PETSc.solve!(x_g, snes);

# Extract solution as well as coordinate vectors
Vx, Vz, P, X, Xc, material = extract_solution_julia(x_g, user_ctx)

# Copy velocity to center points for plotting
Vx_c = (Vx[2:end,:] + Vx[1:end-1,:])/2;
Vz_c = (Vz[:,2:end] + Vz[:,1:end-1])/2;


# plot
fig = Figure()
ax = Axis(fig[1,1], title="Vz field", xlabel="x", ylabel="z")
hm = heatmap!(ax, X[1], Xc[2], Vz)
Colorbar(fig[1,2], hm)

# Subsample for clearer visualization
step = 2  # plot every 2nd point
arrows2d!(ax, Xc[1][1:step:end], Xc[2][1:step:end], 
        Vx_c[1:step:end, 1:step:end], Vz_c[1:step:end, 1:step:end])
        #tipwidth=20, tiplength=150, color=:black)

# overlay a single isocontour at value 0.5 (between rho1=0 and rho2=1)
contour!(ax, Xc[1], X[2], material.ρ_vz; levels=[1.5], color=:white, linewidth=2)
display(fig)
=#