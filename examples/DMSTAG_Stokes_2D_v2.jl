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
using ForwardDiff, SparseArrays, LinearAlgebra

petsclib = first(PETSc.petsclibs);
PETSc.initialized(petsclib) || PETSc.initialize(petsclib)


function residual_mass(x_vec, Δx, Δz, κ)
    Vx = x_vec[1:2]
    Vz = x_vec[3:4]
    P  = x_vec[5]

    rMass = (Vx[2] - Vx[1]) / Δx + (Vz[2] - Vz[1]) / Δz + (1/κ) * P;
    
    return rMass
end

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
    Eii = sqrt(Exx_l[2,2]^2 + Ezz_l[2,2]^2 + 2*sum(Exz_l.^2)/4);
            
    # once we have that we can compute Tii
    Tii = 2*params.ηl * Eii
    η   = Tii / (2*Eii + 1e-16)   # avoid division by zero

    # stresses, assuming linear viscosity
    Txx_l = 2η .* Exx_l;        # Txx = 2*eta*Exx
    Tzz_l = 2η .* Ezz_l;        # Tzz = 2*eta*Ezz
    Txz_l = 2η .* Exz_l;        # Txz = 2*eta*Exz

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
    rVz = -dPdz_l[1] + dTzzdz_l[1] + dTxzdx_l[1] + params.ρ;
            
    # mass conservation
    rP = Exx_l[2,2] + Ezz_l[2,2] + (1/params.κ) .* P_l[2,2];

    return rVx, rVz, rP
end

"""
    res = local_residual_vec(x_in::Vector, params::NamedTuple) 
Residual function that calls `local_residuals` and can be used with ForwardDiff
"""
function local_residual_vec(x_in::Vector, params)
    nVx  = (3,3)
    nVz  = (3,3)
    nP   = (2,2)
    Vx_l = reshape(x_in[1:9],   nVx)
    Vz_l = reshape(x_in[10:18], nVz)
    P_l  = reshape(x_in[19:22], nP )
                
    rVx,rVz,rP = local_residuals(Vx_l, Vz_l, P_l, params)

    return [rVx,rVz,rP]
end

function FormRes!(r_g, snes, x_g, user_ctx)
    #dm = PETSc.getDM(snes)
    dm = user_ctx.dm

    # Copy global to local vectors
    LibPETSc.VecSet(petsclib, user_ctx.r_l, 0.0) # set residual to zero before accumulating contributions
    LibPETSc.VecSet(petsclib, r_g, 0.0) # set residual to zero before accumulating contributions
    
    LibPETSc.VecSet(petsclib, user_ctx.x_l, 0.0) # set solution to zero before accumulating contributions
    PETSc.dm_global_to_local!(x_g, user_ctx.x_l, dm)
    #PETSc.dm_global_to_local!(r_g,user_ctx.r_l, dm)  

    # get coordinates
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    # get the corners 
    corners       = PETSc.getcorners_dmstag(dm)
    ghost_corners = PETSc.getghostcorners_dmstag(dm)
    shift         = corners.lower - ghost_corners.lower;
    #@show corners

    LibPETSc.VecSet(petsclib,user_ctx.r_l, 0.0) # set residual to zero before accumulating contributions
    Xlocal = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.x_l)
    Rlocal = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.r_l)
    
    P  = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);
    
    @show P Vx Vz
    
    #Rlocal[:,:, :] .= 0.0
    rP  = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    rVx = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    rVz = @view(Rlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

   # rVx.= 0.0
   # rVz.= 0.0
   # rP.= 0.0
    
    Δx, Δz = user_ctx.dx, user_ctx.dz;
    
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

    # first momentum balance
    dP_dx = diff(P[2:end-1,2:end-1],dims=1) ./ Δx;
    dTxx_dx = diff(Txx[2:end,2:end-1],dims=1) ./ Δx;
    dTxz_dz = diff(Txz[2:end-1,1:end],dims=2) ./ Δz;
    rM1 = -dP_dx + dTxx_dx + dTxz_dz        

    # 2nd momentum balance
    dP_dz   = diff(P[2:end-1,2:end-1],dims=2) ./ Δz;
    dTzz_dz = diff(Tzz[2:end-1,2:end],dims=2) ./ Δz;
    dTxz_dx = diff(Txz[1:end,2:end-1],dims=1) ./ Δx;
    rM2 = -dP_dz + dTzz_dz + dTxz_dx     # to be added:  gravity term

    # conservation of mass
    rMass = Exx[2:end,2:end-1] + Ezz[2:end-1,2:end] + (1/user_ctx.kappa) .* P[2:end-1,2:end-1];
=#
    # ----

    #=
    # Staggered grid layout (4x4 cells):

    Vx   P   Vx   P   Vx   P   Vx   P   Vx  P
         
    o - Vz - o - Vz - o - Vz - o - Vz - o   Vz
    |        |        |        |        |
    Vx   P   Vx   P   Vx   P   Vx   P   Vx  P
    |        |        |        |        |
    o - Vz - o - Vz - o - Vz - o - Vz - o   Vz
    |        |        |        |        |
    Vx   P   Vx   P   Vx   P   Vx   P   Vx  P
    |        |        |        |        |
    o - Vz - o - Vz - o - Vz - o - Vz - o   Vz
    |        |        |        |        |
    Vx   P   Vx   P   Vx   P   Vx   P   Vx  P
    |        |        |        |        |
    o - Vz - o - Vz - o - Vz - o - Vz - o   Vz
    |        |        |        |        |
    Vx   P   Vx   P   Vx   P   Vx   P   Vx  P  
    |        |        |        |        |
    o - Vz - o - Vz - o - Vz - o - Vz - o   Vz
    
    Note: (1) the grid points outside the mesh are added by default to make the arrays
              regular in size and can be used to set ghost point boundary conditions.
          (2) Adding ghost points adds layers of "P" cells around it, so adding one ghostpoints
              would make the P grid be 4x4 "real" points, but "6x6" including additional points (at left and bottom).
    =#

    # create an anonymous function for residuals
    #r_local = (x) -> local_residual_vec(x, params)

    # collect the residuals
    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        for iy=corners.lower[2]+shift[2] : corners.upper[2] + shift[2]
            @show ix,iy
            Vx_l = [Vx[ix-1, iy-1] Vx[ix-1, iy] Vx[ix-1, iy+1];
                    Vx[ix  , iy-1] Vx[ix  , iy] Vx[ix  , iy+1];
                    Vx[ix+1, iy-1] Vx[ix+1, iy] Vx[ix+1, iy+1]];
            Vz_l = [Vz[ix-1, iy-1] Vz[ix-1, iy] Vz[ix-1, iy+1];
                    Vz[ix  , iy-1] Vz[ix  , iy] Vz[ix  , iy+1];
                    Vz[ix+1, iy-1] Vz[ix+1, iy] Vz[ix+1, iy+1]];
            P_l  = [P[ix-1, iy-1]  P[ix-1, iy];
                    P[ix  , iy-1]  P[ix,   iy] ];

            # local parameters needed in the local residual routine
            # (if parameters )

            x,z = X_coord[ix-shift[1],1], Z_coord[iy-shift[2],2];         # coordinate of Vz points
            if GetPhase(user_ctx,x,1)
                ρ = user_ctx.rho1;
            else
                ρ = user_ctx.rho2;  
            end
            @show x z ρ
            
            params = (Δx=Δx, Δz=Δz, κ=user_ctx.kappa, ηl=user_ctx.eta1, ρ=ρ); 

            rVx_l,rVz_l,rP_l = local_residuals(Vx_l, Vz_l, P_l, params)
            
            #=
            # use AD to compute the coefficients of the local coefficients (note that in this case they are trivial)
            x_in = vcat(Vx_l[:], Vz_l[:], P_l[:])
            Jder = ForwardDiff.jacobian(r_local, x_in)

            #@show Jder

            =#
            # momentum residuals
            if iy>corners.lower[2]+shift[2] 
                rVz[ix, iy] = rVz_l;
            end
            if ix>corners.lower[1]+shift[1]
                rVx[ix, iy] = rVx_l;
            end
            # mass conservation
            rP[ix, iy]=rP_l
            @show rP_l, rVx_l, rVz_l
        end
    end
    #=
    =#

   # @show norm(Rlocal[:])
   # @show rP rVx rVz
    # restore arrays
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, user_ctx.r_l, Rlocal)
    LibPETSc.DMStagVecRestoreArray(petsclib, dm, user_ctx.x_l, Xlocal)
    
   # Base.finalize(Xlocal)
   # Base.finalize(Rlocal)

    #     
    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dm, X_coord,Z_coord,nothing)


    # Copy local into global residual vector
    PETSc.dm_local_to_global!(user_ctx.r_l, r_g, dm)
    
    return 0
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

function  func(out, x, user_ctx)

    ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);          # array with all local x-data (solition array)
    ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, out);        # array with all local residual
    
    ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

    return nothing
end

function FormJacobian!(J, snes, x_g, user_ctx)
    #dm = PETSc.getDMDA(snes)
    dm = user_ctx.dm

    # Extract the local vector
    PETSc.dm_global_to_local!(x_g, user_ctx.x_l, dm, PETSc.INSERT_VALUES)

    # Copy global to local vectors
    #PETSc.dm_global_to_local!(x_g,user_ctx.x_l, dm)
    
    # get the corners 
    corners       = PETSc.getcorners_dmstag(dm)
    ghost_corners = PETSc.getghostcorners_dmstag(dm)
    shift         = corners.lower - ghost_corners.lower;
    
    Xlocal = LibPETSc.DMStagVecGetArray(petsclib, dm, user_ctx.x_l)
    
    P  = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(Xlocal[:,:,PETSc.DMStagDOF_Slot(dm, LibPETSc.DMSTAG_DOWN,   0)]);

    Δx, Δz = user_ctx.dx, user_ctx.dz;

    # collect the residuals 
    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        for iy=corners.lower[2]+shift[2] : corners.upper[2] + shift[2]
            Vx_l = [Vx[ix-1, iy-1] Vx[ix-1, iy] Vx[ix-1, iy+1];
                    Vx[ix  , iy-1] Vx[ix  , iy] Vx[ix  , iy+1];
                    Vx[ix+1, iy-1] Vx[ix+1, iy] Vx[ix+1, iy+1]];
            Vz_l = [Vz[ix-1, iy-1] Vz[ix-1, iy] Vz[ix-1, iy+1];
                    Vz[ix  , iy-1] Vz[ix  , iy] Vz[ix  , iy+1];
                    Vz[ix+1, iy-1] Vz[ix+1, iy] Vz[ix+1, iy+1]];
            P_l  = [P[ix-1, iy-1]  P[ix-1, iy];
                    P[ix  , iy-1]  P[ix,   iy] ];

            # location of all points in stencil format (following the ordering above)
            iix = ix - shift[1] - 1       
            iiy = iy - shift[2] - 1
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
            # (if parameters )
            params = (Δx=Δx, Δz=Δz, κ=user_ctx.kappa, ηl=user_ctx.eta1); 
            
            # use AD to compute the coefficients of the local coefficients (note that in this case they are trivial)
            r_local = (x) -> local_residual_vec(x, params)

            x_in = vcat(Vx_l[:], Vz_l[:], P_l[:])
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
                if Jder[3,i] != 0.0
                    push!(iP_coeff, i_vec[i])
                    push!(iP_val,   Jder[3,i])
                end
                if Jder[1,i] != 0.0
                    push!(iVx_coeff, i_vec[i])
                    push!(iVx_val,   Jder[1,i])
                end
                if Jder[2,i] != 0.0
                    push!(iVz_coeff, i_vec[i])
                    push!(iVz_val,   Jder[2,i])
                end
            end

            if ix-shift[1] == 1 
                # Dirichlet BC on left boundary for Vx
                iVx_coeff = [iVx]
                iVx_val = Float64[1.0]
            end
            if iy-shift[2] == 1 
                # Dirichlet BC on left boundary for Vz
                iVz_coeff = [iVz]
                iVz_val = Float64[1.0]
            end

            # mass conservation   
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iP], length(iP_coeff), iP_coeff, iP_val, PETSc.INSERT_VALUES)
            
            # x-momentum
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], length(iVx_coeff), iVx_coeff, iVx_val, PETSc.INSERT_VALUES)
            
            # z-momentum
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], length(iVz_coeff), iVz_coeff, iVz_val, PETSc.INSERT_VALUES)
            
        end
    end

    # dirichlet BC's on right & top boundaries for Vx,Vz 
    for iy=corners.lower[2]+shift[2] : corners.upper[2] + shift[2]
        ix = corners.upper[1] + shift[1] + 1
        iix = ix - shift[1] -1
        iiy = iy - shift[2] -1
        iVx = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], 1, [iVx], [1.0], PETSc.INSERT_VALUES)

        ix = corners.lower[1] + shift[1] 
        iVx = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVx], 1, [iVx], [1.0], PETSc.INSERT_VALUES)
    end

    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        iy = corners.upper[2] + shift[2] + 1
        iix = ix - shift[1] - 1
        iiy = iy - shift[2] - 1
        iVz = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], 1, [iVz], [1.0], PETSc.INSERT_VALUES)

        iy = corners.lower[2] + shift[2] 
        iVz = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,   iix,iiy,0,0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, J, 1, [iVz], 1, [iVz], [1.0], PETSc.INSERT_VALUES)
    end
    
    PETSc.assemble!(J);
  
    # Store Julia matrix and coloring
    #user_ctx.jac    =   J;
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
    # jacobian and sparsity pattern
    jac
    jac_cache
    colors
end
user_ctx = Data_Stokes2D(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing);  # holds data we need in the local 

# TO BE REMOVED
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
    ind                 =   PETSc.DMStagGetIndices(dm);          # indices of (center/element) points, not including ghost values.
    #sx, sn  = PETSc.DMStagGetCentralNodes(user_ctx.dm); #indices of central points

    ix      =   ind.center[1];           
    iz      =   ind.center[2];

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
    # Note: Do not call Base.finalize on arrays obtained from PETSc DMStagVecGetArrayLocationSlot
    # as they are managed by PETSc. Only restore them with DMStagVecRestoreArray.                                    
end

# TO BE REMOVED
function ComputeStresses!(user_ctx, ArrayLocal_x)
    
    # Compute strain rates
    Exx,Ezz,Exz = ComputeStrainRates!(user_ctx, ArrayLocal_x);

    # Getting Eta at the center and corner of the cells
    EtaE    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_ELEMENT, 0);
    EtaC    = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff,user_ctx.coeff_l,     PETSc.DMSTAG_DOWN_LEFT, 0);
    
    # Getting indices for center nodes (not ghost)
    ind     = PETSc.DMStagGetIndices(user_ctx.dm); #indices of central points
    #ix      =   sx[1]:sn[1];           
    #iz      =   sx[2]:sn[2];
    ix  = ind.center[1];
    iz  = ind.center[2];
    
    # Compute shear stresses
    Txx    = 2 .* EtaE[ix,iz] .* Exx;
    Tzz    = 2 .* EtaE[ix,iz] .* Ezz;
    Txz    = 2 .* EtaC[ix[1]:ix[end]+1,iz[1]:iz[end]+1] .* Exz;

    return Txx,Tzz,Txz

    # Note: Do not finalize arrays obtained from PETSc  
    
end

# TO BE REMOVED
function ComputeStrainRates!(user_ctx, ArrayLocal_x)

    # Getting velocity vectors
    Vx      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_LEFT, 0);
    Vz      = PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,ArrayLocal_x,     PETSc.DMSTAG_DOWN, 0);
    
    # dimensions
    dx,dz   = user_ctx.dx, user_ctx.dz;

    # Getting indices for center nodes (not ghost)
    ind     = PETSc.DMStagGetIndices(user_ctx.dm); #indices of central points
    #ix      =   sx[1]:sn[1];           
    #iz      =   sx[2]:sn[2];
    ix      =   ind.center[1];
    iz      =   ind.center[2];
    

    # Set ghost points (free slip)
    Vx[:,iz[1]-1]   .=  Vx[:,iz[1]];                # ghost points on left and right size (here: free slip)
    Vx[:,iz[end]+1] .=  Vx[:,iz[end]];
    Vz[ix[1]-1,:]   .=  Vz[ix[1],:];
    Vz[ix[end]+1,:] .=  Vz[ix[end],:];

    # Compute deviatoric Strain rates
    DivV  = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .+ (Vz[ix,iz.+1].-Vz[ix,iz])./dz;
    Exx   = (Vx[ix.+1,iz].-Vx[ix,iz])./dx .- 1/3 .* DivV;
    Ezz   = (Vz[ix,iz.+1].-Vz[ix,iz])./dz .- 1/3 .* DivV;
    Exz   = 0.5.*((Vx[ix[1]:ix[end].+1,iz[1]:iz[end].+1].-Vx[ix[1]:ix[end].+1,iz[1].-1:iz[end]])./dz .+
                  (Vz[ix[1]:ix[end].+1,iz[1]:iz[end].+1].-Vz[ix[1].-1:ix[end],iz[1]:iz[end].+1])./dx);
    
    return Exx,Ezz,Exz

end

function PopulateCoefficientData!(user_ctx)
    # NOTE: I don't think that ghost nodes are propely set this way; to be checked (also why the arrays are larger than needed)

    # Create coefficient local vector and array (to store eta and rho)

    user_ctx.coeff_l    =   PETSc.DMLocalVec(user_ctx.dmCoeff); # DO WE NEED TO DO Global2Local before this in parallel?
    coeff_array         =   LibPETSc.DMStagVecGetArray(petsclib,user_ctx.dmCoeff,user_ctx.coeff_l);

    # Get 1D coordinate arrays of the DM that contain the coordinates of the vertex and center points
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dmCoeff)

    # Get the correct entries for each of our variables in local element-wise storage
    iec = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_DOWN_LEFT, 0);  # location eta corner
    iee = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_ELEMENT,   0);  # location eta element
    irc = PETSc.DMStagDOF_Slot(user_ctx.dmCoeff, LibPETSc.DMSTAG_LEFT,      0);  # location rho @ vz points
    
    # Helper views
    η_center = view(coeff_array,:,:,iee);
    η_vertex = view(coeff_array,:,:,iec);
    ρ_vz     = view(coeff_array,:,:,irc);

    nx,nz = size(user_ctx.dmCoeff)[1:2]

    # set properties at centers
    for i=1:nx, j=1:nz
        x,z = X_coord[i,2], Z_coord[j,2];  # coordinate of center points
        if GetPhase(user_ctx,x,1)
            η_center[i,j] = user_ctx.eta1;
        else
            η_center[i,j] = user_ctx.eta2;
        end
    end
    
    # Vertexes
    for i=1:nx+1, j=1:nz+1
        x,z = X_coord[i,1], Z_coord[j,1];  # coordinate of vertex points
        if GetPhase(user_ctx,x,1)
            η_vertex[i,j] = user_ctx.eta1;
#            ρ_vertex[i,j] = user_ctx.rho1;
        else
            η_vertex[i,j] = user_ctx.eta2;
 #           ρ_vertex[i,j] = user_ctx.rho2;  
        end
    end

    # Vz points
    for i=1:nx, j=1:nz+1
        x,z = X_coord[i,2], Z_coord[j,1];  # coordinate of vertex points
        if GetPhase(user_ctx,x,1)
            ρ_vz[i,j] = user_ctx.rho1;
        else
            ρ_vz[i,j] = user_ctx.rho2;  
        end
    end


    # restore arrays
    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dmCoeff, X_coord,Z_coord,nothing)
    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dmCoeff,user_ctx.coeff_l, coeff_array)
    # Note: Do not call Base.finalize on arrays obtained from PETSc; they are managed automatically
        
    return nothing
end

function GetPhase(ctx,x,n)
    # Divide domain in phases according to x-coordinate (similar to SolCx benchmark)
    if x < (ctx.xlim[2]-ctx.xlim[1])/2
        if n == 1 return true else return false end
    else
        if n == 1 return false else return true end
    end
end

#=
# TO BE REMOVED
function SetVecX!(user_ctx,x_g)
    # Set first guess for vector x to 0
    user_ctx.x_l .= 0; 
    PETSc.dm_local_to_global!(user_ctx.x_l,x_g, user_ctx.dm, PETSc.INSERT_VALUES)
    return nothing
end
=#


function set_initial_solution!(x_g, user_ctx, petsclib)
    # Sets 
    LibPETSc.VecSet(petsclib, x_g, 0.0) # set global solution vector to zero (including ghost points)
    
    corners       = PETSc.getcorners_dmstag(user_ctx.dm)
    ghost_corners = PETSc.getghostcorners_dmstag(user_ctx.dm)
    shift         = corners.lower - ghost_corners.lower;
    X_coord,Z_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, user_ctx.dm)

    LibPETSc.VecSet(petsclib,user_ctx.x_l, 0.0)
    X_write = LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dm, user_ctx.x_l)
    
    # add views (makes the code below more readable)
    P  = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_ELEMENT,0)]);
    Vx = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_LEFT,   0)]);
    Vz = @view(X_write[:,:,PETSc.DMStagDOF_Slot(user_ctx.dm, LibPETSc.DMSTAG_DOWN,   0)]);

    Vx .= 0.0
    P .= 0.0
    Vz .= 0.0
    @show "before setting initial" Vx Vz P

    # Loop over P points (center)
    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        for iy=corners.lower[2]+ shift[2] : corners.upper[2] + shift[2]
            x,z = X_coord[ix-shift[1],2], Z_coord[iy-shift[2],2];         # coordinate of center points
            P[ix, iy] = π^2 * cos(π*x)*cos(π*z);  # Set DOF at element center
        end
    end
    # No need to set ghost points for P

    # Vx points
    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        for iy=corners.lower[2]+ shift[2] : corners.upper[2] + shift[2] 
            x,z = X_coord[ix-shift[1],1], Z_coord[iy-shift[2],2];         # coordinate of center points
            Vx[ix, iy] = π*sin(π*x)*cos(π*z);     # Set DOF at Vx points
        end
    end

    # set ghost points for Vx on bottom & top boundary (free slip)
    # NOTE: this may require modfications in parallel as we only want to set ghostpoints on the physical boundary
    for ix=corners.lower[1]+ shift[1] : corners.upper[1] + shift[1]
        # (τxz=0, implies that dVx/dz=0 at boundary)
        Vx[ix, ghost_corners.lower[2]+shift[2]] = Vx[ix, corners.lower[2]+shift[2]]
        Vx[ix, ghost_corners.upper[2]+shift[2]] = Vx[ix, corners.upper[2]+shift[2]]
    end

    # Vz points (not including ghost points)
    for ix=corners.lower[1]+shift[1] : corners.upper[1] + shift[1]
        for iy=corners.lower[2]+ shift[2] : corners.upper[2] + shift[2]
            x,z = X_coord[ix-shift[1],2], Z_coord[iy-shift[2],1];         # coordinate of center points
            
            # Set DOF at Vz points
            Vz[ix, iy] = -π*cos(π*x)*sin(π*z);
        end
    end

    # set ghost points for Vz on left & right boundary (free slip)
    for iy=corners.lower[2]+ shift[2] : corners.upper[2] + shift[2]
        # (τxz=0, implies that dVz/dx=0 at boundary)
        Vz[ghost_corners.lower[1]+shift[1], iy] = Vz[corners.lower[1]+shift[1], iy]
        Vz[ghost_corners.upper[1]+shift[1], iy] = Vz[corners.upper[1]+shift[1], iy]
    end

    @show "after setting initial" Vx Vz P

    LibPETSc.DMStagVecRestoreArray(petsclib, user_ctx.dm, user_ctx.x_l, X_write)


    #Base.finalize(X_write)  # Not needed; PETSc manages the array
    LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, user_ctx.dm, X_coord,Z_coord,nothing)

    PETSc.dm_local_to_global!(user_ctx.x_l,x_g, user_ctx.dm)

    return nothing
end



# Main Solver
#nx, nz           =   32, 32;                  # number of nodes is x and z direction
nx, nz           =   4, 4;                  # number of nodes is x and z direction

user_ctx.xlim    =   [0.0,2.0];                   # x and z dimensions
user_ctx.zlim    =   [0.0,1.0];
xlim             =   user_ctx.xlim;
zlim             =   user_ctx.zlim;
user_ctx.dx      =   (xlim[2]-xlim[1])/nx;   # x-resolution
user_ctx.dz      =   (zlim[2]-zlim[1])/nz;   # z-resolution
user_ctx.eta1    =   1;                      # viscosity phase 1
user_ctx.eta2    =   2;                      # viscosity phase 2
user_ctx.rho1    =   3;                      # density phase 1
user_ctx.rho2    =   4;                      # density phase 2
user_ctx.gz      =   -1;                     # gravity magnitude
user_ctx.kappa   =   1e2;                    # incompressible penalty term

# Create Solution and coefficient DMs
comm        = MPI.COMM_SELF
dofVertex   = 0
dofEdge     = 1
dofCenter   = 1
user_ctx.dm = PETSc.DMStag(petsclib, comm,
                      (PETSc.DM_BOUNDARY_GHOSTED, PETSc.DM_BOUNDARY_GHOSTED),
                      #(PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
                      (nx, nz),
                      (dofVertex, dofEdge, dofCenter),
                      1,
                      PETSc.DMSTAG_STENCIL_BOX)

# set coordinates
PETSc.setuniformcoordinates_stag!( user_ctx.dm, (xlim[1],zlim[1]), (xlim[2],zlim[2]))

# create coefficient DM (for rho and eta)
user_ctx.dmCoeff      =   LibPETSc.DMStagCreateCompatibleDMStag(petsclib, user_ctx.dm,1,1,1,0);   # rho and eta on VERTEX, eta on ELEMENT
PETSc.setuniformcoordinates_stag!(user_ctx.dmCoeff, (xlim[1],zlim[1]), (xlim[2],zlim[2]))

# Populate phases
PopulateCoefficientData!(user_ctx);

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


# Compute jacobian sparsity pattern
#x_l     =   PETSc.unsafe_localarray(Float64, user_ctx.x_l.ptr;  write=false, read=true)

#user_ctx.jac, user_ctx.colors = ComputeSparsityPatternJacobian_automatic(user_ctx.x_l, user_ctx);
#ind     =   PETSc.LocalInGlobalIndices(user_ctx.dm);    # extract indices
#J       =   PETSc.MatSeqAIJ(sparse(user_ctx.jac[ind,ind]));

# Setting up SNES
snes = PETSc.SNES(petsclib,MPI.COMM_SELF; 
        snes_rtol=1e-12, 
        snes_monitor=true, 
        snes_max_it = 500,
        snes_monitor_true_residual=true, 
        snes_view=true,
        snes_linesearch_monitor=true,
        snes_linesearch_view=true,
        snes_linesearch_type="basic",
        snes_maxit=10,
        snes_converged_reason=true, help=false);

snes.user_ctx  =       user_ctx;       # crashes

PETSc.setDM!(snes, user_ctx.dm)

# Set first guess values for solution vector
set_initial_solution!(x_g, user_ctx, petsclib)

PETSc.setfunction!(snes, FormRes!, r_g)
PETSc.setjacobian!(snes, FormJacobian!, J, J)
FormRes!(r_g, snes, x_g, user_ctx)

#Solve
PETSc.solve!(x_g, snes);

#=
# Copy solution to local vector
PETSc.update!(user_ctx.x_l, x_g,  PETSc.INSERT_VALUES);

# Extract solution
Vx  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_LEFT,      0); 
Vz  =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_DOWN,      0); 
P   =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm, user_ctx.x_l,   PETSc.DMSTAG_ELEMENT,   0); 
Eta =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dmCoeff, user_ctx.coeff_l,   PETSc.DMSTAG_ELEMENT,   0); 

# Getting indices for center nodes (not ghost)
ind     = PETSc.DMStagGetIndices(user_ctx.dm);
ix      =   ind.center[1];           
iz      =   ind.center[2];       

# Getting rid of ghost points
Vx  =   Vx[ix[1]:ix[end]+1,iz];
Vz  =   Vz[ix,iz[1]:iz[end]+1];
P   =    P[ix,iz];
Eta =  Eta[ix,iz];

# Compute central average velocities
Vx_cen = (Vx[2:end,:]  + Vx[1:end-1,:])/2;
Vz_cen = (Vz[:,2:end]  + Vz[:,1:end-1])/2;

#get coordinates
dm_coord  = PETSc.getcoordinateDM(user_ctx.dmCoeff);
vec_coord = PETSc.getcoordinateslocal(user_ctx.dmCoeff);

XCoord_e = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_ELEMENT  ,      0); # location coord corner
ZCoord_e = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_ELEMENT  ,      1); # location coord corner
XCoord_c = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_DOWN_LEFT,      0);   # location coord element
ZCoord_c = PETSc.DMStagGetGhostArrayLocationSlot(dm_coord, vec_coord,   PETSc.DMSTAG_DOWN_LEFT,      1);   # location coord element

XCoord_e = XCoord_e[ix,iz];
ZCoord_e = ZCoord_e[ix,iz];
XCoord_c = XCoord_c[ix[1]:ix[end]+1,iz[1]:iz[end]+1];
ZCoord_c = ZCoord_c[ix[1]:ix[end]+1,iz[1]:iz[end]+1];

xe_1D    = XCoord_e[:,1];
ze_1D    = ZCoord_e[1,:];
xc_1D    = XCoord_c[:,1];
zc_1D    = ZCoord_c[1,:];

# Plot
CreatePlots = false;        # set to false by default such that testing does not require Plots
if CreatePlots==true
    using Plots

    p1 = heatmap(xe_1D,ze_1D, P', xlabel="Width", ylabel="Depth", title="Pressure",aspect_ratio = 1)
    p2 = heatmap(xe_1D,zc_1D, Vz', xlabel="Width", ylabel="Depth", title="Vz",aspect_ratio = 1)
    p3 = heatmap(xc_1D,ze_1D, Vx', xlabel="Width", ylabel="Depth", title="Vx",aspect_ratio = 1)
    quiver!(XCoord_e[:],ZCoord_e[:],quiver=(Vx_cen[:]*0.02,Vz_cen[:]*0.02), color=:white,aspect_ratio = 1)
    p4 = heatmap(xe_1D,ze_1D, Eta', xlabel="Width", ylabel="Depth", title="Eta",aspect_ratio = 1)
    plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

end

# check if the solution makes sense
using LinearAlgebra, Test
check =  norm(user_ctx.x_l)
@test check ≈ 31.54688642064483 

=#
