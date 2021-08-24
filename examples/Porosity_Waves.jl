# EXCLUDE FROM TESTING
#=
In this example we solve the 1D viscoelastic porosity wave equations which describe how magma ascends in the Earth's mantle. 
The equations are discussed in [Vasyliev et al. Geophysical Research Letters (25), 17. p. 3239-3242](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/98GL52358):
```math
   {∂ϕ/∂t}  = -(De {∂Pe/∂t} + ϕ^m Pe) \\
De {∂Pe/∂t} =  ∇( ϕ^n (∇Pe + e_z) ) - ϕ^m Pe 
```
with zero Dirichlet boundary conditions for Pe one Dirichlet conditions for ϕ=1 on either side.

To solve the problem we use the standard central, finite difference approximation in 1D, 
but we simultaneously solve for Pe and ϕ on a collocated grid.
=#

using MPI
using PETSc
using OffsetArrays: OffsetArray
using LinearAlgebra: norm
using ForwardDiff
using SparseArrays, SparseDiffTools, ForwardDiff


opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (ksp_monitor = false, 
    ksp_rtol=1e-10,
    ksp_converged_reason = false, 
    snes_monitor = true, 
    snes_converged_reason=true, 
    snes_view=false, 
    snes_max_it=100, 
    snes_max_funcs=10000, 
    snes_max_linear_solve_fail=1000,
    snes_mf_operator=false,
    )
   
end

CreatePlots = false;
if isinteractive()
    CreatePlots=true
    using Plots
end



# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64
PetscInt    = Int64


# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

# dimensionality of the problem
dim = haskey(opts, :dim) ? opts.dim : 1

# Set the total number of grid points in each direction
Nq = ntuple(_ -> 501, dim)

# Set the boundary conditions on each side
bcs = (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_GHOSTED, PETSc.DM_BOUNDARY_GHOSTED)
bcs = bcs[1:dim];

# Set parameters
n = PetscScalar(3)
m = PetscScalar(2)
De = PetscScalar(1e2)
dt = PetscScalar(1e-2)


# Create the PETSC dmda object
da = PETSc.DMDA(
    petsclib,
    comm,
    bcs,                     # boundary conditions
    Nq,                      # Global grid size
    2,                       # Number of DOF per node
    1,                       # Stencil width
    PETSc.DMDA_STENCIL_STAR; # Stencil type
    opts...,
)

# Set coordinates of the DMDA (2D and 3D ignored for 1D DMDA)
PETSc.setuniformcoordinates!(da, (-20, -10, -10), (150, 10, 10))

# Create the PETSC snes object
snes = PETSc.SNES(petsclib, comm; opts...)

# add the da to the snes
PETSc.setDM!(snes, da)

# Sets the initial profiles
g_x = PETSc.DMGlobalVec(da)
PETSc.withlocalarray!(g_x; read = false) do l_x   
    corners = PETSc.getcorners(da)

    Phi = PETSc.getlocalarraydof(da, l_x, dof = 1)
    Pe = PETSc.getlocalarraydof(da, l_x, dof = 2)
    Coord = PETSc.getlocalcoordinatearray(da)     # retrieve arrays with local coordinates

    Phi0 = PetscScalar(1.0)
    dPhi1 = PetscScalar(8.0)
    dPhi2 = PetscScalar(1.0)
    z1 = PetscScalar(0.0)
    z2 = PetscScalar(40.0)
    lambda = PetscScalar(1.0)
    dPe1 = dPhi1 / De
    dPe2 = dPhi2 / De

    # Loop over all the points on the processor and set the initial condition
    for i in ((corners.lower):(corners.upper))
        Phi[i] =
            Phi0 +
            dPhi1 * exp(-((Coord.X[i] - z1)^2.0) / lambda^2) +
            dPhi2 * exp(-((Coord.X[i] - z2)^2.0) / lambda^2)
        Pe[i] =
            -dPe1 * exp(-((Coord.X[i] - z1)^2.0) / lambda^2) -
            dPe2 * exp(-((Coord.X[i] - z2)^2.0) / lambda^2)
    end
end

# Initialize x_old on every processor
g_xold = deepcopy(g_x); # initialize Pe and Phi of last timestep 
l_xold = PETSc.DMLocalVec(da)
PETSc.update!(l_xold, g_xold, PETSc.INSERT_VALUES)
x_old = PETSc.unsafe_localarray(l_xold,read=true, write=false)


# Routine wuth julia-only input/output vectors that computes the local residual
function ComputeLocalResidual!(fx, x, x_old, da)

    # Compute the local residual. The local vectors of x/x_old include ghost points 
    Phi     = PETSc.getlocalarraydof(da, x,     dof = 1) 
    Pe      = PETSc.getlocalarraydof(da, x,     dof = 2) 
    Phi_old = PETSc.getlocalarraydof(da, x_old, dof = 1) 
    Pe_old  = PETSc.getlocalarraydof(da, x_old, dof = 2) 
    
    # The local residual vectors do not include ghost points
    res_Phi = PETSc.getlocalarraydof(da, fx, dof = 1) 
    res_Pe  = PETSc.getlocalarraydof(da, fx, dof = 2) 

    # Coordinates and spacing (assumed constant )
    Coord   = PETSc.getlocalcoordinatearray(da)
    Δx      = Coord.X[CartesianIndex(2)]-Coord.X[CartesianIndex(1)];
    
    # Global grid size
    Nq      = PETSc.getinfo(da).global_size
    
    # Local sizes
    corners = PETSc.getcorners(da)
    
    ix_p1 =  CartesianIndex( 1, 0, 0);  # ix + 1 
    ix_m1 =  CartesianIndex(-1, 0, 0);  # ix - 1
    
    for i in ((corners.lower):(corners.upper))

        if (i[1] == 1 || i[1]==Nq[1])
            res_Phi[i] = Phi[i] - 1.0;  # Dirichlet upper/lower BC's
            res_Pe[i]  = Pe[i] - 0.0;   # Dirichlet upper/lower BC's
            
        else
            res_Phi[i] = (Phi[i] - Phi_old[i])/dt + De*(Pe[i]-Pe_old[i])/dt + (Phi[i]^m)* Pe[i];

            Phi_c_p1   = 0.5*(Phi[i + ix_p1] + Phi[i]);
            Phi_c_m1   = 0.5*(Phi[i + ix_m1] + Phi[i]);
              
            res_Pe[i]  = De*(Pe[i]-Pe_old[i])/dt - ((Phi_c_p1^n) * ( (Pe[i + ix_p1] - Pe[i        ])/Δx + 1.0)  -
                                                    (Phi_c_m1^n) * ( (Pe[i        ] - Pe[i + ix_m1])/Δx + 1.0))/Δx  +  
                                                    (Phi[i]^m) * Pe[i];

        end
    end

end


function  ForwardDiff_res(x, x_old, da)
    fx  = similar(x)               # vector of zeros, of same type as x (local vector)

    ComputeLocalResidual!(fx, x, x_old, da)

    return fx;
end

f_Residual  =   (x -> ForwardDiff_res(x, x_old, da));  

# Set up the nonlinear function
r = similar(g_x)
PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get a local vector and transfer the data from the global->local vector
    l_x     = PETSc.DMLocalVec(da)
    PETSc.update!(l_x,      g_x,  PETSc.INSERT_VALUES)

    # Get local arrays
    PETSc.withlocalarray!(
        (g_fx, l_x);
        read = (false, true),
        write = (true, false),
    ) do fx, x

        # Compute the local residual (using local julia vectors)
        ComputeLocalResidual!(fx, x, x_old, da)

    end

    # Clean up the local vectors
    PETSc.destroy(l_x)
    
    return 0
end

J = PETSc.MatAIJ(da)    # initialize space for the matrix from the dmda
PETSc.setjacobian!(snes, J, J) do J, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get a local vector and transfer the data from the global->local vector
    l_x     = PETSc.DMLocalVec(da)

    PETSc.update!(l_x,      g_x,  PETSc.INSERT_VALUES)

    # Get a local array of the solution vector
    PETSc.withlocalarray!(
        (l_x);
        read = (true),
        write = (false),
        ) do x_julia   
            # Note that we need to extract x_old as well, as it is used inside the residual routine

            # use ForwardDiff to compute (local part of) jacobian
            # NOTE: this initializes a full matrix which is afterwards converted to a sparse matrix.
            #  that is SLOW; using SparsityDetection & SparseDiffTools this can likely be made faster 
            S     =   sparse(ForwardDiff.jacobian(f_Residual, x_julia));
            
            # Get the non-ghosted part of the local matrix
            ind_local = PETSc.LocalInGlobalIndices(da) 
            S = S[ind_local,ind_local];
            
            # Copy sparse to PETSc matrix
            copyto!(J,S);   

    end

    # Assemble the Jacobian matrix
    PETSc.assemble!(J)

    # Clean up the local vectors
    PETSc.destroy(l_x)
    
    return 0
end


# Timestep loop
time, it = 0.0, 1;
while (it<4000 && time<25)

    global time, it, x_old, g_xold, g_x

    # Solve nonlinear system of equations
    PETSc.solve!(g_x, snes)
   
    # Just for checking: compute norm of nonlinear solution
    g = similar(g_x)
    snes.f!(g, snes, g_x)
    nm = norm(g)
    if MPI.Comm_rank(comm) == 0
        @show nm
    end

    # Update time
    time += dt;
    it   += 1;

    # Update the x_old values (Phi_old, Pe_old) on every processor
    PETSc.update!(l_xold, g_x, PETSc.INSERT_VALUES)
    x_old = PETSc.unsafe_localarray(l_xold,read=true, write=false)

    # Visualisation
    if (mod(it,20)==0 && CreatePlots==true)  
        Coord   =   PETSc.getlocalcoordinatearray(da)
        x = Coord.X[:];
        time_vec = ones(size(x))*time;
        #=
        p1 = plot!(x,time_vec,g_x[2:2:end], zlabel="Pe",title="De=$(De)", ylabel="time", xlabel="depth", 
                linecolor=:blue, 
                legend=:none, 
                xlims=(-20, 150),
                ylims=(0, 25),
                zlims=(-1.5, 1.5),
                camera = (20, 70),
                dpi=300)    # Pe
        savefig(p1, "Pe_porositywave") 
        =#

        p2 = plot!(x,time_vec,g_x[1:2:end], zlabel="Phi",title="De=$(De)", ylabel="time", xlabel="depth", 
                linecolor=:blue, 
                legend=:none, 
                xlims=(-20, 150),
                ylims=(0, 25),
                zlims=(0, 2),
                camera = (20, 70),
                dpi=300)    # Pe
        savefig(p2, "Phi_porositywave")     end

    println("Timestep $it, time=$time")
end 

# Do some clean up
PETSc.destroy(J)
PETSc.destroy(g_x)
PETSc.destroy(g_xold)
PETSc.destroy(r)
PETSc.destroy(da)
PETSc.destroy(snes)


#PETSc.finalize(petsclib)
