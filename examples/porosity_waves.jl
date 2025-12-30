# EXCLUDE FROM TESTING
# INCLUDE IN MPI TEST
#=
In this example we solve the 1D viscoelastic porosity wave equations which
describe how magma ascends in the Earth's mantle.  The equations are discussed
in [Vasyliev et al. Geophysical Research Letters (25), 17. p.
3239-3242](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/98GL52358):
```math
        { ∂ϕ/∂t} = -De {∂Pe/∂t} - ϕ^m Pe
     De {∂Pe/∂t} =  ∇( ϕ^n (∇Pe + e_x) ) - ϕ^m Pe
```
with Dirichlet boundary conditions Pe=0 and ϕ=1 on either side.  n,m are
parameters and De is the Deborah number which indicates the importance of
elasticity (De>>1).  e_x is the unit vector in the direction of gravity (here
taken to be the x-direction)

To solve the problem we use a standard central, finite difference approximation
in 1D, and simultaneously solve for Pe and ϕ on a collocated grid.  We employ a
DMDA/SNES and compute the Jacobian by automatic differentiation.

Note: the code is working for 1D problems, but 2D and 3D are currently not functional due to an issue with the lateral boundary conditions
=#

using MPI
using PETSc
using OffsetArrays: OffsetArray
using LinearAlgebra: norm
using ForwardDiff
using SparseArrays
using SparseDiffTools

# when not run interactive options can be passed via commandline arguments:
# julia --project=examples examples/porosity_waves.jl -Nq1 101 -snes_monitor
opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (
        ksp_monitor = true,
        ksp_rtol = 1e-10,
        ksp_converged_reason = false,
        snes_monitor = true,
        snes_converged_reason = true,
        snes_view = false,
        snes_max_it = 100,
        snes_max_funcs = 10000,
        snes_max_linear_solve_fail = 1000,
        snes_mf_operator = false,
        Nq1 = 101,
        max_it = 4000,        # 4000 for a real simulation, 40 for testing
        max_time = 25,
        dim = 1
    )
end

# To use plotting either `]add Plots`
CreatePlots = isinteractive() && try
    using GLMakie
    true
catch
    false
end
CreatePlots=false
@show CreatePlots

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64
PetscInt = Int64

# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar, PetscInt = PetscInt)

# Initialize PETSc
PETSc.initialize(petsclib)

# dimensionality of the problem
dim = PETSc.typedget(opts, :dim, 1)

# Set the total number of grid points in each direction
Nq1 = PETSc.typedget(opts, :Nq1, 11)

Nq = (Nq1, 5, 5)
Nq = Nq[1:dim]

max_it = PETSc.typedget(opts, :max_it, 10)
max_time = PETSc.typedget(opts, :max_time, 25)

if CreatePlots
    fig = Figure(size = (1200, 600))
    if dim==1
        ax1 = Axis3(fig[1, 1], title = "Pe", ylabel= "time", xlabel = "depth", zlabel="Pe", azimuth=-2.191039610705103, elevation=1.0426209566987232)
        ax2 = Axis3(fig[1, 2], title = "ϕ", ylabel= "time", xlabel = "depth", zlabel="ϕ", azimuth=-2.191039610705103, elevation=1.0426209566987232)
    end
end

# Set the boundary conditions on each side (Dirichlet in direction of gravity)
bcs = (
    PETSc.DM_BOUNDARY_NONE,
    PETSc.DM_BOUNDARY_GHOSTED,
    PETSc.DM_BOUNDARY_GHOSTED,
   )[1:dim]

# Set parameters
n = 3
m = 2
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
PETSc.setuniformcoordinates_dmda!(da, (-20, -10, -10), (150, 10, 10))

# Create the PETSC snes object
snes = PETSc.SNES(petsclib, comm; opts...)

# add the da to the snes
PETSc.setDM!(snes, da)

# Sets the initial profiles
g_x = PETSc.DMGlobalVec(da)
PETSc.withlocalarray!(g_x; read = false) do l_x
    corners = PETSc.getcorners(da)

    x = PETSc.reshapelocalarray(l_x, da)
    Phi = @view x[1, :, :, :]
    Pe = @view x[2, :, :, :]

    # retrieve arrays with local coordinates
    coord = PETSc.getlocalcoordinatearray(da)

    Phi0 = 1
    dPhi1 = 8
    dPhi2 = 1
    z1 = 0
    z2 = 40
    lambda = 1
    dPe1 = dPhi1 / De
    dPe2 = dPhi2 / De

    # Loop over all the points on the processor and set the initial condition
    for i in ((corners.lower):(corners.upper))
        Phi[i] =
            Phi0 +
            dPhi1 * exp(-((coord[1, i] - z1)^2) / lambda^2) +
            dPhi2 * exp(-((coord[1, i] - z2)^2) / lambda^2)
        Pe[i] =
            -dPe1 * exp(-((coord[1, i] - z1)^2) / lambda^2) -
            dPe2 * exp(-((coord[1, i] - z2)^2) / lambda^2)
    end
end

# Initialize x_old on every processor
g_xold = deepcopy(g_x); # initialize Pe and Phi of last timestep
l_xold = PETSc.DMLocalVec(da)
PETSc.dm_global_to_local!(g_xold, l_xold, da, PETSc.INSERT_VALUES)

x_old = PETSc.unsafe_localarray(l_xold, read = true, write = false)

# Routine wuth julia-only input/output vectors that computes the local residual
function ComputeLocalResidual!(l_fx, l_x)

    # Compute the local residual.
    # The local vectors of l_x/x_old include ghost points
    x = PETSc.reshapelocalarray(l_x, da)
    Phi = @view x[1, :, :, :]
    Pe = @view x[2, :, :, :]
    x_old_reshaped = PETSc.reshapelocalarray(x_old, da)
    Phi_old = @view x_old_reshaped[1, :, :, :]
    Pe_old = @view x_old_reshaped[2, :, :, :]

    # The local residual vectors do not include ghost points
    fx = PETSc.reshapelocalarray(l_fx, da)
    res_Phi = @view fx[1, :, :, :]
    res_Pe = @view fx[2, :, :, :]

    # Global grid size
    Nq = PETSc.getinfo(da).global_size

    # Local sizes
    corners = PETSc.getcorners(da)

    # set ghost boundaries (flux free conditions)
    if Nq[2] > 1
        if corners.lower[2] == 1
            Phi[:, 0, :] = Phi[:, 1, :]
            Pe[:, 0, :] = Pe[:, 1, :]
        end
        if corners.upper[2] == Nq[2]
            Phi[:, Nq[2] + 1, :] = Phi[:, Nq[2], :]
            Pe[:, Nq[2] + 1, :] = Pe[:, Nq[2], :]
        end
    end
    if Nq[3] > 1
        if corners.lower[3] == 1
            Phi[:, :, 0] = Phi[:, :, 1]
            Pe[:, :, 0] = Pe[:, :, 1]
        end
        if corners.upper[3] == Nq[3]
            Phi[:, :, Nq[3] + 1] = Phi[:, :, Nq[3]]
            Pe[:, :, Nq[3] + 1] = Pe[:, :, Nq[3]]
        end
    end

    # Stencil
    ix_p1 = CartesianIndex(1, 0, 0)   # ix + 1
    ix_m1 = CartesianIndex(-1, 0, 0)  # ix - 1
    iy_p1 = CartesianIndex(0, 1, 0)   # iy + 1
    iy_m1 = CartesianIndex(0, -1, 0)  # iy - 1
    iz_p1 = CartesianIndex(0, 0, 1)   # iz + 1
    iz_m1 = CartesianIndex(0, 0, -1)  # iz - 1

    # Coordinates and spacing (assumed constant)
    coord = PETSc.getlocalcoordinatearray(da)
    Δx = coord[1, corners.lower + ix_p1] - coord[1, corners.lower]
    if dim > 1
        Δy = coord[2, corners.lower + iy_p1] - coord[2, corners.lower]
    end
    if dim == 3
        Δz = coord[1, corners.lower + iz_p1] - coord[1, corners.lower]
    end

    # corners.lower
    for i in ((corners.lower):(corners.upper))
        if (i[1] == 1 || i[1] == Nq[1])   # Dirichlet upper/lower BC's
            res_Phi[i] = Phi[i] - 1
            res_Pe[i] = Pe[i] - 0

        else                            # Central points
            res_Phi[i] =
                (Phi[i] - Phi_old[i]) / dt +
                De * (Pe[i] - Pe_old[i]) / dt +
                (Phi[i]^m) * Pe[i]

            Phi_c_p1 = (Phi[i + ix_p1] + Phi[i]) / 2
            Phi_c_m1 = (Phi[i + ix_m1] + Phi[i]) / 2

            res_Pe[i] =
                De * (Pe[i] - Pe_old[i]) / dt -
                (
                    (Phi_c_p1^n) * ((Pe[i + ix_p1] - Pe[i]) / Δx + 1) -
                    (Phi_c_m1^n) * ((Pe[i] - Pe[i + ix_m1]) / Δx + 1)
                ) / Δx + (Phi[i]^m) * Pe[i]
            if dim > 1
                # Derivatives in y-direction (for 2D and 3D)
                Phi_c_p1 = (Phi[i + iy_p1] + Phi[i]) / 2
                Phi_c_m1 = (Phi[i + iy_m1] + Phi[i]) / 2
                res_Pe[i] =
                    res_Pe[i] +
                    (
                        (Phi_c_p1^n) * ((Pe[i + iy_p1] - Pe[i]) / Δy) -
                        (Phi_c_m1^n) * ((Pe[i] - Pe[i + iy_m1]) / Δy)
                    ) / Δy
            end
            if dim == 3
                # Derivatives in z-direction (for 3D)
                Phi_c_p1 = (Phi[i + iz_p1] + Phi[i]) / 2
                Phi_c_m1 = (Phi[i + iz_m1] + Phi[i]) / 2
                res_Pe[i] =
                    res_Pe[i] +
                    (
                        (Phi_c_p1^n) * ((Pe[i + iz_p1] - Pe[i]) / Δz) -
                        (Phi_c_m1^n) * ((Pe[i] - Pe[i + iz_m1]) / Δz)
                    ) / Δz
            end
        end
    end
end

# Helper function, for autmatic differentiation
function ForwardDiff_res(x)
    fx = similar(x) # vector of zeros, of same type as x (local vector)

    ComputeLocalResidual!(fx, x)

    return fx
end

# Set up the nonlinear function
r = similar(g_x)
PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get a local vector and transfer the data from the global->local vector
    l_x = PETSc.DMLocalVec(da)
    PETSc.dm_global_to_local!(g_x, l_x, da, PETSc.INSERT_VALUES)

    # Get local arrays
    PETSc.withlocalarray!(
        (g_fx, l_x);
        read = (false, true),
        write = (true, false),
    ) do fx, x

        # Compute the local residual (using local julia vectors)
        ComputeLocalResidual!(fx, x)
    end

    # Clean up the local vectors
    #PETSc.destroy(l_x)

    return 0
end

# Compute the sparsity pattern and coloring of the jacobian @ every processor
# Ideally this should be done with an automatic sparsity detection algorithm;
# yet this is currently not working in combination with OffsetArrays (see
# https://github.com/JuliaDiff/SparseDiffTools.jl/issues/154) Therefore, we use
# a less efficient approach here.
l_x = PETSc.DMLocalVec(da)

input = rand(length(l_x))
sparsity_pattern =
    sparse(abs.(ForwardDiff.jacobian(ForwardDiff_res, input)) .> 0);
jac = Float64.(sparsity_pattern)
colors = matrix_colors(jac)

# Compute the Jacobian, using automatic differentiation
J = PETSc.MatAIJ(da)        # initialize space for the matrix from the dmda
PETSc.setjacobian!(snes, J) do J, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get a local vector and transfer the data from the global->local vector
    l_x = PETSc.DMLocalVec(da)
    PETSc.dm_global_to_local!(g_x, l_x, da, PETSc.INSERT_VALUES)

    # Get a local array of the solution vector
    PETSc.withlocalarray!((l_x); read = (true), write = (false)) do x_julia
        # Note that we need to extract x_old as well, as it is used inside the
        # residual routine

        # use ForwardDiff to compute (local part of) jacobian. This is slow, as
        # it forms a dense matrix (and therefor commented)
        # S =   sparse(ForwardDiff.jacobian(ForwardDiff_res, x_julia));

        # Using SparseDiffTools this is more efficient, but requires
        # initializing the sparsity pattern & coloring:
        S = forwarddiff_color_jacobian(
            ForwardDiff_res,
            x_julia,
            colorvec = colors,
            sparsity = sparsity_pattern,
            jac_prototype = jac,
        )

        ind_local = PETSc.localinteriorlinearindex(da)

        S = S[ind_local, ind_local]

        # Copy local sparse matrix to parallel PETSc matrix
        copyto!(J, S)
    end

    # Assemble the Jacobian matrix
    PETSc.assemble!(J)

    # Clean up the local vectors
    #PETSc.destroy(l_x)
    #PETSc.destroy(da)
    return 0
end

# Timestep loop
@show max_it
time, it = PetscScalar(0), 1;
#while (it < max_it && time < max_time)
    @show it
    global time, it, x_old, g_xold, g_x, ax1, fig

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
    time += dt
    it += 1

    # Update the x_old values (Phi_old, Pe_old) on every processor
    #finalize(x_old) # Free the previous x_old vector
    PETSc.dm_global_to_local!(g_x, l_xold, da, PETSc.INSERT_VALUES)

    x_old = PETSc.unsafe_localarray(l_xold, read = true, write = false)

    # Visualisation
    if (mod(it, 40) == 0 && CreatePlots == true)
        coord = PETSc.getlocalcoordinatearray(da)
        if dim==1
            x = coord[1, :, :, :]
            x = x[:,1,1]
            time_vec = ones(size(x)) * time
            
            Pe = g_x[2:2:end]
            ϕ  = g_x[1:2:end]
            
            # Plot Pe
            lines!(ax1,x,time_vec,Pe, color=:black)
            ax1.title = "De=$(De)"
            xlims!(ax1,(-20, 150))
            ylims!(ax1,(0, 25))

            # Plot ϕ   
            lines!(ax2,x,time_vec,ϕ, color=:black)
            ax2.title = "De=$(De)"
            xlims!(ax2,(-20, 150))
            ylims!(ax2,(0, 25))
            zlims!(ax2,(0, 2))
        end
        display(fig)

    end
    
    if MPI.Comm_rank(comm) == 0
        println("Timestep $it, time=$time")
    end
#end

if CreatePlots
    save("porositywave.png",fig)
end

# Do some clean up
#=
PETSc.destroy(J)
PETSc.destroy(g_x)
finalize(x_old);
PETSc.destroy(r)
PETSc.destroy(snes)
PETSc.destroy(da)

PETSc.finalize(petsclib)

=#