
# # Finite Difference Laplacian
#
# This example solves a simple 2-D Laplacian with zero Dirichlet boundary 
# conditions using the PETSc DMDA interface.
#
# The approach taken is that the 5-point Laplacian is applied in the interior
# and boundary points are included in the solve, but the forcing for the
# boundary data is the Dirichlet boundary condition and the matrix entry for
# these points is set to the identity.
#
# The exact solution is u(x,y) = sin(2πx)sin(2πy) on [0,1]².
#
# ## Command Line Options:
#
# Grid configuration:
#   --N <value>          Grid size (default: 65, creates N×N grid)
#
# PETSc options (use single dash):
#   -ksp_type <type>     KSP solver type (preonly, gmres, cg, etc.)
#   -pc_type <type>      Preconditioner type (lu, mg, sor, ilu, etc.)
#   -ksp_rtol <value>    Relative convergence tolerance
#   -ksp_monitor         Monitor convergence
#
# ## Usage Examples:
#
# 1. Run with default settings (65×65 grid, LU solver):
#
#     julia --project=.. dmda_laplacian.jl
#
# 2. Solve on a 129×129 grid with direct LU solver:
#
#     julia --project=.. dmda_laplacian.jl --N 129 -ksp_type preonly -pc_type lu
#
# 3. Large grid with iterative SOR solver:
#
#     julia --project=.. dmda_laplacian.jl --N 257 -pc_type sor -ksp_rtol 1e-10 -ksp_monitor
#
# 4. Small grid with ILU preconditioner:
#
#     julia --project=.. dmda_laplacian.jl --N 33 -pc_type ilu -ksp_rtol 1e-10 -ksp_monitor
#
# 5. Geometric multigrid solver on 513×513 grid:
#
#     julia --project=.. dmda_laplacian.jl --N 513 -pc_type mg -ksp_rtol 1e-12 -mg_coarse_pc_type lu
#
# 6. Run in parallel with MPI (4 processes):
#
#     julia> using MPI
#     julia> mpiexec(cmd -> run(`$cmd -n 4 julia --project dmda_laplacian.jl --N 129 -pc_type mg`))
#
# For convergence analysis across multiple grid resolutions, see dmda_laplacian_convergence.jl

using MPI
using PETSc
using UnicodePlots: heatmap
using ForwardDiff: derivative
using Printf

# boundary data and forcing  
exact(x, y) = sin(2π * x) * sin(2π * y)
# The discrete operator is -∇², so forcing = -∇²(exact)
# For exact(x,y) = sin(2πx)sin(2πy), we have ∇²u = -2(2π)²sin(2πx)sin(2πy)
# So forcing = -∇²u = 8π²sin(2πx)sin(2πy)
forcing(x, y) = 8 * π^2 * sin(2π * x) * sin(2π * y)

# Parse options
opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (ksp_monitor = false, ksp_type = "preonly", pc_type = "lu", ksp_rtol = 1e-12)
end

# Parse grid size from command line (default: 65)
N = parse(Int, get(opts, Symbol("-N"), "65"))

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64

# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

# Grid size (set via --N option, default: 65)
Nq = (N, N)

# Create the DMDA
da = PETSc.DMDA(
    petsclib,
    comm,
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
    Nq,
    1,
    1,
    PETSc.DMDA_STENCIL_STAR;
    opts...,
)

# Setup the Krylov solver
ksp = PETSc.KSP(da; opts...)

# Define the operator assembly function
function assemble_operator!(A, _, ksp)
    da = PETSc.getDM(ksp)
    corners = PETSc.getcorners(da)
    Nq = PETSc.getinfo(da).global_size[1:2]
    
    Δx = PetscScalar(1 / (Nq[1] - 1))
    Δy = PetscScalar(1 / (Nq[2] - 1))
    
    interior = (CartesianIndex(2, 2, 1)):(CartesianIndex(Nq[1] - 1, Nq[2] - 1, 1))
    
    sten = (
        CartesianIndex(-1, 0, 0),
        CartesianIndex(1, 0, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, 0, 0),
    )
    vals = (-1 / Δx^2, -1 / Δx^2, -1 / Δy^2, -1 / Δy^2, 2 / Δx^2 + 2 / Δy^2)
    
    for i in ((corners.lower):(corners.upper))
        if i ∈ interior
            for (s, v) in zip(sten, vals)
                A[i, i + s] = v
            end
        else
            A[i, i] = 1
        end
    end
    
    PETSc.assemble!(A)
    return 0
end

PETSc.setcomputeoperators!(ksp, assemble_operator!)

# Set the right-hand side
PETSc.setcomputerhs!(ksp) do petsc_b, ksp
    da = PETSc.getDM(ksp)
    corners = PETSc.getcorners(da)
    Nq = PETSc.getinfo(da).global_size[1:2]
    
    g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
    g_y = range(PetscScalar(0), length = Nq[2], stop = 1)
    
    l_x = g_x[(corners.lower[1]):(corners.upper[1])]
    l_y = g_y[(corners.lower[2]):(corners.upper[2])]
    
    PETSc.withlocalarray!(petsc_b; read = false) do b
        b = reshape(b, Int64(corners.size[1]), Int64(corners.size[2]))
        b .= forcing.(l_x, l_y')
        
        if corners.lower[1] == 1
            b[1, :] .= exact.(l_x[1], l_y)
        end
        if corners.lower[2] == 1
            b[:, 1] .= exact.(l_x, l_y[1])
        end
        if corners.upper[1] == Nq[1]
            b[end, :] .= exact.(l_x[end], l_y)
        end
        if corners.upper[2] == Nq[2]
            b[:, end] .= exact.(l_x, l_y[end])
        end
    end
    
    return 0
end

# Solve the problem
PETSc.solve!(ksp)

# Check convergence
reason = PETSc.LibPETSc.KSPGetConvergedReason(petsclib, ksp)

# Get the solution and compute error
sol = PETSc.get_solution(ksp)
corners = PETSc.getcorners(da)
Nq = PETSc.getinfo(da).global_size[1:2]

g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
g_y = range(PetscScalar(0), length = Nq[2], stop = 1)

l_x = g_x[(corners.lower[1]):(corners.upper[1])]
l_y = g_y[(corners.lower[2]):(corners.upper[2])]

Δx = PetscScalar(l_x.step)

x = sol[:]
u = reshape(x, corners.size[1], corners.size[2])

# Compute error only in interior (excluding boundaries)
loc_err2_sum = 0.0
loc_max_err = 0.0
loc_npts = 0

for i in 1:length(l_x), j in 1:length(l_y)
    global loc_err2_sum, loc_max_err, loc_npts
    gi = corners.lower[1] + i - 1
    gj = corners.lower[2] + j - 1
    
    if gi > 1 && gi < Nq[1] && gj > 1 && gj < Nq[2]
        err = abs(u[i,j] - exact(l_x[i], l_y[j]))
        loc_err2_sum += err^2
        loc_max_err = max(loc_max_err, err)
        loc_npts += 1
    end
end

root = 0
err2_sum = MPI.Reduce(loc_err2_sum, +, root, comm)
npts_total = MPI.Reduce(loc_npts, +, root, comm)
max_err = MPI.Reduce(loc_max_err, max, root, comm)

# Report results on root processor
if root == MPI.Comm_rank(comm)
    L2_error = sqrt(err2_sum / npts_total)
    
    # Report convergence status
    if Integer(reason) > 0
        println("KSP converged with reason: $reason")
    elseif Integer(reason) < 0
        println("WARNING: KSP diverged with reason: $reason")
    else
        println("KSP still iterating: $reason")
    end
    
    println("L2-error (interior): $(Printf.@sprintf("%.6e", L2_error))")
    println("Max-error (interior): $(Printf.@sprintf("%.6e", max_err))")
    println("Mesh spacing h: $(Printf.@sprintf("%.6e", Δx))")
    
    # Compute error field for plotting
    u_exact = exact.(l_x, l_y')
    error_field = abs.(u .- u_exact)
    
    # Plot solution
    println("\nSolution heatmap:")
    println(heatmap(u'))
    
    # Plot error
    println("\nError heatmap:")
    println(heatmap(error_field'))
end

# Clean up
PETSc.destroy(ksp)
PETSc.destroy(da)

PETSc.finalize(petsclib)

