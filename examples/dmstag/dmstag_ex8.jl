# INCLUDE IN MPI TEST
# PETSc DMStag tutorial ex8.c
# Solves the 1D Poisson equation on a vertex-centered DMStag grid,
# demonstrating the simplest use of geometric multigrid.
#
#   -u''(x) = f(x)  on [0, 1]
#   u(0) = 0,  u(1) = 0   (Dirichlet BCs)
#   f(x) = x - 0.5          (source term, scaled by h²)
#
# The second-order FD stencil on vertices is:
#   (u[i-1] - 2*u[i] + u[i+1]) = h² * f(x_i)
#
# Usage:
#   julia --project=. examples/dmstag/dmstag_ex8.jl
#   julia --project=. examples/dmstag/dmstag_ex8.jl -pc_type mg -pc_mg_galerkin -pc_mg_levels 4
#   mpiexec -n 2 julia --project=. examples/dmstag/dmstag_ex8.jl -pc_type mg -pc_mg_galerkin -pc_mg_levels 7 -stag_grid_x 512 -ksp_converged_reason

using PETSc, MPI

if !MPI.Initialized()
    MPI.Init()
end

petsclib = first(PETSc.petsclibs)
PETSc.isinitialized(petsclib) || PETSc.initialize(petsclib)

PetscScalar = PETSc.scalartype(petsclib)
PetscInt    = PETSc.inttype(petsclib)

# ── Parse command-line options ──────────────────────────────────────────────
cli = PETSc.parse_options(ARGS)
opts = (; [(k => cli[k]) for k in keys(cli)]...)

# ── Create 1D DMStag ────────────────────────────────────────────────────────
# 1 DOF per vertex (0-cells), 0 DOF per element (1-cells)
comm = MPI.COMM_WORLD

dm = PETSc.DMStag(petsclib, comm,
    (PETSc.DM_BOUNDARY_NONE,),       # boundary type
    (8,),                              # global element count (overridable via -stag_grid_x)
    (1, 0);                            # (dof_vertex, dof_element)
    stencil_width = 1,
    stencil_type = PETSc.DMSTAG_STENCIL_BOX,
    opts...)

# ── Assemble the linear system  Ax = b ──────────────────────────────────────

function assemble_system(dm, PetscScalar)
    A = PETSc.PetscMat(dm)
    b = PETSc.global_vec(dm)

    corners = PETSc.corners(dm)
    N       = size(dm)[1]       # global element count  (= number of intervals)
    n_extra = corners.nextra[1] # 1 on the rightmost rank, 0 elsewhere
    vertex  = PETSc.face_location(dm, 1)    # in 1D the vertices are the faces

    # Loop over owned elements, including the extra partial element on the
    # rightmost rank (which holds the right-boundary vertex N). The left
    # vertex of element i is vertex number i - 1.
    for i in corners.lower[1]:(corners.upper[1] + n_extra)
        row = PETSc.stencil(dm, vertex, (i,))

        if i == 1 || i == N + 1
            # ── Boundaries: Dirichlet u(0) = u(1) = 0 ──
            PETSc.set_values!(A, dm, [row], [row], PetscScalar[1.0])
        else
            # ── Interior: standard 3-point Laplacian stencil ──
            col = [PETSc.stencil(dm, vertex, (i - 1,)), row, PETSc.stencil(dm, vertex, (i + 1,))]
            PETSc.set_values!(A, dm, [row], col, PetscScalar[1.0, -2.0, 1.0])
        end

        # ── Right-hand side (forcing) ──
        h = PetscScalar(1.0 / N)
        x = PetscScalar(i - 1) / PetscScalar(N)
        f = (x - PetscScalar(0.5)) * h * h        # scaled by h²
        PETSc.set_values!(b, dm, [row], PetscScalar[f])
    end

    PETSc.assemble!(A)
    PETSc.assemble!(b)

    return A, b
end

A, b = assemble_system(dm, PetscScalar)

# ── Solve ────────────────────────────────────────────────────────────────────
# Create KSP from the DM (enables geometric multigrid via DM hierarchy), then
# supply our own A and keep the DM for its geometry only. The options are
# applied at the first solve, after the operators and the DM mode are set, so
# command-line choices such as -pc_type mg take effect.
ksp = PETSc.KSP(dm;
    ksp_type = "fgmres",
    ksp_monitor = true,
    ksp_converged_reason = true,
    options_at = :first_solve,
    opts...)
PETSc.set_operators!(ksp, A)
PETSc.set_dm_active!(ksp, false)

x = similar(b)
PETSc.solve!(x, ksp, b; check = true)    # throws ConvergenceError if it diverges

# ── Optional: print solution on rank 0 ───────────────────────────────────────
if MPI.Comm_rank(comm) == 0
    N = size(dm)[1]
    println("Solved 1D Poisson on $N elements ($(N+1) vertices)")
end

# ── Cleanup ──────────────────────────────────────────────────────────────────
PETSc.destroy!(ksp)
PETSc.destroy!(A)
PETSc.destroy!(b)
PETSc.destroy!(x)
PETSc.destroy!(dm)
