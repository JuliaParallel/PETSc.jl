# Translation of PETSc DMStag tutorial ex8.c
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
#   julia --project=. examples/dmstag_ex8.jl
#   julia --project=. examples/dmstag_ex8.jl -pc_type mg -pc_mg_galerkin -pc_mg_levels 4
#   mpiexec -n 2 julia --project=. examples/dmstag_ex8.jl -pc_type mg -pc_mg_galerkin -pc_mg_levels 7 -stag_grid_x 512 -ksp_converged_reason

using PETSc, MPI

if !MPI.Initialized()
    MPI.Init()
end

petsclib = first(PETSc.petsclibs)
PETSc.initialized(petsclib) || PETSc.initialize(petsclib)

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
    (1, 0),                            # (dof_vertex, dof_element)
    1,                                 # stencil width
    PETSc.DMSTAG_STENCIL_BOX;
    opts...)

# ── Assemble the linear system  Ax = b ──────────────────────────────────────

function assemble_system!(dm, petsclib, PetscScalar, PetscInt)
    A   = LibPETSc.DMCreateMatrix(petsclib, dm)
    b   = PETSc.DMGlobalVec(dm)

    corners = PETSc.getcorners_dmstag(dm)
    N       = size(dm)[1]       # global element count  (= number of intervals)
    n_extra = corners.nextra[1] # 1 on the rightmost rank, 0 elsewhere

    # Loop over owned elements, including the extra partial element on the
    # rightmost rank (which holds the right-boundary vertex N).
    for e in (corners.lower[1] - 1) : (corners.upper[1] - 1 + n_extra)
        # Stencil row for the LEFT vertex of element e  (vertex index = e)
        row = LibPETSc.DMStagStencil(
            LibPETSc.DMSTAG_LEFT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0))

        if e == 0
            # ── Left boundary: Dirichlet u(0) = 0 ──
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A,
                PetscInt(1), [row],
                PetscInt(1), [row],
                PetscScalar[1.0],
                PETSc.INSERT_VALUES)

        elseif e == N
            # ── Right boundary: Dirichlet u(1) = 0 ──
            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A,
                PetscInt(1), [row],
                PetscInt(1), [row],
                PetscScalar[1.0],
                PETSc.INSERT_VALUES)

        else
            # ── Interior: standard 3-point Laplacian stencil ──
            col = [
                LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e - 1), PetscInt(0), PetscInt(0), PetscInt(0)),
                LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e),     PetscInt(0), PetscInt(0), PetscInt(0)),
                LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e + 1), PetscInt(0), PetscInt(0), PetscInt(0)),
            ]
            val = PetscScalar[1.0, -2.0, 1.0]

            LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A,
                PetscInt(1), [row],
                PetscInt(3), col,
                val,
                PETSc.INSERT_VALUES)
        end

        # ── Right-hand side (forcing) ──
        h = PetscScalar(1.0 / N)
        x = PetscScalar(e) / PetscScalar(N)
        f = (x - PetscScalar(0.5)) * h * h        # scaled by h²

        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm, b,
            PetscInt(1), [row],
            PetscScalar[f],
            PETSc.INSERT_VALUES)
    end

    PETSc.assemble!(A)
    PETSc.assemble!(b)

    return A, b
end

A, b = assemble_system!(dm, petsclib, PetscScalar, PetscInt)

# ── Solve ────────────────────────────────────────────────────────────────────
# Create KSP from the DM (enables geometric multigrid via DM hierarchy),
# then set operators and disable DM-active mode (we supply our own A).
ksp = PETSc.KSP(dm;
    ksp_type = "fgmres",
    ksp_monitor = true,
    ksp_converged_reason = true,
    opts...)

LibPETSc.KSPSetOperators(petsclib, ksp, A, A)
LibPETSc.KSPSetDMActive(petsclib, ksp, LibPETSc.PETSC_FALSE)

# Re-apply options so that command-line overrides (e.g. -pc_type mg) take effect
# after KSPSetOperators / KSPSetDMActive
if !isempty(opts)
    popts = PETSc.Options(petsclib; opts...)
    push!(popts)
    LibPETSc.KSPSetFromOptions(petsclib, ksp)
    pop!(popts)
else
    LibPETSc.KSPSetFromOptions(petsclib, ksp)
end

x = similar(b)
PETSc.solve!(x, ksp, b)

# ── Check convergence ────────────────────────────────────────────────────────
reason = LibPETSc.KSPGetConvergedReason(petsclib, ksp)
if Integer(reason) < 0
    error("Linear solve failed with reason $reason")
end

# ── Optional: print solution on rank 0 ───────────────────────────────────────
if MPI.Comm_rank(comm) == 0
    N = size(dm)[1]
    println("Solved 1D Poisson on $N elements ($(N+1) vertices)")
end

# ── Cleanup ──────────────────────────────────────────────────────────────────
PETSc.destroy(ksp)
PETSc.destroy(A)
PETSc.destroy(b)
PETSc.destroy(x)
PETSc.destroy(dm)
