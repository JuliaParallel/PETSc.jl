# INCLUDE IN MPI TEST
# 2D incompressible Stokes flow around a dense, stiff sinker, on a staggered grid (DMStag).
#
#   −∇·(2ηε̇) + ∇P = −ρg ŷ      momentum
#          −∇·v   = 0          mass
#
# on the unit square with free-slip walls and g = 1. 
# The sinker is a disc in the middle, ρ = 1 and η = 10⁴ inside, ρ = 0 and η = 1 outside. 
# Walkthrough: docs/src/man/dmstag_stokes.md.
#
#   julia --project=examples examples/dmstag_stokes.jl
#   mpiexec -n 4 julia --project=examples examples/dmstag_stokes.jl -ksp_monitor

using PETSc, MPI, ForwardDiff
using PETSc: LibPETSc

# Where the unknowns of cell (i, j) live:
#
#              Vy(i,j+1)
#         +-------↑-------+
#         |               |
#   Vx(i,j) →   P(i,j)    → Vx(i+1,j)
#         |               |
#         +-------↑-------+
#    (i,j)     Vy(i,j)
#
# Vx sits on the left face, Vy on the bottom face, P in the centre, and the shear stress τxy
# on the lower-left vertex. With n cells per side, Vx runs to i = n + 1 and Vy to j = n + 1.

# Utils
# ----------------------------------------------------------------------------
const VX = LibPETSc.DMSTAG_LEFT
const VY = LibPETSc.DMSTAG_DOWN
const PC = LibPETSc.DMSTAG_ELEMENT

on_wall(g, i, j) = i in (1, g.n + 1) || j in (1, g.n + 1)

# Geometry
# ----------------------------------------------------------------------------
viscosity(x, y, contrast) = (x - 0.5)^2 + (y - 0.5)^2 < 0.15^2 ? contrast : 1.0
density(x, y) = (x - 0.5)^2 + (y - 0.5)^2 < 0.15^2 ? 1.0 : 0.0

# The equations
# ----------------------------------------------------------------------------
#
# Each equation is the residual of one row, written with the unknowns as functions
# vx(i, j), vy(i, j) and p(i, j). `g` holds the grid: n cells per side of width h.

# Zero on the walls, so τxy vanishes there: the free-slip condition
η_cell(g, i, j) = viscosity((i - 0.5) * g.h, (j - 0.5) * g.h, g.contrast)
η_vertex(g, i, j) = on_wall(g, i, j) ? 0.0 : viscosity((i - 1) * g.h, (j - 1) * g.h, g.contrast)

# Stresses: τxx and τyy in the cell centres, τxy on the vertices
τxx(g, vx, i, j) = 2η_cell(g, i, j) * (vx(i + 1, j) - vx(i, j)) / g.h
τyy(g, vy, i, j) = 2η_cell(g, i, j) * (vy(i, j + 1) - vy(i, j)) / g.h
τxy(g, vx, vy, i, j) = η_vertex(g, i, j) * ((vx(i, j) - vx(i, j - 1)) / g.h + (vy(i, j) - vy(i - 1, j)) / g.h)

# x momentum at Vx(i, j): −∂τxx/∂x − ∂τxy/∂y + ∂P/∂x = 0
x_momentum(g, vx, vy, p, i, j) =
    -(τxx(g, vx, i, j) - τxx(g, vx, i - 1, j)) / g.h -
     (τxy(g, vx, vy, i, j + 1) - τxy(g, vx, vy, i, j)) / g.h +
     (p(i, j) - p(i - 1, j)) / g.h

# y momentum at Vy(i, j): −∂τxy/∂x − ∂τyy/∂y + ∂P/∂y + ρg = 0
y_momentum(g, vx, vy, p, i, j) =
    -(τxy(g, vx, vy, i + 1, j) - τxy(g, vx, vy, i, j)) / g.h -
     (τyy(g, vy, i, j) - τyy(g, vy, i, j - 1)) / g.h +
     (p(i, j) - p(i, j - 1)) / g.h + density((i - 0.5) * g.h, (j - 1) * g.h)

# Mass in cell (i, j): −∂Vx/∂x − ∂Vy/∂y = 0
mass(g, vx, vy, p, i, j) = -(vx(i + 1, j) - vx(i, j)) / g.h - (vy(i, j + 1) - vy(i, j)) / g.h

# The wall-normal velocities, which are zero
is_wall(g, loc, i, j) = (loc == VX && i in (1, g.n + 1)) || (loc == VY && j in (1, g.n + 1))

# From an equation to a matrix row
# ----------------------------------------------------------------------------

# The unknowns a row can reach: each field in the 3×3 box of cells around the row's cell
const BOX = [(loc, di, dj) for loc in (VX, VY, PC) for dj in -1:1 for di in -1:1]

"""
    linearize(equation, g, dm, i, j) -> entries, constant

Row (i, j) of a linear `equation`, as `stencil => coefficient` pairs, and its constant
part: the equation's value is `Σ coefficient * unknown + constant`. 
The coefficients are the gradient with respect to the unknowns in `BOX`, from ForwardDiff. 
A wall velocity is read as the constant 0, so its column is left out, which keeps the matrix symmetric.
"""
function linearize(equation, g, dm, i, j)
    function residual(u)
        unknown(loc, k, l) = is_wall(g, loc, k, l) ? 0.0 : u[findfirst(==((loc, k - i, l - j)), BOX)]
        vx(k, l) = unknown(VX, k, l)
        vy(k, l) = unknown(VY, k, l)
        p( k, l) = unknown(PC, k, l)
        return equation(g, vx, vy, p, i, j)
    end
    u0 = zeros(length(BOX))
    coefficients = ForwardDiff.gradient(residual, u0)
    entries = [PETSc.stencil(dm, loc, (i + di, j + dj)) => c
               for ((loc, di, dj), c) in zip(BOX, coefficients) if c != 0]
    return entries, residual(u0)
end

# The grid and the system
# ----------------------------------------------------------------------------

# The cells and faces this rank owns. The last rank along an axis also owns the faces on
# the upper wall: `nextra` is 1 there.
function owned(dm)
    c = PETSc.corners(dm)
    cells_i = c.lower[1]:c.upper[1]
    cells_j = c.lower[2]:c.upper[2]
    faces_i = c.lower[1]:(c.upper[1] + c.nextra[1])
    faces_j = c.lower[2]:(c.upper[2] + c.nextra[2])
    return (; cells_i, cells_j, faces_i, faces_j)
end

"""
    stokes_system(petsclib, comm, n, contrast) -> dm, A, P, b

Assemble the Stokes equations on an `n × n` grid: the matrix `A`, the right-hand side `b`,
and `P`, the matrix the preconditioner is built from. `P` is `A` plus −1/(2η) in the
pressure block, the approximation of the Schur complement. The caller destroys all four.
"""
function stokes_system(petsclib, comm, n, contrast)
    ghosted = PETSc.DM_BOUNDARY_GHOSTED
    dm = PETSc.DMStag(petsclib, comm, 
                      (ghosted, ghosted),           # Boundary type 
                      (n, n),                       # Global dimensions
                      (0, 1, 1),                    # nothing on vertices, 1 per face, 1 per cell
                      1, PETSc.DMSTAG_STENCIL_BOX)  # Stencil type
    A = PETSc.PetscMat(dm)
    P = PETSc.PetscMat(dm)
    b = PETSc.global_vec(dm)
    g = (; n, h = 1 / n, contrast)

    o = owned(dm)
    rows = ((VX, x_momentum, o.faces_i, o.cells_j),
            (VY, y_momentum, o.cells_i, o.faces_j),
            (PC, mass, o.cells_i, o.cells_j))
    for (loc, equation, irange, jrange) in rows, j in jrange, i in irange
        row = PETSc.stencil(dm, loc, (i, j))
        if is_wall(g, loc, i, j)
            entries, constant = [row => 1.0], 0.0
        else
            entries, constant = linearize(equation, g, dm, i, j)
        end
        PETSc.set_values!(A, dm, [row], first.(entries), last.(entries))
        loc == PC && push!(entries, row => -1 / (2η_cell(g, i, j)))
        PETSc.set_values!(P, dm, [row], first.(entries), last.(entries))
        PETSc.set_values!(b, dm, [row], [-constant])
    end

    foreach(PETSc.assemble!, (A, P, b))
    return dm, A, P, b
end

# The solver
# ----------------------------------------------------------------------------

"""
    stokes_solver(petsclib, comm, dm, A, P; levels = 3, options...) -> ksp

FGMRES with a Schur-complement fieldsplit: a geometric multigrid V-cycle with `levels`
levels for the velocity, the diagonal of `P`'s pressure block for the Schur complement.
`options` go to the KSP and override these.
"""
function stokes_solver(petsclib, comm, dm, A, P; levels = 3, options...)
    # Set solver options
    ksp = PETSc.KSP(A, P;
        ksp_type = "fgmres", 
        ksp_rtol = 1e-8,
        pc_type  = "fieldsplit", 
        pc_fieldsplit_type = "schur",
        pc_fieldsplit_schur_fact_type = "upper",
        pc_fieldsplit_schur_precondition = "a11",   # the Schur approximation is P's pressure block
        fieldsplit_u_ksp_type = "preonly",
        fieldsplit_p_ksp_type = "preonly", 
        fieldsplit_p_pc_type = "jacobi",
        options...)

    velocity = LibPETSc.IS(dm, VX => 0, VY => 0)
    pressure = LibPETSc.IS(dm, PC => 0)

    # In a closed box only pressure differences are determined: tell the pressure split
    PETSc.compose_constant_nullspace!(petsclib, comm, pressure)
    PETSc.set_fieldsplit_is!(PETSc.pc(ksp), "u", velocity)
    PETSc.set_fieldsplit_is!(PETSc.pc(ksp), "p", pressure)
    PETSc.destroy!(velocity)                # the fieldsplit keeps its own reference
    PETSc.destroy!(pressure)

    # The velocity split's solver exists only once the fieldsplit has been set up
    PETSc.set_post_setup!(PETSc.pc(ksp)) do fieldsplit
        velocity_multigrid!(fieldsplit, petsclib, comm, dm, levels)
    end
    return ksp
end

# Geometric multigrid on the velocity block. Each level's interpolation is DMStag's, 
# cut down to the velocity unknowns; the coarse operators are Galerkin, RᵀAR.
function velocity_multigrid!(fieldsplit, petsclib, comm, dm, levels)
    _, split_ksps = LibPETSc.PCFieldSplitGetSubKSP(petsclib, fieldsplit)
    mg = PETSc.pc(split_ksps[1])
    PETSc.type_name(mg) === :mg && return nothing      # if done already at an earlier setup
    PETSc.set_type!(mg, :mg)
    LibPETSc.PCMGSetLevels(petsclib, mg, levels, C_NULL)
    LibPETSc.PCMGSetGalerkin(petsclib, mg, LibPETSc.PC_MG_GALERKIN_BOTH)

    fine = dm
    for level in (levels - 1):-1:1         # level 0 is the coarsest
        coarse = LibPETSc.DMCoarsen(petsclib, fine, comm)
        interp, scaling = LibPETSc.DMCreateInterpolation(petsclib, coarse, fine)

        # Get fine and coarsed velocity dm
        fine_u   = LibPETSc.IS(PETSc.narrow(fine)  , VX => 0, VY => 0)
        coarse_u = LibPETSc.IS(PETSc.narrow(coarse), VX => 0, VY => 0)
       
        # Interpolate
        interp_u = LibPETSc.MatCreateSubMatrix(petsclib, interp, fine_u, coarse_u, LibPETSc.MAT_INITIAL_MATRIX)
        LibPETSc.PCMGSetInterpolation(petsclib, mg, level, interp_u)

        # Destroy objects
        foreach(PETSc.destroy!, (interp_u, coarse_u, fine_u, scaling, interp))
        
        # Rotate
        fine === dm || PETSc.destroy!(fine)
        fine = coarse
    end
    fine === dm || PETSc.destroy!(fine)

    # Let -fieldsplit_u_mg_* options reach the levels just made
    LibPETSc.PCSetFromOptions(petsclib, mg)
    return nothing
end

# The solution
# ----------------------------------------------------------------------------

# The owned Vx, Vy and P of the global vector x, as matrices
function owned_fields(dm, x)
    xl = PETSc.local_vec(dm)
    PETSc.global_to_local!(xl, dm, x, PETSc.INSERT_VALUES)
    o = owned(dm)
    fields = PETSc.with_field_views!(dm, xl; fields = (VX => 0, VY => 0, PC => 0), write = false) do (Vx, Vy, P)
        (Vx = [Vx[i, j] for i in o.faces_i, j in o.cells_j],
         Vy = [Vy[i, j] for i in o.cells_i, j in o.faces_j],
         P  = [P[ i, j] for i in o.cells_i, j in o.cells_j])
    end
    PETSc.destroy!(xl)
    return fields
end

# Run
# ----------------------------------------------------------------------------

# Initialize PETSc and MPI
MPI.Initialized() || MPI.Init()
petsclib = first(PETSc.petsclibs)
PETSc.isinitialized(petsclib) || PETSc.initialize(petsclib)
comm = MPI.COMM_WORLD

dm, A, P, b = stokes_system(petsclib, comm, 32, 1e4)
ksp         = stokes_solver(petsclib, comm, dm, A, P; PETSc.parse_options(ARGS)...)
x           = PETSc.global_vec(dm)
PETSc.solve!(x, ksp, b; check = true)

# Gather solution from the ranks
sol = owned_fields(dm, x)
vmax = MPI.Allreduce(maximum(abs, sol.Vy), max, comm)
MPI.Comm_rank(comm) == 0 &&
    println("$(PETSc.iteration_number(ksp)) FGMRES iterations, max |Vy| = $vmax")
foreach(PETSc.destroy!, (ksp, x, b, P, A, dm))
