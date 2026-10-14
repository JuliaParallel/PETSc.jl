# Stokes flow on a DMStag

[`examples/dmstag/dmstag_stokes.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/dmstag/dmstag_stokes.jl) solves 2D incompressible Stokes flow around a dense, stiff sinker: a staggered grid, the equations written as residuals and turned into matrix rows by ForwardDiff, free-slip walls, and a Schur-complement fieldsplit with a geometric multigrid V-cycle on the velocity. This page walks through it. The functions it uses are documented on the [DMStag](dmstag.md) and [KSP](ksp.md) pages.

```
julia --project=examples examples/dmstag/dmstag_stokes.jl
mpiexec -n 4 julia --project=examples examples/dmstag/dmstag_stokes.jl -ksp_monitor
```

## The problem

On the unit square, with viscosity η, density ρ and g = 1:

−∇·(2ηε̇) + ∇P = −ρg ŷ,   −∇·v = 0,   ε̇ = ½(∇v + ∇vᵀ)

The walls are free slip: the normal velocity and the shear stress τxy vanish there. The sinker is a disc of radius 0.15 in the middle, ρ = 1 and η = 10⁴ inside, ρ = 0 and η = 1 outside.

The signs are chosen so that the velocity block A of the matrix is symmetric positive definite, and the whole matrix [A B; Bᵀ 0] symmetric. Multigrid smoothers such as Chebyshev and SOR need a definite A.

## The grid

The DMStag has one unknown on each face and one in each cell:

```
             Vy(i,j+1)
        +-------↑-------+
        |               |
  Vx(i,j) →   P(i,j)    → Vx(i+1,j)
        |               |
        +-------↑-------+
  vertex (i,j)  Vy(i,j)
```

`Vx` sits on the left face of cell (i, j), `Vy` on its bottom face, `P` in its centre, and τxy on its lower-left vertex. DMStag names these locations `DMSTAG_LEFT`, `DMSTAG_DOWN` and `DMSTAG_ELEMENT`. With n cells per side, `Vx` runs to i = n + 1 and `Vy` to j = n + 1: the last rank along an axis owns those faces on the upper wall, and `corners(dm).nextra` counts them.

```julia
ghosted = PETSc.DM_BOUNDARY_GHOSTED
dm = PETSc.DMStag(petsclib, comm, (ghosted, ghosted), (n, n), (0, 1, 1), 1, PETSc.DMSTAG_STENCIL_BOX)
```

`DM_BOUNDARY_GHOSTED` adds ghost points outside the walls. A stencil next to a wall then still has a valid index, and PETSc drops what is written to a point that has no global row. With `DM_BOUNDARY_NONE` the same stencil is out of range.

## The equations

Each equation is a Julia function that returns the residual of one row, with the unknowns read through `vx(i, j)`, `vy(i, j)` and `p(i, j)`. It reads like the finite-difference form of the PDE:

```julia
τxx(g, vx, i, j) = 2η_cell(g, i, j) * (vx(i + 1, j) - vx(i, j)) / g.h
τxy(g, vx, vy, i, j) = η_vertex(g, i, j) * ((vx(i, j) - vx(i, j - 1)) / g.h + (vy(i, j) - vy(i - 1, j)) / g.h)

# x momentum at Vx(i, j): −∂τxx/∂x − ∂τxy/∂y + ∂P/∂x = 0
x_momentum(g, vx, vy, p, i, j) =
    -(τxx(g, vx, i, j) - τxx(g, vx, i - 1, j)) / g.h -
    (τxy(g, vx, vy, i, j + 1) - τxy(g, vx, vy, i, j)) / g.h +
    (p(i, j) - p(i - 1, j)) / g.h
```

`g` holds the grid: n cells per side of width h. The y momentum and mass equations are written the same way.

The matrix row comes from ForwardDiff. The unknowns a row can reach are the three fields in the 3×3 box of cells around (i, j): 27 values. `linearize` puts them in a vector `u`, evaluates the equation at u = 0 and takes the gradient. Since the equations are linear, the gradient is the matrix row and the value at zero is the constant part, here the −ρg of the y momentum:

```julia
coefficients = ForwardDiff.gradient(residual, zeros(27))
entries = [PETSc.stencil(dm, loc, (i + di, j + dj)) => c
           for ((loc, di, dj), c) in zip(BOX, coefficients) if c != 0]
```

[`stencil`](@ref PETSc.stencil) takes the location and the 1-based cell index, and [`set_values!`](@ref PETSc.set_values!) writes the row:

```julia
PETSc.set_values!(A, dm, [row], first.(entries), last.(entries))
PETSc.set_values!(b, dm, [row], [-constant])
```

Assembling a 128² grid this way takes 0.4 s. With a viscosity that depends on the strain rate, the same functions give the residual of a nonlinear solve, and ForwardDiff its exact Jacobian.

## The walls

Free slip has two parts:

- **Zero shear stress.** The viscosity of a vertex on a wall is set to zero, so τxy vanishes there and the equations need no special case.
- **Zero normal velocity.** These unknowns get an identity row, with zero on the right-hand side. Inside the other equations a wall velocity is read as the constant 0, not as an unknown, so its column drops out of the gradient. Keeping those columns would give the same solution but an unsymmetric A, and the Galerkin coarse operators RᵀAR of the multigrid would mix them into the coarse wall rows: here that took the iteration count from 12 to 84 at 32², growing with the grid.

## The solver

### Why a Schur fieldsplit

A multigrid on the whole velocity-pressure matrix needs a smoother for the pressure rows, and those have a zero diagonal, since −∇·v = 0 has no P in it. Point smoothers (Jacobi, SOR) divide by that diagonal. A small penalty on it gives them something to divide by, but they then relax each pressure with the penalty's stiffness, while its real stiffness, through the velocity, is of order 1/η.

A block preconditioner treats the two fields separately. With K = [A B; Bᵀ 0] and the Schur complement S = −Bᵀ A⁻¹ B, the upper block factorization

[A B; 0 S]⁻¹

is the exact inverse when its two solves are exact, and FGMRES then converges in two iterations. The preconditioner approximates each solve.

### The Schur complement: −1/(2η)

For uniform η and incompressible flow, ∇·(2ηε̇(v)) = η∇²v + η∇(∇·v). Applied to a gradient, A∇φ = −2η∇(∇²φ), so A⁻¹∇ψ = −∇(∇⁻²ψ)/(2η), and

S p = ∇·(A⁻¹∇p) = −p/(2η)

for every pressure but the constant. S is a multiple of the identity, and the cell viscosity gives the approximation Ŝ = −diag(1/(2η)) when η varies. The example builds a second matrix P, the same as A plus −1/(2η) in the pressure block, and lets PETSc build the preconditioner from it:

```julia
ksp = PETSc.KSP(A, P;
    ksp_type = "fgmres", ksp_rtol = 1e-8,
    pc_type = "fieldsplit", pc_fieldsplit_type = "schur",
    pc_fieldsplit_schur_fact_type = "upper",
    pc_fieldsplit_schur_precondition = "a11",   # the Schur approximation is P's pressure block
    fieldsplit_u_ksp_type = "preonly",
    fieldsplit_p_ksp_type = "preonly", fieldsplit_p_pc_type = "jacobi")
```

The two splits are given as index sets of whole fields:

```julia
velocity = LibPETSc.IS(dm, DMSTAG_LEFT => 0, DMSTAG_DOWN => 0)
pressure = LibPETSc.IS(dm, DMSTAG_ELEMENT => 0)
PETSc.set_fieldsplit_is!(PETSc.pc(ksp), "u", velocity)
PETSc.set_fieldsplit_is!(PETSc.pc(ksp), "p", pressure)
```

PETSc's `selfp`, which forms Bᵀ diag(A)⁻¹ B, does as well on this grid (13 and 25 iterations at 128², against 12 and 30): the pressure solve applies only its diagonal, and with uniform spacing that diagonal is also a multiple of 1/η.

### The velocity block: geometric multigrid

One V-cycle of Galerkin multigrid on DMStag's grid hierarchy. The velocity split doesn't know about the DM, so the example builds the levels itself: it coarsens the DM, takes DMStag's interpolation between each pair of grids, and cuts it down to the velocity unknowns with `MatCreateSubMatrix`. The split's solver exists only once the fieldsplit has been set up, so this goes in a [`set_post_setup!`](@ref PETSc.set_post_setup!) hook:

```julia
PETSc.set_post_setup!(PETSc.pc(ksp)) do fieldsplit
    velocity_multigrid!(fieldsplit, petsclib, comm, dm, levels)
end
```

The hook runs at every setup, so `velocity_multigrid!` returns early once the split has its multigrid. It ends with `PCSetFromOptions`, so options such as `-fieldsplit_u_mg_levels_ksp_max_it 4` reach the levels it made.

### The constant pressure

In a closed box only pressure differences are determined. A constant null space composed on the pressure index set is handed by the fieldsplit to the Schur solve, which removes the mean:

```julia
PETSc.compose_constant_nullspace!(petsclib, comm, pressure)
```

On this problem it changes little: without it the iteration counts are the same and the mean pressure is 10⁻¹² instead of 10⁻¹⁷. It makes the pressure well defined.

## Reading the solution

[`with_field_views!`](@ref PETSc.with_field_views!) hands out one array per field of a local vector, indexed by the same 1-based (i, j) as the stencils:

```julia
PETSc.with_field_views!(dm, xl; fields = (DMSTAG_LEFT => 0, DMSTAG_DOWN => 0, DMSTAG_ELEMENT => 0), write = false) do (Vx, Vy, P)
    Vy[i, j]
end
```

## Results

FGMRES iterations to a relative residual of 10⁻⁸, serial, three multigrid levels (Linux, PETSc_jll 3.25.4):

| grid | η = 1 everywhere | sinker, η = 10⁴ |
|---|---|---|
| 32² | 12 | 24 |
| 64² | 12 | 29 |
| 128² | 12 | 30 |

The count is flat in h. At 32² the sinker's outline is only a few cells wide, which is why that count is lower. On 2 and 4 ranks the 32² sinker takes 27 and 29 iterations, since SOR is applied per rank.

`examples/dmstag/dmstag_stokes_test.jl` checks the 16² solution against LU to 10⁻⁷, and that the count stays flat from 64² to 128². LU needs a nonsingular matrix, so the test fixes the pressure of one cell at zero. It factors A itself: `ksp_type = "preonly"` with `pc_type = "lu"` factors the matrix the preconditioner is built from, which in the example is P.

## Taking it further

- Command-line options reach the solver: `-ksp_monitor`, `-fieldsplit_u_mg_levels_ksp_max_it 4`, or `-ksp_view` to see the whole configuration.
- Many stiff inclusions slow the V-cycle down: the rigid motions of each inclusion are near the null space of A, the smoother doesn't damp them, and the coarse grids can't represent them.
- A viscosity that depends on the strain rate makes the problem nonlinear. A SNES then evaluates the same equation functions over `with_field_views!` for the residual, and the same gradient, taken at the current iterate instead of at zero, gives the Jacobian.
