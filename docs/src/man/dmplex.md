# DMPlex (Unstructured Meshes + FEM)

DMPlex manages unstructured meshes and provides a complete finite element
workflow: mesh creation, field discretisation, residual/Jacobian assembly via
pointwise callbacks, boundary conditions, VTK output, and Lagrangian mesh
advection.  

## Overview

The typical workflow is:

1. **Create a mesh** — from a Gmsh file or a built-in box (or triangular) mesh.
2. **Set up finite elements** — velocity, pressure, auxiliary fields.
3. **Define the PetscDS** — attach pointwise residual/Jacobian callbacks.
4. **Set boundary conditions** — Dirichlet via `add_boundary!`.
5. **Configure SNES** — attach the FEM assembly and null spaces to solve nonlinear problems.
6. **Time loop** — solve, update history fields, write VTK.

A complete, worked example is in [`examples/ex62b.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex62b.jl)
(2D/3D Stokes with Maxwell viscoelasticity, power-law rheology, and Lagrangian
advection).

---

## Creating a Mesh

### From a Gmsh file

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib(; PetscScalar = Float64)
PETSc.initialize(petsclib)

dm = PETSc.DMPlex(petsclib, MPI.COMM_WORLD; dm_plex_filename = "mesh.msh")
```

Any option accepted by `DMSetFromOptions` can be passed as a keyword argument.
Useful options include `-dm_refine N` (uniform refinement) and
`-dm_plex_dim D` to force the mesh dimension.

### Box mesh (built-in)

```julia
dm = PETSc.DMPlex(petsclib, MPI.COMM_WORLD, 2, true, (10, 10))
#                                            ^dim  ^simplex  ^faces
```

`simplex = true` produces triangles (2D) / tetrahedra (3D); `false` gives
quads / hexahedra.

---

## Finite Element Setup

### Creating FE spaces

```julia
dim     = LibPETSc.DMGetDimension(petsclib, dm)
simplex = PETSc.isplexsimplex(dm)

# P2 velocity (dim components, degree 2)
fe_vel = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, dim, simplex;
                                  degree = 2, prefix = "vel_")

# Discontinuous P1 pressure (1 component, degree 1)
let opts = PETSc.Options(petsclib; pres_petscdualspace_lagrange_continuity = 0)
    push!(opts)
    global fe_pres = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, simplex;
                                              degree = 1, prefix = "pres_")
    pop!(opts)
end

# Lagrange P1 (explicit degree, no PetscOptions lookup)
fe_p1 = PETSc.fe_create_lagrange(petsclib, MPI.COMM_SELF, dim, dim, simplex, 1)

# Copy quadrature rule from velocity to pressure for consistency
PETSc.fe_copy_quadrature!(petsclib, fe_vel, fe_pres)
```

### Attaching FE spaces to the DM

```julia
PETSc.setfield!(dm, 0, fe_vel)    # field 0 = velocity
PETSc.setfield!(dm, 1, fe_pres)   # field 1 = pressure
PETSc.createds!(dm)               # finalise the PetscDS
```

### Naming objects

```julia
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe_vel), "velocity")
PETSc.petsc_setname!(petsclib, dm, "stokes")
```

---

## Pointwise Callbacks

PETSc's FEM assembly calls user-supplied C function pointers at each
quadrature point.  Three macros generate the required `@cfunction` wrappers
from plain Julia functions.

### `@petsc_simple_fn` — boundary / exact-solution functions

Signature: `f(t, x, u, ctx)` where `x` is the coordinate vector and `u` is
the output vector.

```julia
function bg_vel(t, x, u, ctx)
    u[1] = x[1] * 0.1    # Vx = exx * x
    u[2] = x[2] * (-0.1) # Vy = eyy * y
end
const bg_vel_ptr = PETSc.@petsc_simple_fn(bg_vel)
```

### `@petsc_residual_fn` — f0 / f1 residual and post-processing functions

Signature: full PETSc pointwise signature with `dim_`, `Nf`, `NfAux`,
`uOff`, `uOff_x`, `u`, `u_t`, `u_x`, `aOff`, `aOff_x`, `a`, `a_t`, `a_x`,
`t`, `x`, `nC`, `cst`, `out`.  The second argument is the number of output
components.

```julia
function f0_p(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
              aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    f0[1] = 0.0
    for d in 0:dim_-1; f0[1] -= u_x[d*dim_+d+1]; end
end
const f0_p_ptr = PETSc.@petsc_residual_fn(f0_p, 1)
```

### `@petsc_jacobian_fn` — g0/g1/g2/g3 Jacobian functions

Same signature as `@petsc_residual_fn` but with an extra `utShift` argument
before `x`.

```julia
function g1_pu(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g1)
    for d in 0:dim_-1; g1[d*dim_+d+1] = -1.0; end
end
const g1_pu_ptr = PETSc.@petsc_jacobian_fn(g1_pu, dim_*dim_)
```

!!! note "Automatic differentiation"
    Jacobian callbacks can be computed automatically using ForwardDiff:
    ```julia
    function g3_uu(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x, ...)
        deviatoric_stress_tangent!(g3, dim_, u_x, ...)  # uses ForwardDiff internally
    end
    ```
    See `examples/ex62b.jl` for a complete example.

---

## PetscDS: Residuals, Jacobians, and Constants

```julia
ds = PETSc.getds(dm)

# Residual callbacks (f0 = volume source, f1 = flux)
PETSc.set_residual!(ds, 0, f0_u_ptr, f1_u_ptr)   # momentum
PETSc.set_residual!(ds, 1, f0_p_ptr, C_NULL)       # continuity

# Jacobian blocks (fieldI, fieldJ, g0, g1, g2, g3)
PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_uu_ptr)
PETSc.set_jacobian!(ds, 0, 1, C_NULL, C_NULL, g2_up_ptr, C_NULL)
PETSc.set_jacobian!(ds, 1, 0, C_NULL, g1_pu_ptr, C_NULL, C_NULL)

# Preconditioner Jacobian (can differ from the true Jacobian)
PETSc.set_jacobian_preconditioner!(ds, 1, 1, g0_pp_ptr, C_NULL, C_NULL, C_NULL)

# Exact solution (for manufactured solution tests / Dirichlet BCs)
PETSc.set_exact_solution!(ds, 0, exact_vel_ptr)
PETSc.set_exact_solution!(ds, 1, exact_pres_ptr)

# Constants (read by all callbacks as the `cst` argument)
PETSc.set_constants!(ds, [mu, rho, gravity, dt])
```

---

## Boundary Conditions

```julia
# Get the boundary label (created by Gmsh or DMPlexCreateBoxMesh)
label = PETSc.getlabel(dm, "Face Sets")

# Dirichlet (essential) BC on velocity (field 0), boundary tag 1
PETSc.add_boundary!(petsclib, dm,
    LibPETSc.DM_BC_ESSENTIAL, "wall", label,
    PetscInt[1],    # boundary tag values
    0,              # field index
    PetscInt[],     # components (empty = all)
    exact_vel_ptr)

# Split a Gmsh "boundary" label into per-face labels (for box meshes)
PETSc.create_split_boundary_labels!(petsclib, dm)
```

---

## Auxiliary Fields

Auxiliary fields carry per-cell data (phase, history stress, etc.) that callbacks read from the `a` / `aOff` arguments.

```julia
# Clone the primary DM and set up DG-P0 auxiliary fields
dm_aux = PETSc.dmclone(dm)

fe_phase = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, simplex;
                                    degree = 0, prefix = "phase_")
fe_tau   = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, dim*dim, simplex;
                                    degree = 0, prefix = "tau_")
PETSc.setfield!(dm_aux, 0, fe_phase)
PETSc.setfield!(dm_aux, 1, fe_tau)
PETSc.createds!(dm_aux)

aux_vec = PETSc.DMGlobalVec(dm_aux)

# Attach aux_vec so callbacks see it via `a`
LibPETSc.DMSetAuxiliaryVec(petsclib, dm,
    LibPETSc.DMLabel(C_NULL), PetscInt(0), PetscInt(0), aux_vec)
```

---

## Projection and L² Norms

```julia
# Project exact functions onto a DM vector (initialisation / BCs)
PETSc.dm_project_function!(petsclib, dm, 0.0,
    [vel_fn_ptr, pres_fn_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, u)

# Project residual-style callbacks onto a DM vector (post-processing)
PETSc.dm_project_field!(petsclib, dm_out, t, u,
    [copy_vel_ptr, copy_pres_ptr, compute_tau_3x3_ptr],
    LibPETSc.INSERT_ALL_VALUES, out_vec)

# L² error vs. exact solution
err = PETSc.dm_compute_l2diff(petsclib, dm, t,
    [exact_vel_ptr, exact_pres_ptr], nothing, u)
```

---

## SNES Integration and Null Spaces

```julia
# Wire FEM residual/Jacobian assembly into the SNES
PETSc.plex_set_snes_local_fem!(petsclib, dm)

# Constant-pressure null space (for incompressible flow)
null_vec = PETSc.DMGlobalVec(dm)
PETSc.dm_project_function!(petsclib, dm, 0.0,
    [zero_vel_ptr, one_pres_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, null_vec)
LibPETSc.VecNormalize(petsclib, null_vec)
nullspace = GC.@preserve null_vec PETSc.mat_null_space_create(petsclib, comm, (null_vec,))
PETSc.snes_set_jacobian_null_space!(snes, nullspace)

# Attach constant null space directly to the pressure FE
PETSc.fe_compose_constant_null_space!(petsclib, comm, fe_pres)

# Propagate discretisation to coarser MG levels
let cdm = dm
    while convert(Ptr{Cvoid}, cdm) != C_NULL
        PETSc.dm_copy_disc!(dm, cdm)
        cdm = PETSc.dm_get_coarse(cdm)
    end
end

# Cleanup
PETSc.mat_null_space_destroy!(petsclib, nullspace)
```

---

## VTK Output

```julia
# Write a named vector to a VTU file (merges parallel pieces automatically)
PETSc.vtk_save!(petsclib, comm, "solution.vtu", out_vec)

# Mark tensor fields so ParaView shows them as tensors
PETSc.vtk_merge_tensor!(fname, "strainrate", "tau")
```

`vtk_save!` calls `PetscViewerVTKOpen` + `DMView` / `VecView` internally and works in parallel (each rank writes its own piece; PETSc merges the XML).
`vtk_merge_tensor!` post-processes the XML header to annotate multiple tensor fields so ParaView's tensor glyph filter recognises them.

### Writing a ParaView PVD animation file

```julia
pvd_entries = Tuple{Float64,String}[]
for step in 1:nsteps
    # ... solve ...
    fname = "out_$(lpad(step, 4, '0')).vtu"
    PETSc.vtk_save!(petsclib, comm, fname, out_vec)
    push!(pvd_entries, (t, abspath(fname)))
    # rewrite PVD after every step so it is always playable
    open("sim.pvd", "w") do io
        println(io, """<?xml version="1.0"?>""")
        println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
        println(io, "  <Collection>")
        for (ts, vtu) in pvd_entries
            rel = relpath(vtu, dirname(abspath("sim.pvd")))
            println(io, """    <DataSet timestep="$ts" group="" part="0" file="$rel"/>""")
        end
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
end
```

---

## Lagrangian Mesh Advection

For moving-mesh simulations, update the coordinate vector directly:

```julia
# Project P2 velocity onto P1 (vertex-based) to get nodal velocities
dm_p1 = PETSc.dmclone(dm)
PETSc.setfield!(dm_p1, 0, PETSc.fe_create_lagrange(petsclib, MPI.COMM_SELF, dim, dim, simplex, 1))
PETSc.createds!(dm_p1)
vel_p1 = PETSc.DMGlobalVec(dm_p1)
PETSc.dm_project_field!(petsclib, dm_p1, 0.0, u, [copy_vel_ptr],
    LibPETSc.INSERT_ALL_VALUES, vel_p1)

# Move each mesh node by dt * v
PL = typeof(petsclib)
coords = LibPETSc.PetscVec{PL}(C_NULL)
LibPETSc.DMGetCoordinates(petsclib, dm, coords)
LibPETSc.VecAXPY(petsclib, coords, PetscScalar(dt), vel_p1)
```

!!! warning "Mesh distortion"
    Pure Lagrangian advection distorts the mesh over time, degrading solver convergence.  For background-deformation (`bg`) setups the mesh does not need to move at all.  For large-displacement `grav` problems, periodic remeshing is needed for long runs.

---

## Examples

### ex17.jl — Linear Elasticity

[`examples/ex17.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex17.jl)
is a Julia port of PETSc's `snes/tutorials/ex17.c`.  It solves linear isotropic elasticity −∇·σ(u) = f on (0,1)^dim with manufactured exact solutions, demonstrating:

- 2D/3D box meshes (simplex triangles/tetrahedra or quad/hex elements)
- Multiple solution types including uniform strain, axial displacement, and geological shear
- Direct, GAMG, and geometric multigrid solvers
- Near-null space (rigid-body modes) for AMG

```bash
# 2D Q1 quad mesh, quadratic MMS, direct solver
julia --project examples/ex17.jl \
  -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
  -sol_type elas_quad -displacement_petscspace_degree 1 -pc_type lu

# GAMG with rigid-body near-null space (recommended for large problems)
julia --project examples/ex17.jl \
  -dm_plex_simplex 0 -dm_plex_box_faces 8,8 \
  -sol_type elas_quad -pc_type gamg -ksp_type cg

# Geometric multigrid — 3 refinement levels, 3D
julia --project examples/ex17.jl \
  -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 2,2,2 \
  -dm_refine_hierarchy 3 -sol_type elas_quad \
  -pc_type mg -pc_mg_type full \
  -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -ksp_type cg

# MPI — 4 ranks
mpiexec -n 4 julia --project examples/ex17.jl \
  -dm_plex_simplex 0 -dm_plex_box_faces 8,8 \
  -sol_type elas_quad -pc_type gamg -ksp_type cg
```

---

### ex62.jl — Isoviscous Stokes (MMS verification)

[`examples/ex62.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex62.jl)
is a Julia port of PETSc's `snes/tutorials/ex62.c`.  It solves constant-viscosity Stokes on the unit square/cube with polynomial or trigonometric manufactured exact solutions, used for convergence verification of the P2/P1 and Q2/Q1 element pairs.

```bash
# 2D P2/P1 simplex, quadratic MMS — L² ≈ machine precision (exact for P2)
julia --project examples/ex62.jl -sol quadratic \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9

# 2D Q2/Q1 quad, trig MMS, refinement convergence study
julia --project examples/ex62.jl -sol trig \
  -dm_plex_simplex 0 -dm_refine 2 \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9

# 3D Q2/Q1 hex, quadratic MMS
julia --project examples/ex62.jl -sol quadratic \
  -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 3,3,3 \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9
```

!!! note "Boundary marker"
    `ex62.jl` uses a single boundary label (id=1) covering the entire box
    boundary — do **not** pass `-dm_plex_separate_marker`.

---

### ex69.jl — Variable-Viscosity Stokes

[`examples/ex69.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex69.jl)
is a Julia port of PETSc's `snes/tutorials/ex69.c`.  It solves variable-viscosity Stokes on (0,1)² with free-slip walls, using analytical benchmark solutions (`solkx` with exponential viscosity variation, `solcx` with a piecewise-constant jump), for convergence rate verification.

```bash
# Q2/Q1 quad, SolKx (exponential viscosity)
julia --project examples/ex69.jl \
  -dm_plex_simplex 0 -dm_plex_separate_marker \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9

# SolCx with viscosity jump etaB=1000
julia --project examples/ex69.jl \
  -dm_plex_simplex 0 -dm_plex_separate_marker \
  -sol_type solcx -etaB 1e3 \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9

# Refinement study — 2 uniform refinements, Q2/Q1
julia --project examples/ex69.jl \
  -dm_plex_simplex 0 -dm_plex_separate_marker -dm_refine 2 \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9

# MPI — 4 ranks, SolCx, viscosity jump 1000
mpiexec -n 4 julia --project examples/ex69.jl \
  -dm_plex_simplex 0 -dm_plex_separate_marker -dm_refine 1 \
  -sol_type solcx -etaB 1e3 \
  -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu -ksp_rtol 1e-9
```

Expected convergence rates (L² velocity+pressure combined):
- **Q2/Q1 quad**: ~4× per refinement (O(h²), pressure-dominated)
- **Q1/P0 quad**: ~2× per refinement (O(h), velocity-dominated)

---

## Examples: Stokes Flow (ex62b.jl)

[`examples/ex62b.jl`](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex62b.jl) solves the incompressible Stokes equations on a square/cube domain containing a circular/spherical inclusion meshed with Gmsh.  It supports four solution types (`-sol quadratic|trig|bg|grav`), power-law and Maxwell viscoelastic rheology, time stepping, and VTK/PVD output.

### Solver options

Three solver families are available, trading exactness for scalability:

| Label | Velocity block | Pressure block | Best for |
|-------|---------------|----------------|----------|
| **(A) Direct** | LU (`lu`) | LU | Small 2D problems, exact reference |
| **(B) GAMG** | CG + GAMG | Jacobi | Large 2D, 3D, moderate viscosity contrast |
| **(C) GMG** | FGMRES + geometric MG | Jacobi | 3D with coarse base mesh + refinement hierarchy |

All three use a Schur-complement fieldsplit preconditioner:
```
-pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur
-pc_fieldsplit_schur_factorization_type full
-pc_fieldsplit_schur_precondition a11
```

### 2D examples

**Pure shear, direct solver** (quick test, serial):
```bash
julia --project=examples examples/ex62b.jl \
  -sol bg -exx_bg 1.0 -eyy_bg -1.0 -mu_inclusion 0.1 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
  -ksp_rtol 1e-9 -snes_rtol 1e-9 -vtk_output out.vtu
```

**Gravity-driven sinking cylinder, viscous, direct solver** (20 time steps):
```bash
julia --project=examples examples/ex62b.jl \
  -sol grav -mesh_h 0.05 -mesh_h_inclusion 0.02 \
  -rho_bg 1.0 -rho_inclusion 2.0 -gravity 1.0 -mu_inclusion 0.1 \
  -dt 0.02 -nsteps 20 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
  -ksp_rtol 1e-9 -snes_rtol 1e-9 -vtk_output grav_viscous.vtu
```

**Gravity-driven sinking cylinder, Maxwell viscoelastic** (τ_rel = η/G = 0.2):
```bash
julia --project=examples examples/ex62b.jl \
  -sol grav -mesh_h 0.05 -mesh_h_inclusion 0.02 \
  -rho_bg 1.0 -rho_inclusion 2.0 -gravity 1.0 -mu_inclusion 0.1 \
  -G_bg 5.0 -G_inclusion 0.5 -dt 0.01 -nsteps 50 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
  -ksp_rtol 1e-9 -snes_rtol 1e-9 -vtk_output grav_viscoelastic.vtu
```

**Pure shear, power-law inclusion (n=3), viscoelastic matrix** (100 steps, GAMG):
```bash
julia --project=examples examples/ex62b.jl \
  -sol bg -exx_bg 1.0 -eyy_bg -1.0 \
  -mesh_h 0.04 -mesh_h_inclusion 0.015 \
  -mu_inclusion 0.1 -n_inclusion 3.0 -eps0_inclusion 0.5 \
  -G_bg 10.0 -G_inclusion 1.0 -dt 0.005 -nsteps 100 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_ksp_type cg -fieldsplit_velocity_pc_type gamg \
  -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi \
  -ksp_type fgmres -ksp_rtol 1e-6 -snes_rtol 1e-6 -snes_max_it 20 \
  -vtk_output shear_powerlaw.vtu
```

**Analytical viscoelastic verification** — Maxwell stress build-up, τ_rel = η/G = 1:
```bash
julia --project=examples examples/ex62b.jl \
  -sol bg -exx_bg 0.1 -eyy_bg -0.1 \
  -mesh_h 0.5 -mu_inclusion 1.0 -G_bg 1.0 -G_inclusion 1.0 \
  -dt 0.1 -nsteps 30 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_pc_type lu \
  -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
  -ksp_rtol 1e-9 -snes_rtol 1e-9
```
The code prints the numerical τ_II alongside the continuous and discrete-BE analytical solutions at each step so convergence can be verified.

### 3D examples

Add `-dm_plex_dim 3` to switch to a sphere-in-cube geometry.  The z strain
rate (`ezz`) is set automatically to `-(exx+eyy)` to enforce incompressibility.

**3D pure shear, algebraic multigrid (GAMG)**:
```bash
julia --project=examples examples/ex62b.jl \
  -dm_plex_dim 3 -sol bg -exx_bg 1.0 -eyy_bg -0.5 -mu_inclusion 0.01 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_ksp_type cg -fieldsplit_velocity_pc_type gamg \
  -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi \
  -ksp_type fgmres -ksp_rtol 1e-6 -snes_rtol 1e-6 -vtk_output out3d.vtu
```

**3D pure shear, geometric multigrid** (coarse mesh h=0.4, 2 uniform refinements → 3 MG levels):
```bash
julia --project=examples examples/ex62b.jl \
  -dm_plex_dim 3 -sol bg -exx_bg 1.0 -eyy_bg -0.5 -mu_inclusion 0.01 \
  -mesh_h 0.4 -dm_refine 2 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_ksp_type fgmres -fieldsplit_velocity_pc_type mg \
  -fieldsplit_velocity_mg_levels_ksp_type chebyshev \
  -fieldsplit_velocity_mg_levels_pc_type sor \
  -fieldsplit_velocity_mg_coarse_pc_type lu \
  -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi \
  -ksp_type fgmres -ksp_rtol 1e-6 -snes_rtol 1e-6 -vtk_output out3d_mg.vtu
```

!!! note "Geometric vs algebraic multigrid"
    With `-dm_refine N`, PETSc builds a refinement hierarchy *inside the preconditioner* only.  The SNES still assembles on the coarse Gmsh mesh (`-mesh_h`), so solution accuracy is at `mesh_h` resolution.  Use a
    smaller `mesh_h` for finer physics; `-dm_refine` only controls the number of MG levels.  GAMG (algebraic MG) requires no mesh hierarchy and is generally preferred for moderate viscosity contrasts.

**3D gravity-driven sinking sphere, viscoelastic, GAMG** (30 time steps):
```bash
julia --project=examples examples/ex62b.jl \
  -dm_plex_dim 3 -sol grav \
  -mesh_h 0.3 -mesh_h_inclusion 0.12 \
  -rho_bg 1.0 -rho_inclusion 2.0 -gravity 1.0 -mu_inclusion 0.1 \
  -G_bg 5.0 -G_inclusion 0.5 -dt 0.02 -nsteps 30 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_ksp_type cg -fieldsplit_velocity_pc_type gamg \
  -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi \
  -ksp_type fgmres -ksp_rtol 1e-6 -snes_rtol 1e-6 -vtk_output sphere_sinking.vtu
```

### Running in parallel with MPI

Every example above runs in parallel without code changes — just prepend
`mpiexec -n N`:

```bash
mpiexec -n 4 julia --project=examples examples/ex62b.jl \
  -dm_plex_dim 3 -sol grav \
  -mesh_h 0.2 -mesh_h_inclusion 0.08 \
  -rho_bg 1.0 -rho_inclusion 2.0 -gravity 1.0 -mu_inclusion 0.1 \
  -G_bg 5.0 -G_inclusion 0.5 -dt 0.02 -nsteps 30 \
  -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
  -pc_fieldsplit_schur_factorization_type full -pc_fieldsplit_schur_precondition a11 \
  -fieldsplit_velocity_ksp_type cg -fieldsplit_velocity_pc_type gamg \
  -fieldsplit_pressure_ksp_type preonly -fieldsplit_pressure_pc_type jacobi \
  -ksp_type fgmres -ksp_rtol 1e-6 -snes_rtol 1e-6 -vtk_output sphere_sinking.vtu
```

The mesh is partitioned automatically by PETSc's DMPlex.  The direct solver
(`-fieldsplit_velocity_pc_type lu`) works in serial only; use GAMG or GMG for
parallel runs.  VTK output is written per-rank and merged into a single XML
collection automatically by `vtk_save!`.

For large HPC runs, see [Running on HPC Systems](@ref) for MPI launch syntax
on different cluster schedulers.

---

## API Reference

```@autodocs
Modules = [PETSc]
Pages   = ["dmplex.jl"]
```
