# EXCLUDE FROM TESTING
#=
  ex19.jl — 2D Driven Cavity: velocity–vorticity–temperature formulation
  Port of PETSc snes/tutorials/ex19.c

  Solves on the unit square [0,1]²:
    −∇²u  − ∂ω/∂y               = 0
    −∇²v  + ∂ω/∂x               = 0
    −∇²ω  + ∇·(uω, vω) − Gr·∂T/∂x = 0
    −∇²T  + Pr·∇·(uT, vT)       = 0

  4 DOFs per node: (u, v, ω, T).  Convective terms are upwinded.

  Boundary conditions:
    All walls:   u = v = 0  (no-slip),  except top lid: u = lidvelocity
    Left:        T = 0  (cold, Dirichlet)
    Right:       T = 1  (hot, Dirichlet, only when Gr > 0)
    Top/bottom:  ∂T/∂n = 0  (insulated, Neumann)
    ω:           derived from the no-slip condition at each wall

  Usage:
    julia --project ex19.jl
    julia --project ex19.jl -snes_monitor -ksp_monitor -da_grid_x 32 -da_grid_y 32
    mpiexec -n 4 julia --project ex19.jl -snes_monitor -pc_type mg -da_grid_x 64 -da_grid_y 64

  To switch to GPU (once withlocalarray! supports device arrays):
    Replace  `KernelAbstractions.CPU()`  with  `CUDABackend()`  and adapt the
    array-wrapping inside the residual callback accordingly.
=#

using MPI
using PETSc
using KernelAbstractions

backend = KernelAbstractions.CPU()


# ── Physical parameters ──────────────────────────────────────────────────────
Base.@kwdef mutable struct AppCtx{T}
    lidvelocity :: T = one(T)
    prandtl     :: T = one(T)
    grashof     :: T = one(T)
end

# ── Residual kernel ───────────────────────────────────────────────────────────
#
# Iterates over the locally-owned grid: (li, lj) ∈ 1..nx_own × 1..ny_own.
#
# Array layout (plain 1-based, no OffsetArray — GPU-compatible):
#   x_par[dof, xi, xj]   ghost input,  xi = li + ox, xj = lj + oy
#   f_par[dof, li, lj]   owned output
#
# ox = xs - xsg, oy = ys - ysg  (ghost offsets; 0 at a domain boundary, 1 otherwise)
#
# Corner priority matches the C code: left/right walls processed last, so they
# take precedence over bottom/top at corners.  Achieved here by checking i==1
# and i==mx before j==1 and j==my.
#
@kernel function cavity_residual_kernel!(
    f_par,
    x_par,
    dhx, dhy,           # mx-1, my-1
    hx,  hy,            # 1/dhx, 1/dhy
    hydhx, hxdhy,       # metric factors for Laplacian scaling
    grashof, prandtl, lid,
    mx :: Int, my :: Int,
    xs :: Int, ys :: Int,
    ox :: Int, oy :: Int,
)
    li, lj = @index(Global, NTuple)

    i  = xs + li - 1   # global 1-based x-coordinate
    j  = ys + lj - 1   # global 1-based y-coordinate
    xi = li + ox        # ghost-array x-index for this point
    xj = lj + oy        # ghost-array y-index for this point

    # ── Boundary conditions ───────────────────────────────────────────────────

    if i == 1                           # left wall — cold (T = 0)
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] -
            (x_par[2, xi+1, xj] - x_par[2, xi, xj]) * dhx
        f_par[4, li, lj] = x_par[4, xi,   xj  ]

    elseif i == mx                      # right wall — hot (T = 1 when Gr > 0)
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] -
            (x_par[2, xi, xj] - x_par[2, xi-1, xj]) * dhx
        f_par[4, li, lj] = x_par[4, xi,   xj  ] -
            (grashof > zero(grashof) ? one(grashof) : zero(grashof))

    elseif j == 1                       # bottom wall — no-slip, insulated
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] +
            (x_par[1, xi, xj+1] - x_par[1, xi, xj]) * dhy
        f_par[4, li, lj] = x_par[4, xi,   xj  ] - x_par[4, xi, xj+1]

    elseif j == my                      # top wall — moving lid, insulated
        f_par[1, li, lj] = x_par[1, xi,   xj  ] - lid
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] +
            (x_par[1, xi, xj] - x_par[1, xi, xj-1]) * dhy
        f_par[4, li, lj] = x_par[4, xi,   xj  ] - x_par[4, xi, xj-1]

    else                                # ── interior point ───────────────────

        # Upwind split of advecting velocities
        vx  = x_par[1, xi, xj];  avx = abs(vx)
        vxp = oftype(vx, 0.5) * (vx + avx)   # max(vx, 0)
        vxm = oftype(vx, 0.5) * (vx - avx)   # min(vx, 0)

        vy  = x_par[2, xi, xj];  avy = abs(vy)
        vyp = oftype(vy, 0.5) * (vy + avy)
        vym = oftype(vy, 0.5) * (vy - avy)

        # u-equation:  −∇²u − ∂ω/∂y = 0
        cu  = x_par[1, xi, xj]
        uxx = (2cu - x_par[1, xi-1, xj] - x_par[1, xi+1, xj]) * hydhx
        uyy = (2cu - x_par[1, xi, xj-1] - x_par[1, xi, xj+1]) * hxdhy
        f_par[1, li, lj] = uxx + uyy -
            oftype(cu, 0.5) * (x_par[3, xi, xj+1] - x_par[3, xi, xj-1]) * hx

        # v-equation:  −∇²v + ∂ω/∂x = 0
        cv  = x_par[2, xi, xj]
        vxx = (2cv - x_par[2, xi-1, xj] - x_par[2, xi+1, xj]) * hydhx
        vyy = (2cv - x_par[2, xi, xj-1] - x_par[2, xi, xj+1]) * hxdhy
        f_par[2, li, lj] = vxx + vyy +
            oftype(cv, 0.5) * (x_par[3, xi+1, xj] - x_par[3, xi-1, xj]) * hy

        # ω-equation:  −∇²ω + ∇·(uω, vω) − Gr·∂T/∂x = 0
        cω  = x_par[3, xi, xj]
        wxx = (2cω - x_par[3, xi-1, xj] - x_par[3, xi+1, xj]) * hydhx
        wyy = (2cω - x_par[3, xi, xj-1] - x_par[3, xi, xj+1]) * hxdhy
        f_par[3, li, lj] = wxx + wyy +
            (vxp * (cω - x_par[3, xi-1, xj]) + vxm * (x_par[3, xi+1, xj] - cω)) * hy +
            (vyp * (cω - x_par[3, xi, xj-1]) + vym * (x_par[3, xi, xj+1] - cω)) * hx -
            oftype(cω, 0.5) * grashof * (x_par[4, xi+1, xj] - x_par[4, xi-1, xj]) * hy

        # T-equation:  −∇²T + Pr·∇·(uT, vT) = 0
        cT  = x_par[4, xi, xj]
        txx = (2cT - x_par[4, xi-1, xj] - x_par[4, xi+1, xj]) * hydhx
        tyy = (2cT - x_par[4, xi, xj-1] - x_par[4, xi, xj+1]) * hxdhy
        f_par[4, li, lj] = txx + tyy + prandtl * (
            (vxp * (cT - x_par[4, xi-1, xj]) + vxm * (x_par[4, xi+1, xj] - cT)) * hy +
            (vyp * (cT - x_par[4, xi, xj-1]) + vym * (x_par[4, xi, xj+1] - cT)) * hx)
    end
end

# ── Setup ─────────────────────────────────────────────────────────────────────
opts = isinteractive() ? NamedTuple() : PETSc.parse_options(ARGS)

petsclib = PETSc.getlib(; PetscScalar = Float64)
PETSc.initialize(petsclib)

T       = Float64
PetscInt = petsclib.PetscInt
comm    = MPI.COMM_WORLD

# DMDA: 4×4 default (matches ex19.c); override via -da_grid_x / -da_grid_y
da = PETSc.DMDA(
    petsclib, comm,
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
    (4, 4),
    4,   # DOFs per node: (u, v, ω, T)
    1,   # stencil width
    PETSc.DMDA_STENCIL_STAR;
    opts...,
)

snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, da)

# Actual grid size after setfromoptions (may differ from the 4×4 default)
info = PETSc.getinfo(da)
mx   = info.global_size[1]
my   = info.global_size[2]

user = AppCtx{T}(
    lidvelocity = T(1) / (mx - 1),
    prandtl     = T(1),
    grashof     = T(1),
)

# Precomputed grid metrics
dhx   = T(mx - 1);   dhy   = T(my - 1)
hx    = one(T) / dhx; hy    = one(T) / dhy
hydhx = hy * dhx;    hxdhy = hx * dhy

# ── Initial condition: u = v = ω = 0, T linear in x ─────────────────────────
x = PETSc.DMGlobalVec(da)

PETSc.withlocalarray!(x; read = false) do x_arr
    corners = PETSc.getcorners(da)
    xs = corners.lower[1];  ys = corners.lower[2]
    xe = corners.upper[1];  ye = corners.upper[2]
    nx_own = xe - xs + 1;   ny_own = ye - ys + 1
    dx = one(T) / (mx - 1)
    x_par = reshape(x_arr, 4, nx_own, ny_own)
    for lj in 1:ny_own, li in 1:nx_own
        ig = xs + li - 1
        x_par[1, li, lj] = zero(T)
        x_par[2, li, lj] = zero(T)
        x_par[3, li, lj] = zero(T)
        x_par[4, li, lj] = user.grashof > 0 ? T(ig - 1) * dx : zero(T)
    end
end

# ── Residual callback ─────────────────────────────────────────────────────────
#
# To run on GPU, replace CPU() with CUDABackend() (or ROCBackend()) and adapt
# the array wrapping once withlocalarray! supports device arrays (see vec.jl).
#
r       = similar(x)

PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    da = PETSc.getDM(snes)

    l_x = PETSc.DMLocalVec(da)
    PETSc.dm_global_to_local!(g_x, l_x, da, PETSc.INSERT_VALUES)

    PETSc.withlocalarray!(
        (g_fx, l_x);
        read  = (false, true),
        write = (true,  false),
    ) do fx, lx
        corners       = PETSc.getcorners(da)
        ghost_corners = PETSc.getghostcorners(da)

        xs  = corners.lower[1];        ys  = corners.lower[2]
        xe  = corners.upper[1];        ye  = corners.upper[2]
        xsg = ghost_corners.lower[1];  ysg = ghost_corners.lower[2]
        xeg = ghost_corners.upper[1];  yeg = ghost_corners.upper[2]

        nx_own = xe  - xs  + 1;  ny_own = ye  - ys  + 1
        nx_g   = xeg - xsg + 1;  ny_g   = yeg - ysg + 1

        # Plain [dof, x, y] arrays — no OffsetArray, safe for KA on GPU
        x_par = reshape(lx, 4, nx_g,   ny_g)
        f_par = reshape(fx, 4, nx_own, ny_own)

        # Ghost offset: ghost-array index for owned start = 1 + ox (0 at domain wall)
        ox = xs - xsg
        oy = ys - ysg

        cavity_residual_kernel!(backend, 64)(
            f_par, x_par,
            dhx, dhy, hx, hy, hydhx, hxdhy,
            user.grashof, user.prandtl, user.lidvelocity,
            mx, my, xs, ys, ox, oy;
            ndrange = (nx_own, ny_own),
        )
        KernelAbstractions.synchronize(backend)
    end

    PETSc.destroy(l_x)
    return PetscInt(0)
end

# ── Jacobian (finite differences via PETSc's built-in column-by-column FD) ───
#
# Pass SNESComputeJacobianDefault directly as the C function pointer, exactly
# like the C code does:
#   SNESSetJacobian(snes, J, J, SNESComputeJacobianDefault, NULL)
# This avoids a nested Julia→C→Julia callback chain and is more robust
# in parallel.  For production, replace with coloring-based FD by swapping
# SNESComputeJacobianDefault → SNESComputeJacobianDefaultColor.
#
J = LibPETSc.DMCreateMatrix(petsclib, da)
LibPETSc.SNESSetJacobian(petsclib, snes, J, J,
    cglobal((:SNESComputeJacobianDefault, petsclib.petsc_library)), C_NULL)

# ── Solve ─────────────────────────────────────────────────────────────────────
PETSc.solve!(x, snes)

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iterations.")
end

# ── Cleanup ───────────────────────────────────────────────────────────────────
# Explicitly destroy the PetscOptions stored on the SNES before finalization.
# Its GC finalizer calls PetscOptionsDestroy, which can use MPI internally.
# If GC runs it after MPI is alive but in a different collective-sync state
# across ranks, it triggers intermittent crashes.  Destroying it explicitly
# here (while all ranks are synchronized and PETSc/MPI are still fully active)
# is safe and prevents any later GC-driven call.
if !isnothing(snes.opts)
    PETSc.destroy(snes.opts)
    snes.opts = nothing
end

# Run a full GC now so any lingering VecRestoreArray finalizers from
# withlocalarray! run while PETSc is still valid, then barrier all ranks.
GC.gc(true)
MPI.Barrier(comm)

# SNES holds internal PETSc references to J, da, and r — destroy it first so
# those reference counts are decremented before we explicitly free the objects.
PETSc.destroy(snes)
PETSc.destroy(J)
PETSc.destroy(x)
PETSc.destroy(r)
PETSc.destroy(da)
PETSc.finalize(petsclib)

# On macOS ARM64 with MPICH ch4:ofi, MPICH's C atexit handler crashes during
# process teardown (SIGSEGV in libfabric/OFI cleanup).  Using quick_exit(0)
# after explicitly finalizing PETSc and MPI bypasses all C atexit() handlers
# (while still running at_quick_exit() handlers) and avoids the crash.
# All MPI communication is already complete at this point.
MPI.Barrier(comm)
MPI.Finalize()
ccall(:quick_exit, Cvoid, (Cint,), 0)
