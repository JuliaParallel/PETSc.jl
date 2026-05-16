# test/dmplex.jl
# Tests for DMPlex (unstructured mesh) functionality.
# Covers both high-level PETSc.* wrappers and low-level LibPETSc.* calls.
# Runs across all available petsclibs (different PetscScalar/PetscInt combos).

using Test
using PETSc
using PETSc: LibPETSc
using MPI

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

# ── PETSc library and scalar types ───────────────────────────────────────────
# @petsc_simple_fn / @petsc_residual_fn / @petsc_jacobian_fn resolve
# PetscInt / PetscScalar / PetscReal from the calling module at macro-expansion
# time, so they MUST be module-level names before the macros are invoked.
# We define them once for the Float64/Int64 lib that owns the FEM callbacks.

const _dmplex_petsclib = PETSc.getlib(PetscScalar = Float64)
PETSc.initialize(_dmplex_petsclib)

const PetscInt    = _dmplex_petsclib.PetscInt
const PetscScalar = Float64
const PetscReal   = Float64

# Serial-safe communicator (shared across all lib loops)
const _TC = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF

# ── Minimal FEM callbacks for testing ────────────────────────────────────────
# Scalar Poisson on the unit square: u = x₁ + x₂ (linear, harmonic).

# f0 body force: −Δu = 0 (u = x₁+x₂ is harmonic)
function _dm_f0(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    f0[1] = 0.0
end
const _dm_f0_ptr = PETSc.@petsc_residual_fn(_dm_f0, 1)

# f1: diffusion flux ∇u
function _dm_f1(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f1)
    for d in 1:dim_
        f1[d] = u_x[d]
    end
end
const _dm_f1_ptr = PETSc.@petsc_residual_fn(_dm_f1, dim_)

# g3: identity stiffness tensor (Laplacian Jacobian)
function _dm_g3(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                aOff, aOff_x, a, a_t, a_x, t, utS, x, nC, cst, g3)
    fill!(g3, 0.0)
    for d in 0:dim_-1
        g3[d * dim_ + d + 1] = 1.0
    end
end
const _dm_g3_ptr = PETSc.@petsc_jacobian_fn(_dm_g3, dim_ * dim_)

# exact solution: u(t, x) = x₁ + x₂  (linear → exact in P1/Q1)
function _dm_exact(t, x, u, ctx)
    u[1] = 0.0
    for d in 1:length(x)
        u[1] += x[d]
    end
end
const _dm_exact_ptr = PETSc.@petsc_simple_fn(_dm_exact)

# zero BC function
function _dm_zero(t, x, u, ctx)
    fill!(u, 0.0)
end
const _dm_zero_ptr = PETSc.@petsc_simple_fn(_dm_zero)

# ─────────────────────────────────────────────────────────────────────────────
@testset "DMPlex" begin

for petsclib in PETSc.petsclibs
    PETSc.initialize(petsclib)
    PetscInt_t    = petsclib.PetscInt
    PetscScalar_t = petsclib.PetscScalar

    @testset "$(PetscScalar_t)/$(PetscInt_t)" begin

    # ── Low-level: DMCreate / DMSetType / DMSetFromOptions ───────────────────
    @testset "Low-level: DMCreate + DMSetType + DMSetFromOptions" begin
        dm = LibPETSc.DMCreate(petsclib, _TC)
        @test dm isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, dm) != C_NULL

        LibPETSc.DMSetType(petsclib, dm, "plex")

        opts = PETSc.Options(petsclib;
                             dm_plex_dim=2, dm_plex_simplex=0,
                             dm_plex_box_faces="4,4")
        push!(opts)
        try
            LibPETSc.DMSetFromOptions(petsclib, dm)
        finally
            pop!(opts)
        end

        @test LibPETSc.DMGetType(petsclib, dm) == "plex"
        @test LibPETSc.DMGetDimension(petsclib, dm) == 2
        @test Bool(LibPETSc.DMPlexIsSimplex(petsclib, dm)) == false

        LibPETSc.DMDestroy(petsclib, dm)
        @test convert(Ptr{Cvoid}, dm) == C_NULL
    end

    # ── Low-level: DMPlexCreateBoxMesh ───────────────────────────────────────
    @testset "Low-level: DMPlexCreateBoxMesh" begin
        PetscReal_t = real(PetscScalar_t)
        faces       = PetscInt_t[4, 4]
        lower       = PetscReal_t[0, 0]
        upper       = PetscReal_t[1, 1]
        periodicity = [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE]

        dm = LibPETSc.DMPlexCreateBoxMesh(
            petsclib, _TC,
            PetscInt_t(2),
            LibPETSc.PetscBool(false),
            faces, lower, upper, periodicity,
            LibPETSc.PetscBool(true),   # interpolate
            PetscInt_t(0),              # localize_height
            LibPETSc.PetscBool(false),  # sparse_localize
        )
        @test dm isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, dm) != C_NULL
        @test LibPETSc.DMGetDimension(petsclib, dm) == 2
        @test Bool(LibPETSc.DMPlexIsSimplex(petsclib, dm)) == false

        # Low-level vector creation
        gvec = LibPETSc.DMCreateGlobalVector(petsclib, dm)
        @test convert(Ptr{Cvoid}, gvec) != C_NULL
        lvec = LibPETSc.DMCreateLocalVector(petsclib, dm)
        @test convert(Ptr{Cvoid}, lvec) != C_NULL

        # Low-level scatter: global → local
        @test_nowarn LibPETSc.DMGlobalToLocal(
            petsclib, dm, gvec, LibPETSc.INSERT_VALUES, lvec)

        # Low-level clone
        dm2 = LibPETSc.PetscDM(petsclib)
        LibPETSc.DMClone(petsclib, dm, dm2)
        @test convert(Ptr{Cvoid}, dm2) != C_NULL
        @test LibPETSc.DMGetDimension(petsclib, dm2) == 2

        # Low-level block size
        @test LibPETSc.DMGetBlockSize(petsclib, dm) == 1

        LibPETSc.VecDestroy(petsclib, gvec)
        LibPETSc.VecDestroy(petsclib, lvec)
        LibPETSc.DMDestroy(petsclib, dm)
        LibPETSc.DMDestroy(petsclib, dm2)
    end

    # ── Low-level: simplex (triangle) box mesh ───────────────────────────────
    # DMPlexCreateBoxMesh with simplex=true requires double-precision geometry;
    # Float32 builds (PetscReal=Float32) fail with PetscError(63).
    if real(PetscScalar_t) != Float32
    @testset "Low-level: DMPlexCreateBoxMesh (simplex)" begin
        PetscReal_t = real(PetscScalar_t)
        faces       = PetscInt_t[4, 4]
        lower       = PetscReal_t[0, 0]
        upper       = PetscReal_t[1, 1]
        periodicity = [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE]

        dm = LibPETSc.DMPlexCreateBoxMesh(
            petsclib, _TC,
            PetscInt_t(2),
            LibPETSc.PetscBool(true),
            faces, lower, upper, periodicity,
            LibPETSc.PetscBool(true),
            PetscInt_t(0),
            LibPETSc.PetscBool(false),
        )
        @test Bool(LibPETSc.DMPlexIsSimplex(petsclib, dm)) == true
        LibPETSc.DMDestroy(petsclib, dm)
    end
    end # real(PetscScalar_t) != Float32

    # ── Low-level: PetscFE creation ──────────────────────────────────────────
    @testset "Low-level: PetscFECreateDefault" begin
        opts = PETSc.Options(petsclib; petscspace_degree=1)
        push!(opts)
        fe = try
            LibPETSc.PetscFECreateDefault(
                petsclib, _TC,
                PetscInt_t(2), PetscInt_t(1),
                LibPETSc.PetscBool(false),    # not simplex
                "",
                PetscInt_t(-1),
            )
        finally
            pop!(opts)
        end
        @test convert(Ptr{Cvoid}, fe) != C_NULL
        # Note: PetscFEDestroy has a broken auto-wrapper; let PETSc manage FE lifetime.
    end

    @testset "Low-level: PetscFECreateLagrange" begin
        fe1 = LibPETSc.PetscFECreateLagrange(
            petsclib, _TC,
            PetscInt_t(2), PetscInt_t(1),
            LibPETSc.PetscBool(false),
            PetscInt_t(1),
            PetscInt_t(-1),
        )
        @test convert(Ptr{Cvoid}, fe1) != C_NULL
        fe2 = LibPETSc.PetscFECreateLagrange(
            petsclib, _TC,
            PetscInt_t(2), PetscInt_t(1),
            LibPETSc.PetscBool(false),
            PetscInt_t(2),
            PetscInt_t(-1),
        )
        @test convert(Ptr{Cvoid}, fe2) != C_NULL
        # PetscFEDestroy auto-wrapper is broken (ReadOnlyMemoryError); FEs are
        # reference-counted by PETSc and freed with the parent DM or at finalize.
    end

    # ── High-level: options-based constructor ────────────────────────────────
    @testset "DMPlex options constructor" begin
        dm = PETSc.DMPlex(petsclib, _TC;
                          dm_plex_dim=2, dm_plex_simplex=0,
                          dm_plex_box_faces="4,4")
        @test dm isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, dm) != C_NULL
        @test LibPETSc.DMGetType(petsclib, dm) == "plex"
        @test LibPETSc.DMGetDimension(petsclib, dm) == 2
    end

    # ── High-level: explicit box constructor ─────────────────────────────────
    @testset "DMPlex explicit box constructor" begin
        dm_hex = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        @test dm_hex isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, dm_hex) != C_NULL
        @test LibPETSc.DMGetDimension(petsclib, dm_hex) == 2

        if real(PetscScalar_t) != Float32
            dm_tri = PETSc.DMPlex(petsclib, _TC, 2, true, [4, 4])
            @test dm_tri isa LibPETSc.PetscDM
            @test LibPETSc.DMGetDimension(petsclib, dm_tri) == 2
        end

        dm_3d = PETSc.DMPlex(petsclib, _TC, 3, false, [2, 2, 2])
        @test dm_3d isa LibPETSc.PetscDM
        @test LibPETSc.DMGetDimension(petsclib, dm_3d) == 3
    end

    # ── isplexsimplex ────────────────────────────────────────────────────────
    @testset "isplexsimplex" begin
        dm_hex = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        @test PETSc.isplexsimplex(dm_hex) == false
        if real(PetscScalar_t) != Float32
            dm_tri = PETSc.DMPlex(petsclib, _TC, 2, true,  [4, 4])
            @test PETSc.isplexsimplex(dm_tri) == true
        end
    end

    # ── plexdistribute! (serial: mesh stays on one rank) ────────────────────
    @testset "plexdistribute! serial" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        # On COMM_SELF (1 rank), DMPlexDistribute returns a DM with NULL ptr
        # (no redistribution needed).  We just check it does not throw.
        dm_par = @test_nowarn PETSc.plexdistribute!(dm; overlap = 0)
        @test dm_par isa LibPETSc.PetscDM
    end

    # ── petsc_setname! ───────────────────────────────────────────────────────
    @testset "petsc_setname!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        @test_nowarn PETSc.petsc_setname!(petsclib, dm, "testmesh")
        @test_nowarn PETSc.petsc_setname!(petsclib, fe, "temperature")
    end

    # ── getlabel ─────────────────────────────────────────────────────────────
    @testset "getlabel" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        # DMPlexCreateBoxMesh always creates a "marker" label for the boundary
        lbl = PETSc.getlabel(dm, "marker")
        @test lbl isa Ptr{Cvoid}
        @test lbl != C_NULL
        # Non-existent label → C_NULL
        lbl_none = PETSc.getlabel(dm, "no_such_label_xyz")
        @test lbl_none == C_NULL
    end

    # ── fe_create_default ────────────────────────────────────────────────────
    @testset "fe_create_default" begin
        fe = PETSc.fe_create_default(petsclib, _TC, 2, 1, false; degree = 1)
        @test convert(Ptr{Cvoid}, fe) != C_NULL
        if real(PetscScalar_t) != Float32
            fe2 = PETSc.fe_create_default(petsclib, _TC, 2, 2, true; degree = 2)
            @test convert(Ptr{Cvoid}, fe2) != C_NULL
        end
    end

    # ── fe_create_lagrange ───────────────────────────────────────────────────
    @testset "fe_create_lagrange" begin
        fe_q1 = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        @test convert(Ptr{Cvoid}, fe_q1) != C_NULL
        fe_q2 = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 2)
        @test convert(Ptr{Cvoid}, fe_q2) != C_NULL
        if real(PetscScalar_t) != Float32
            fe_p1 = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, true, 1)
            @test convert(Ptr{Cvoid}, fe_p1) != C_NULL
        end
        fe_3d = PETSc.fe_create_lagrange(petsclib, _TC, 3, 3, false, 1)
        @test convert(Ptr{Cvoid}, fe_3d) != C_NULL
    end

    # ── setfield! / createds! / getds ────────────────────────────────────────
    @testset "setfield! / createds! / getds" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe, "u")

        @test_nowarn PETSc.setfield!(dm, 0, fe)
        @test_nowarn PETSc.createds!(dm)

        ds = PETSc.getds(dm)
        @test ds isa PETSc.PetscDS
    end

    # ── PetscDS: set_constants! ──────────────────────────────────────────────
    @testset "PetscDS: set_constants!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)
        ds = PETSc.getds(dm)
        @test_nowarn PETSc.set_constants!(ds, [1.0, 2.0, 3.0])
    end

    # ── PetscDS: set_residual! / set_jacobian! / set_exact_solution! ─────────
    @testset "PetscDS: set_residual! / set_jacobian! / set_exact_solution!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)
        ds = PETSc.getds(dm)

        @test_nowarn PETSc.set_residual!(ds, 0, _dm_f0_ptr, _dm_f1_ptr)
        @test_nowarn PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, _dm_g3_ptr)
        @test_nowarn PETSc.set_exact_solution!(ds, 0, _dm_exact_ptr)
    end

    # ── PetscDS: set_jacobian_preconditioner! ────────────────────────────────
    @testset "PetscDS: set_jacobian_preconditioner!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)
        ds = PETSc.getds(dm)
        @test_nowarn PETSc.set_jacobian_preconditioner!(
            ds, 0, 0, C_NULL, C_NULL, C_NULL, _dm_g3_ptr)
    end

    # ── dm_create_global_vec / dm_create_local_vec ───────────────────────────
    @testset "dm_create_global_vec / dm_create_local_vec" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        gvec = PETSc.dm_create_global_vec(dm)
        @test convert(Ptr{Cvoid}, gvec) != C_NULL

        lvec = PETSc.dm_create_local_vec(dm)
        @test convert(Ptr{Cvoid}, lvec) != C_NULL
    end

    # ── dm_global_to_local! ──────────────────────────────────────────────────
    @testset "dm_global_to_local!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        gvec = PETSc.dm_create_global_vec(dm)
        lvec = PETSc.dm_create_local_vec(dm)
        @test_nowarn PETSc.dm_global_to_local!(dm, gvec, lvec)
    end

    # ── dmclone ──────────────────────────────────────────────────────────────
    @testset "dmclone" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        dm2 = PETSc.dmclone(dm)
        @test dm2 isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, dm2) != C_NULL
        # Clone has distinct pointer but same topology
        @test convert(Ptr{Cvoid}, dm2) != convert(Ptr{Cvoid}, dm)
        @test LibPETSc.DMGetDimension(petsclib, dm2) == 2
        @test PETSc.isplexsimplex(dm2) == false
    end

    # ── dm_get_coarse (non-hierarchical mesh → NULL coarse) ──────────────────
    @testset "dm_get_coarse (no hierarchy)" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        cdm = PETSc.dm_get_coarse(dm)
        @test cdm isa LibPETSc.PetscDM
        @test convert(Ptr{Cvoid}, cdm) == C_NULL
    end

    # ── dm_copy_disc! ────────────────────────────────────────────────────────
    @testset "dm_copy_disc!" begin
        dm_src = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe, "u")
        PETSc.setfield!(dm_src, 0, fe)
        PETSc.createds!(dm_src)
        ds_src = PETSc.getds(dm_src)
        PETSc.set_exact_solution!(ds_src, 0, _dm_exact_ptr)

        dm_dst = PETSc.dmclone(dm_src)
        @test_nowarn PETSc.dm_copy_disc!(dm_src, dm_dst)

        # After copy_disc, the destination should also have a DS
        ds_dst = PETSc.getds(dm_dst)
        @test ds_dst isa PETSc.PetscDS
    end

    # ── Multi-field DS setup ──────────────────────────────────────────────────
    @testset "Multi-field DS (2 fields)" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        # field 0: scalar u; field 1: scalar p
        fe_u = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 2)
        fe_p = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe_u, "velocity")
        PETSc.petsc_setname!(petsclib, fe_p, "pressure")
        @test_nowarn PETSc.setfield!(dm, 0, fe_u)
        @test_nowarn PETSc.setfield!(dm, 1, fe_p)
        @test_nowarn PETSc.createds!(dm)

        ds = PETSc.getds(dm)
        @test ds isa PETSc.PetscDS

        gvec = PETSc.dm_create_global_vec(dm)
        @test convert(Ptr{Cvoid}, gvec) != C_NULL
    end

    # ── add_boundary! ────────────────────────────────────────────────────────
    @testset "add_boundary!" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe, "u")
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        label = PETSc.getlabel(dm, "marker")
        @test label != C_NULL

        bd = PETSc.add_boundary!(
            petsclib, dm,
            LibPETSc.DM_BC_ESSENTIAL, "dirichlet", label,
            PetscInt_t[1], 0, PetscInt_t[], _dm_zero_ptr,
        )
        @test bd isa PetscInt_t
        @test bd >= 0
    end

    # ── dm_project_function! + dm_compute_l2diff ─────────────────────────────
    # Only for real Float64 scalar type: callbacks are typed for Float64.
    if PetscScalar_t == Float64
    @testset "dm_project_function! + dm_compute_l2diff" begin
        # 2D hex mesh, 8×8: project u = x₁+x₂ (linear, exact in Q1).
        # The L² error of the projected function vs itself must be ≈ 0.
        dm = PETSc.DMPlex(petsclib, _TC, 2, false, [8, 8])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe, "u")
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        u = PETSc.dm_create_global_vec(dm)
        PETSc.dm_project_function!(petsclib, dm, 0.0,
                                    [_dm_exact_ptr], nothing,
                                    LibPETSc.INSERT_ALL_VALUES, u)

        l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0,
                                         [_dm_exact_ptr], nothing, u)
        @test l2err ≈ 0.0 atol = 1e-12
    end

    # ── dm_project_function! (simplex mesh) ──────────────────────────────────
    @testset "dm_project_function! simplex" begin
        dm = PETSc.DMPlex(petsclib, _TC, 2, true, [8, 8])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 2, 1, true, 1)
        PETSc.petsc_setname!(petsclib, fe, "u")
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        u = PETSc.dm_create_global_vec(dm)
        PETSc.dm_project_function!(petsclib, dm, 0.0,
                                    [_dm_exact_ptr], nothing,
                                    LibPETSc.INSERT_ALL_VALUES, u)

        l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0,
                                         [_dm_exact_ptr], nothing, u)
        @test l2err ≈ 0.0 atol = 1e-12
    end

    # ── 3D hex mesh basic sanity ──────────────────────────────────────────────
    @testset "3D hex mesh + FE + project" begin
        dm = PETSc.DMPlex(petsclib, _TC, 3, false, [4, 4, 4])
        fe = PETSc.fe_create_lagrange(petsclib, _TC, 3, 1, false, 1)
        PETSc.petsc_setname!(petsclib, fe, "u3d")
        PETSc.setfield!(dm, 0, fe)
        PETSc.createds!(dm)

        # _dm_exact sums all coordinate components, works in any dimension.
        u = PETSc.dm_create_global_vec(dm)
        PETSc.dm_project_function!(petsclib, dm, 0.0,
                                    [_dm_exact_ptr], nothing,
                                    LibPETSc.INSERT_ALL_VALUES, u)
        l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0,
                                         [_dm_exact_ptr], nothing, u)
        @test l2err ≈ 0.0 atol = 1e-12
    end
    end # PetscScalar_t == Float64

    end # @testset "$(PetscScalar_t)/$(PetscInt_t)"
end # for petsclib

end # @testset "DMPlex"
