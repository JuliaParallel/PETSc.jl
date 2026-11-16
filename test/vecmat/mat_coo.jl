# test/vecmat/mat_coo.jl
# COO assembly of a matrix from a DM: set_matrix_preallocate_skip!,
# set_preallocation_coo_local!, set_preallocation_coo! and set_values_coo! give the
# same matrix as set_values! with stencils. Runs on any number of ranks.

using Test
using PETSc, MPI
MPI.Initialized() || MPI.Init()

# The 5-point Laplacian on the points this rank owns, as COO triplets in 0-based
# local (ghosted) indices, and as stencil rows for set_values!
function laplacian_coo(da)
    c = PETSc.corners(da)
    gc = PETSc.ghost_corners(da)
    nx, ny = size(da)
    # 0-based local index of the 1-based global grid point (i, j)
    local_0b(i, j) = (i - gc.lower[1]) + (j - gc.lower[2]) * gc.size[1]
    rows, cols, vals = Int[], Int[], Float64[]
    stencil_rows = Tuple{LibPETSc.MatStencil, Vector{LibPETSc.MatStencil}, Vector{Float64}}[]
    for j in c.lower[2]:c.upper[2], i in c.lower[1]:c.upper[1]
        nbrs = [(i, j, 4.0)]
        for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
            p, q = i + di, j + dj
            1 <= p <= nx && 1 <= q <= ny && push!(nbrs, (p, q, -1.0))
        end
        for (p, q, v) in nbrs
            push!(rows, local_0b(i, j)); push!(cols, local_0b(p, q)); push!(vals, v)
        end
        push!(stencil_rows, (
            LibPETSc.MatStencil(0, j - 1, i - 1, 0),
            [LibPETSc.MatStencil(0, q - 1, p - 1, 0) for (p, q, _) in nbrs],
            [v for (_, _, v) in nbrs],
        ))
    end
    return rows, cols, vals, stencil_rows
end

@testset "COO assembly on a DMDA" begin
    comm = MPI.COMM_WORLD
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    none2 = (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE)
    da = PETSc.DMDA(petsclib, comm, none2, (7, 6), 1; stencil_width = 1, stencil_type = :star)
    rows, cols, vals, stencil_rows = laplacian_coo(da)

    # the reference, before the DM skips preallocation for every later matrix
    A_ref = PETSc.PetscMat(da)
    for (row, scols, svals) in stencil_rows
        PETSc.set_values!(A_ref, [row], scols, svals)
    end
    PETSc.assemble!(A_ref)

    @test PETSc.set_matrix_preallocate_skip!(da, true) === da
    A = PETSc.PetscMat(da)
    @test PETSc.set_preallocation_coo_local!(A, rows, cols) === A
    @test PETSc.set_values_coo!(A, vals) === A
    @test PETSc.isassembled(A)
    @test A == A_ref

    # the values of a repeated pattern add up; ADD_VALUES adds to the entries
    B = PETSc.PetscMat(da)
    PETSc.set_preallocation_coo_local!(B, vcat(rows, rows), vcat(cols, cols))
    PETSc.set_values_coo!(B, vcat(vals, vals))
    PETSc.set_values_coo!(B, vcat(vals, vals); mode = PETSc.ADD_VALUES)
    r = first(PETSc.ownership_range(B))
    @test B[r, r] == 16

    # on one rank the local and the global numbering agree
    if MPI.Comm_size(comm) == 1
        G = PETSc.PetscMat(da)
        PETSc.set_preallocation_coo!(G, rows, cols)
        PETSc.set_values_coo!(G, vals)
        @test G == A_ref
        PETSc.destroy!(G)
    end

    @test_throws "as many column indices as row indices" PETSc.set_preallocation_coo!(
        B, [0, 1], [0])

    # the COO matrix solves as the reference does
    b = PETSc.global_vec(da)
    fill!(b, 1)
    ksp = PETSc.KSP(A; ksp_rtol = 1e-10)
    ksp_ref = PETSc.KSP(A_ref; ksp_rtol = 1e-10)
    x = ksp \ b
    x_ref = ksp_ref \ b
    @test x ≈ x_ref

    foreach(PETSc.destroy!, (x_ref, x, ksp_ref, ksp, b, B, A, A_ref, da))
    PETSc.finalize(petsclib)
end
