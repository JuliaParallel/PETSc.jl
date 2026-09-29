# test/core/performance.jl
# The performance contract for the calls that go inside a loop: each one infers a
# concrete return type, and its allocation count is the same at two grid sizes, so
# the cost per call does not grow with the problem. The budgets are upper bounds
# with some room, not the measured numbers; a call that drops below one can have
# its budget lowered, but a call that grows past one is a regression.

using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

@testset "allocations do not grow with the problem" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscScalar = petsclib.PetscScalar
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc

    # One case per call: a name, a budget, and a closure built from the objects of
    # a problem of size `n`. `setup` returns the closure and everything to destroy.
    function cases(n)
        m = round(Int, sqrt(n))
        none2 = (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE)
        dm = PETSc.DMStag(petsclib, comm, none2, (m, m), (0, 1, 1), 1)
        x = PETSc.PetscVec(petsclib, comm, rand(PetscScalar, n))
        y = PETSc.PetscVec(petsclib, comm, rand(PetscScalar, n))
        S = spdiagm(-1 => -ones(PetscScalar, n - 1), 0 => 2ones(PetscScalar, n),
                    1 => -ones(PetscScalar, n - 1))
        A = PETSc.PetscMat(petsclib, comm, S)
        B = PETSc.PetscMat(petsclib, comm, S)   # written to, so left unassembled
        ksp = PETSc.KSP(A; ksp_rtol = 1e-8, pc_type = "jacobi")
        gl = PETSc.global_vec(dm)
        lo = PETSc.local_vec(dm)
        J = PETSc.PetscMat(dm)
        e = PETSc.element_location(dm)
        I = CartesianIndex(1, 1)
        rows = [PETSc.stencil(dm, e, (1, 1)), PETSc.stencil(dm, e, (2, 1))]
        vals = PetscScalar[1, 2, 3, 4]
        fields = (PETSc.face_location(dm, 1) => 0, PETSc.element_location(dm) => 0)

        calls = (
            ("fill!",                  0, () -> fill!(x, 2)),
            ("copyto!",                0, () -> copyto!(x, y)),
            ("norm",                   0, () -> norm(y)),
            ("length",                 0, () -> length(x)),
            ("assemble! a Vec",        0, () -> PETSc.assemble!(x)),
            ("ownership_range",        0, () -> PETSc.ownership_range(x)),
            ("v[i] = a",               0, () -> (x[2] = one(PetscScalar))),
            ("v[i]",                   4, () -> x[2]),
            ("with_local_array!",      4, () -> PETSc.with_local_array!(a -> a[1], x)),
            ("with_local_array! read", 4, () -> PETSc.with_local_array!(a -> a[1], x; write = false)),
            ("with_local_array! x2",   6, () -> PETSc.with_local_array!((a, b) -> a[1] + b[1], x, y)),
            ("broadcast into a Vec",   8, () -> (x .= 2 .* y .+ 1; nothing)),
            ("size(A)",                0, () -> size(A)),
            ("A[i, j] = a",            0, () -> (B[2, 2] = one(PetscScalar))),
            ("isassembled",            0, () -> PETSc.isassembled(A)),
            ("mul!",                   2, () -> mul!(x, A, y)),
            ("KSP solve!",             2, () -> PETSc.solve!(x, ksp, y)),
            ("KSP iteration_number",   0, () -> PETSc.iteration_number(ksp)),
            ("KSP converged_reason",   0, () -> PETSc.converged_reason(ksp)),
            ("corners",                0, () -> PETSc.corners(dm)),
            ("ghost_corners",          0, () -> PETSc.ghost_corners(dm)),
            ("stencil",                0, () -> PETSc.stencil(dm, e, I)),
            ("dof_slot",               0, () -> PETSc.dof_slot(dm, e, 0)),
            ("set_values! a Mat",      0, () -> PETSc.set_values!(J, dm, rows, rows, vals, LibPETSc.ADD_VALUES)),
            ("set_values! a Vec",      0, () -> PETSc.set_values!(gl, dm, rows, view(vals, 1:2))),
            ("global_to_local!",       0, () -> PETSc.global_to_local!(lo, dm, gl)),
            ("local_to_local!",        0, () -> PETSc.local_to_local!(lo, dm)),
            ("with_field_views!",     12, () -> PETSc.with_field_views!(v -> v[1][I], dm, lo; fields)),
            ("with_field_views! all", 12, () -> PETSc.with_field_views!(a -> a[I, 1], dm, lo)),
        )
        return calls, (J, lo, gl, ksp, B, A, y, x, dm)
    end

    # `f` twice to compile and settle, then its allocation count
    function allocations(f)
        f()
        f()
        return @allocations f()
    end

    small, small_objs = cases(64)
    large, large_objs = cases(4096)
    @testset "$name" for ((name, budget, f), (_, _, g)) in zip(small, large)
        @test (@inferred f()) === (@inferred f())   # concrete, and a second call agrees
        a, b = allocations(f), allocations(g)
        @test a == b            # the cost does not grow with the problem
        @test a <= budget
    end
    foreach(PETSc.destroy!, small_objs)
    foreach(PETSc.destroy!, large_objs)
    PETSc.finalize(petsclib)
end
