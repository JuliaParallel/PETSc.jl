using Test
using PETSc
using MPI
using LinearAlgebra: diag
using SparseArrays: sparse

isdefined(Main, :PETScTestUtils) || include(joinpath(@__DIR__, "..", "testutils.jl"))
using .PETScTestUtils: without_petsc_traceback

@testset "MatShell" begin
    for petsclib in PETSc.petsclibs
        petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        test_comm = MPI.COMM_SELF

        local_rows = 10
        local_cols = 5
        function f!(p_x, p_y)
            PETSc.with_local_array!((p_x, p_y), write=(true,false)) do x, y
                x .= [2y; 3y]
            end
        end

        matshell =  PETSc.MatShell(f!, petsclib, test_comm, local_rows, local_cols)
        x = PetscScalar.(collect(1:5))
        petsc_x = LibPETSc.VecCreateSeqWithArray(petsclib, test_comm, PetscInt(1), PetscInt(local_cols), x)
        petsc_y = LibPETSc.VecCreateSeqWithArray(petsclib, test_comm, PetscInt(1), PetscInt(local_rows), zeros(PetscScalar, local_rows))
        
        LibPETSc.MatMult(petsclib, matshell, petsc_x, petsc_y)
        @test petsc_y[:] == [2x; 3x]
        @test matshell * x == [2x; 3x]

        PETSc.destroy!(petsc_x)
        PETSc.destroy!(petsc_y)
        PETSc.destroy!(matshell)
        PETSc.finalize(petsclib)
    end
end

@testset "MatShell operations" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        comm = MPI.COMM_SELF
        n = 6
        z = PetscScalar <: Complex ? im : 0
        # Nonsymmetric, diagonally dominant, with a varying diagonal so Jacobi matters
        Ad = PetscScalar[i == j ? 4i : (j == i + 1 ? 1 + z : (j == i - 1 ? -2 : 0)) for i in 1:n, j in 1:n]
        apply!(M) = (y, x) -> PETSc.with_local_array!((y, x); write = (true, false)) do ya, xa
            ya .= M * xa
        end
        calls = Ref(0)
        mult_add = function (zv, xv, yv)
            calls[] += 1
            PETSc.with_local_array!((zv, xv, yv); write = (true, false, false)) do za, xa, ya
                za .= ya .+ Ad * xa
            end
        end
        diagonal = d -> PETSc.with_local_array!(da -> (da .= diag(Ad)), d; write = true)
        A = PETSc.MatShell(apply!(Ad), petsclib, comm, n, n;
                           mult_transpose = apply!(transpose(Ad)), diagonal, mult_add)

        xs = PetscScalar[i + z for i in 1:n]
        ys = PetscScalar[2i for i in 1:n]
        x = PETSc.PetscVec(petsclib, copy(xs))
        y = PETSc.PetscVec(petsclib, copy(ys))

        @test (A * x)[:] ≈ Ad * xs
        @test (transpose(A) * x)[:] ≈ transpose(Ad) * xs
        # Through MatMultTranspose with x conjugated, in a complex build
        @test (A' * x)[:] ≈ Ad' * xs

        d = PETSc.PetscVec(petsclib, n)
        @test PETSc.diagonal!(d, A)[:] == diag(Ad)

        w = PETSc.PetscVec(petsclib, n)
        LibPETSc.MatMultAdd(petsclib, A, x, y, w)
        @test w[:] ≈ ys + Ad * xs
        LibPETSc.MatMultAdd(petsclib, A, x, y, y)  # z may be y
        @test y[:] ≈ ys + Ad * xs
        @test calls[] == 2

        # Jacobi on the shell matches Jacobi on the assembled matrix
        B = PETSc.PetscMat(petsclib, comm, sparse(Ad))
        b = PETSc.PetscVec(petsclib, copy(ys))
        opts = (ksp_type = "gmres", pc_type = "jacobi", ksp_rtol = 1e-10)
        ksp_shell = PETSc.KSP(A; opts...)
        ksp_mat = PETSc.KSP(B; opts...)
        u_shell = PETSc.PetscVec(petsclib, n)
        u_mat = PETSc.PetscVec(petsclib, n)
        PETSc.solve!(u_shell, ksp_shell, b; check = true)
        PETSc.solve!(u_mat, ksp_mat, b; check = true)
        @test PETSc.iteration_number(ksp_shell) == PETSc.iteration_number(ksp_mat)
        @test u_shell[:] ≈ u_mat[:]
        @test u_shell[:] ≈ Ad \ ys

        # An operation left out fails when PETSc needs it
        bare = PETSc.MatShell(apply!(Ad), petsclib, comm, n, n)
        @test_throws LibPETSc.PetscError without_petsc_traceback(() -> transpose(bare) * x, petsclib)
        @test_throws LibPETSc.PetscError without_petsc_traceback(() -> PETSc.diagonal!(d, bare), petsclib)

        # An exception in an operation reaches the caller
        failing = PETSc.MatShell(apply!(Ad), petsclib, comm, n, n;
                                 mult_transpose = (y, x) -> error("transpose failed"))
        @test_throws "transpose failed" without_petsc_traceback(() -> transpose(failing) * x, petsclib)

        foreach(PETSc.destroy!, (ksp_shell, ksp_mat, A, B, bare, failing,
                                 x, y, d, w, b, u_shell, u_mat))
    end
end
