using Test
using PETSc
using MPI
using SparseArrays
using LinearAlgebra
MPI.Initialized() || MPI.Init()

@testset "PC" begin
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar

        # 1D Laplacian, split below into its first and last three unknowns
        n = 6
        S = spdiagm(
            -1 => -ones(PetscScalar, n - 1),
            0 => 2ones(PetscScalar, n),
            1 => -ones(PetscScalar, n - 1),
        )
        b = PetscScalar.(1:n)

        @testset "pc, set_type!, type_name ($PetscScalar)" begin
            ksp = PETSc.KSP(petsclib, comm, S)
            p = PETSc.pc(ksp)
            @test p isa LibPETSc.PC
            @test PETSc.pc(ksp).ptr == p.ptr   # every read is the same PETSc object
            PETSc.set_type!(p, :jacobi)
            @test PETSc.type_name(p) === :jacobi
            @test PETSc.type_name(PETSc.pc(ksp)) === :jacobi
            @test ksp \ b ≈ Matrix(S) \ b
            @test PETSc.type_name(PETSc.pc(ksp)) === :jacobi   # the solve kept it

            # borrowed from ksp (naming.md §3.3): destroy! leaves it, and ksp, usable
            @test !PETSc.owns(p)
            @test PETSc.destroy!(p) === nothing
            @test p.ptr != C_NULL
            @test PETSc.type_name(p) === :jacobi
            @test ksp \ b ≈ Matrix(S) \ b
            PETSc.destroy!(ksp)
        end

        @testset "set_fieldsplit_is! ($PetscScalar)" begin
            ksp = PETSc.KSP(petsclib, comm, S)
            p = PETSc.pc(ksp)
            is_a = LibPETSc.ISCreateStride(petsclib, comm, 3, 0, 1)
            is_b = LibPETSc.ISCreateStride(petsclib, comm, 3, 3, 1)

            # PETSc would silently ignore this on a PC that is not a fieldsplit
            @test_throws ArgumentError PETSc.set_fieldsplit_is!(p, "a", is_a)

            PETSc.set_type!(p, :fieldsplit)
            PETSc.set_fieldsplit_is!(p, "a", is_a)
            PETSc.set_fieldsplit_is!(p, "b", is_b)
            @test ksp \ b ≈ Matrix(S) \ b

            nsplits, subksps = LibPETSc.PCFieldSplitGetSubKSP(petsclib, p)
            @test nsplits == 2
            @test LibPETSc.KSPGetOptionsPrefix(petsclib, subksps[1]) == "fieldsplit_a_"
            @test LibPETSc.KSPGetOptionsPrefix(petsclib, subksps[2]) == "fieldsplit_b_"

            LibPETSc.ISDestroy(petsclib, is_a)
            LibPETSc.ISDestroy(petsclib, is_b)
            PETSc.destroy!(ksp)
        end

        PETSc.finalize(petsclib)
    end
end
