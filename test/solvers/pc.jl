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

        @testset "shell ($PetscScalar)" begin
            ksp = PETSc.KSP(petsclib, comm, S; ksp_rtol = 1e-4)
            @test_throws ArgumentError PETSc.set_shell_apply!((y, p, x) -> nothing, PETSc.pc(ksp))
            @test_throws ArgumentError PETSc.set_shell_setup!(p -> nothing, PETSc.pc(ksp))
            PETSc.set_type!(PETSc.pc(ksp), :shell)

            # Jacobi by hand: the diagonal of S is 2. Registered through handles
            # that are dropped at once, so only the PETSc object keeps the closures.
            napply = Ref(0)
            nsetup = Ref(0)
            PETSc.set_shell_apply!(PETSc.pc(ksp)) do y, p, x
                @test p isa LibPETSc.PC
                napply[] += 1
                PETSc.with_local_array!(y, x; read = (false, true), write = (true, false)) do ya, xa
                    ya .= xa ./ 2
                end
                return nothing
            end
            q = PETSc.set_shell_setup!(PETSc.pc(ksp)) do p
                nsetup[] += 1
                return nothing
            end
            @test q isa LibPETSc.PC && q.ptr == PETSc.pc(ksp).ptr
            GC.gc()
            @test ksp \ b ≈ Matrix(S) \ b rtol = 1e-3
            @test napply[] > 0
            @test nsetup[] == 1

            # a second apply! replaces the first, and the setup! stays
            nreplaced = Ref(0)
            PETSc.set_shell_apply!(PETSc.pc(ksp)) do y, p, x
                nreplaced[] += 1
                PETSc.with_local_array!(y, x; read = (false, true), write = (true, false)) do ya, xa
                    ya .= xa ./ 2
                end
                return nothing
            end
            napply_before = napply[]
            @test ksp \ b ≈ Matrix(S) \ b rtol = 1e-3
            @test nreplaced[] > 0
            @test napply[] == napply_before

            # an exception in the callback comes out of the solve as itself (naming.md §18.4)
            PETSc.set_shell_apply!(PETSc.pc(ksp)) do y, p, x
                throw(DomainError(-1.0, "shell apply! failed on purpose"))
            end
            @test_throws DomainError ksp \ b
            # started through LibPETSc there is no high-level call to rethrow it
            petsc_b = LibPETSc.VecCreateSeqWithArray(petsclib, comm, 1, n, b)
            petsc_x = similar(petsc_b)
            @test_logs (:error, r"shell apply!") match_mode = :any begin
                @test_throws LibPETSc.PetscError LibPETSc.KSPSolve(petsclib, ksp, petsc_b, petsc_x)
            end
            PETSc.destroy!(petsc_b)
            PETSc.destroy!(petsc_x)

            # the closures die with the PETSc object, not with a wrapper
            state = PETSc.object_state(PETSc.pc(ksp))
            @test state isa PETSc.PCState && state.alive
            PETSc.destroy!(ksp)
            @test !state.alive
            PETSc.finalize(petsclib)
            @test !haskey(PETSc.object_states, typeof(petsclib))
            PETSc.initialize(petsclib)
        end

        @testset "ownership ($PetscScalar)" begin
            # made through LibPETSc, the PC is owned and destroy! frees it
            p = LibPETSc.PCCreate(petsclib, comm)
            @test PETSc.owns(p)
            PETSc.destroy!(p)
            @test p.ptr == C_NULL
        end

        PETSc.finalize(petsclib)
    end
end
