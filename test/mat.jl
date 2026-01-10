using Test
using PETSc, MPI
using LinearAlgebra: norm, mul!, Adjoint, Transpose, issymmetric, ishermitian
using SparseArrays: sprand, spdiagm
using Random

MPI.Initialized() || MPI.Init()
comm = MPI.COMM_WORLD


@testset "MatSeqAIJ" begin
    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[7]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt    

        nz_int = PetscInt(2)
        num_rows, num_cols = PetscInt(5), PetscInt(7)
        nz_vec = petsclib.PetscInt.([0, 3, 2, 1, 0])

        A = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)

        @test LibPETSc.MatGetSize(petsclib, A) == (num_rows, num_cols)
        @test size(A) == (num_rows, num_cols)
        
        B = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)

        #PETSc.assemble!(A)
        LibPETSc.MatAssemblyBegin(petsclib, A, PETSc.MAT_FINAL_ASSEMBLY)
        LibPETSc.MatAssemblyEnd(petsclib, A, PETSc.MAT_FINAL_ASSEMBLY)

        PETSc.assemble!(B)
        

        C = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, PetscInt(0), nz_vec)
        #@test size(A) == (num_rows, num_cols)
        @test size(C) == (num_rows, num_cols)

        PETSc.assemble!(C)
        
        # both empty
        @test A == C

        if PetscScalar <: Real
            D = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)

            # NOTE: there appears to be an issue in displaying this in the REPL
            D[1, [1, 2]] = PetscScalar.([1, 2]);
            D[2, [3, 4]] = PetscScalar.([3, 4]);
            D[5, [3, 4]] = PetscScalar.([3, 4]);
            PETSc.assemble!(D);

            DJ = zeros(PetscScalar, num_rows, num_cols);
            DJ[1, [1, 2]] .= PetscScalar.([1, 2]);
            DJ[2, [3, 4]] .= PetscScalar.([3, 4]);
            DJ[5, [3, 4]] .= PetscScalar.([3, 4]);
        else
            D = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)

            D[1, [1, 2]] = PetscScalar.([1, 2im]);
            D[2, [3, 4]] = PetscScalar.([3, 4im]);
            D[5, [3, 4]] = PetscScalar.([3, 4im]);
            PETSc.assemble!(D);

            DJ = zeros(PetscScalar, num_rows, num_cols);
            DJ[1, [1, 2]] .= PetscScalar.([1, 2im]);
            DJ[2, [3, 4]] .= PetscScalar.([3, 4im]);
            DJ[5, [3, 4]] .= PetscScalar.([3, 4im]);
        end

        E = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, PetscInt(0), nz_vec)
        E[2, [1, 3, 4]] = PetscScalar.([1, 3, 4]);
        E[3, [3, 4]] = PetscScalar.([3, 4]);
        E[4, [5]] = PetscScalar.([6]);
        PETSc.assemble!(E);

     
        # Test norm
        @test norm(A) == 0
        @test norm(D) ≈ norm(DJ)

        x = PetscScalar.(Array(1:num_cols))
        y = zeros(PetscScalar, num_rows)
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), copy(x))
        vec_y = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(y)), copy(y))

        # Test mul!
        mul!(vec_y, D, vec_x)
        y = DJ * x
        @test all(vec_y[1:num_rows] .≈ y)
        #@test all(DJ * x .≈ D * x)

        # Zero out vec_x before the transpose multiplication to avoid stale data
        fill!(vec_x, zero(PetscScalar))
        mul!(vec_x, Transpose(D), vec_y)
        x = Transpose(DJ) * y
        @test all(vec_x[1:end] .≈ x)

        # test issymmetric and ishermitian
        if PetscScalar <: Real
            A = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            A[1, 1] = 1;
            A[2, 1] = -2;
            A[1, 2] = -2;
            PETSc.assemble!(A);

            B = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            B[1, 1] = 1;
            B[2, 1] = 2;
            B[1, 2] = -2;
            PETSc.assemble!(B);

            @test issymmetric(A)
            @test ishermitian(A)
            @test !issymmetric(B)
            @test !ishermitian(B)
        else
            A = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            A[1, 1] = PetscScalar(1)
            A[2, 1] = PetscScalar(-2 + im)
            A[1, 2] = PetscScalar(-2 - im)
            PETSc.assemble!(A)

            B = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            B[1, 1] = PetscScalar(1)
            B[2, 1] = PetscScalar(-2 + im)
            B[1, 2] = PetscScalar(-2 - im)
            PETSc.assemble!(B)
            @test !issymmetric(A)

            # TODO: fix hermitian for complex matrix
          #  @test ishermitian(A)
            @test !issymmetric(B)
          #  @test !ishermitian(B)
        end

        Random.seed!(777)
        A1 = sprand(PetscScalar, 10, 11, 0.2)
        B1 = PETSc.MatSeqAIJWithArrays(petsclib, comm, A1)
        sleep(0.1) # seems to help the tests
        @test sum(B1[:,:] - Matrix(A1)) == 0.0

        PETSc.destroy(A)
        PETSc.destroy(B)
        PETSc.destroy(B1)
        PETSc.destroy(C)
        PETSc.destroy(D)
        PETSc.destroy(E) 
        PETSc.destroy(vec_x)
        PETSc.destroy(vec_y)   
        
        PETSc.finalize(petsclib)
    end
end


@testset "MatSeqDense" begin
    for petsclib in PETSc.petsclibs
    #    petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        Ajl = PetscScalar.([
            1 2 3 4
            5 6 7 8
            9 10 11 12
            13 14 15 16
        ])

        x = PetscScalar.(collect(1:4))

        A = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(size(Ajl,1)), PetscInt(size(Ajl,2)), Ajl[:])
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), copy(x))

        @test sum(A[:,:] - Ajl) == 0.0
        #@test all(A * x .≈ Ajl * x)

        PETSc.destroy(A)
        PETSc.destroy(vec_x)
        PETSc.finalize(petsclib)
    end
end

@testset "MatSeqAIJ_Sparse" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        n = 10

        A_sp = spdiagm(
            -1 => ones(PetscScalar, n - 1),
            0 => -2ones(PetscScalar, n),
            1 => ones(PetscScalar, n - 1),
        )

        A = PETSc.MatCreateSeqAIJ(petsclib, comm, A_sp)
       
        @test sum(A[1:10,1:10] - Matrix(A_sp)) == 0.0 == 0.0
        PETSc.destroy(A)
        PETSc.finalize(petsclib)
    end
end