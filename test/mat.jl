using Test
using PETSc
using LinearAlgebra: norm, mul!, Adjoint, Transpose, issymmetric, ishermitian
using SparseArrays: sprand
using Random

@testset "MatSeqAIJ" begin
    num_rows, num_cols = 5, 7
    nz_int = 2
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        nz_vec = petsclib.PetscInt.([0, 3, 2, 1, 0])

        A = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_int)

        @test size(A) == (num_rows, num_cols)
        B = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_int)
        PETSc.assemble!(A)
        PETSc.assemble!(B)
        # both empty
        @test A == B

        C = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_vec)
        @test size(A) == (num_rows, num_cols)
        PETSc.assemble!(C)
        # both empty
        @test A == C

        if PetscScalar <: Real
            D = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_int)
            D[1, [1, 2]] .= [1, 2]
            D[2, [3, 4]] .= [3, 4]
            D[5, [3, 4]] .= [3, 4]
            PETSc.assemble!(D)

            DJ = zeros(PetscScalar, num_rows, num_cols)
            DJ[1, [1, 2]] .= [1, 2]
            DJ[2, [3, 4]] .= [3, 4]
            DJ[5, [3, 4]] .= [3, 4]
        else
            D = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_int)
            D[1, [1, 2]] .= [1, 2im]
            D[2, [3, 4]] .= [3, 4im]
            D[5, [3, 4]] .= [3, 4im]
            PETSc.assemble!(D)

            DJ = zeros(PetscScalar, num_rows, num_cols)
            DJ[1, [1, 2]] .= [1, 2im]
            DJ[2, [3, 4]] .= [3, 4im]
            DJ[5, [3, 4]] .= [3, 4im]
        end

        E = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_vec)
        E[2, [1, 3, 4]] .= [1, 3, 4]
        E[3, [3, 4]] .= [3, 4]
        E[4, [5]] .= [6]
        PETSc.assemble!(E)

        @test A != E
        @test D != E

        # Test norm
        @test norm(A) == 0
        @test norm(D) ≈ norm(DJ)

        x = PetscScalar.(Array(1:num_cols))
        y = zeros(PetscScalar, num_rows)
        vec_x = PETSc.VecSeqWithArray(petsclib, copy(x))
        vec_y = PETSc.VecSeqWithArray(petsclib, copy(y))

        # Test mul!
        mul!(vec_y, D, vec_x)
        y = DJ * x
        @test all(vec_y.array .≈ y)
        @test all(DJ * x .≈ D * x)

        mul!(vec_x, Adjoint(D), vec_y)
        x = Adjoint(DJ) * y
        @test all(vec_x.array .≈ x)

        mul!(vec_x, Transpose(D), vec_y)
        x = Transpose(DJ) * y
        @test all(vec_x.array .≈ x)

        # test issymmetric and ishermitian
        if PetscScalar <: Real
            A = PETSc.MatSeqAIJ(petsclib, 5, 5, 2)
            A[1, 1] = 1
            A[2, 1] = -2
            A[1, 2] = -2
            PETSc.assemble!(A)

            B = PETSc.MatSeqAIJ(petsclib, 5, 5, 2)
            B[1, 1] = 1
            B[2, 1] = 2
            B[1, 2] = -2
            PETSc.assemble!(B)

            @test issymmetric(A)
            @test ishermitian(A)
            @test !issymmetric(B)
            @test !ishermitian(B)
        else
            A = PETSc.MatSeqAIJ(petsclib, 5, 5, 2)
            A[1, 1] = 1
            A[2, 1] = -2 + im
            A[1, 2] = -2 - im
            PETSc.assemble!(A)

            B = PETSc.MatSeqAIJ(petsclib, 5, 5, 2)
            B[1, 1] = 1
            B[2, 1] = -2 + im
            B[1, 2] = -2 + im
            PETSc.assemble!(B)
            @test !issymmetric(A)
            @test ishermitian(A)
            @test issymmetric(B)
            @test !ishermitian(B)
        end

        Random.seed!(777)
        A = sprand(PetscScalar, 10, 10, 0.2)
        @test A == PETSc.MatSeqAIJ(petsclib, A)

        PETSc.finalize(petsclib)
    end
end

@testset "MatSeqDense" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        Ajl = PetscScalar.([
            1 2 3 4
            5 6 7 8
            9 10 11 12
            13 14 15 16
        ])

        x = PetscScalar.(collect(1:4))

        A = PETSc.MatSeqDense(petsclib, copy(Ajl))

        @test all(A * x .≈ Ajl * x)

        PETSc.finalize(petsclib)
    end
end
