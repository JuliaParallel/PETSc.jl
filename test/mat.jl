using Test
using PETSc, MPI
using LinearAlgebra: norm, mul!, Adjoint, Transpose, issymmetric, ishermitian
using SparseArrays: sprand, spdiagm
using Random

# Windows PETSc binaries are built without MPI support, skip MPI initialization
if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end
# Windows PETSc binaries are built without MPI support, use PETSC_COMM_SELF
comm = LibPETSc.PETSC_COMM_SELF
# Intel Mac has sporadic issues with complex numbers
isintelmac = Sys.isapple() && Sys.ARCH == :x86_64


@testset "MatSeqAIJ" begin
    for petsclib in PETSc.petsclibs
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
        LibPETSc.MatAssemblyBegin(petsclib, A, PETSc.MAT_FINAL_ASSEMBLY)
        LibPETSc.MatAssemblyEnd(petsclib, A, PETSc.MAT_FINAL_ASSEMBLY)
        PETSc.assemble!(B)

        C = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, PetscInt(0), nz_vec)
        @test size(C) == (num_rows, num_cols)
        PETSc.assemble!(C)
        @test A == C

        if PetscScalar <: Real
            D = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)
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

        @test norm(A) == 0
        @test norm(D) ≈ norm(DJ)

        x = PetscScalar.(Array(1:num_cols))
        y = zeros(PetscScalar, num_rows)
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), copy(x))
        vec_y = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(y)), copy(y))

        # Zero vectors before multiplication to avoid stale data
        fill!(vec_x, zero(PetscScalar))
        fill!(vec_y, zero(PetscScalar))

        # Test mul!
        vec_x .= x
        mul!(vec_y, D, vec_x)
        y = DJ * x
        @test vec_y[1:num_rows] ≈ y rtol=1e-5

        if !(isintelmac)
            fill!(vec_x, zero(PetscScalar))
            mul!(vec_x, Transpose(D), vec_y)
            x = Transpose(DJ) * y
            @test vec_x[1:end] ≈ x rtol=1e-5
        end

        PETSc.destroy(vec_x)
        PETSc.destroy(vec_y)

        if PetscScalar <: Real
            Asym = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            Asym[1, 1] = 1;
            Asym[2, 1] = -2;
            Asym[1, 2] = -2;
            PETSc.assemble!(Asym);
            Bsym = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            Bsym[1, 1] = 1;
            Bsym[2, 1] = 2;
            Bsym[1, 2] = -2;
            PETSc.assemble!(Bsym);
            @test issymmetric(Asym)
            @test ishermitian(Asym)
            @test !issymmetric(Bsym)
            @test !ishermitian(Bsym)
            PETSc.destroy(Asym)
            PETSc.destroy(Bsym)
        else
            Asym = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            Asym[1, 1] = PetscScalar(1)
            Asym[2, 1] = PetscScalar(-2 + im)
            Asym[1, 2] = PetscScalar(-2 - im)
            PETSc.assemble!(Asym)
            Bsym = LibPETSc.MatCreateSeqAIJ(petsclib, comm, PetscInt(5), PetscInt(5), PetscInt(2), C_NULL)
            Bsym[1, 1] = PetscScalar(1)
            Bsym[2, 1] = PetscScalar(-2 + im)
            Bsym[1, 2] = PetscScalar(-2 - im)
            PETSc.assemble!(Bsym)
            @test !issymmetric(Asym)
            @test !issymmetric(Bsym)
            PETSc.destroy(Asym)
            PETSc.destroy(Bsym)
        end

        Random.seed!(777)
        A1 = sprand(PetscScalar, 10, 11, 0.2)
        B1 = PETSc.MatSeqAIJWithArrays(petsclib, comm, A1)
        sleep(0.1)
        @test sum(B1[:,:] - Matrix(A1)) == 0.0

        PETSc.destroy(A)
        PETSc.destroy(B)
        PETSc.destroy(B1)
        PETSc.destroy(C)
        PETSc.destroy(D)
        PETSc.destroy(E)
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
        fill!(vec_x, zero(PetscScalar))

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

@testset "MatSeqAIJ constructor" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Test with integer nonzeros (same for all rows)
        num_rows, num_cols = 5, 7
        nonzeros = 3
        A = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nonzeros)
        @test size(A) == (num_rows, num_cols)
        
        # Set some values
        A[1, 1] = PetscScalar(1.0)
        A[2, 3] = PetscScalar(2.0)
        A[3, 5] = PetscScalar(3.0)
        PETSc.assemble!(A)
        
        @test A[1, 1] == PetscScalar(1.0)
        @test A[2, 3] == PetscScalar(2.0)
        @test A[3, 5] == PetscScalar(3.0)
        
        PETSc.destroy(A)
        
        # Test with vector of nonzeros (one per row)
        nz_vec = PetscInt.([2, 3, 1, 4, 2])
        B = PETSc.MatSeqAIJ(petsclib, num_rows, num_cols, nz_vec)
        @test size(B) == (num_rows, num_cols)
        
        # Set values matching the nonzero pattern
        B[1, 1] = PetscScalar(10.0)
        B[1, 2] = PetscScalar(11.0)
        B[2, 1] = PetscScalar(20.0)
        B[2, 2] = PetscScalar(21.0)
        B[2, 3] = PetscScalar(22.0)
        PETSc.assemble!(B)
        
        @test B[1, 1] == PetscScalar(10.0)
        @test B[1, 2] == PetscScalar(11.0)
        @test B[2, 1] == PetscScalar(20.0)
        
        PETSc.destroy(B)
        
        PETSc.finalize(petsclib)
    end
end

@testset "MatSeqDense constructor" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Create a Julia matrix
        if PetscScalar <: Real
            Ajl = PetscScalar.([
                1 2 3 4
                5 6 7 8
                9 10 11 12
            ])
        else
            Ajl = PetscScalar.([
                1 2+im 3 4-im
                5 6+2im 7 8
                9-im 10 11+im 12
            ])
        end

        # Create PETSc dense matrix from Julia matrix
        A = PETSc.MatSeqDense(petsclib, Ajl)
        
        # Test that the matrix was created successfully
        @test A !== nothing
        @test A.ptr != C_NULL
        
        # Test size
        @test size(A) == size(Ajl)
        @test size(A) == (3, 4)
        
        # Test individual value retrieval
        @test A[1, 1] == Ajl[1, 1]
        @test A[2, 3] == Ajl[2, 3]
        @test A[3, 4] == Ajl[3, 4]
        
        # Test retrieving all values
        A_retrieved = A[:, :]
        @test A_retrieved == Ajl
        
        # Test matrix-vector multiplication
        x = PetscScalar.(collect(1:4))
        y_expected = Ajl * x
        
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), copy(x))
        vec_y = LibPETSc.VecCreateSeq(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(3))
        
        mul!(vec_y, A, vec_x)
        @test vec_y[:] ≈ y_expected rtol=1e-5
        
        # Test matrix modification
        A[1, 1] = PetscScalar(100.0)
        PETSc.assemble!(A)
        @test A[1, 1] == PetscScalar(100.0)
        
        # Explicitly destroy objects
        PETSc.destroy(vec_x)
        PETSc.destroy(vec_y)
        PETSc.destroy(A)
        
        PETSc.finalize(petsclib)
    end
end