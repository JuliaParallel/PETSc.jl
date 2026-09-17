using Test
using PETSc, MPI
using LinearAlgebra: norm, mul!, Adjoint, Transpose, issymmetric, ishermitian
using SparseArrays: sprand, spdiagm
using Random

MPI.Initialized() || MPI.Init()
comm = LibPETSc.PETSC_COMM_SELF
# Intel Mac has sporadic issues with complex numbers
isintelmac = Sys.isapple() && Sys.ARCH == :x86_64


@testset "PetscMat seqaij" begin
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
            D[2, [2, 3]] = PetscScalar.([3, 4]);
            D[3, [3, 4]] = PetscScalar.([0, 0]);
            D[4, [4, 5]] = PetscScalar.([0, 0]);
            D[5, [4, 5]] = PetscScalar.([3, 4]);
            PETSc.assemble!(D);
            DJ = zeros(PetscScalar, num_rows, num_cols);
            DJ[1, [1, 2]] .= PetscScalar.([1, 2]);
            DJ[2, [2, 3]] .= PetscScalar.([3, 4]);
            DJ[3, [3, 4]] .= PetscScalar.([0, 0]);
            DJ[4, [4, 5]] .= PetscScalar.([0, 0]);
            DJ[5, [4, 5]] .= PetscScalar.([3, 4]);
        else
            D = LibPETSc.MatCreateSeqAIJ(petsclib, comm, num_rows, num_cols, nz_int, C_NULL)
            D[1, [1, 2]] = PetscScalar.([1 + 0im, 0 + 2im]);
            D[2, [2, 3]] = PetscScalar.([3 + 0im, 0 + 4im]);
            D[3, [3, 4]] = PetscScalar.([0 + 0im, 0 + 0im]);
            D[4, [4, 5]] = PetscScalar.([0 + 0im, 0 + 0im]);
            D[5, [4, 5]] = PetscScalar.([3 + 0im, 0 + 4im]);
            PETSc.assemble!(D);
            DJ = zeros(PetscScalar, num_rows, num_cols);
            DJ[1, [1, 2]] .= PetscScalar.([1 + 0im, 0 + 2im]);
            DJ[2, [2, 3]] .= PetscScalar.([3 + 0im, 0 + 4im]);
            DJ[3, [3, 4]] .= PetscScalar.([0 + 0im, 0 + 0im]);
            DJ[4, [4, 5]] .= PetscScalar.([0 + 0im, 0 + 0im]);
            DJ[5, [4, 5]] .= PetscScalar.([3 + 0im, 0 + 4im]);
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

        xc = copy(x)
        yc = copy(y)
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), xc)
        vec_y = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(y)), yc)
        
        # Test mul!
        @test vec_x[:] ≈ x rtol=1e-5
        
        mul!(vec_y, D, vec_x)
        
        y = DJ * x
        @test vec_y[:] ≈ y rtol=1e-5
        
        PETSc.destroy!(D)
        PETSc.destroy!(vec_x)
        PETSc.destroy!(vec_y)
        finalize(xc)
        finalize(yc)

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
            PETSc.destroy!(Asym)
            PETSc.destroy!(Bsym)
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
            PETSc.destroy!(Asym)
            PETSc.destroy!(Bsym)
        end

        Random.seed!(777)
        A1 = sprand(PetscScalar, 10, 11, 0.2)
        B1 = PETSc.PetscMat(petsclib, comm, A1; with_arrays = true)
        sleep(0.1)
        @test sum(B1[:,:] - Matrix(A1)) == 0.0

        PETSc.destroy!(A)
        PETSc.destroy!(B)
        PETSc.destroy!(B1)
        PETSc.destroy!(C)
        PETSc.destroy!(E)
        PETSc.finalize(petsclib)
    end
end


@testset "PetscMat dense" begin
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

        PETSc.destroy!(A)
        PETSc.destroy!(vec_x)
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

        A = PETSc.PetscMat(petsclib, comm, A_sp)
       
        @test sum(A[1:10,1:10] - Matrix(A_sp)) == 0.0 == 0.0
        PETSc.destroy!(A)
        PETSc.finalize(petsclib)
    end
end

@testset "PetscMat seqaij constructor" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Test with integer nonzeros (same for all rows)
        num_rows, num_cols = 5, 7
        nonzeros = 3
        A = PETSc.PetscMat(petsclib, num_rows, num_cols, nonzeros)
        @test size(A) == (num_rows, num_cols)
        
        # Set some values
        A[1, 1] = PetscScalar(1.0)
        A[2, 3] = PetscScalar(2.0)
        A[3, 5] = PetscScalar(3.0)
        PETSc.assemble!(A)
        
        @test A[1, 1] == PetscScalar(1.0)
        @test A[2, 3] == PetscScalar(2.0)
        @test A[3, 5] == PetscScalar(3.0)
        
        PETSc.destroy!(A)
        
        # Test with vector of nonzeros (one per row)
        nz_vec = PetscInt.([2, 3, 1, 4, 2])
        B = PETSc.PetscMat(petsclib, num_rows, num_cols, nz_vec)
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
        
        PETSc.destroy!(B)
        
        PETSc.finalize(petsclib)
    end
end

@testset "PetscMat dense constructor" begin
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
        A = PETSc.PetscMat(petsclib, Ajl)
        
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
        
        xc = copy(x)
        vec_x = LibPETSc.VecCreateSeqWithArray(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(1), PetscInt(length(x)), xc)
        vec_y = LibPETSc.VecCreateSeq(petsclib, LibPETSc.PETSC_COMM_SELF, PetscInt(3))
        
        mul!(vec_y, A, vec_x)
        @test vec_y[:] ≈ y_expected rtol=1e-5
        
        finalize(xc)
        
        # Test matrix modification
        A[1, 1] = PetscScalar(100.0)
        PETSc.assemble!(A)
        @test A[1, 1] == PetscScalar(100.0)
        
        # Explicitly destroy the objects
        PETSc.destroy!(vec_x)
        PETSc.destroy!(vec_y)
        PETSc.destroy!(A)
        
        PETSc.finalize(petsclib)
    end
end

@testset "MatDenseGetArray and RestoreArray" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar

        n, m = 5, 3
        # 1. Create a Dense Matrix
        Ajl = zeros(PetscScalar, n, m)
        A = PETSc.PetscMat(petsclib, Ajl)

        # 2. Test GetArray
        # This should return a Julia Matrix view pointing to PETSc's memory
        jl_view = LibPETSc.MatDenseGetArray(petsclib, A)

        @test size(jl_view) == (n, m)
        @test jl_view isa Matrix{PetscScalar}

        # 3. Test Mutability: Writing to the Julia view should update PETSc
        test_val = PetscScalar(42.0)
        jl_view[2, 2] = test_val

        # 4. Test RestoreArray
        # Passing the array back to PETSc to unlock it
        LibPETSc.MatDenseRestoreArray(petsclib, A, jl_view)

        # 5. Verify the value is actually in PETSc now
        # assembly! is good practice after direct memory modification
        PETSc.assemble!(A)
        @test A[2, 2] == test_val

        # 6. Test copyto!
        A2_jl = PetscScalar.(rand(n, m))
        jl_view2 = LibPETSc.MatDenseGetArray(petsclib, A)
        try
            copyto!(jl_view2, A2_jl)
        finally
            LibPETSc.MatDenseRestoreArray(petsclib, A, jl_view2)
        end
        PETSc.assemble!(A)

        @test A[:, :] ≈ A2_jl

        PETSc.destroy!(A)
        PETSc.finalize(petsclib)
    end
end

@testset "MatDenseGetArray/RestoreArray variants" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        n, m = 4, 3
        Ajl = reshape(PetscScalar.(1:(n*m)), n, m)
        A = PETSc.PetscMat(petsclib, Ajl)

        # MatDenseGetArray
        arr = LibPETSc.MatDenseGetArray(petsclib, A)
        @test size(arr) == (n, m)
        @test arr isa Matrix{PetscScalar}
        arr[1,1] = PetscScalar(123)
        LibPETSc.MatDenseRestoreArray(petsclib, A, arr)
        PETSc.assemble!(A)
        @test A[1,1] == PetscScalar(123)

        # MatDenseGetArrayRead
        arr_read = LibPETSc.MatDenseGetArrayRead(petsclib, A)
        @test size(arr_read) == (n, m)
        @test arr_read[2,2] == A[2,2]
        LibPETSc.MatDenseRestoreArrayRead(petsclib, A, arr_read)

        # MatDenseGetArrayWrite
        arr_write = LibPETSc.MatDenseGetArrayWrite(petsclib, A)
        @test size(arr_write) == (n, m)
        arr_write[3,1] = PetscScalar(456)
        LibPETSc.MatDenseRestoreArrayWrite(petsclib, A, arr_write)
        PETSc.assemble!(A)
        @test A[3,1] == PetscScalar(456)

        # MatDenseGetArrayAndMemType
        arr_and_mem, mtype = LibPETSc.MatDenseGetArrayAndMemType(petsclib, A)
        @test size(arr_and_mem) == (n, m)
        arr_and_mem[4,2] = PetscScalar(789)
        LibPETSc.MatDenseRestoreArrayAndMemType(petsclib, A, arr_and_mem)
        PETSc.assemble!(A)
        @test A[4,2] == PetscScalar(789)

        # MatDenseGetArrayReadAndMemType
        arr_read_mem, mtype2 = LibPETSc.MatDenseGetArrayReadAndMemType(petsclib, A)
        @test size(arr_read_mem) == (n, m)
        @test arr_read_mem[1,2] == A[1,2]
        # No explicit restore needed for read test here

        # MatDenseGetArrayWriteAndMemType
        arr_write_mem, mtype3 = LibPETSc.MatDenseGetArrayWriteAndMemType(petsclib, A)
        @test size(arr_write_mem) == (n, m)
        arr_write_mem[2,3] = PetscScalar(321)
        LibPETSc.MatDenseRestoreArrayWriteAndMemType(petsclib, A, arr_write_mem)
        PETSc.assemble!(A)
        @test A[2,3] == PetscScalar(321)

        PETSc.destroy!(A)
        PETSc.finalize(petsclib)
    end
end

@testset "MatGetGhosts, MatGetRow, MatRestoreRow, MatGetOwnershipRanges" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        # Use a small matrix for ghost/row tests
        n, m = 4, 4
        Ajl = PetscScalar.([1 2 0 0; 0 3 4 0; 0 0 5 6; 7 0 0 8])
        A = PETSc.PetscMat(petsclib, n, m, 2)
        for i in 1:n, j in 1:m
            if Ajl[i,j] != 0
                A[i,j] = Ajl[i,j]
            end
        end
        PETSc.assemble!(A)

        # # MatGetRow and MatRestoreRow
        for row in 1:n
            ncols, cols, vals = LibPETSc.MatGetRow(petsclib, A, PetscInt(row - 1))
            @test ncols == count(!iszero, Ajl[row,:])
            @test length(cols) == ncols
            @test length(vals) == ncols
            # PETSc uses 0-based indices, Julia uses 1-based
            @test all(Ajl[row, cols .+ 1] .== vals)
            # MatRestoreRow hands the row data back to PETSc
            LibPETSc.MatRestoreRow(petsclib, A, PetscInt(row - 1), ncols, cols, vals)
        end

        # MatGetOwnershipRanges and MatGetOwnershipRangesColumn
        # These are meaningful for parallel, but should return [0, n] for sequential
        ranges = LibPETSc.MatGetOwnershipRanges(petsclib, A)
        @test length(ranges) == 2
        @test ranges[1] == 0
        @test ranges[2] == n
        ranges_col = LibPETSc.MatGetOwnershipRangesColumn(petsclib, A)
        @test length(ranges_col) == 2
        @test ranges_col[1] == 0
        @test ranges_col[2] == m

        # MatGetGhosts (should work, but for non-parallel, expect empty or zero ghosts)
        nghosts, ghosts = LibPETSc.MatGetGhosts(petsclib, A)
        @test Int(nghosts) == length(ghosts)
        # For sequential, usually zero ghosts
        @test Int(nghosts) == 0 || all(isa.(ghosts, Int))

        PETSc.destroy!(A)
        PETSc.finalize(petsclib)
    end
end

@testset "MatSeqAIJGetArray and RestoreArray" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # 1. Create a simple Sparse AIJ matrix (Diagonal)
        n = 5
        A = PETSc.PetscMat(petsclib, n, n, 1)
        for i in 1:n
            A[i, i] = PetscScalar(i)
        end
        PETSc.assemble!(A)

        # 2. Test GetArray (Fixed implementation)
        # For AIJ, this returns a 1D view of the non-zero values
        val_view = LibPETSc.MatSeqAIJGetArray(petsclib, A)

        # In a diagonal 5x5 matrix, there are 5 non-zeros.
        # Note: Depending on PETSc internals, the returned length might be
        # the full local buffer size.
        @test length(val_view) >= n
        @test val_view[1] == PetscScalar(1.0)

        # 3. Test Mutability: Update values through the pointer
        new_val = PetscScalar(99.0)
        val_view[1] = new_val

        # 4. Test RestoreArray
        LibPETSc.MatSeqAIJRestoreArray(petsclib, A, val_view)
        PETSc.assemble!(A)

        # 5. Verify PETSc sees the change
        @test A[1, 1] == new_val

        PETSc.destroy!(A)
        PETSc.finalize(petsclib)
    end
end
@testset "CSR access: MatGetRowIJ/MatRestoreRowIJ, MatSeqAIJGetCSRAndMemType, MatMPIAIJGetSeqAIJ, local-to-global mapping" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        comm = MPI.COMM_SELF

        # tridiagonal 5x5 seqaij
        n = 5
        A = PETSc.PetscMat(petsclib, n, n, 3)
        for i in 1:n
            A[i, i] = PetscScalar(2i)
            i > 1 && (A[i, i - 1] = PetscScalar(-1))
            i < n && (A[i, i + 1] = PetscScalar(-1))
        end
        PETSc.assemble!(A)
        nnz_expected = 3n - 2

        # MatGetRowIJ (0-based) / MatRestoreRowIJ
        nr, ia, ja, done = LibPETSc.MatGetRowIJ(petsclib, A, PetscInt(0), LibPETSc.PETSC_FALSE, LibPETSc.PETSC_FALSE)
        @test done == LibPETSc.PETSC_TRUE
        @test nr == n
        @test length(ia) == n + 1
        @test ia[1] == 0 && ia[end] == nnz_expected
        @test length(ja) == nnz_expected
        @test ja[1:2] == PetscInt[0, 1]              # row 1: columns 0, 1
        @test ja[ia[3]+1:ia[4]] == PetscInt[1, 2, 3] # row 3: columns 1, 2, 3
        done2 = LibPETSc.MatRestoreRowIJ(petsclib, A, PetscInt(0), LibPETSc.PETSC_FALSE, LibPETSc.PETSC_FALSE, ia, ja)
        @test done2 == LibPETSc.PETSC_TRUE

        # 1-based variant
        _, ia1, ja1, done = LibPETSc.MatGetRowIJ(petsclib, A, PetscInt(1), LibPETSc.PETSC_FALSE, LibPETSc.PETSC_FALSE)
        @test done == LibPETSc.PETSC_TRUE
        @test ia1[1] == 1 && ia1[end] == nnz_expected + 1 && ja1[1] == 1
        LibPETSc.MatRestoreRowIJ(petsclib, A, PetscInt(1), LibPETSc.PETSC_FALSE, LibPETSc.PETSC_FALSE, ia1, ja1)

        # MatSeqAIJGetCSRAndMemType: same structure plus the values
        i, j, a, mtype = LibPETSc.MatSeqAIJGetCSRAndMemType(petsclib, A)
        @test length(i) == n + 1 && i[end] == nnz_expected
        @test length(j) == nnz_expected && length(a) == nnz_expected
        @test a[1] == PetscScalar(2) && a[2] == PetscScalar(-1)
        @test mtype == LibPETSc.PETSC_MEMTYPE_HOST

        # local-to-global mapping: set a permutation on A and read it back / apply it
        perm = PetscInt[4, 3, 2, 1, 0]
        ltog = LibPETSc.ISLocalToGlobalMappingCreate(petsclib, comm, PetscInt(1), PetscInt(n), perm, LibPETSc.PETSC_COPY_VALUES)
        LibPETSc.MatSetLocalToGlobalMapping(petsclib, A, ltog, ltog)
        rmap, cmap = LibPETSc.MatGetLocalToGlobalMapping(petsclib, A)
        @test rmap != C_NULL && cmap != C_NULL
        out = LibPETSc.ISLocalToGlobalMappingApply(petsclib, rmap, PetscInt(3), PetscInt[0, 2, 4])
        @test out == PetscInt[4, 2, 0]
        outb = LibPETSc.ISLocalToGlobalMappingApplyBlock(petsclib, rmap, PetscInt(2), PetscInt[1, 3])
        @test outb == PetscInt[3, 1]
        LibPETSc.ISLocalToGlobalMappingDestroy(petsclib, ltog)
        PETSc.destroy!(A)

        # MatMPIAIJGetSeqAIJ on a one-process mpiaij matrix: the diagonal block is the whole matrix
        B = LibPETSc.MatCreate(petsclib, comm)
        LibPETSc.MatSetSizes(petsclib, B, PetscInt(n), PetscInt(n), PetscInt(n), PetscInt(n))
        LibPETSc.MatSetType(petsclib, B, "mpiaij")
        LibPETSc.MatMPIAIJSetPreallocation(petsclib, B, PetscInt(3), C_NULL, PetscInt(0), C_NULL)
        for i in 1:n
            LibPETSc.MatSetValues(petsclib, B, PetscInt(1), PetscInt[i - 1], PetscInt(1), PetscInt[i - 1], PetscScalar[3i], LibPETSc.INSERT_VALUES)
        end
        LibPETSc.MatAssemblyBegin(petsclib, B, LibPETSc.MAT_FINAL_ASSEMBLY)
        LibPETSc.MatAssemblyEnd(petsclib, B, LibPETSc.MAT_FINAL_ASSEMBLY)
        @test LibPETSc.MatGetType(petsclib, B) == "mpiaij"
        let
            Ad, Ao, colmap = LibPETSc.MatMPIAIJGetSeqAIJ(petsclib, B)
            @test LibPETSc.MatGetType(petsclib, Ad) == "seqaij"
            md, nd = LibPETSc.MatGetLocalSize(petsclib, Ad)
            @test (md, nd) == (n, n)
            _, no = LibPETSc.MatGetLocalSize(petsclib, Ao)
            @test no == 0 && length(colmap) == 0
            ncols, cols, vals = LibPETSc.MatGetRow(petsclib, Ad, PetscInt(2))
            @test ncols == 1 && vals[1] == PetscScalar(9)
            LibPETSc.MatRestoreRow(petsclib, Ad, PetscInt(2), ncols, cols, vals)
        end
        LibPETSc.MatDestroy(petsclib, B)
        PETSc.finalize(petsclib)
    end
end
