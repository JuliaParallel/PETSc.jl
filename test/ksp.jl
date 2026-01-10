using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: mul!
using SparseArrays: spdiagm


@testset "KSP" begin
    comm = MPI.COMM_SELF
    mpisize = MPI.Comm_size(comm)
    mpirank = MPI.Comm_rank(comm)

    for (ilib, petsclib) in enumerate(PETSc.petsclibs)
        #petsclib = PETSc.petsclibs[8]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        
        loc_num_rows = PetscInt(10)
        loc_num_cols = PetscInt(10)
        diag_nonzeros = PetscInt(3)
        off_diag_non_zeros = PetscInt(3)

        A = LibPETSc.MatCreateAIJ(petsclib, 
                                comm, 
                                loc_num_rows, 
                                loc_num_cols, 
                                PetscInt(LibPETSc.PETSC_DETERMINE),
                                PetscInt(LibPETSc.PETSC_DETERMINE),
                                diag_nonzeros, 
                                C_NULL,
                                off_diag_non_zeros, 
                                C_NULL)

        # Get compatible vectors
        x,b  = LibPETSc.MatCreateVecs(petsclib,A)

        # set coefficients in matrix
        for i=1:size(A)[1]
            if i==1
                A[i, [i,i+1]] = [-2, 1]
            elseif i == mpisize * loc_num_rows 
                A[i, [i-1,i]] = [1, -2]
            else
                A[i, [i-1,i,i+1]] = [1, -2, 1]
            end
            x[i] = PetscScalar(i)^3
        end
        PETSc.assemble!(A)
        PETSc.assemble!(x)

     
        mul!(b, A, x)
        y = similar(x)

        # Create KSP solver manually ---
        ksp = LibPETSc.KSPCreate(petsclib,comm)
        opts = PETSc.Options(petsclib; 
                                ksp_rtol = 1e-16, pc_type = "jacobi", 
                                ksp_monitor = false, ksp_view = false);

        LibPETSc.KSPSetOperators(petsclib, ksp,A,A)
        
        # Push options to PETSc options database
        push!(opts)
        LibPETSc.KSPSetFromOptions(petsclib, ksp)
        pop!(opts)
        # ------------------------------
        
        # Solve system of equations
        LibPETSc.KSPSolve(petsclib,ksp, b, y)

        # Verify solution
        @test y[1:10] ≈ x[1:10]
        PETSc.destroy(x)

        # do the same with backslash
        x1 = ksp \ b
        @test y ≈ x1
        #PETSc.destroy(x)

        # Create a KSP solver in a simpler way --
        ksp1 = PETSc.KSP(A; ksp_rtol = 1e-16, pc_type = "jacobi", ksp_monitor=false)
        x2 = ksp1 \ b
        @test y ≈ x2
        #PETSc.destroy(ksp1)
        # ------------------------------

        # test some of the get functions:
        b1 = LibPETSc.KSPGetRhs(petsclib, ksp) 
        @test b1 ≈ b

        # this segfaults:
        x3 = LibPETSc.KSPGetSolution(petsclib, ksp) 
        @test x3 ≈ x2
        PETSc.destroy(x2)
        PETSc.destroy(x3)
        

        A1, P1 = LibPETSc.KSPGetOperators(petsclib, ksp) 
        @test A1[1:3,1:3] ≈ A[1:3,1:3]
        @test P1[1:end,1:end] ≈ A[1:end,1:end]

        
        if petsclib== PETSc.petsclibs[1]
            it = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
            @show it
            @test (it == 34) || (it == 33) || (it == 36)   # depending on PETSc version

            it1 = LibPETSc.KSPGetTotalIterations(petsclib, ksp)
            @show it1
            @test (it1 == 68) || (it1 == 66) || (it1 == 72)   # depending on PETSc version

        end

        nrm = LibPETSc.KSPGetResidualNorm(petsclib, ksp)
        @test nrm < 1e-10

        type = LibPETSc.KSPGetType(petsclib,ksp)
        @test type == "gmres" 
        
        if petsclib== PETSc.petsclibs[1]
            rtol1, abstol1, dtol1, maxits1 = LibPETSc.KSPGetTolerances(petsclib, ksp)
            @test rtol1 == 1e-16
            @test abstol1 ==  1.0e-50
            @test dtol1 ==  10000.0f0
            @test maxits1 == 10000

            rtol,abstol,dtol,maxits = 1e-7, 1e-10, 1e-13, 1000;
            LibPETSc.KSPSetTolerances(petsclib, ksp, rtol, abstol, dtol, maxits)

            rtol1, abstol1, dtol1, maxits1 = LibPETSc.KSPGetTolerances(petsclib, ksp)
            @test rtol1 == rtol
            @test abstol1 == abstol
            @test dtol1 == dtol
            @test maxits1 == maxits
        end

        PETSc.destroy(y)
        PETSc.destroy(b)
        PETSc.destroy(b1)
        PETSc.destroy(A)
        PETSc.destroy(A1)
        PETSc.destroy(ksp)
        PETSc.finalize(petsclib)
        

    end
end


if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @testset "KSP with SparseMatrixCSC" begin
        comm = MPI.COMM_SELF
        n = 10

        for petsclib in PETSc.petsclibs 
            #petsclib = PETSc.petsclibs[5]
            PETSc.initialize(petsclib)
            PetscScalar = petsclib.PetscScalar
            PetscInt = petsclib.PetscInt

            A = spdiagm(
                -1 => ones(PetscScalar, n - 1),
                0 => -2ones(PetscScalar, n),
                1 => ones(PetscScalar, n - 1),
            )

            ksp = PETSc.KSP(petsclib, comm, A)
            b = rand(PetscScalar, 10)

            petsc_b = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(length(b)), PetscScalar.(b))
            petsc_x = ksp \ petsc_b
            # test when we supply a petsc vector:
            @test petsc_x[:] ≈ Matrix(A) \ b

            # test with julia b vecror
            @test ksp \ b ≈ Matrix(A) \ b

            PETSc.destroy(ksp)
            PETSc.destroy(petsc_x)
            PETSc.destroy(petsc_b)
            PETSc.finalize(petsclib)
        end
    end
else
    # Even though only rank 0 is running the test all ranks need to initialize
    # PETSc
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PETSc.finalize(petsclib)
    end
end
