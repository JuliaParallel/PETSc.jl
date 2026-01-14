using Test
using PETSc, MPI
using LinearAlgebra: norm

MPI.Initialized() || MPI.Init()
comm = MPI.COMM_WORLD

@testset "VecBasics" begin
    for petsclib in PETSc.petsclibs
        #petsclib        = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar     = petsclib.PetscScalar
        PetscInt        = petsclib.PetscInt
        N               = PetscInt(10)
        v1 =  LibPETSc.VecCreateSeq(petsclib,comm, N)
        v2 =  LibPETSc.VecCreateMPI(petsclib,comm, PetscInt(LibPETSc.PETSC_DECIDE), N)

        LibPETSc.VecZeroEntries(petsclib,v1)

        indices = PetscInt.([0, 1, 2, 8])  # Note that PETSc uses 0-based indexing
        values  = PetscScalar.([1.0, 2.0, 3.0, 4.0])
        LibPETSc.VecSetValues(petsclib, v1, PetscInt(4), indices, values, LibPETSc.INSERT_VALUES)
        LibPETSc.VecAssemblyBegin(petsclib,v1)
        LibPETSc.VecAssemblyEnd(petsclib,v1)

        x = LibPETSc.VecGetArray(petsclib,v1)
        @test x[9] == 4.0

        x1 = rand(PetscScalar,10)
        bs = PetscInt(1)
        
        # Use GC.@preserve to ensure x1 stays alive while PETSc is using it
        # VecCreateSeqWithArray wraps the Julia array, so the array must not be moved/collected
        GC.@preserve x1 begin
            v3 = LibPETSc.VecCreateSeqWithArray(petsclib,comm, bs, PetscInt(length(x1)), x1)

            # Check get/restore array
            x3 = LibPETSc.VecGetArray(petsclib,v3)
            x3[10] = PetscScalar(42.0)
            LibPETSc.VecRestoreArray(petsclib,v3, x3)
            
            # After restore, x1 contains the modified data (v3 wraps x1)
            expected_sum = sum(x1)

            # get values from PETSc vector (note 0-based indexing)
            indices = PetscInt.([8,9]) # in 0-based indexing! 
            vals = zeros(PetscScalar, length(indices))
            vals =LibPETSc.VecGetValues(petsclib,v3,PetscInt(length(indices)), indices)
            expected_vals = [x1[9], PetscScalar(42.0)]  # Use original array x1, with modified last value
            @test vals == expected_vals

            # create a duplicate vector 
            v4 = LibPETSc.VecDuplicate(petsclib,v3)

            # copy content (note that this function is not correctly parsed automatically)
            LibPETSc.VecCopy(petsclib,v3, v4)
            
            # Verify the sum using the underlying array x1 which should have the data
            @test LibPETSc.VecSum(petsclib,v3) ≈ expected_sum rtol=1e-10
            @test LibPETSc.VecSum(petsclib,v4) ≈ expected_sum rtol=1e-10
            
            PETSc.destroy(v3)
            PETSc.destroy(v4)
        end

        # Julia candy:
        v5      =   LibPETSc.VecCreateSeq(petsclib,comm, N)
        
        v5[1]   =   PetscScalar(3.14)
        v5[2:3] =   PetscScalar.([2.71, 1.61])
        @test v5[1:4] == PetscScalar.([ 3.14, 2.71,1.61,0.0])

        fill!(v5, PetscScalar(1.11))
        @test v5[1] == PetscScalar(1.11)

        PETSc.destroy(v1)
        PETSc.destroy(v2)
        PETSc.destroy(v5)
        PETSc.finalize(petsclib)
    end
end


@testset "VecCreateSeqWithArray" begin
    N = 10
    for petsclib in PETSc.petsclibs[1:2]
        #petsclib = PETSc.petsclibs[5]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt    = petsclib.PetscInt
        x           = rand(PetscScalar, N)
        
        # Use GC.@preserve to ensure x stays alive while PETSc is using it
        GC.@preserve x begin
            petsc_x     = LibPETSc.VecCreateSeqWithArray(petsclib,comm, PetscInt(1), PetscInt(length(x)), x)

            @test LibPETSc.VecGetSize(petsclib, petsc_x) == N
            val = LibPETSc.VecNorm(petsclib,petsc_x, PETSc.NORM_2)
            @test val ≈ norm(x)

            # make sure the viewer works
            #=
            _stdout = stdout
            (rd, wr) = redirect_stdout()
            @show petsc_x
            @test readline(rd) == "petsc_x = Vec Object: 1 MPI process"
            @test readline(rd) == "  type: seq"
            redirect_stdout(_stdout)

            _stdout = stdout
            (rd, wr) = redirect_stdout()
            show(stdout, "text/plain", petsc_x)
            @test readline(rd) == "Vec Object: 1 MPI process"
            @test readline(rd) == "  type: seq"
            redirect_stdout(_stdout)
            =#
            
            @test LibPETSc.VecGetOwnershipRange(petsclib,petsc_x) == (0, N)
            
            x2 = LibPETSc.VecGetArray(petsclib,petsc_x)
            @test x2 == x
            
            PETSc.destroy(petsc_x)
        end
        
        PETSc.finalize(petsclib)
    end
end

@testset "VecSeq" begin
    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt    = petsclib.PetscInt
        N           = PetscInt(10)
        petsc_x = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, N)
        @test LibPETSc.VecGetSize(petsclib, petsc_x) == N

        @test LibPETSc.VecGetOwnershipRange(petsclib,petsc_x) == (0, N)

        x = rand(PetscScalar, N)
        x2 = LibPETSc.VecGetArray(petsclib,petsc_x)
        x2 .= x
        LibPETSc.VecRestoreArray(petsclib,petsc_x, x2)
        @test LibPETSc.VecNorm(petsclib, petsc_x, PETSc.NORM_2) ≈ norm(x)
    
        @test LibPETSc.VecGetType(petsclib, petsc_x) == "seq"

        PETSc.destroy(petsc_x)
        PETSc.finalize(petsclib)
    end
end

@testset "VecSeq constructor with array" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        
        # Test with simple array
        x = ones(PetscScalar, 3)
        # Use GC.@preserve since VecSeq wraps the Julia array
        GC.@preserve x begin
            v = PETSc.VecSeq(petsclib, x)
            
            # Verify the vector was created successfully
            @test v !== nothing
            @test v.ptr != C_NULL
            @test LibPETSc.VecGetSize(petsclib, v) == 3
            
            # Verify values
            @test v[1:3] == ones(PetscScalar, 3)
            @test LibPETSc.VecSum(petsclib, v) == PetscScalar(3.0)
            
            # Test that modifications to the vector affect the underlying array
            # (since VecSeq uses VecCreateSeqWithArray)
            v[1] = PetscScalar(42.0)
            @test x[1] == PetscScalar(42.0)
            
            PETSc.destroy(v)
        end
        
        # Test with different array values
        if PetscScalar <: Real
            x2 = PetscScalar.([1.0, 2.0, 3.0, 4.0, 5.0])
        else
            x2 = PetscScalar.([1.0, 2.0+im, 3.0, 4.0-im, 5.0])
        end
        GC.@preserve x2 begin
            v2 = PETSc.VecSeq(petsclib, x2)
            
            @test LibPETSc.VecGetSize(petsclib, v2) == 5
            @test v2[1:5] == x2
            @test LibPETSc.VecNorm(petsclib, v2, PETSc.NORM_2) ≈ norm(x2)
            
            PETSc.destroy(v2)
        end
        
        # Test with blocksize parameter
        x3 = PetscScalar.([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        GC.@preserve x3 begin
            v3 = PETSc.VecSeq(petsclib, x3; blocksize=2)
            @test LibPETSc.VecGetSize(petsclib, v3) == 6
            @test LibPETSc.VecGetBlockSize(petsclib, v3) == 2
            
            PETSc.destroy(v3)
        end
        
        PETSc.finalize(petsclib)
    end
end


@testset "withlocalarray!" begin
    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt    = petsclib.PetscInt
        N           = PetscInt(10)
        petsc_x     = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, N)
        petsc_y     = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, N)

        # extract one array, write
        PETSc.withlocalarray!(
            petsc_x;
            read = false, write = true,
        ) do x
            @test all(x .== 0)
            for i in eachindex(x)
                x[i] = PetscScalar(i)
            end
        end
        @test petsc_x[1:N] == PetscScalar.(1:N)
        
        # with tuple as input
        PETSc.withlocalarray!(
            (petsc_x, );
            read = (false,), write = (true,),
        ) do x
            for i in eachindex(x)
                x[i] = PetscScalar(i)
            end
        end
        @test petsc_x[1:N] == PetscScalar.(1:N)

        # with 2-variable tuple as input
        PETSc.withlocalarray!(
            (petsc_x, petsc_y);
            read = (false,false), write = (true,true),
        ) do x, y
            for i in eachindex(x)
                x[i] = PetscScalar(i)
                y[i] = PetscScalar(2i)
            end
        end
        @test petsc_x[1:N] == PetscScalar.(1:N)
        @test petsc_y[1:N] == PetscScalar.(2:2:2N)

        # with 2-variable tuple as input
        PETSc.withlocalarray!(
            petsc_x, petsc_y;
            read = (false,false), write = (true,true),
        ) do x, y
            for i in eachindex(x)
                x[i] = PetscScalar(i)
                y[i] = PetscScalar(2i)
            end
        end
        @test petsc_x[1:N] == PetscScalar.(1:N)
        @test petsc_y[1:N] == PetscScalar.(2:2:2N)


        PETSc.destroy(petsc_x)
        PETSc.destroy(petsc_y)

        PETSc.finalize(petsclib)
    end
end