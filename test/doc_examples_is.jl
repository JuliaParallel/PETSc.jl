using Test
using PETSc
using MPI

@testset "Documentation examples for IS" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF
    PetscInt = PETSc.LibPETSc.PetscInt
    
    @testset "Basic Usage" begin
        # Create an index set from an array of indices (0-based)
        # Note: IS expects indices to be PetscInt (usually Int64), not Int32
        indices = PetscInt[0, 2, 4, 6, 8]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, test_comm, length(indices), indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        @test is isa PETSc.LibPETSc.IS
        @test is.ptr != C_NULL
        
        # Create a stride index set: indices = first:step:(first + step*(n-1))
        is_stride = PETSc.LibPETSc.ISCreateStride(petsclib, test_comm, 10, 0, 2)  # 0, 2, 4, ..., 18
        @test is_stride isa PETSc.LibPETSc.IS
        @test is_stride.ptr != C_NULL
        
        # Get the size of an index set (returns value directly)
        n = PETSc.LibPETSc.ISGetSize(petsclib, is)
        @test n == 5
        
        # Get local size (returns value directly)
        local_n = PETSc.LibPETSc.ISGetLocalSize(petsclib, is)
        @test local_n == 5
        
        # Cleanup
        PETSc.LibPETSc.ISDestroy(petsclib, is)
        @test is.ptr == C_NULL
        
        PETSc.LibPETSc.ISDestroy(petsclib, is_stride)
        @test is_stride.ptr == C_NULL
    end
    
    @testset "Creating Index Sets" begin
        # General index set from array
        indices = PetscInt[1, 3, 5, 7, 9]
        is_general = PETSc.LibPETSc.ISCreateGeneral(petsclib, test_comm, 5, indices, 
                                                    PETSc.LibPETSc.PETSC_COPY_VALUES)
        @test is_general.ptr != C_NULL
        
        # Stride index set: first, first+step, first+2*step, ...
        is_stride = PETSc.LibPETSc.ISCreateStride(petsclib, test_comm, 5, 10, 3)
        @test is_stride.ptr != C_NULL
        size = PETSc.LibPETSc.ISGetSize(petsclib, is_stride)
        @test size == 5
        
        # Block index set: block-structured indices
        block_indices = PetscInt[0, 2, 4]  # Blocks starting at these indices
        is_block = PETSc.LibPETSc.ISCreateBlock(petsclib, test_comm, 2, 3, block_indices, 
                                                PETSc.LibPETSc.PETSC_COPY_VALUES)
        @test is_block.ptr != C_NULL
        
        PETSc.LibPETSc.ISDestroy(petsclib, is_general)
        PETSc.LibPETSc.ISDestroy(petsclib, is_stride)
        PETSc.LibPETSc.ISDestroy(petsclib, is_block)
    end
    
    @testset "Set Operations - skipped due to complex API" begin
        # Note: IS set operations (ISSum, ISDifference, ISExpand) have complex
        # usage patterns that require careful setup. Skip for now.
        @test_skip true
    end
    
    @testset "Querying Properties" begin
        # Create sorted index set
        indices = PetscInt[0, 2, 4, 6, 8]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, test_comm, 5, indices, 
                                            PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Check if index set is sorted (returns Bool directly)
        is_sorted = PETSc.LibPETSc.ISSorted(petsclib, is)
        @test is_sorted == PETSc.LibPETSc.PETSC_TRUE
        
        # Check if identity permutation (returns Bool directly)
        is_identity = PETSc.LibPETSc.ISIdentity(petsclib, is)
        @test typeof(is_identity) == Bool
        
        # Check if a permutation (returns Bool directly)
        is_perm = PETSc.LibPETSc.ISPermutation(petsclib, is)
        @test typeof(is_perm) == Bool
        
        # Get min/max values
        min_val, max_val = PETSc.LibPETSc.ISGetMinMax(petsclib, is)
        @test min_val == 0
        @test max_val == 8
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    @testset "Duplicate and Copy" begin
        indices = PetscInt[1, 3, 5, 7]
        is_original = PETSc.LibPETSc.ISCreateGeneral(petsclib, test_comm, 4, indices, 
                                                      PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Duplicate creates a new IS
        is_dup = PETSc.LibPETSc.ISDuplicate(petsclib, is_original)
        @test is_dup.ptr != C_NULL
        @test is_dup.ptr != is_original.ptr
        
        # Verify duplicated IS has same size
        size_orig = PETSc.LibPETSc.ISGetSize(petsclib, is_original)
        size_dup = PETSc.LibPETSc.ISGetSize(petsclib, is_dup)
        @test size_orig == size_dup
        
        PETSc.LibPETSc.ISDestroy(petsclib, is_original)
        PETSc.LibPETSc.ISDestroy(petsclib, is_dup)
    end
    
    PETSc.finalize(petsclib)
end
