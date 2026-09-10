using Test
using PETSc
using MPI

# Initialize PETSc
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF

@testset "PetscSection Low-Level API" begin
    
    @testset "Basic Section Creation and Query" begin
        # Create a section (wrapper has wrong signature, use ccall)
        section = Ref{LibPETSc.PetscSection}()
        err = ccall(
            (:PetscSectionCreate, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (MPI.MPI_Comm, Ptr{LibPETSc.PetscSection}),
            test_comm, section
        )
        @test err == 0
        @test section[] != C_NULL
        
        # Set chart: range of valid point indices [pStart, pEnd)
        LibPETSc.PetscSectionSetChart(petsclib, section[], 0, 10)
        
        # Verify chart
        pStart, pEnd = LibPETSc.PetscSectionGetChart(petsclib, section[])
        @test pStart == 0
        @test pEnd == 10
        
        # Set DOF count for each point
        for p in 0:9
            num_dofs = (p < 4) ? 1 : 2  # Different DOFs per point
            LibPETSc.PetscSectionSetDof(petsclib, section[], p, num_dofs)
        end
        
        # Setup: compute offsets
        LibPETSc.PetscSectionSetUp(petsclib, section[])
        
        # Query DOF counts
        for p in 0:9
            dof = LibPETSc.PetscSectionGetDof(petsclib, section[], p)
            expected_dof = (p < 4) ? 1 : 2
            @test dof == expected_dof
        end
        
        # Query offsets
        offset5 = LibPETSc.PetscSectionGetOffset(petsclib, section[], 5)
        @test offset5 == 4 * 1 + 1 * 2  # Points 0-3 have 1 DOF each (4), point 4 has 2 DOFs (2)
        
        # Get total storage size
        storage_size = LibPETSc.PetscSectionGetStorageSize(petsclib, section[])
        expected_size = 4 * 1 + 6 * 2  # 4 points with 1 DOF + 6 points with 2 DOFs
        @test storage_size == expected_size
        
        # Cleanup (wrapper has wrong signature, use ccall)
        err = ccall(
            (:PetscSectionDestroy, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (Ptr{LibPETSc.PetscSection},),
            section
        )
        @test err == 0
    end
    
    @testset "Multi-Field Section" begin
        # Create section with 2 fields
        section = Ref{LibPETSc.PetscSection}()
        err = ccall(
            (:PetscSectionCreate, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (MPI.MPI_Comm, Ptr{LibPETSc.PetscSection}),
            test_comm, section
        )
        @test err == 0
        
        LibPETSc.PetscSectionSetNumFields(petsclib, section[], 2)
        
        # Verify number of fields
        num_fields = LibPETSc.PetscSectionGetNumFields(petsclib, section[])
        @test num_fields == 2
        
        # Set field names
        LibPETSc.PetscSectionSetFieldName(petsclib, section[], 0, "velocity")
        LibPETSc.PetscSectionSetFieldName(petsclib, section[], 1, "pressure")
        
        # Note: PetscSectionGetFieldName wrapper has signature issues, skip testing retrieval
        
        # Set chart
        LibPETSc.PetscSectionSetChart(petsclib, section[], 0, 5)
        
        # Set field components
        LibPETSc.PetscSectionSetFieldComponents(petsclib, section[], 0, 3)
        LibPETSc.PetscSectionSetFieldComponents(petsclib, section[], 1, 1)
        
        # Verify field components
        vel_comp = LibPETSc.PetscSectionGetFieldComponents(petsclib, section[], 0)
        pres_comp = LibPETSc.PetscSectionGetFieldComponents(petsclib, section[], 1)
        @test vel_comp == 3
        @test pres_comp == 1
        
        # Set DOFs per field per point
        for p in 0:4
            LibPETSc.PetscSectionSetFieldDof(petsclib, section[], p, 0, 3)  # 3 velocity DOFs
            LibPETSc.PetscSectionSetFieldDof(petsclib, section[], p, 1, 1)  # 1 pressure DOF
            LibPETSc.PetscSectionSetDof(petsclib, section[], p, 4)          # Total: 4 DOFs
        end
        
        LibPETSc.PetscSectionSetUp(petsclib, section[])
        
        # Verify field DOFs
        for p in 0:4
            vel_dof = LibPETSc.PetscSectionGetFieldDof(petsclib, section[], p, 0)
            pres_dof = LibPETSc.PetscSectionGetFieldDof(petsclib, section[], p, 1)
            total_dof = LibPETSc.PetscSectionGetDof(petsclib, section[], p)
            @test vel_dof == 3
            @test pres_dof == 1
            @test total_dof == 4
        end
        
        # Get total storage
        storage_size = LibPETSc.PetscSectionGetStorageSize(petsclib, section[])
        @test storage_size == 5 * 4  # 5 points with 4 DOFs each
        
        # Cleanup
        err = ccall(
            (:PetscSectionDestroy, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (Ptr{LibPETSc.PetscSection},),
            section
        )
        @test err == 0
    end
    
    @testset "Section with Uniform DOFs" begin
        # Create section
        section = Ref{LibPETSc.PetscSection}()
        err = ccall(
            (:PetscSectionCreate, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (MPI.MPI_Comm, Ptr{LibPETSc.PetscSection}),
            test_comm, section
        )
        @test err == 0
        
        # Set chart
        npoints = 100
        LibPETSc.PetscSectionSetChart(petsclib, section[], 0, npoints)
        
        # Set uniform DOF count
        dofs_per_point = 5
        for p in 0:npoints-1
            LibPETSc.PetscSectionSetDof(petsclib, section[], p, dofs_per_point)
        end
        
        LibPETSc.PetscSectionSetUp(petsclib, section[])
        
        # Verify storage size
        storage_size = LibPETSc.PetscSectionGetStorageSize(petsclib, section[])
        @test storage_size == npoints * dofs_per_point
        
        # Verify max DOF
        max_dof = LibPETSc.PetscSectionGetMaxDof(petsclib, section[])
        @test max_dof == dofs_per_point
        
        # Verify offsets are sequential
        for p in 0:npoints-1
            expected_offset = p * dofs_per_point
            actual_offset = LibPETSc.PetscSectionGetOffset(petsclib, section[], p)
            @test actual_offset == expected_offset
        end
        
        # Cleanup
        err = ccall(
            (:PetscSectionDestroy, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (Ptr{LibPETSc.PetscSection},),
            section
        )
        @test err == 0
    end
    
    @testset "Section Copy" begin
        # Create original section
        section = Ref{LibPETSc.PetscSection}()
        err = ccall(
            (:PetscSectionCreate, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (MPI.MPI_Comm, Ptr{LibPETSc.PetscSection}),
            test_comm, section
        )
        @test err == 0
        
        LibPETSc.PetscSectionSetChart(petsclib, section[], 0, 5)
        for p in 0:4
            LibPETSc.PetscSectionSetDof(petsclib, section[], p, p + 1)
        end
        LibPETSc.PetscSectionSetUp(petsclib, section[])
        
        # Create new section for copy
        section_copy = Ref{LibPETSc.PetscSection}()
        err = ccall(
            (:PetscSectionCreate, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (MPI.MPI_Comm, Ptr{LibPETSc.PetscSection}),
            test_comm, section_copy
        )
        @test err == 0
        
        # Copy section
        LibPETSc.PetscSectionCopy(petsclib, section[], section_copy[])
        
        # Verify copy has same properties
        pStart_orig, pEnd_orig = LibPETSc.PetscSectionGetChart(petsclib, section[])
        pStart_copy, pEnd_copy = LibPETSc.PetscSectionGetChart(petsclib, section_copy[])
        @test pStart_orig == pStart_copy
        @test pEnd_orig == pEnd_copy
        
        storage_orig = LibPETSc.PetscSectionGetStorageSize(petsclib, section[])
        storage_copy = LibPETSc.PetscSectionGetStorageSize(petsclib, section_copy[])
        @test storage_orig == storage_copy
        
        # Cleanup
        err = ccall(
            (:PetscSectionDestroy, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (Ptr{LibPETSc.PetscSection},),
            section
        )
        @test err == 0
        
        err = ccall(
            (:PetscSectionDestroy, petsclib.petsc_library),
            PETSc.LibPETSc.PetscErrorCode,
            (Ptr{LibPETSc.PetscSection},),
            section_copy
        )
        @test err == 0
    end
    
end

PETSc.finalize(petsclib)
