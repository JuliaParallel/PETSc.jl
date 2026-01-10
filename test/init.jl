using Test
using PETSc

@testset "init" begin
    for petsclib in PETSc.petsclibs
        # The first time through finalize should be false since we have never
        # initialized petsclib yet...
        initial_finalized_value = false

        # since we haven't called anything these should be false!
        @test !(PETSc.initialized(petsclib))
        @test PETSc.finalized(petsclib) == initial_finalized_value

        # The second time through  time through finalize should be true since
        # we have initialized petsclib yet...
        initial_finalized_value = true

        # initialize PETSc
        PETSc.initialize(petsclib)

        # Check values again
        @test PETSc.initialized(petsclib)
        @test !(PETSc.finalized(petsclib))

        PETSc.finalize(petsclib)

        # Check values again
        @test !(PETSc.initialized(petsclib))
        @test PETSc.finalized(petsclib)
    end

    # Test initialize with log_view enabled
    for petsclib in PETSc.petsclibs
        # Test basic initialize with log_view
        @test !(PETSc.initialized(petsclib))
        
        # Store original PETSC_OPTIONS if it exists
        original_opts = get(ENV, "PETSC_OPTIONS", nothing)
        
        # Initialize with log view and a test file
        test_log_file = tempname() * ".txt"
        PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file"])
        
        # Check initialization succeeded
        @test PETSc.initialized(petsclib)
        @test !(PETSc.finalized(petsclib))
        
        # Check that PETSC_OPTIONS was cleaned up
        current_opts = get(ENV, "PETSC_OPTIONS", nothing)
        @test current_opts == original_opts
        
        # Finalize
        PETSc.finalize(petsclib)
        @test !(PETSc.initialized(petsclib))
        @test PETSc.finalized(petsclib)
        
        # Check that log file was created
        @test isfile(test_log_file)
        
        # Verify log file contains expected content
        log_content = read(test_log_file, String)
        @test occursin("PETSc Performance Summary", log_content)
        
        # Clean up test log file
        rm(test_log_file; force = true)
    end
    
    # Test initialize with options but without log_view
    for petsclib in PETSc.petsclibs
        @test !(PETSc.initialized(petsclib))
        
        # Store original PETSC_OPTIONS
        original_opts = get(ENV, "PETSC_OPTIONS", nothing)
        
        # Initialize with custom options only (no log_view)
        PETSc.initialize(petsclib; options = ["-malloc_debug"])
        
        @test PETSc.initialized(petsclib)
        
        # Check that PETSC_OPTIONS was cleaned up
        current_opts = get(ENV, "PETSC_OPTIONS", nothing)
        @test current_opts == original_opts
        
        PETSc.finalize(petsclib)
    end
    
    # Test initialize with multiple options and log_view
    for petsclib in PETSc.petsclibs
        @test !(PETSc.initialized(petsclib))
        
        # Initialize with multiple options
        test_log_file = tempname() * ".txt"
        PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file", "-log_view_memory"])
        
        @test PETSc.initialized(petsclib)
        
        PETSc.finalize(petsclib)
        
        # Check log file was created and contains memory info
        @test isfile(test_log_file)
        log_content = read(test_log_file, String)
        @test occursin("Memory usage", log_content) || occursin("memory", lowercase(log_content))
        
        rm(test_log_file; force = true)
    end
    
    # Test that PETSC_OPTIONS is preserved if it already exists
    for petsclib in PETSc.petsclibs
        @test !(PETSc.initialized(petsclib))
        
        # Set a custom PETSC_OPTIONS with a valid but harmless option
        ENV["PETSC_OPTIONS"] = "-malloc_debug 0"
        
        test_log_file = tempname() * ".txt"
        PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file"])
        
        @test PETSc.initialized(petsclib)
        
        # After initialization, original option should be restored
        @test get(ENV, "PETSC_OPTIONS", "") == "-malloc_debug 0"
        
        PETSc.finalize(petsclib)
        
        # Clean up
        delete!(ENV, "PETSC_OPTIONS")
        @test isfile(test_log_file)
        rm(test_log_file; force = true)
    end
end

