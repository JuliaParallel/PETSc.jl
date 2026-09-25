using Test
using PETSc

@testset "init" begin
    for petsclib in PETSc.petsclibs
        # The first time through finalize should be false since we have never
        # initialized petsclib yet...
        initial_finalized_value = false

        # since we haven't called anything these should be false!
        @test !(PETSc.isinitialized(petsclib))
        @test PETSc.isfinalized(petsclib) == initial_finalized_value

        # The second time through  time through finalize should be true since
        # we have initialized petsclib yet...
        initial_finalized_value = true

        # initialize PETSc
        PETSc.initialize(petsclib)

        # Check values again
        @test PETSc.isinitialized(petsclib)
        @test !(PETSc.isfinalized(petsclib))

        PETSc.finalize(petsclib)

        # Check values again
        @test !(PETSc.isinitialized(petsclib))
        @test PETSc.isfinalized(petsclib)
    end

    # Test initialize with log_view enabled
    for petsclib in PETSc.petsclibs
        # Test basic initialize with log_view
        @test !(PETSc.isinitialized(petsclib))
        
        # Store original PETSC_OPTIONS if it exists
        original_opts = get(ENV, "PETSC_OPTIONS", nothing)
        
        # Initialize with log view and a test file
        # Use current directory to avoid Windows temp directory issues with short paths
        test_log_file = joinpath(pwd(), "test_log_$(hash(petsclib)).txt")
        
        try
            PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file"])
            
            # Check initialization succeeded
            @test PETSc.isinitialized(petsclib)
            @test !(PETSc.isfinalized(petsclib))
            
            # Check that PETSC_OPTIONS was cleaned up
            current_opts = get(ENV, "PETSC_OPTIONS", nothing)
            @test current_opts == original_opts
            
            # Finalize
            PETSc.finalize(petsclib)
            @test !(PETSc.isinitialized(petsclib))
            @test PETSc.isfinalized(petsclib)
            
            # Check that log file was created
            # Note: On some platforms/configurations, PETSc may not create the log file
            if isfile(test_log_file)
                # Verify log file contains expected content
                log_content = read(test_log_file, String)
                @test occursin("PETSc Performance Summary", log_content)
            else
                @warn "Log file not created on this platform: $test_log_file"
            end
        finally
            # Clean up test log file even if tests fail
            rm(test_log_file; force = true)
        end
    end
    
    # Test initialize with options but without log_view
    for petsclib in PETSc.petsclibs
        @test !(PETSc.isinitialized(petsclib))
        
        # Store original PETSC_OPTIONS
        original_opts = get(ENV, "PETSC_OPTIONS", nothing)
        
        # Initialize with custom options only (no log_view)
        PETSc.initialize(petsclib; options = ["-malloc_debug"])
        
        @test PETSc.isinitialized(petsclib)
        
        # Check that PETSC_OPTIONS was cleaned up
        current_opts = get(ENV, "PETSC_OPTIONS", nothing)
        @test current_opts == original_opts
        
        PETSc.finalize(petsclib)
                @testset "wrappers_version" begin
                    # Ensure the wrappers version check runs and returns the expected fields
                    result = PETSc.check_wrappers_version()
                    @test isa(result, NamedTuple)
                    @test :wrappers_version in keys(result)
                    @test :installed_version in keys(result)
                    @test :match in keys(result)
                    @test result.match === nothing || isa(result.match, Bool)
                end
    end
    
    # Test initialize with multiple options and log_view
    for petsclib in PETSc.petsclibs
        @test !(PETSc.isinitialized(petsclib))
        
        # Initialize with multiple options
        # Use current directory to avoid Windows temp directory issues
        test_log_file = joinpath(pwd(), "test_log_multi_$(hash(petsclib)).txt")
        
        try
            PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file", "-log_view_memory"])
            
            @test PETSc.isinitialized(petsclib)
            
            PETSc.finalize(petsclib)
            
            # Check log file was created and contains memory info
            # Note: On some platforms/configurations, PETSc may not create the log file
            if isfile(test_log_file)
                log_content = read(test_log_file, String)
                @test occursin("PETSc Performance Summary", log_content)
                # Memory summary may not always be present depending on PETSc build
            else
                @warn "Log file not created on this platform: $test_log_file"
            end
        finally
            # Clean up even if tests fail
            rm(test_log_file; force = true)
        end
    end
    
    # Test that PETSC_OPTIONS is preserved if it already exists
    for petsclib in PETSc.petsclibs
        @test !(PETSc.isinitialized(petsclib))
        
        # Set a custom PETSC_OPTIONS with a valid but harmless option
        ENV["PETSC_OPTIONS"] = "-malloc_debug 0"
        
        # Use current directory to avoid Windows temp directory issues
        test_log_file = joinpath(pwd(), "test_log_preserved_$(hash(petsclib)).txt")
        
        try
            PETSc.initialize(petsclib; log_view = true, options = [":$test_log_file"])
            
            @test PETSc.isinitialized(petsclib)
            
            # After initialization, original option should be restored
            @test get(ENV, "PETSC_OPTIONS", "") == "-malloc_debug 0"
            
            PETSc.finalize(petsclib)
        finally
            # Clean up
            delete!(ENV, "PETSC_OPTIONS")
            # Note: On some platforms/configurations, PETSc may not create the log file
            rm(test_log_file; force = true)
        end
    end

    # Test the PetscLibType constructor with a precompiled library
    @testset "PetscLibType constructor" begin
        # Get path from one of the precompiled libraries
        precompiled_lib = PETSc.petsclibs[1]
        lib_path = precompiled_lib.petsc_library
        
        # Create a custom library instance using the PetscLibType constructor
        custom_lib = PETSc.LibPETSc.PetscLibType(lib_path; PetscScalar=Float64, PetscInt=Int64)
        
        # Verify type parameters
        @test custom_lib isa PETSc.LibPETSc.PetscLibType{Float64, Int64, String}
        @test custom_lib.petsc_library == lib_path
        @test PETSc.scalartype(custom_lib) == Float64
        @test PETSc.inttype(custom_lib) == Int64
        
        # Test that it can be initialized and works
        @test !PETSc.isinitialized(custom_lib)
        PETSc.initialize(custom_lib)
        @test PETSc.isinitialized(custom_lib)
        
        # Test basic functionality
        version = PETSc.LibPETSc.PetscGetVersionNumber(custom_lib)
        @test version[1] == 3  # PETSc version 3.x
        @test version[2] >= 0  # Minor version
        
        PETSc.finalize(custom_lib)
        @test !PETSc.isinitialized(custom_lib)
        
        # Test with different scalar/int types (only when more than one library is
        # configured; with set_library! there is just one)
        if length(PETSc.petsclibs) >= 2
            lib_path_f32 = PETSc.petsclibs[2].petsc_library  # Float32 library
            custom_lib_f32 = PETSc.LibPETSc.PetscLibType(lib_path_f32; PetscScalar=Float32, PetscInt=Int64)
            @test custom_lib_f32 isa PETSc.LibPETSc.PetscLibType{Float32, Int64, String}
            @test PETSc.scalartype(custom_lib_f32) == Float32
            @test PETSc.inttype(custom_lib_f32) == Int64
        end
    end

    # naming.md §12: `library_info` returns the values it used to print, and the
    # report is its `show` method.
    @testset "library_info" begin
        info = PETSc.library_info()
        @test info isa NamedTuple
        @test keys(info) === (:source, :path, :scalar, :int, :real)
        @test info.source in (:preferences, :jll)
        @test info.scalar === PETSc.petsclibs[1].PetscScalar
        @test info.int === PETSc.petsclibs[1].PetscInt
        @test info.real === real(PETSc.petsclibs[1].PetscScalar)
        @test info.path isa AbstractString

        report = sprint(show, MIME"text/plain"(), info)
        @test occursin("Source  : ", report)
        @test occursin("Loaded libraries (this session):", report)
        if info.source === :preferences
            @test occursin("Path    : $(info.path)", report)
        else
            @test occursin("PETSc_jll", report)
        end
    end

    # `options` reach PETSc on every platform: they are handed to PetscInitialize
    # as its command line, where PETSC_OPTIONS set from Julia is invisible on Windows
    @testset "options reach PETSc" begin
        petsclib = PETSc.petsclibs[1]
        PETSc.isinitialized(petsclib) && PETSc.finalize(petsclib)
        options = ["-petscjl_test_int", "7", "-petscjl_test_pair 8"]
        PETSc.initialize(petsclib; options)
        global_options = PETSc.LibPETSc.PetscOptions{typeof(petsclib)}(C_NULL, petsclib.age; own = false)
        value(name) = PETSc.LibPETSc.PetscOptionsGetInt(petsclib, global_options, "", name)
        @test value("-petscjl_test_int") == (7, PETSc.LibPETSc.PETSC_TRUE)
        @test value("-petscjl_test_pair") == (8, PETSc.LibPETSc.PETSC_TRUE)   # one entry, two words
        @test options == ["-petscjl_test_int", "7", "-petscjl_test_pair 8"]  # left as given
        PETSc.finalize(petsclib)
    end

    # naming.md §14: user input raises `ArgumentError`, not a bare `error`.
    @testset "ArgumentError on user input" begin
        @test_throws ArgumentError PETSc.set_library!(
            joinpath(mktempdir(), "no-such-libpetsc.so"),
        )
    end
end

