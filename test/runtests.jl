using Test
using MPI: MPI, mpiexec
using PETSc, Pkg

# Make sure that all dependencies are installed also on a clean system
Pkg.instantiate()

# When set_library! has been used, petsc_library is a path string and PETSc_jll
# is not loaded.  Only import JLL-specific symbols when using the default binaries.
const USING_CUSTOM_LIB = PETSc.petsclibs[1].petsc_library isa AbstractString

import MPIPreferences  # always import to ensure local MPI API is configured

if USING_CUSTOM_LIB
    @info "Testing PETSc.jl with custom library" path=PETSc.petsclibs[1].petsc_library MPIPreferences.binary MPIPreferences.abi
else
    using PETSc_jll
    @info "Testing PETSc.jl with" MPIPreferences.binary MPIPreferences.abi PETSc_jll.host_platform
end

include("testutils.jl")     # helpers the test files share; each includes it when run alone

# The test groups, one folder each, in the order they run. `Pkg.test(test_args = ["dm"])`
# runs only the named groups; without arguments every group runs.
const GROUPS = ("core", "vecmat", "solvers", "dm", "lowlevel", "regression", "examples", "mpi")
const SELECTED = isempty(ARGS) ? GROUPS : Tuple(ARGS)
for group in SELECTED
    group in GROUPS || error("unknown test group \"$group\"; the groups are $(join(GROUPS, ", "))")
end
selected(group) = group in SELECTED

# Runs `file` in a Julia process of its own, on `nranks` MPI ranks if given
function run_separately(file; nranks = nothing)
    julia = `$(Base.julia_cmd()) --project=. $(joinpath(@__DIR__, file))`
    cmd = isnothing(nranks) ? julia : `$(mpiexec()) -n $nranks $julia`
    return success(pipeline(cmd, stderr = stderr))
end

if selected("core")
    include("core/init.jl")             # must run first: it checks a library never initialized
    include("core/blas_threads.jl")     # -blas_num_threads sizes the BLAS pool PETSc runs in
    include("core/lib.jl")
    include("core/options.jl")
    include("core/bang_returns.jl")     # ! functions return the object they mutate
    include("core/audit.jl")            # leak auditor
    include("core/errors.jl")           # argument validation
    include("core/smoke.jl")            # Vec, Mat and KSP end to end
    include("core/petscbool.jl")        # PetscBool is one byte (PETSc >= 3.24)
    include("core/destroy.jl")          # destroy! guards: stale cycle, double destroy
    include("core/handles.jl")          # destroy! on IS, AO, PF, Tao; ISColoringGetIS ownership
    include("core/handle_arrays.jl")    # C arrays of handles go back to their release function
    include("core/lifetimes.jl")        # arrays a Vec or Mat wraps live as long as it does
    include("core/deprecations.jl")     # the v0.4 names still work and warn (the only file calling them)
    include("core/api_surface.jl")      # scripts/api_surface.jl --check: the register covers the surface
    @testset "method ambiguities do not grow" begin
        # 130 with the regenerated wrappers (all in the high-level layer); regenerate or rename
        # without adding new ones
        @test length(detect_ambiguities(PETSc; recursive = true)) <= 130
    end
end

if selected("vecmat")
    include("vecmat/vec.jl")
    include("vecmat/mat.jl")
    include("vecmat/mat_vec_methods.jl")  # copyto!, fill!, zero_rows!, set_option!, diagonal!, isassembled
    include("vecmat/matshell.jl")
end

if selected("solvers")
    include("solvers/ksp.jl")
    include("solvers/pc.jl")              # high-level PC
    include("solvers/snes.jl")
    include("solvers/solver_api.jl")      # KSP/SNES/TS/PC readers, operators, prefixes, set_from_options!
    include("solvers/ts.jl")              # high-level TS interface
    include("solvers/ts_long_runs.jl")    # equation type, step number, pre/post-step hooks, failure limits
    include("solvers/ts_scenarios.jl")    # ex16 equivalence and a damping sweep
end

if selected("dm")
    include("dm/dmda.jl")
    include("dm/dmstag.jl")
    include("dm/dmstag_serial.jl")
    include("dm/dmstag_stencils.jl")      # DMStag locations by axis, stencils, stencil assembly
    include("dm/dmstag_halos.jl")         # local_to_local!, with_product_coordinates, preallocate only
    include("dm/dmstag_views.jl")         # with_field_views!: types, allocations, hand-back
    include("dm/dmplex.jl")
    include("dm/dmnetwork.jl")
    include("dm/dmshell.jl")
    include("dm/dmproduct.jl")
    include("dm/dm_dimension.jl")         # dimension-correct DM returns: shape, inference, allocations
end

if selected("lowlevel")
    include("lowlevel/readers.jl")             # LibPETSc readers return plain values
    include("lowlevel/wrapper_signatures.jl")  # wrapper arguments take the abstract types
    include("lowlevel/wrapper_quality.jl")     # every generated method infers a concrete return type; no allocations
    include("lowlevel/wrapper_leaks.jl")       # repeated create/destroy and Get/Restore pairs do not grow memory
    include("lowlevel/viewer.jl")
    include("lowlevel/ts.jl")
    include("lowlevel/is.jl")
    include("lowlevel/petscsection.jl")
    include("lowlevel/petscsf.jl")             # graph and communication functions
    include("lowlevel/tao.jl")
end

if selected("regression")
    include("regression/ts_ex51.jl")           # repeated ex51 solves
    include("regression/ts_ex51_implicit.jl")  # repeated implicit Gauss solves
    include("regression/ts_ex16.jl")           # the van der Pol IMEX example
end

if selected("examples")
    include("examples/manual_is.jl")      # the IS examples in the manual
    # The lmvm Tao types work only in the first PETSc cycle of a process (see
    # _reset_stale_register_flags in src/init.jl), so the Tao manual examples run in a fresh one
    @testset "Tao manual examples, fresh process" begin
        @test run_separately("examples/manual_tao.jl")
    end
    include("examples/manual_ts.jl")      # the TS examples in the manual
    include("examples/manual_viewer.jl")  # the PetscViewer examples in the manual

    include("examples/examples.jl")       # every script in examples/
    include("examples/mpi_examples.jl")   # the examples marked `# INCLUDE IN MPI TEST`
end

# Last, so that no MPI run happens inside another
if selected("mpi")
    @testset "MPI Tests" begin
        for file in ("mpi/vec.jl", "mpi/mat.jl", "solvers/ksp.jl", "dm/dmstag.jl",
                     "mpi/dmstag_halos.jl", "mpi/blas_threads.jl")
            @test run_separately(file; nranks = 4)
        end
    end
end
