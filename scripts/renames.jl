# The rename register for the v0.5 high-level API (docs/src/man/naming.md §1.1).
#
# This file is plain data and is edited by hand. Five things are generated from it, so
# that they cannot drift apart:
#
#   src/deprecations.jl      the forwarding shims (§17)
#   src/public_names.jl      the `public` declaration (§13.1)
#   src/audit_names.jl       the creator/destroyer name sets used by src/audit.jl
#   test/test_deprecations.jl  one test per shim (§17.3)
#   (the rename table at the end of docs/src/man/naming.md)
#
# Run `julia --project=. scripts/generate_renames.jl` after editing.
#
# Every high-level binding must appear in exactly one of
#
#   RENAMES            renamed, with a shim
#   INTERNAL           renamed freely or kept, no shim, not part of the API
#   UNCHANGED_PUBLIC   kept, part of the API
#
# Macros appear as `Symbol("@old") => Symbol("@new")`.

# ---------------------------------------------------------------------------
# Renamed bindings: old => new. A shim is generated for every entry.
# ---------------------------------------------------------------------------

const RENAMES = Pair{Symbol, Symbol}[
    # dm.jl
    :destroy => :destroy!,
    :getinfo => :info,
    :getcorners => :corners,
    :getcorners_dmda => :corners,
    :getghostcorners => :ghost_corners,
    :getghostcorners_dmda => :ghost_corners,
    :dm_local_to_global! => :local_to_global!,
    :dm_global_to_local! => :global_to_local!,
    :setuniformcoordinates_dmda! => :set_uniform_coordinates!,
    :coordinatesDMLocalVec => :local_coordinates,
    :getlocalcoordinatearray => :local_coordinate_array,
    :DMGlobalVec => :global_vec,
    :DMLocalVec => :local_vec,
    :getdimension => :ndims,
    :setfromoptions! => :set_from_options!,

    # dmstag.jl
    :getcorners_dmstag => :corners,
    :getghostcorners_dmstag => :ghost_corners,
    :local_indices_dmstag => :local_indices,
    :global_indices_dmstag => :global_indices,
    :setuniformcoordinates_stag! => :set_uniform_coordinates!,
    :DMStagDOF_Slot => :dof_slot,

    # dmda.jl
    :reshapelocalarray => :reshape_local_array,
    :localinteriorlinearindex => :local_interior_linear_index,
    :dmda_star_fd_coloring => :star_fd_coloring,

    # dmplex.jl
    :isplexsimplex => :issimplex,
    :plexdistribute! => :distribute!,
    :petsc_setname! => :set_name!,
    :getds => :ds,
    :createds! => :create_ds!,
    :getlabel => :label,
    :dm_project_function! => :project_function!,
    :dm_project_field! => :project_field!,
    :dm_compute_l2diff => :l2diff,
    :dm_create_global_vec => :global_vec,
    :dm_create_local_vec => :local_vec,
    :dm_set_auxiliary_vec! => :set_auxiliary_vec!,
    :dm_coarsen_hook_add! => :add_coarsen_hook!,
    :dm_copy_disc! => :copy_disc!,
    :dm_get_coarse => :coarse_dm,
    :fe_copy_quadrature! => :copy_quadrature!,
    :mat_null_space_create => :mat_nullspace_create,
    :mat_set_null_space! => :set_nullspace!,
    :mat_null_space_destroy! => :destroy!,
    :vtk_save! => :save_vtk!,
    :vtk_save_fields! => :save_vtk!,
    :setfield! => :set_field!,
    :dmclone => :clone,
    :plex_set_snes_local_fem! => :set_snes_local_fem!,
    :snes_set_jacobian_null_space! => :set_jacobian_nullspace!,
    :fe_compose_constant_null_space! => :compose_constant_nullspace!,
    Symbol("@petsc_residual_fn") => Symbol("@residual_fn"),
    Symbol("@petsc_jacobian_fn") => Symbol("@jacobian_fn"),
    Symbol("@petsc_bd_fn") => Symbol("@bd_fn"),
    Symbol("@petsc_simple_fn") => Symbol("@simple_fn"),

    # ksp.jl, snes.jl
    :getDM => :dm,
    :setDM! => :set_dm!,
    :get_solution => :solution,
    :gettype => :type_name,
    :type => :type_name,
    :setcomputeoperators! => :set_compute_operators!,
    :setcomputerhs! => :set_compute_rhs!,
    :setfunction! => :set_function!,
    :setjacobian! => :set_snes_jacobian!,
    :setconvergencetest! => :set_convergence_test!,

    # vec.jl, mat.jl: construction goes through the type (§5.1, §6)
    :VecSeq => :PetscVec,
    :MatSeqAIJ => :PetscMat,
    :MatSeqDense => :PetscMat,
    :MatCreateSeqAIJ => :PetscMat,
    :MatSeqAIJWithArrays => :PetscMat,
    :MatAIJ => :PetscMat,
    :unsafe_localarray => :unsafe_local_array,
    :wrap_localarray => :wrap_local_array,
    :acquire_petsc_local_array => :acquire_local_array,
    :release_petsc_local_array => :release_local_array,
    :get_petsc_arrays => :local_arrays,
    :restore_petsc_arrays => :restore_local_arrays!,
    :withlocalarray! => :with_local_array!,
    :ghostupdate! => :ghost_update!,
    :ghostupdatebegin! => :ghost_update_begin!,
    :ghostupdateend! => :ghost_update_end!,
    :ownershiprange => :ownership_range,
    :setvalues! => :set_values!,
    :addindex! => :add_index!,
    :determine_memtype => :memtype,

    # init.jl, options.jl, sys.jl, audit.jl
    :Options => :PetscOptions,
    :set_petsclib => :PetscLibType,
    :initialized => :isinitialized,
    :finalized => :isfinalized,
    :check_petsc_wrappers_version => :check_wrappers_version,
    :getcomm => :comm,
    :typedget => :parse_option,
    :audit_petsc_file => :audit_file,

    # types (§5.5): a `const` alias, no warning
    :Fn_KSPComputeRHS => :KSPComputeRHSFn,
    :Fn_KSPComputeOperators => :KSPComputeOperatorsFn,
    :Fn_SNESSetFunction => :SNESSetFunctionFn,
    :Fn_SNESSetJacobian => :SNESSetJacobianFn,
    :Fn_SNESSetConvergenceTest => :SNESSetConvergenceTestFn,
    :AbstractPETScMemBackend => :AbstractPetscMemBackend,
]

# Old names that name a type: the shim is a `const` alias rather than a forwarding
# method, so no warning is emitted and `x isa Fn_SNESSetFunction` keeps working.
const TYPE_RENAMES = Set{Symbol}([
    :Fn_KSPComputeRHS,
    :Fn_KSPComputeOperators,
    :Fn_SNESSetFunction,
    :Fn_SNESSetJacobian,
    :Fn_SNESSetConvergenceTest,
    :AbstractPETScMemBackend,
])

# Type aliases that are already written by hand in the source rather than generated.
const TYPE_RENAMES_INSOURCE = Set{Symbol}([:AbstractPETScMemBackend])

# New names that are methods of a function owned by Base: not declared `public` here.
const BASE_TARGETS = Set{Symbol}([:ndims])

# Shims the generated `old(args...; kwargs...) = new(args...; kwargs...)` gets wrong,
# because the argument list changed. Written out in full instead: the reorders and
# the calls that dropped `petsclib` (§8, §17.2).
const CUSTOM_SHIMS = Dict{Symbol, String}(
    # `MatSeqAIJWithArrays(petsclib, comm, A::SparseMatrixCSC)` and
    # `MatCreateSeqAIJ(petsclib, comm, S)` had the same argument list and
    # different meanings, so only one of them can keep it. §6 gives the CSR
    # arrays to `PetscMat(petsclib, rowptr, colval, nzval)`; the shim converts.
    :MatSeqAIJWithArrays => """
    function MatSeqAIJWithArrays(args...; kwargs...)
        @warn "MatSeqAIJWithArrays is deprecated, use PetscMat" maxlog = 1
        return mat_seqaij_with_arrays(args...; kwargs...)
    end
    """,

    # ── Step 3: the reorders and the dropped `petsclib` (§8, §17.2) ──────────
    #
    # Each shim warns first and then forwards in the new order, so the warning
    # is emitted even when the forwarded call cannot succeed. The old argument
    # lists are spelled out rather than slurped so that a call written against
    # the v0.4 order lands on the right parameters.

    :dm_project_function! => """
    function dm_project_function!(args...; kwargs...)
        @warn "dm_project_function! is deprecated, use project_function!" maxlog = 1
        petsclib, dm, time, funcs, ctxs, mode, X = args
        return project_function!(X, dm, time, funcs, ctxs, mode; kwargs...)
    end
    """,

    :dm_project_field! => """
    function dm_project_field!(args...; kwargs...)
        @warn "dm_project_field! is deprecated, use project_field!" maxlog = 1
        petsclib, dm, time, U, funcs, mode, X = args
        return project_field!(X, dm, time, U, funcs, mode; kwargs...)
    end
    """,

    :dm_compute_l2diff => """
    function dm_compute_l2diff(args...; kwargs...)
        @warn "dm_compute_l2diff is deprecated, use l2diff" maxlog = 1
        petsclib, dm, time, funcs, ctxs, X = args
        return l2diff(dm, time, funcs, ctxs, X; kwargs...)
    end
    """,

    # v0.4 spelled this one `dm_global_to_local!(gvec, lvec, dm, mode)`.
    :dm_global_to_local! => """
    function dm_global_to_local!(args...; kwargs...)
        @warn "dm_global_to_local! is deprecated, use global_to_local!" maxlog = 1
        gvec, lvec, dm, rest... = args
        return global_to_local!(lvec, dm, gvec, rest...; kwargs...)
    end
    """,

    # v0.4 spelled this one `dm_local_to_global!(lvec, gvec, dm, mode)`.
    :dm_local_to_global! => """
    function dm_local_to_global!(args...; kwargs...)
        @warn "dm_local_to_global! is deprecated, use local_to_global!" maxlog = 1
        lvec, gvec, dm, rest... = args
        return local_to_global!(gvec, dm, lvec, rest...; kwargs...)
    end
    """,

    :plex_set_snes_local_fem! => """
    function plex_set_snes_local_fem!(args...; kwargs...)
        @warn "plex_set_snes_local_fem! is deprecated, use set_snes_local_fem!" maxlog = 1
        petsclib, dm = args
        return set_snes_local_fem!(dm; kwargs...)
    end
    """,

    :dmda_star_fd_coloring => """
    function dmda_star_fd_coloring(args...; kwargs...)
        @warn "dmda_star_fd_coloring is deprecated, use star_fd_coloring" maxlog = 1
        petsclib, da = args
        return star_fd_coloring(da; kwargs...)
    end
    """,

    # `vtk_save!(petsclib, comm, filename, vec)` loses both `petsclib` and
    # `comm`: the vector carries the library and PETSc answers for the
    # communicator (§8).
    :vtk_save! => """
    function vtk_save!(args...; kwargs...)
        @warn "vtk_save! is deprecated, use save_vtk!" maxlog = 1
        petsclib, comm, filename, vec = args
        return save_vtk!(vec, filename; kwargs...)
    end
    """,

    :vtk_save_fields! => """
    function vtk_save_fields!(args...; kwargs...)
        @warn "vtk_save_fields! is deprecated, use save_vtk!" maxlog = 1
        petsclib, comm, filename, vecs = args
        return save_vtk!(vecs, filename; kwargs...)
    end
    """,

    # `dm_coarsen_hook_add!(dm, hook, restrict)` puts the callback first (§8.1).
    :dm_coarsen_hook_add! => """
    function dm_coarsen_hook_add!(args...; kwargs...)
        @warn "dm_coarsen_hook_add! is deprecated, use add_coarsen_hook!" maxlog = 1
        dm, hook, rest... = args
        return add_coarsen_hook!(hook, dm, rest...; kwargs...)
    end
    """,
)

# ---------------------------------------------------------------------------
# Deprecated methods of names that did not change (§17.2, "argument dropped").
#
# These keep an argument the v0.5 signature no longer takes and warn when it is
# passed, so they cannot be expressed as `old => new`. They are emitted verbatim
# into src/deprecations.jl and removed with it in v0.6. `test` is a test body
# appended to test/test_deprecations.jl.
# ---------------------------------------------------------------------------

const EXTRA_SHIMS = NamedTuple{(:name, :code, :test), Tuple{String, String, String}}[
    (
        name = "ownership_range(A, base_one)",
        code = """
        # §12.1: `ownership_range` is 1-based only. The positional form is kept for
        # one release because that is the form v0.4 has; a keyword shim would compile
        # and never fire.
        function ownership_range(
            obj::Union{LibPETSc.AbstractPetscVec, LibPETSc.AbstractPetscMat},
            base_one::Bool,
        )
            @warn "ownership_range(A, base_one) is deprecated, use ownership_range(A), " *
                  "which is 1-based" maxlog = 1
            r = ownership_range(obj)
            return base_one ? r : ((first(r) - 1):(last(r) - 1))
        end
        """,
        test = """
        @testset "ownership_range(A, base_one)" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            v = PETSc.PetscVec(petsclib, 5)
            one_based = PETSc.ownership_range(v)
            zero_based = warns() do
                PETSc.ownership_range(v, false)
            end
            @test zero_based == ((first(one_based) - 1):(last(one_based) - 1))
            PETSc.destroy!(v)
        end
        """,
    ),
    (
        name = "set_type!(obj, ::AbstractString)",
        code = """
        # §3.1: type names are `Symbol` at the Julia API. The `String` spelling is
        # accepted for one release and warns (§17.2, "argument dropped").
        function set_type!(obj, type::AbstractString)
            @warn "set_type!(obj, \\"\$type\\") is deprecated, use set_type!(obj, :\$type)" maxlog = 1
            return set_type!(obj, Symbol(type))
        end
        """,
        test = """
        @testset "set_type!(obj, ::AbstractString)" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            # The low-level creator, because the `KSP` constructor wants the
            # operators and this test only needs an object with a type.
            ksp = PETSc.LibPETSc.KSPCreate(petsclib, PETSc.MPI.COMM_SELF)
            warns() do
                PETSc.set_type!(ksp, "cg")
            end
            @test PETSc.type_name(ksp) === :cg
            PETSc.destroy!(ksp)
        end
        """,
    ),
    (
        name = "add_boundary!(petsclib, dm, ...)",
        code = """
        # §8: `petsclib` never leads a high-level call. `add_boundary!` and
        # `add_natural_boundary!` keep their names, so the dropped argument is a
        # deprecated method rather than a rename.
        function add_boundary!(petsclib::LibPETSc.PetscLibType, dm, args...; kwargs...)
            @warn "add_boundary!(petsclib, dm, ...) is deprecated, use add_boundary!(dm, ...)" maxlog = 1
            return add_boundary!(dm, args...; kwargs...)
        end

        function add_natural_boundary!(petsclib::LibPETSc.PetscLibType, dm, args...; kwargs...)
            @warn "add_natural_boundary!(petsclib, dm, ...) is deprecated, use " *
                  "add_natural_boundary!(dm, ...)" maxlog = 1
            return add_natural_boundary!(dm, args...; kwargs...)
        end
        """,
        test = """
        @testset "add_boundary! and add_natural_boundary! without petsclib" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            for f in (PETSc.add_boundary!, PETSc.add_natural_boundary!)
                warns() do
                    try
                        f(petsclib, DeprecationProbe())
                    catch
                    end
                end
            end
        end
        """,
    ),
]

# ---------------------------------------------------------------------------
# Internal helpers: renamed freely, no shim, not part of the API (§1.1).
# ---------------------------------------------------------------------------

const INTERNAL = Set{Symbol}([
    # leading underscore dropped (naming.md "Internals")
    :petsc_link,
    :petsc_subst,
    :vtk_merge_one_tensor!,
    :build_petsc_options,
    :ensure_library_handle,
    :ensure_mpi_initialized,
    :library_ptr,
    :post_initialize,
    :release_library_handle,
    :doc_external,
    # the alias src/autowrapped/ interpolates, and the docstring helpers
    :_doc_external,
    :doc_borrowed,
    Symbol("@renamed"),
    :library_path_string,
    :lib_handles,
    :petsc_program_name,
    :errorcode,
    :run_callback,
    :with_options,
    # the underscore stays: inner worker beside a wrapper of the same name
    :_mul!,
    :_unsafe_local_array,
    :_local_arrays,
    :_restore_local_arrays!,
    # kept as they are
    :_taoterm_resettable,
    :_reset_stale_register_flags,
    :_MATSEQAIJ_WITHARRAYS_STORAGE,
    :_PETSC_ERR_LIB,
    # defined into `PETSc` at runtime, the first time `check_wrappers_version`
    # includes src/autowrapped/petsc_wrappers_version.jl
    :PETSC_WRAPPERS_VERSION,
    # other internals
    :check_initialized,
    :isdestroyable,
    :as_petsc_vec,
    :csr_from_csc,
    :mat_seqaij_with_arrays,
    :own_dm!,
    :axis_names,
    :axis_ranges,
    :type_name_symbol,
    :stencil_type_enum,
    :petsclib_of,
    :make_local_array,
    :to_petscint_tuple,
    :audit_walk,
    :audit_targets,
    :audit_report,
    :audit_creator,
    :audit_destroyer,
    :audit_callee,
    :audit_argnames,
    :audit_isbroadcast,
    :audit_hasparseerror,
    :AUDIT_TYPE_CREATORS,
    :AUDIT_NAMED_CREATORS,
    :AUDIT_DESTROYER_NAMES,
])

# ---------------------------------------------------------------------------
# Kept names that are part of the API.
# ---------------------------------------------------------------------------

const UNCHANGED_PUBLIC = Symbol[
    # dm.jl / dmda.jl / dmstag.jl / dmplex.jl
    :setup!,
    :DMDA,
    :DMStag,
    :DMPlex,
    :narrow,
    :ndofs,
    :PetscDS,
    :AbstractPetscDS,
    :add_boundary!,
    :add_natural_boundary!,
    :create_split_boundary_labels!,
    :fe_create_default,
    :fe_create_lagrange,
    :set_constants!,
    :set_exact_solution!,
    :set_residual!,
    :set_jacobian!,
    :set_jacobian_preconditioner!,
    :vtk_merge_tensor!,
    # vec.jl / mat.jl
    :VecPtr,
    :MatPtr,
    :MatShell,
    :MatOp,
    :assemble!,
    :array_type,
    :memtype_backend,
    :owns,
    :AbstractPetscMemBackend,
    # ksp.jl / snes.jl / ts.jl
    :KSP,
    :SNES,
    :TS,
    :solve!,
    :step!,
    :reset!,
    :interpolate!,
    :set_type!,
    :set_adapt_type!,
    :set_problem_type!,
    :set_exact_final_time!,
    :set_dm!,
    :set_solution!,
    :set_time!,
    :set_timestep!,
    :set_max_time!,
    :set_max_steps!,
    :set_tolerances!,
    :set_user_ctx!,
    :set_from_options!,
    :set_rhs_function!,
    :set_rhs_jacobian!,
    :set_ifunction!,
    :set_ijacobian!,
    :set_monitor!,
    :tolerances,
    :user_ctx,
    :current_time,
    :timestep,
    :max_time,
    :max_steps,
    :step_number,
    :step_rejections,
    :converged_reason,
    :solve_time,
    :snes_iterations,
    :ksp_iterations,
    :snes_failures,
    :snes,
    :ksp,
    :pc,                     # new in v0.5.1
    :set_fieldsplit_is!,     # new in v0.5.1
    :TSSetRHSFunctionFn,
    :TSSetRHSJacobianFn,
    :TSSetIFunctionFn,
    :TSSetIJacobianFn,
    :TSMonitorSetFn,
    # init.jl / options.jl
    :PetscNotInitialized,
    :initialize,
    :finalize,
    :scalartype,
    :inttype,
    :library_info,
    :petsclibs,
    :set_library!,
    :unset_library!,
    :set_petscint!,          # new in v0.5 (#241): one PetscInt width per process
    :tao_usable_after_reinitialize,
    :parse_options,
]

# Names `PETSc` exports (§13): the nine types
# and construction entry points, `petsclibs`, and the `LibPETSc` submodule.
# A name cannot be both exported and `public`, so these are subtracted from the
# generated `public` declaration. The twelve v0.4 exports are gone; the ones that
# were renamed stay reachable through their shims, qualified.
const EXPORTED = Symbol[
    :LibPETSc,
    :DMDA,
    :DMStag,
    :DMPlex,
    :PetscVec,
    :PetscMat,
    :PetscOptions,
    :KSP,
    :SNES,
    :TS,
    :petsclibs,
]

# ---------------------------------------------------------------------------
# Names src/audit.jl matches against (naming.md "Two consequences worth knowing").
# ---------------------------------------------------------------------------

# Type constructors: the name is also the kind of object.
const AUDIT_TYPE_CREATORS = Dict{Symbol, String}(
    :KSP => "KSP",
    :SNES => "SNES",
    :TS => "TS",
    :DMDA => "DM",
    :DMStag => "DM",
    :DMPlex => "DM",
    :PetscVec => "Vec",
    :PetscMat => "Mat",
    :PetscOptions => "Options",
)

# Creators whose name carries no `Create`/`Duplicate` marker.
const AUDIT_NAMED_CREATORS = Dict{Symbol, String}(
    :global_vec => "Vec",
    :local_vec => "Vec",
    :DMGetCoordinateDM => "DM",
    :DMStagCreateCompatibleDMStag => "DM",
    :DMCreateMatrix => "Mat",
    :MatCreateVecs => "Vec",
    :MatShell => "Mat",
    :clone => "DM",
    # v0.4 spellings, still reachable through the shims
    :VecSeq => "Vec",
    :MatAIJ => "Mat",
    :MatSeqAIJ => "Mat",
    :MatSeqDense => "Mat",
    :DMGlobalVec => "Vec",
    :DMLocalVec => "Vec",
    :dm_create_global_vec => "Vec",
    :dm_create_local_vec => "Vec",
    :dmclone => "DM",
)

# Releases. `finalizer` hands the release to the garbage collector, which still
# accounts for the object.
const AUDIT_DESTROYERS = Symbol[:destroy!, :destroy, :finalizer]
