# The rename register for the v0.5 high-level API (docs/src/man/naming.md §1.1).
#
# This file is plain data and is edited by hand. src/public_names.jl, the `public`
# declaration (§13.1), is generated from it so that the two cannot drift apart.
# scripts/audit.jl reads its creator and destroyer name sets.
#
# Run `julia --project=. scripts/generate_renames.jl` after editing.
#
# Every high-level binding must appear in exactly one of
#
#   RENAMES            renamed from v0.4 to v0.5 (the old names are gone in v0.6)
#   INTERNAL           renamed freely or kept, not part of the API
#   UNCHANGED_PUBLIC   kept, part of the API
#
# Macros appear as `Symbol("@old") => Symbol("@new")`.

# ---------------------------------------------------------------------------
# Renamed bindings, v0.4 => v0.5. The old names have no shims from v0.6; the
# new names are public. The rename table in docs/src/man/naming.md follows it.
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
    :plexdistribute! => :distribute,
    :distribute! => :distribute,   # returns a new DM and leaves `dm` as it is
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
    :vtk_save! => :save_vtk,
    :vtk_save_fields! => :save_vtk,
    :save_vtk! => :save_vtk,       # v0.5.1: it writes a file and mutates no argument (§7)
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
    :acquire_petsc_local_array => :acquire_local_array,
    :release_petsc_local_array => :release_local_array,
    :withlocalarray! => :with_local_array!,
    :ghostupdate! => :ghost_update!,
    :ghostupdatebegin! => :ghost_update_begin!,
    :ghostupdateend! => :ghost_update_end!,
    :ownershiprange => :ownership_range,
    :setvalues! => :set_values!,
    :addindex! => :add_index!,
    :determine_memtype => :memtype,

    # init.jl, options.jl, sys.jl
    :Options => :PetscOptions,
    :set_petsclib => :PetscLibType,
    :initialized => :isinitialized,
    :finalized => :isfinalized,
    :check_petsc_wrappers_version => :check_wrappers_version,
    :getcomm => :comm,
    :typedget => :parse_option,
]

# New names that are methods of a function owned by Base: not declared `public` here.
const BASE_TARGETS = Set{Symbol}([:ndims])

# ---------------------------------------------------------------------------
# Internal helpers: renamed freely, not part of the API (§1.1).
# ---------------------------------------------------------------------------

const INTERNAL = Set{Symbol}([
    # the checkout behind with_local_arrays!
    :checkout_local_arrays,
    # the per-library body of the MatShell constructor
    :create_matshell,
    # the finalizer queue behind objects the garbage collector destroys
    :destroy_queue,
    :destroy_queue_lock,
    :destroy_queue_length,
    :queue_for_destroy!,
    :destroy_at_gc!,
    :own_serial!,
    :drain_destroy_queue!,
    # exception-safe setup of a new solver, and the one-at-a-time Vec checkout
    :own_solver!,
    :acquire_each,
    # the expansions behind the pointwise macros and their type lookup
    :simple_fn_expr,
    :residual_fn_expr,
    :jacobian_fn_expr,
    :bd_fn_expr,
    :pointwise_types,
    # leading underscore dropped (naming.md "Internals")
    :petsc_link,
    :petsc_subst,
    :build_petsc_options,
    :ensure_library_handle,
    :ensure_mpi_initialized,
    :library_ptr,
    :post_initialize,
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
    # the underscore stays: inner worker beside a wrapper of the same name
    :_mul!,
    :_zero_rows_vecs,
    :check_axis,
    :stencil_vector,
    :value_vector,
    :field_slots,
    :call_with_views,
    :field_views,
    :TSStepHookFn,
    # callback trampolines, internal from v0.6
    :KSPComputeRHSFn,
    :KSPComputeOperatorsFn,
    :SNESSetFunctionFn,
    :SNESSetJacobianFn,
    :SNESSetConvergenceTestFn,
    :TSSetRHSFunctionFn,
    :TSSetRHSJacobianFn,
    :TSSetIFunctionFn,
    :TSSetIJacobianFn,
    :TSMonitorSetFn,
    :apply_blas_num_threads,
    :finalize_at_exit,
    :register_finalize_at_exit,
    :BLAS_THREAD_VARIABLES,
    :ranks_on_node,
    :_local_arrays,
    :_restore_local_arrays!,
    # kept as they are
    # state kept with the PETSc object, and running callbacks (naming.md §18.3, §18.4)
    :ObjectState,
    :state_type,
    :object_states,
    :object_states_lock,
    :object_state_destroy,
    :object_state,
    :object_state!,
    :state_pointer,
    :drop_object_states!,
    :library_type,
    :wrapper_fields,
    :forward_getproperty,
    :forward_setproperty!,
    :forward_propertynames,
    :CALLBACK_ERROR,
    :capture_callback_errors,
    :rethrow_callback_error,
    :doc_callback,
    :SNESState,
    :KSPState,
    :TSState,
    :PCState,
    :MatShellState,
    :WrappedArrays,
    :keep_alive!,
    :shell_state,
    :DM_ACTIVE_PARTS,
    :PCShellSetApplyFn,      # a trampoline type (§18.6)
    :PCShellSetSetUpFn,
    :PCSetPostSetUpFn,
    :KSPSetPreSolveFn,
    :_taoterm_resettable,
    :_reset_stale_register_flags,
    :_PETSC_ERR_LIB,
    # defined into `PETSc` at runtime, the first time `check_wrappers_version`
    # includes src/autowrapped/petsc_wrappers_version.jl
    :PETSC_WRAPPERS_VERSION,
    # other internals
    :check_initialized,
    :isdestroyable,
    :destroy_handle!,
    :finish_dm!,
    :destroy_options!,
    :as_petsc_vec,
    :csr_from_csc,
    :check_mul_sizes,
    :check_renamed_dm_keywords,
    :destroy_each!,
    :RENAMED_DM_KEYWORDS,
    :mat_seqaij_with_arrays,
    :own_dm!,
    :axis_names,
    :axis_ranges,
    :type_name_symbol,
    :stencil_type_enum,
    :petsclib_of,
    :make_local_array,
    :synchronize_device,
    :device_memtypes,
    :host_access,
    :HostCopy,
    :acquire_host_array,
    :acquire_native_array,
    :array_getter,
    :with_broadcast_arrays,
    :aliases_dest,
    :first_petscvec,
    :check_scalar_result,
    :default_rtol,
    :entrywise_isapprox,
    :comm_of,
    :reduce_owned,
    :compare_owned,
    :check_nonempty,
    :compare_rows,
    :rows_match,
    :with_row,
    :sparse_rows_match,
    :reduction_type,
    :gather_partials,
    :reduce_empty_result,
    :with_init,
    :reduce_zeros,
    :PetscVecStyle,
    :TransposedPetscMat,
    :seq_vec_with_array,
    :with_convergence_check,
    :convergence_error,
    :convergence_hint,
    :reason_text,
    :caller_location,
    :warn_divergence,
    :SOURCE_DIR,
    :DestroyList,
    :options_stacks,
    :options_stack,
    :layered_keys,
    :set_by_enclosing_block,
    :with_layered_options,
    :option_entries,
    :apply_options!,
    :apply_global_options!,
    :OPTIONS_AT,
    :check_options_at,
    :options_due,
    :with_solver_options,
    :store_options!,
    :KSP_REASON_TEXT,
    :SNES_REASON_TEXT,
    :TS_REASON_TEXT,
    :block_values,
    :c_vector,
    # 1-based Julia indices to PETSc's 0-based ones, and the axes of the result
    :indices_0b,
    :index_axes,
    :index_shaped,
    # whether a command-line argument is an option or a value
    :isoption,
    :check_block,
    :checkout_arrays,
    :checked_array,
    :flags,
    :to_petscint_tuple,
    :finalize_if_serial!,
    :check_same_length,
    :waxpy,
    :scaled_copy,
    :add_scaled,
    :vec_from_ptr,
    :mat_from_ptr,
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
    :vertex_location,        # new in v0.5.2
    :edge_location,          # new in v0.5.2
    :face_location,          # new in v0.5.2
    :element_location,       # new in v0.5.2
    :on_lower_side,          # new in v0.5.2
    :stencil,                # new in v0.5.2
    :local_to_local!,        # new in v0.5.2
    :set_matrix_preallocate_only!,  # new in v0.5.2
    :with_product_coordinates,      # new in v0.5.2
    :with_field_views!,             # new in v0.5.2
    :PetscDS,
    :AbstractPetscDS,
    :add_boundary!,
    :add_natural_boundary!,
    :fe_create_default,
    :fe_create_lagrange,
    :set_constants!,
    :set_exact_solution!,
    :set_residual!,
    :set_jacobian!,
    :set_jacobian_preconditioner!,
    # handle types, defined in LibPETSc (naming.md §5.2, §5.5); new in v0.5.4
    :AbstractPetscVec,
    :AbstractPetscMat,
    :AbstractPetscDM,
    :AbstractKSP,
    :AbstractPC,
    :AbstractSNES,
    :AbstractTS,
    :AbstractPetscOptions,
    :AbstractIS,
    :AbstractAO,
    :AbstractPF,
    :AbstractTao,
    :PetscDM,
    :PC,
    :IS,
    :AO,
    :PF,
    :Tao,
    # vec.jl / mat.jl
    :MatShell,
    :MatOp,
    :assemble!,
    :array_type,
    :memtype_backend,
    :owns,
    :AbstractPetscMemBackend,
    :isassembled,            # new in v0.5.2
    :zero_rows!,             # new in v0.5.2
    :zero_rows_local!,       # new in v0.5.2
    :set_option!,            # new in v0.5.2
    :diagonal!,              # new in v0.5.2
    # ksp.jl / snes.jl / ts.jl
    :KSP,
    :SNES,
    :TS,
    :solve!,
    :step!,
    :reset!,
    :interpolate!,
    :interpolate,
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
    :set_from_options!,
    :set_rhs_function!,
    :set_rhs_jacobian!,
    :set_ifunction!,
    :set_ijacobian!,
    :set_monitor!,
    :tolerances,
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
    :set_shell_apply!,       # new in v0.5.1
    :set_shell_setup!,       # new in v0.5.1
    :set_post_setup!,        # new in v0.5.3
    :set_pre_solve!,         # new in v0.5.3
    :iteration_number,       # new in v0.5.1
    :function_norm,          # new in v0.5.1
    :set_function_domain_error!,  # new in v0.5.1
    :prev_time,              # new in v0.5.1
    :set_operators!,         # new in v0.5.1
    :set_dm_active!,         # new in v0.5.1
    :set_options_prefix!,    # new in v0.5.1
    :options_prefix,         # new in v0.5.1
    :equation_type,          # new in v0.5.2
    :set_equation_type!,     # new in v0.5.2
    :set_step_number!,       # new in v0.5.2
    :set_max_snes_failures!, # new in v0.5.2
    :set_max_step_rejections!,  # new in v0.5.2
    :set_error_if_step_fails!,  # new in v0.5.2
    :set_pre_step!,          # new in v0.5.2
    :set_post_step!,         # new in v0.5.2
    # init.jl / options.jl
    :PetscNotInitialized,
    :ConvergenceError,       # new in v0.5.3
    :with_options,           # new in v0.5.3
    :with_local_arrays!,     # new in v0.5.5
    :with_local_coordinate_array,  # new in v0.5.5
    :with_destroy,           # new in v0.5.3
    :track!,                 # new in v0.5.3
    :untrack!,               # new in v0.5.3
    :equal,                  # new in v0.6.0
    :initialize,
    :finalize,
    :scalartype,
    :inttype,
    :realtype,               # new in v0.5.4
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
# generated `public` declaration. The twelve v0.4 exports are gone.
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
# Names scripts/audit.jl matches against (naming.md "Two consequences worth knowing").
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
)

# Releases. `finalizer` hands the release to the garbage collector, which still
# accounts for the object.
const AUDIT_DESTROYERS = Symbol[:destroy!, :finalizer, :destroy_at_gc!]
