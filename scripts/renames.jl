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
# because the argument list changed. Written out in full instead. Empty for step 1:
# every rename so far keeps its argument list, and the four reorders land in step 3.
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
)

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
    :SNESConvergenceTestBox,
    :_MATSEQAIJ_WITHARRAYS_STORAGE,
    :_PETSC_ERR_LIB,
    # other internals
    :check_initialized,
    :isdestroyable,
    :as_petsc_vec,
    :csr_from_csc,
    :mat_seqaij_with_arrays,
    :own_dm!,
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
    :set_library!,
    :unset_library!,
    :tao_usable_after_reinitialize,
    :parse_options,
]

# Names `PETSc` exports. Step 1 leaves the export list alone (§13 lands in step 4), so
# these are the v0.4 exports minus `HostBackend`, which was never defined.
# A name cannot be both exported and `public`, so these are subtracted from the
# generated `public` declaration.
const EXPORTED = Symbol[
    :LibPETSc,
    :audit_petsc_file,
    :set_petsclib,
    :set_library!,
    :unset_library!,
    :library_info,
    :AbstractPetscMemBackend,
    :AbstractPETScMemBackend,
    :determine_memtype,
    :get_petsc_arrays,
    :restore_petsc_arrays,
    :dmda_star_fd_coloring,
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
