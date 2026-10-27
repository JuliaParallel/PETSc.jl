# The rename register for the v0.5 high-level API (docs/src/man/naming.md §1.1).
#
# This file is plain data and is edited by hand. Five things are generated from it, so
# that they cannot drift apart:
#
#   src/deprecations.jl      the forwarding shims (§17)
#   src/public_names.jl      the `public` declaration (§13.1)
#   src/audit_names.jl       the creator/destroyer name sets used by src/audit.jl
#   test/core/deprecations.jl  one test per shim (§17.3)
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
    :vtk_merge_tensor! => :vtk_merge_tensor,   # v0.5.1: likewise
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
        return project_function!(X, dm, time, funcs; ctxs, mode, kwargs...)
    end
    """,

    :dm_project_field! => """
    function dm_project_field!(args...; kwargs...)
        @warn "dm_project_field! is deprecated, use project_field!" maxlog = 1
        petsclib, dm, time, U, funcs, mode, X = args
        return project_field!(X, dm, time, U, funcs; mode, kwargs...)
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
        mode = isempty(rest) ? INSERT_VALUES : only(rest)
        return global_to_local!(lvec, dm, gvec; mode, kwargs...)
    end
    """,

    # v0.4 spelled this one `dm_local_to_global!(lvec, gvec, dm, mode)`.
    :dm_local_to_global! => """
    function dm_local_to_global!(args...; kwargs...)
        @warn "dm_local_to_global! is deprecated, use local_to_global!" maxlog = 1
        lvec, gvec, dm, rest... = args
        mode = isempty(rest) ? INSERT_VALUES : only(rest)
        return local_to_global!(gvec, dm, lvec; mode, kwargs...)
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
        @warn "vtk_save! is deprecated, use save_vtk" maxlog = 1
        petsclib, comm, filename, vec = args
        return save_vtk(vec, filename; kwargs...)
    end
    """,

    :vtk_save_fields! => """
    function vtk_save_fields!(args...; kwargs...)
        @warn "vtk_save_fields! is deprecated, use save_vtk" maxlog = 1
        petsclib, comm, filename, vecs = args
        return save_vtk(vecs, filename; kwargs...)
    end
    """,
)

# ---------------------------------------------------------------------------
# Deprecated methods of names that did not change (§17.2, "argument dropped").
#
# These keep an argument the v0.5 signature no longer takes and warn when it is
# passed, so they cannot be expressed as `old => new`. They are emitted verbatim
# into src/deprecations.jl and removed with it in v0.6. `test` is a test body
# appended to test/core/deprecations.jl.
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
    (
        name = "positional insert mode, scatter mode and diag",
        code = """
        # §9: an argument with a sensible default is a keyword. The insert mode,
        # scatter mode and `diag` of these calls were positional; those forms
        # warn for one release (§17.2, "argument dropped").
        for f in (:ghost_update!, :ghost_update_begin!, :ghost_update_end!)
            @eval function \$f(
                vec::LibPETSc.AbstractPetscVec,
                mode::InsertMode,
                scatter::LibPETSc.ScatterMode = SCATTER_FORWARD,
            )
                @warn string(\$(string(f)), "(v, mode, scatter) is deprecated, use ",
                             \$(string(f)), "(v; mode, scatter)") maxlog = 1
                return \$f(vec; mode, scatter)
            end
        end

        for f in (:local_to_global!, :global_to_local!, :local_to_local!)
            @eval function \$f(
                dst::LibPETSc.AbstractPetscVec,
                dm::LibPETSc.AbstractPetscDM,
                src::LibPETSc.AbstractPetscVec,
                mode::InsertMode,
            )
                @warn string(\$(string(f)), "(dst, dm, src, mode) is deprecated, use ",
                             \$(string(f)), "(dst, dm, src; mode)") maxlog = 1
                return \$f(dst, dm, src; mode)
            end
        end

        function local_to_local!(v::LibPETSc.AbstractPetscVec, dm::LibPETSc.AbstractPetscDM, mode::InsertMode)
            @warn "local_to_local!(v, dm, mode) is deprecated, use local_to_local!(v, dm; mode)" maxlog = 1
            return local_to_local!(v, dm; mode)
        end

        function set_values!(
            M::LibPETSc.AbstractPetscMat,
            rows_0b::AbstractVector,
            cols_0b::AbstractVector,
            rowvals::AbstractArray,
            mode::InsertMode;
            kwargs...,
        )
            @warn "set_values!(M, rows, cols, vals, mode) is deprecated, use " *
                  "set_values!(M, rows, cols, vals; mode)" maxlog = 1
            return set_values!(M, rows_0b, cols_0b, rowvals; mode, kwargs...)
        end

        function set_values!(J::LibPETSc.AbstractPetscMat, dm::DMStag, rows, cols, vals, mode::InsertMode)
            @warn "set_values!(J, dm, rows, cols, vals, mode) is deprecated, use " *
                  "set_values!(J, dm, rows, cols, vals; mode)" maxlog = 1
            return set_values!(J, dm, rows, cols, vals; mode)
        end

        function set_values!(v::LibPETSc.AbstractPetscVec, dm::DMStag, positions, vals, mode::InsertMode)
            @warn "set_values!(v, dm, positions, vals, mode) is deprecated, use " *
                  "set_values!(v, dm, positions, vals; mode)" maxlog = 1
            return set_values!(v, dm, positions, vals; mode)
        end

        function zero_rows!(A::LibPETSc.AbstractPetscMat, rows_0b::AbstractVector{<:Integer}, diag; kwargs...)
            @warn "zero_rows!(A, rows, diag) is deprecated, use zero_rows!(A, rows; diag)" maxlog = 1
            return zero_rows!(A, rows_0b; diag, kwargs...)
        end

        function zero_rows_local!(A::LibPETSc.AbstractPetscMat, rows_0b::AbstractVector{<:Integer}, diag; kwargs...)
            @warn "zero_rows_local!(A, rows, diag) is deprecated, use zero_rows_local!(A, rows; diag)" maxlog = 1
            return zero_rows_local!(A, rows_0b; diag, kwargs...)
        end

        function zero_rows_local!(J::LibPETSc.AbstractPetscMat, dm::DMStag, rows, diag; kwargs...)
            @warn "zero_rows_local!(J, dm, rows, diag) is deprecated, use " *
                  "zero_rows_local!(J, dm, rows; diag)" maxlog = 1
            return zero_rows_local!(J, dm, rows; diag, kwargs...)
        end

        function project_function!(X, dm, time, funcs, ctxs, mode)
            @warn "project_function!(X, dm, time, funcs, ctxs, mode) is deprecated, use " *
                  "project_function!(X, dm, time, funcs; ctxs, mode)" maxlog = 1
            return project_function!(X, dm, time, funcs; ctxs, mode)
        end

        function project_field!(X, dm, time, U, funcs, mode)
            @warn "project_field!(X, dm, time, U, funcs, mode) is deprecated, use " *
                  "project_field!(X, dm, time, U, funcs; mode)" maxlog = 1
            return project_field!(X, dm, time, U, funcs; mode)
        end
        """,
        test = """
        @testset "positional insert mode, scatter mode and diag" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            comm = PETSc.MPI.COMM_SELF
            LibPETSc = PETSc.LibPETSc
            PetscInt = petsclib.PetscInt
            B = PETSc.DM_BOUNDARY_NONE

            # A serial ghosted vector with no ghosts: the updates have nothing to move
            gv = LibPETSc.VecCreateGhost(petsclib, comm, PetscInt(3), PetscInt(LibPETSc.PETSC_DECIDE),
                                         PetscInt(0), PetscInt[])
            for f in (PETSc.ghost_update!, PETSc.ghost_update_begin!, PETSc.ghost_update_end!)
                @test (warns() do
                    f(gv, PETSc.INSERT_VALUES, PETSc.SCATTER_FORWARD)
                end) === gv
            end

            da = PETSc.DMDA(petsclib, comm, (B,), (4,), 1)
            g = PETSc.global_vec(da)
            l = PETSc.local_vec(da)
            fill!(g, 1)
            warns(() -> PETSc.global_to_local!(l, da, g, PETSc.INSERT_VALUES))
            warns(() -> PETSc.local_to_local!(l, da, l, PETSc.INSERT_VALUES))
            warns(() -> PETSc.local_to_local!(l, da, PETSc.INSERT_VALUES))
            warns(() -> PETSc.local_to_global!(g, da, l, PETSc.ADD_VALUES))
            @test g[:] == fill(2, 4)

            A = PETSc.PetscMat(petsclib, 2, 2; type = :dense)
            warns(() -> PETSc.set_values!(A, [0], [0], [3.0], PETSc.INSERT_VALUES))
            PETSc.assemble!(A)
            @test A[1, 1] == 3
            warns(() -> PETSc.zero_rows!(A, [0], 5))
            @test A[1, 1] == 5

            J = PETSc.PetscMat(da)
            PETSc.assemble!(J)
            warns(() -> PETSc.zero_rows_local!(J, [0], 7))
            @test J[1, 1] == 7

            dm = PETSc.DMStag(petsclib, comm, (B,), (3,), (1, 1))
            loc = PETSc.element_location(dm)
            at = PETSc.stencil(dm, loc, (1,))
            K = PETSc.PetscMat(dm)
            warns(() -> PETSc.set_values!(K, dm, [at], [at], [2.0], PETSc.ADD_VALUES))
            PETSc.assemble!(K)
            warns(() -> PETSc.zero_rows_local!(K, dm, [at], 4))
            v = PETSc.global_vec(dm)
            warns(() -> PETSc.set_values!(v, dm, [at], [6.0], PETSc.INSERT_VALUES))

            # The projections need a DMPlex with a finite-element setup; the probe
            # reaches the warning and then fails in the forwarded call
            for call in (() -> PETSc.project_function!(ntuple(_ -> DeprecationProbe(), 6)...),
                         () -> PETSc.project_field!(ntuple(_ -> DeprecationProbe(), 6)...))
                warns() do
                    try
                        call()
                    catch
                    end
                end
            end

            foreach(PETSc.destroy!, (v, K, dm, J, A, l, g, da, gv))
        end
        """,
    ),
    (
        name = "DMDA and DMStag with a positional stencil width",
        code = """
        # §9: the stencil width and type have defaults, so they are keywords. The
        # positional forms warn for one release (§17.2, "argument dropped"). The
        # renamed keywords `setfromoptions`/`dmsetup` are caught in `finish_dm!`.
        function DMDA(
            petsclib::LibPETSc.PetscLibType,
            comm::MPI.Comm,
            boundary_type::Tuple,
            global_dim::Tuple,
            dof_per_node::Integer,
            stencil_width::Integer,
            stencil_type = :box;
            kwargs...,
        )
            @warn "DMDA(petsclib, comm, boundary, dims, dof, stencil_width, stencil_type) is " *
                  "deprecated, use DMDA(petsclib, comm, boundary, dims, dof; stencil_width, " *
                  "stencil_type)" maxlog = 1
            return DMDA(petsclib, comm, boundary_type, global_dim, dof_per_node;
                        stencil_width, stencil_type = something(stencil_type, :box), kwargs...)
        end

        function DMStag(
            petsclib::LibPETSc.PetscLibType,
            comm::MPI.Comm,
            boundary_type::Tuple,
            global_dim::Tuple,
            dof_per_node::Tuple,
            stencil_width::Integer,
            stencil_type = :box;
            kwargs...,
        )
            @warn "DMStag(petsclib, comm, boundary, dims, dof, stencil_width, stencil_type) is " *
                  "deprecated, use DMStag(petsclib, comm, boundary, dims, dof; stencil_width, " *
                  "stencil_type)" maxlog = 1
            return DMStag(petsclib, comm, boundary_type, global_dim, dof_per_node;
                          stencil_width, stencil_type, kwargs...)
        end
        """,
        test = """
        @testset "DMDA and DMStag with a positional stencil width" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            comm = PETSc.MPI.COMM_SELF
            LibPETSc = PETSc.LibPETSc
            B = PETSc.DM_BOUNDARY_NONE

            da = warns(() -> PETSc.DMDA(petsclib, comm, (B, B), (5, 5), 1, 2, :star))
            @test PETSc.info(da).stencil_width == 2
            @test PETSc.info(da).stencil_type == LibPETSc.DMDA_STENCIL_STAR
            dm = warns(() -> PETSc.DMStag(petsclib, comm, (B,), (5,), (1, 1), 2))
            @test LibPETSc.DMStagGetStencilWidth(petsclib, dm) == 2
            @test LibPETSc.DMStagGetStencilType(petsclib, dm) == LibPETSc.DMSTAG_STENCIL_BOX

            # The renamed keywords keep their meaning instead of becoming options
            for kw in (:dmsetup, :setfromoptions)
                d = warns(() -> PETSc.DMDA(petsclib, comm, (B,), (5,), 1; (kw => false,)...))
                PETSc.destroy!(d)
            end
            lazy = warns(() -> PETSc.DMStag(dm, (1, 1); dmsetfromoptions = false, dmsetup = false))
            PETSc.setup!(lazy)
            foreach(PETSc.destroy!, (lazy, dm, da))
        end
        """,
    ),
    (
        name = "add_coarsen_hook! with the hook first",
        code = """
        # The hooks are C pointers, which `do` cannot supply, so naming.md §8.1 does
        # not apply and the DM comes first. The hook-first order warns for one release.
        function add_coarsen_hook!(
            coarsenhook::Ptr{Cvoid},
            dm::AbstractPetscDM,
            restricthook::Ptr{Cvoid} = C_NULL,
        )
            @warn "add_coarsen_hook!(coarsenhook, dm, restricthook) is deprecated, use " *
                  "add_coarsen_hook!(dm, coarsenhook, restricthook)" maxlog = 1
            return add_coarsen_hook!(dm, coarsenhook, restricthook)
        end
        """,
        test = """
        @testset "add_coarsen_hook! with the hook first" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            B = PETSc.DM_BOUNDARY_NONE
            da = PETSc.DMDA(petsclib, PETSc.MPI.COMM_SELF, (B,), (9,), 1)
            # PETSc skips a null hook, so the shim's forwarding is all this exercises
            @test warns(() -> PETSc.add_coarsen_hook!(Ptr{Cvoid}(C_NULL), da)) === da
            PETSc.destroy!(da)
        end
        """,
    ),
    (
        name = "positional own",
        code = """
        # §9: `own` has a default, so it is a keyword, as in `narrow(dm; own)`.
        function VecPtr(petsclib::LibPETSc.PetscLibType, ptr::CVec, own::Bool)
            @warn "VecPtr(petsclib, ptr, own) is deprecated, use VecPtr(petsclib, ptr; own)" maxlog = 1
            return VecPtr(petsclib, ptr; own)
        end

        function MatPtr(petsclib::LibPETSc.PetscLibType, ptr::CMat, own::Bool)
            @warn "MatPtr(petsclib, ptr, own) is deprecated, use MatPtr(petsclib, ptr; own)" maxlog = 1
            return MatPtr(petsclib, ptr; own)
        end

        function LibPETSc.PetscVec(petsclib::LibPETSc.PetscLibType, ptr::CVec, own::Bool)
            @warn "PetscVec(petsclib, ptr, own) is deprecated, use PetscVec(petsclib, ptr; own)" maxlog = 1
            return LibPETSc.PetscVec(petsclib, ptr; own)
        end
        """,
        test = """
        @testset "positional own" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            v = PETSc.PetscVec(petsclib, 3)
            m = PETSc.PetscMat(petsclib, 3, 3, 1)
            @test !warns(() -> PETSc.VecPtr(petsclib, v.ptr, false)).own
            @test !warns(() -> PETSc.MatPtr(petsclib, m.ptr, false)).own
            @test !warns(() -> PETSc.PetscVec(petsclib, v.ptr, false)).own
            foreach(PETSc.destroy!, (m, v))
        end
        """,
    ),
    (
        name = "CSR PetscMat with comm as a keyword",
        code = """
        # §9: `comm` is positional second, as in every other constructor. The
        # keyword form is caught in the CSR `PetscMat` method in src/mat.jl.
        """,
        test = """
        @testset "CSR PetscMat with comm as a keyword" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            PetscInt, PetscScalar = petsclib.PetscInt, petsclib.PetscScalar
            A = warns(() -> PETSc.PetscMat(petsclib, PetscInt[0, 1], PetscInt[0], PetscScalar[2];
                                           comm = PETSc.MPI.COMM_SELF))
            @test A[1, 1] == PetscScalar(2)
            PETSc.destroy!(A)
        end
        """,
    ),    (
        name = "MatShell with the library first",
        code = """
        # §8.1: the operator comes first so `do` works. The old order dispatches on
        # the library in argument 1.
        function MatShell(
            petsclib::LibPETSc.PetscLibType,
            obj,
            comm::MPI.Comm,
            local_rows::Integer,
            local_cols::Integer,
            global_rows::Integer = LibPETSc.PETSC_DECIDE,
            global_cols::Integer = LibPETSc.PETSC_DECIDE,
        )
            @warn "MatShell(petsclib, obj, comm, local_rows, local_cols, global_rows, global_cols) " *
                  "is deprecated, use MatShell(obj, petsclib, comm, local_rows, local_cols; " *
                  "global_rows, global_cols)" maxlog = 1
            return MatShell(obj, petsclib, comm, local_rows, local_cols; global_rows, global_cols)
        end

        # Two libraries match both orders. A library is never an operator.
        function MatShell(
            ::LibPETSc.PetscLibType,
            ::LibPETSc.PetscLibType,
            ::MPI.Comm,
            ::Integer,
            ::Integer,
        )
            throw(ArgumentError("MatShell takes the operator first, then the library"))
        end
        """,
        test = """
        @testset "MatShell with the library first" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            comm = PETSc.MPI.COMM_SELF
            A = warns(() -> PETSc.MatShell(petsclib, (y, x) -> (y .= 3 .* x), comm, 2, 2))
            x = PETSc.PetscVec(petsclib, petsclib.PetscScalar[1, 2])
            y = PETSc.PetscVec(petsclib, 2)
            PETSc.LinearAlgebra.mul!(y, A, x)
            @test y[:] == petsclib.PetscScalar[3, 6]
            @test_throws "operator first" PETSc.MatShell(petsclib, petsclib, comm, 2, 2)
            foreach(PETSc.destroy!, (y, x, A))
        end
        """,
    ),    (
        name = "local_arrays and restore_local_arrays!",
        code = """
        # The pair warns in src/vec.jl: `with_local_arrays!` replaces it and restores
        # the arrays also when the kernel throws.
        """,
        test = """
        @testset "local_arrays and restore_local_arrays!" begin
            petsclib = PETSc.petsclibs[1]
            PETSc.initialize(petsclib)
            T = petsclib.PetscScalar
            x = PETSc.PetscVec(petsclib, T[1, 2])
            fx = PETSc.PetscVec(petsclib, 2)
            arrays = warns(() -> PETSc.local_arrays(petsclib, fx, x))
            arrays[1] .= 3 .* arrays[2]
            @test warns(() -> PETSc.restore_local_arrays!(petsclib, fx, x, arrays...)) === nothing
            @test fx[:] == T[3, 6]
            foreach(PETSc.destroy!, (fx, x))
        end
        """,
    ),
]

# ---------------------------------------------------------------------------
# Internal helpers: renamed freely, no shim, not part of the API (§1.1).
# ---------------------------------------------------------------------------

const INTERNAL = Set{Symbol}([
    # the checkout behind with_local_arrays! and the deprecated local_arrays
    :checkout_local_arrays,
    # the per-library body of the MatShell constructor
    :create_matshell,
    # leading underscore dropped (naming.md "Internals")
    :petsc_link,
    :petsc_subst,
    :vtk_merge_one_tensor,
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
    :apply_blas_num_threads,
    :BLAS_THREAD_VARIABLES,
    :ranks_on_node,
    :_unsafe_local_array,
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
    :take_renamed_dm_keywords,
    :RENAMED_DM_KEYWORDS,
    :destroy_options!,
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
    :with_broadcast_arrays,
    :aliases_dest,
    :check_undistributed,
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
    :c_vector,
    :check_block,
    :checkout_arrays,
    :checked_array,
    :flags,
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
    :create_split_boundary_labels!,
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
    :VecPtr,
    :MatPtr,
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
    :TSSetRHSFunctionFn,
    :TSSetRHSJacobianFn,
    :TSSetIFunctionFn,
    :TSSetIJacobianFn,
    :TSMonitorSetFn,
    # init.jl / options.jl
    :PetscNotInitialized,
    :ConvergenceError,       # new in v0.5.3
    :with_options,           # new in v0.5.3
    :with_local_arrays!,     # new in v0.5.5
    :with_destroy,           # new in v0.5.3
    :track!,                 # new in v0.5.3
    :untrack!,               # new in v0.5.3
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
