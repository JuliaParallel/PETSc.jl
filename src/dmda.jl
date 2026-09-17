import .LibPETSc: AbstractPetscDM, PetscDM, CDM

# `stencil_type` is a PETSc enum (§3.1's enum rule), and naming.md also spells
# it as a `Symbol` in the constructor examples. Both work: the Symbol is mapped
# here, the enum passes through.
stencil_type_enum(::Type{DMDAStencilType}, t::DMDAStencilType) = t
stencil_type_enum(::Type{DMDAStencilType}, ::Nothing) = nothing
function stencil_type_enum(::Type{DMDAStencilType}, t::Symbol)
    t === :star && return DMDA_STENCIL_STAR
    t === :box && return DMDA_STENCIL_BOX
    throw(ArgumentError("unknown DMDA stencil type :$t, expected :star or :box"))
end

stencil_type_enum(::Type{DMStagStencilType}, t::DMStagStencilType) = t
function stencil_type_enum(::Type{DMStagStencilType}, t::Symbol)
    t === :none && return DMSTAG_STENCIL_NONE
    t === :star && return DMSTAG_STENCIL_STAR
    t === :box && return DMSTAG_STENCIL_BOX
    throw(ArgumentError(
        "unknown DMStag stencil type :$t, expected :none, :star or :box",
    ))
end


"""
    DMDA(
        petsclib::PetscLib
        comm::MPI.Comm,
        boundary_type::NTuple{D, DMBoundaryType},
        global_dim::NTuple{D, Integer},
        dof_per_node::Integer,
        stencil_width::Integer,
        stencil_type;
        points_per_proc::Tuple,
        processors::Tuple,
        setfromoptions = true,
        dmsetup = true,
        prefix = "",
        options...
    )

Creates a `D`-dimensional distributed array with the options specified using
keyword arguments.

If keyword argument `points_per_proc[k] isa Vector{petsclib.PetscInt}` then this
specifies the points per processor in dimension `k`.

If keyword argument `processors[k] isa Integer` then this specifies the number of
processors used in dimension `k`; ignored when `D == 1`.

If keyword argument `setfromoptions == true` then `set_from_options!` called.

If keyword argument `dmsetup == true` then `setup!` is called.

When `D == 1` the `stencil_type` argument is not required and ignored if specified.

# External Links
$(doc_external("DMDA/DMDACreate1d"))
$(doc_external("DMDA/DMDACreate2d"))
$(doc_external("DMDA/DMDACreate3d"))
"""
function DMDA(
    petsclib::PetscLib,
    comm::MPI.Comm,
    boundary_type::NTuple{N, DMBoundaryType},
    global_dim::NTuple{N, Integer},
    dof_per_node::Integer,
    stencil_width::Integer,
    stencil_type = nothing;
    points_per_proc::Union{Tuple, Nothing} = nothing,
    processors = nothing,
    setfromoptions = true,
    dmsetup = true,
    prefix = "",
    options...,
) where {PetscLib, N}
    PetscInt = inttype(PetscLib)
    stencil_type = stencil_type_enum(DMDAStencilType, stencil_type)

    if isnothing(points_per_proc)
        points_per_proc = ntuple(_ -> nothing, N)
    end
    if isnothing(processors)
        processors = ntuple(_ -> PETSC_DECIDE, N)
    end

    ref_points_per_proc = ntuple(N) do d
        if isnothing(points_per_proc[d]) || points_per_proc[d] == PETSC_DECIDE
            C_NULL
        else
            points_per_proc[d] isa Array{PetscLib.PetscInt} || throw(
                ArgumentError(
                    "points_per_proc[$d] must be an Array{$(PetscLib.PetscInt)}, " *
                    "got $(typeof(points_per_proc[d]))",
                ),
            )
            length(points_per_proc[d]) == MPI.Comm_size(comm) || throw(
                DimensionMismatch(
                    "points_per_proc[$d] has $(length(points_per_proc[d])) entries, " *
                    "but the communicator has $(MPI.Comm_size(comm)) ranks",
                ),
            )
            points_per_proc[d]
        end
    end
    
    if N==1  
        da = LibPETSc.DMDACreate1d(petsclib,
                                   comm, 
                                   boundary_type[1], 
                                   PetscInt(global_dim[1]), 
                                   PetscInt(dof_per_node), 
                                   PetscInt(stencil_width), 
                                   ref_points_per_proc[1]
                                   )
    elseif N==2
        da =   LibPETSc.DMDACreate2d(
                                    petsclib,
                                    comm,
                                    boundary_type[1], boundary_type[2],
                                    stencil_type,
                                    PetscInt(global_dim[1]), PetscInt(global_dim[2]),
                                    PetscInt(processors[1]), PetscInt(processors[2]),
                                    PetscInt(dof_per_node),
                                    PetscInt(stencil_width),
                                    ref_points_per_proc[1], ref_points_per_proc[2]
                                )                                    
     elseif N==3
        da =   LibPETSc.DMDACreate3d(
                                    petsclib,
                                    comm,
                                    boundary_type[1], boundary_type[2], boundary_type[3],
                                    stencil_type,
                                    PetscInt(global_dim[1]), PetscInt(global_dim[2]), PetscInt(global_dim[3]),
                                     PetscInt(processors[1]), PetscInt(processors[2]), PetscInt(processors[3]),
                                    PetscInt(dof_per_node),
                                    PetscInt(stencil_width),
                                    ref_points_per_proc[1], ref_points_per_proc[2], ref_points_per_proc[3]
                                )                    
    end

    # Take ownership of the handle the creator returned: from here on `da` is a
    # `DMDA{PetscLib, N}` and every DMDA method dispatches on it.
    da = DMDA{PetscLib, N}(da.ptr, petsclib.age, true)

    if !isempty(prefix)
        # options prefix
        LibPETSc.DMSetOptionsPrefix(petsclib, da, prefix)
    end

    if setfromoptions
        # set options (if any)
        opts = PetscOptions(petsclib; options...)
        push!(opts)
        LibPETSc.DMSetFromOptions(PetscLib, da)
        pop!(opts)
    end
    
    if dmsetup
        # initialize dmda
        setup!(da)            
    end

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    return own_dm!(da, comm)
end

"""
    ndofs(da::DMDA)

Return the number of dofs in for `da`

# External Links
$(doc_external("DMDA/DMDAGetDof"))
"""
function ndofs(da::DMDA{PetscLib}) where PetscLib
    PetscInt = PetscLib.PetscInt
    ndof = [PetscInt(0)]

    # number of DOF that the DMDA has
    ndof = LibPETSc.DMDAGetDof(PetscLib, da)
    return @inbounds ndof[1]
end


"""
    reshape_local_array(Arr, da::Union{DMDA, DMStag}, ndof = ndofs(da))

Returns an array with the same data as `Arr` but reshaped as an array that can
be addressed with global indexing.
"""
function reshape_local_array(
    Arr,
    da::Union{DMDA{PetscLib, N}, DMStag{PetscLib, N}},
    ndof::Integer = ndofs(da),
) where {PetscLib, N}

    # First we try to use a ghosted size
    gc = ghost_corners(da)
    c = length(Arr) < prod(gc.size) * ndof ? corners(da) : gc
    length(Arr) == prod(c.size) * ndof || throw(
        DimensionMismatch(
            "array has $(length(Arr)) entries, but the local domain needs " *
            "$(prod(c.size) * ndof) ($(c.size) points x $ndof dofs)",
        ),
    )

    # `N + 1` axes: the dof axis, then one axis per dimension of the DM. v0.4
    # always built four axes, which only worked because everything was padded
    # to three dimensions (§12).
    return OffsetArray(
        reshape(Arr, Int64(ndof), ntuple(i -> Int64(c.size[i]), Val(N))...),
        1:Int64(ndof),
        ntuple(i -> (c.lower[i]):(c.upper[i]), Val(N))...,
    )
end

"""
    ind = local_interior_linear_index(dmda::DMDA)

Returns the linear indices associated with the degrees of freedom own by this MPI rank embedded in the ghost index space for the `dmda`
"""
function local_interior_linear_index(da::DMDA{PetscLib, N}) where {PetscLib, N}
    # Determine the indices of the linear indices of the local part of the
    # matrix we own
    gc = ghost_corners(da)
    c = corners(da)

    # First compute the Cartesian indices for the local portion we own
    offset = gc.lower - oneunit(CartesianIndex{N})
    l_inds = ((c.lower):(c.upper)) .- offset

    # Create a grid of indices with ghost then extract only the local part
    lower = CartesianIndex(1, gc.lower)
    upper = CartesianIndex(ndofs(da), gc.upper)
    ind_local = LinearIndices(lower:upper)[:, l_inds][:]
    return ind_local
end

"""
    star_fd_coloring(da::DMDA{PetscLib, 2})

Build all data needed for manual FD coloring of a **2-D** DMDA with a STAR
stencil, using `IS_COLORING_LOCAL` and ghost-local COO indexing.

`petsclib` is not an argument: the `DMDA` carries it as a type parameter
(docs/src/man/naming.md §8). The dimension is a type parameter too, so a DMDA
of any other dimension is a `MethodError` rather than a silently wrong answer.

Specifically, this function:
1. Creates an `IS_COLORING_LOCAL` `ISColoring` via `DMCreateColoring` and
   extracts the per-DOF color vector (ghost-local layout).
2. Enumerates all STAR-stencil (row, col) pairs for every owned node and
   records their ghost-local 0-based indices and colors.
3. Builds per-color index arrays (`perturb_cols_1b`, `coo_idxs_1b`, `local_rows_1b`)
   ready for use in an FD coloring Newton loop.

Returns a `NamedTuple`. The index vectors carry the base they use in their name
(§12.1), because this function returns both: the COO arrays go straight to
`LibPETSc.MatSetPreallocationCOOLocal` and keep PETSc's base, while the rest
index Julia arrays and are 1-based.
- `n_colors`         — number of colors
- `n_local_dofs`     — number of locally owned DOFs (owned nodes × dof/node)
- `nnz_coo`          — total number of COO entries
- `row_coo_local_0b` — ghost-local 0-based row indices (`Vector{PetscInt}`)
- `col_coo_local_0b` — ghost-local 0-based column indices (`Vector{PetscInt}`)
- `perturb_cols_1b`  — `perturb_cols_1b[c]`: 1-based owned-local column indices
                       with color `c-1`; used to scatter `+h` perturbations.
- `coo_idxs_1b`      — `coo_idxs_1b[c]`: 1-based COO entry indices for color `c-1`
- `local_rows_1b`    — `local_rows_1b[c]`: corresponding 1-based owned-local
                       residual-row indices; used to read `(f1-f0)/h`.

!!! note "2-D DMDA STAR stencil only"
    Neighbor enumeration covers only `±x` and `±y` directions.  The ghost-local
    flat-index formula is `d + ix*dof + iy*dof*nx_g`.  For a 3-D DMDA:
    - add `±z` neighbors guarded by `kk > 1` / `kk < mz`,
    - extend the flat-index formula with `+ iz*dof*nx_g*ny_g`,
    - reshape `col_colors_mat` to `(dof, nx_g, ny_g, nz_g)`,
    - decode `z_owned` in the `perturb_cols_1b` loop.
"""
function star_fd_coloring(da::DMDA{PetscLib, 2}) where {PetscLib}
    petsclib = getlib(PetscLib)
    CPetscInt = petsclib.PetscInt

    # ── ISColoring ────────────────────────────────────────────────────────────
    # IS_COLORING_LOCAL covers owned + ghost DOFs; ghost colors are consistent
    # with the owning rank so no extra MPI communication is needed here.
    iscoloring = LibPETSc.DMCreateColoring(petsclib, da, LibPETSc.IS_COLORING_LOCAL)
    _, nc_pi, col_colors_local = LibPETSc.ISColoringGetColors(petsclib, iscoloring)
    n_colors = Int(nc_pi)
    LibPETSc.ISColoringDestroy(petsclib, iscoloring)

    # ── DMDA geometry ─────────────────────────────────────────────────────────
    da_info = info(da)
    mx            = Int(da_info.global_size[1])
    my            = Int(da_info.global_size[2])
    dof_per_node  = Int(da_info.ndofs)
    c             = corners(da)
    gc            = ghost_corners(da)
    xs_da  = c.lower[1];       ys_da  = c.lower[2]
    xe_da  = c.upper[1];       ye_da  = c.upper[2]
    xsg_da = gc.lower[1]; ysg_da = gc.lower[2]
    xeg_da = gc.upper[1]; yeg_da = gc.upper[2]
    nx_g_da = xeg_da - xsg_da + 1
    ny_g_da = yeg_da - ysg_da + 1
    nx_own  = xe_da  - xs_da  + 1
    ny_own  = ye_da  - ys_da  + 1
    ox_coo  = xs_da  - xsg_da   # ghost offset in x (grid nodes)
    oy_coo  = ys_da  - ysg_da   # ghost offset in y (grid nodes)
    n_local_dofs = nx_own * ny_own * dof_per_node

    # col_colors_mat[d+1, ix_ghost+1, iy_ghost+1] → color of that ghost DOF
    col_colors_mat = reshape(col_colors_local, dof_per_node, nx_g_da, ny_g_da)

    # ── COO triplets from 2-D STAR stencil ────────────────────────────────────
    # Ghost-local 0-based row/col indices (for MatSetPreallocationCOOLocal).
    # Both use DMDA ghost-local numbering: idx = d + ix_g*dof + iy_g*dof*nx_g.
    #
    # Pre-compute exact nnz_coo analytically so we can allocate once and fill
    # by index rather than growing via push! (avoids O(nnz) realloc+copy work).
    # Each owned node contributes dof² entries per STAR neighbor it has.
    # Interior nodes have 5 neighbors (self + 4); boundary nodes have fewer.
    nbr_left   = xs_da == 1  ? nx_own - 1 : nx_own   # columns with a left neighbor
    nbr_right  = xe_da == mx ? nx_own - 1 : nx_own   # columns with a right neighbor
    nbr_bottom = ys_da == 1  ? ny_own - 1 : ny_own   # rows with a bottom neighbor
    nbr_top    = ye_da == my ? ny_own - 1 : ny_own   # rows with a top neighbor
    total_nbr_pairs = nx_own * ny_own +               # self
                      nbr_left   * ny_own +
                      nbr_right  * ny_own +
                      nbr_bottom * nx_own +
                      nbr_top    * nx_own
    nnz_coo = total_nbr_pairs * dof_per_node^2

    row_coo_local_0b     = Vector{CPetscInt}(undef, nnz_coo)
    col_coo_local_0b     = Vector{CPetscInt}(undef, nnz_coo)
    local_row_per_coo = Vector{CPetscInt}(undef, nnz_coo)  # 0-based owned-local row
    color_per_coo     = Vector{CPetscInt}(undef, nnz_coo)  # 0-based color of column

    # Inner helper: write dof² entries for a (row-node, col-node) pair.
    # Captures ix_gh, iy_gh, ix_ow, iy_ow from the outer scope.
    @inline function fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, nix_gh, njy_gh)
        r_base = ix_gh  * dof_per_node + iy_gh  * dof_per_node * nx_g_da
        c_base = nix_gh * dof_per_node + njy_gh * dof_per_node * nx_g_da
        p_base = ix_ow  * dof_per_node + iy_ow  * nx_own       * dof_per_node
        for d_row in 0:dof_per_node-1
            r_local = CPetscInt(d_row + r_base)
            p_owned = CPetscInt(d_row + p_base)
            for d_col in 0:dof_per_node-1
                @inbounds begin
                    row_coo_local_0b[k]     = r_local
                    col_coo_local_0b[k]     = CPetscInt(d_col + c_base)
                    color_per_coo[k]     = CPetscInt(col_colors_mat[d_col+1, nix_gh+1, njy_gh+1])
                    local_row_per_coo[k] = p_owned
                end
                k += 1
            end
        end
        return k
    end

    k = 1
    for jj in ys_da:ye_da, ii in xs_da:xe_da
        ix_gh = ii - xsg_da   # 0-based ghost-x of this owned node
        iy_gh = jj - ysg_da   # 0-based ghost-y of this owned node
        ix_ow = ii - xs_da    # 0-based owned-x
        iy_ow = jj - ys_da    # 0-based owned-y

        # self
        k = fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, ix_gh, iy_gh)
        # left
        ii > 1  && (k = fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, ix_gh - 1, iy_gh))
        # right
        ii < mx && (k = fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, ix_gh + 1, iy_gh))
        # bottom
        jj > 1  && (k = fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, ix_gh, iy_gh - 1))
        # top
        jj < my && (k = fill_nbr!(k, ix_gh, iy_gh, ix_ow, iy_ow, ix_gh, iy_gh + 1))
    end

    # ── Per-color index arrays ────────────────────────────────────────────────
    # perturb_cols_1b[c]: 1-based owned-local column indices with color c-1.
    # col_colors_local uses the ghost-local layout, but VecGetArray returns only
    # owned DOFs re-indexed 1..n_local_dofs.  We convert owned-local → ghost-local
    # before looking up the color.
    hint_cols = max(1, n_local_dofs ÷ n_colors)
    perturb_cols_1b = [sizehint!(Int32[], hint_cols) for _ in 1:n_colors]
    for p_local in 1:n_local_dofs
        p0      = p_local - 1
        d       = p0 % dof_per_node
        x_owned = (p0 ÷ dof_per_node) % nx_own
        y_owned = (p0 ÷ dof_per_node) ÷ nx_own
        k_ghost = d + (x_owned + ox_coo) * dof_per_node +
                      (y_owned + oy_coo) * dof_per_node * nx_g_da + 1
        c = Int(col_colors_local[k_ghost]) + 1
        push!(perturb_cols_1b[c], Int32(p_local))
    end

    hint_coo = max(1, nnz_coo ÷ n_colors)
    coo_idxs_1b   = [sizehint!(Int32[], hint_coo) for _ in 1:n_colors]
    local_rows_1b = [sizehint!(Int32[], hint_coo) for _ in 1:n_colors]
    for k in 1:nnz_coo
        c = Int(color_per_coo[k]) + 1
        push!(coo_idxs_1b[c],   Int32(k))
        push!(local_rows_1b[c], Int32(local_row_per_coo[k] + 1))
    end

    return (;
        n_colors,
        n_local_dofs,
        nnz_coo,
        row_coo_local_0b,
        col_coo_local_0b,
        perturb_cols_1b,
        coo_idxs_1b,
        local_rows_1b,
    )
end