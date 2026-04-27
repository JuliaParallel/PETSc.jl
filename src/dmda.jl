import .LibPETSc: AbstractPetscDM, PetscDM, CDM


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

If keyword argument `setfromoptions == true` then `setfromoptions!` called.

If keyword argument `dmsetup == true` then `setup!` is called.

When `D == 1` the `stencil_type` argument is not required and ignored if specified.

# External Links
$(_doc_external("DMDA/DMDACreate1d"))
$(_doc_external("DMDA/DMDACreate2d"))
$(_doc_external("DMDA/DMDACreate3d"))
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
    opts = Options(petsclib; options...)

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
            @assert points_per_proc[d] isa Array{PetscLib.PetscInt}
            @assert length(points_per_proc[d]) == MPI.Comm_size(comm)
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

    if !isempty(prefix)
        # options prefix
        LibPETSc.DMSetOptionsPrefix(petsclib, da, prefix)
    end

    if setfromoptions
        # set options (if any)
        opts = PETSc.Options(petsclib; options...);
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
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, da)
    end
    return da
end

"""
    ndofs(da::AbstractPetscDM)

Return the number of dofs in for `da`

# External Links
$(_doc_external("DMDA/DMDAGetDof"))
"""
function ndofs(da::AbstractPetscDM{PetscLib}) where PetscLib
    PetscInt = PetscLib.PetscInt
    ndof = [PetscInt(0)]

    # number of DOF that the DMDA has
    ndof = LibPETSc.DMDAGetDof(PetscLib, da)
    return @inbounds ndof[1]
end


"""
    reshapelocalarray(Arr, da::AbstractPetscDM{PetscLib}, ndof = ndofs(da))

Returns an array with the same data as `Arr` but reshaped as an array that can
be addressed with global indexing.
"""
function reshapelocalarray(
    Arr,
    da::AbstractPetscDM{PetscLib},
    ndof::Integer = ndofs(da),
) where {PetscLib}

    # First we try to use a ghosted size
    corners = getghostcorners(da)
    # If this is two big for the array use non-ghosted
    if length(Arr) < prod(corners.size) * ndof
        corners = getcorners(da)
    end
    @assert length(Arr) == prod(corners.size) * ndof

    oArr = OffsetArray(
        reshape(Arr, Int64(ndof), Int64.(corners.size)...),
        1:ndof,
        (corners.lower[1]):(corners.upper[1]),
        (corners.lower[2]):(corners.upper[2]),
        (corners.lower[3]):(corners.upper[3]),
    )

    return oArr
end

"""
    ind = localinteriorlinearindex(dmda::AbstractPetscDM)

Returns the linear indices associated with the degrees of freedom own by this MPI rank embedded in the ghost index space for the `dmda`
"""
function localinteriorlinearindex(da::AbstractPetscDM{PetscLib}) where PetscLib
    # Determine the indices of the linear indices of the local part of the
    # matrix we own
    @assert gettype(da) == "da" 
    ghost_corners = PETSc.getghostcorners(da)
    corners = PETSc.getcorners(da)

    # First compute the Cartesian indices for the local portion we own
    offset = ghost_corners.lower - CartesianIndex(1, 1, 1)
    l_inds = ((corners.lower):(corners.upper)) .- offset

    # Create a grid of indices with ghost then extract only the local part
    lower = CartesianIndex(1, ghost_corners.lower)
    upper = CartesianIndex(ndofs(da), ghost_corners.upper)
    ind_local = LinearIndices(lower:upper)[:, l_inds][:]
    return ind_local
end
"""
    dmda_star_fd_coloring(petsclib, da)

Build all data needed for manual FD coloring of a **2-D** DMDA with a STAR
stencil, using `IS_COLORING_LOCAL` and ghost-local COO indexing.

Specifically, this function:
1. Creates an `IS_COLORING_LOCAL` `ISColoring` via `DMCreateColoring` and
   extracts the per-DOF color vector (ghost-local layout).
2. Enumerates all STAR-stencil (row, col) pairs for every owned node and
   records their ghost-local 0-based indices and colors.
3. Builds per-color index arrays (`perturb_cols`, `coo_idxs`, `local_rows`)
   ready for use in an FD coloring Newton loop.

Returns a `NamedTuple`:
- `n_colors`      — number of colors
- `n_local_dofs`  — number of locally owned DOFs (owned nodes × dof/node)
- `nnz_coo`       — total number of COO entries
- `row_coo_local` — ghost-local 0-based row indices (`Vector{PetscInt}`)
- `col_coo_local` — ghost-local 0-based column indices (`Vector{PetscInt}`)
- `perturb_cols`  — `perturb_cols[c]`: 1-based owned-local column indices
                    with color `c-1`; used to scatter `+h` perturbations.
- `coo_idxs`      — `coo_idxs[c]`: 1-based COO entry indices for color `c-1`
- `local_rows`    — `local_rows[c]`: corresponding 1-based owned-local
                    residual-row indices; used to read `(f1-f0)/h`.

!!! note "2-D DMDA STAR stencil only"
    Neighbor enumeration covers only `±x` and `±y` directions.  The ghost-local
    flat-index formula is `d + ix*dof + iy*dof*nx_g`.  For a 3-D DMDA:
    - add `±z` neighbors guarded by `kk > 1` / `kk < mz`,
    - extend the flat-index formula with `+ iz*dof*nx_g*ny_g`,
    - reshape `col_colors_mat` to `(dof, nx_g, ny_g, nz_g)`,
    - decode `z_owned` in the `perturb_cols` loop.
"""
function dmda_star_fd_coloring(petsclib::PetscLib, da::AbstractPetscDM{PetscLib}) where PetscLib
    CPetscInt = petsclib.PetscInt

    # ── ISColoring ────────────────────────────────────────────────────────────
    # IS_COLORING_LOCAL covers owned + ghost DOFs; ghost colors are consistent
    # with the owning rank so no extra MPI communication is needed here.
    iscoloring = LibPETSc.DMCreateColoring(petsclib, da, LibPETSc.IS_COLORING_LOCAL)
    _, nc_pi, col_colors_local = LibPETSc.ISColoringGetColors(petsclib, iscoloring)
    n_colors = Int(nc_pi)
    LibPETSc.ISColoringDestroy(petsclib, iscoloring)

    # ── DMDA geometry ─────────────────────────────────────────────────────────
    info          = getinfo(da)
    mx            = Int(info.global_size[1])
    my            = Int(info.global_size[2])
    dof_per_node  = Int(info.dof)
    corners       = getcorners(da)
    ghost_corners = getghostcorners(da)
    xs_da  = corners.lower[1];       ys_da  = corners.lower[2]
    xe_da  = corners.upper[1];       ye_da  = corners.upper[2]
    xsg_da = ghost_corners.lower[1]; ysg_da = ghost_corners.lower[2]
    xeg_da = ghost_corners.upper[1]; yeg_da = ghost_corners.upper[2]
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
    row_coo_local     = CPetscInt[]
    col_coo_local     = CPetscInt[]
    local_row_per_coo = CPetscInt[]  # 0-based owned-local row (VecGetArray index)
    color_per_coo     = CPetscInt[]  # 0-based color of each COO entry's column

    for jj in ys_da:ye_da, ii in xs_da:xe_da
        ix_gh = ii - xsg_da   # 0-based ghost-x of this owned node
        iy_gh = jj - ysg_da   # 0-based ghost-y of this owned node
        ix_ow = ii - xs_da    # 0-based owned-x
        iy_ow = jj - ys_da    # 0-based owned-y

        neighbors = Tuple{Int,Int}[(ii, jj)]
        ii > 1  && push!(neighbors, (ii-1, jj))
        ii < mx && push!(neighbors, (ii+1, jj))
        jj > 1  && push!(neighbors, (ii, jj-1))
        jj < my && push!(neighbors, (ii, jj+1))

        for (ni, nj) in neighbors
            nix_gh = ni - xsg_da
            njy_gh = nj - ysg_da
            for d_row in 0:dof_per_node-1, d_col in 0:dof_per_node-1
                r_local = d_row + ix_gh  * dof_per_node + iy_gh  * dof_per_node * nx_g_da
                c_local = d_col + nix_gh * dof_per_node + njy_gh * dof_per_node * nx_g_da
                p_owned = d_row + ix_ow  * dof_per_node + iy_ow  * nx_own       * dof_per_node
                push!(row_coo_local,     CPetscInt(r_local))
                push!(col_coo_local,     CPetscInt(c_local))
                push!(color_per_coo,     CPetscInt(col_colors_mat[d_col+1, nix_gh+1, njy_gh+1]))
                push!(local_row_per_coo, CPetscInt(p_owned))
            end
        end
    end
    nnz_coo = length(row_coo_local)

    # ── Per-color index arrays ────────────────────────────────────────────────
    # perturb_cols[c]: 1-based owned-local column indices with color c-1.
    # col_colors_local uses the ghost-local layout, but VecGetArray returns only
    # owned DOFs re-indexed 1..n_local_dofs.  We convert owned-local → ghost-local
    # before looking up the color.
    perturb_cols = [Int32[] for _ in 1:n_colors]
    for p_local in 1:n_local_dofs
        p0      = p_local - 1
        d       = p0 % dof_per_node
        x_owned = (p0 ÷ dof_per_node) % nx_own
        y_owned = (p0 ÷ dof_per_node) ÷ nx_own
        k_ghost = d + (x_owned + ox_coo) * dof_per_node +
                      (y_owned + oy_coo) * dof_per_node * nx_g_da + 1
        c = Int(col_colors_local[k_ghost]) + 1
        push!(perturb_cols[c], Int32(p_local))
    end

    coo_idxs   = [Int32[] for _ in 1:n_colors]
    local_rows = [Int32[] for _ in 1:n_colors]
    for k in 1:nnz_coo
        c = Int(color_per_coo[k]) + 1
        push!(coo_idxs[c],   Int32(k))
        push!(local_rows[c], Int32(local_row_per_coo[k] + 1))
    end

    return (;
        n_colors,
        n_local_dofs,
        nnz_coo,
        row_coo_local,
        col_coo_local,
        perturb_cols,
        coo_idxs,
        local_rows,
    )
end