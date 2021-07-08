mutable struct DMDALocalInfo{IT}
    dim::IT
    dof_per_node::IT
    stencil_width::IT
    global_size::NTuple{3, IT}
    local_start::NTuple{3, IT}
    local_size::NTuple{3, IT}
    ghosted_local_start::NTuple{3, IT}
    ghosted_local_size::NTuple{3, IT}
    boundary_type::NTuple{3, DMBoundaryType}
    stencil_type::DMDAStencilType
    ___padding___::NTuple{5, IT}
    DMDALocalInfo{IT}() where {IT} = new{IT}()
end

"""
    DMDACreate1d(
        ::PetscLib
        comm::MPI.Comm,
        boundary_type::DMBoundaryType,
        global_dim,
        dof_per_node,
        stencil_width,
        points_per_proc::Union{Nothing, Vector{PetscInt}};
        options...
    )

Creates a 1-D distributed array with the options specified using keyword
arguments.

# External Links
$(_doc_external("DMDA/DMDACreate1d"))
"""
function DMDACreate1d end

@for_petsc function DMDACreate1d(
    ::$UnionPetscLib,
    comm::MPI.Comm,
    boundary_type::DMBoundaryType,
    global_dim,
    dof_per_node,
    stencil_width,
    points_per_proc::Union{Nothing, Vector{$PetscInt}};
    options...,
)
    opts = Options($petsclib, options...)
    ref_points_per_proc = if isnothing(points_per_proc)
        C_NULL
    else
        @assert length(points_per_proc) == MPI.Comm_size(comm)
        points_per_proc
    end
    da = DM{$PetscLib}(C_NULL, opts)
    with(da.opts) do
        @chk ccall(
            (:DMDACreate1d, $petsc_library),
            PetscErrorCode,
            (
                MPI.MPI_Comm,
                DMBoundaryType,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                Ptr{$PetscInt},
                Ptr{CDM},
            ),
            comm,
            boundary_type,
            global_dim,
            dof_per_node,
            stencil_width,
            ref_points_per_proc,
            da,
        )
    end
    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if comm == MPI.COMM_SELF
        finalizer(destroy, da)
    end
    return da
end

"""
    DMDACreate2d(
        ::PetscLib
        comm::MPI.Comm,
        boundary_type_x::DMBoundaryType,
        boundary_type_y::DMBoundaryType,
        stencil_type::DMDAStencilType,
        global_dim_x,
        global_dim_y,
        procs_x,
        procs_y,
        dof_per_node,
        stencil_width,
        points_per_proc_x::Union{Nothing, Vector{PetscInt}};
        points_per_proc_y::Union{Nothing, Vector{PetscInt}};
        options...
    )

Creates a 2-D distributed array with the options specified using keyword
arguments.

# External Links
$(_doc_external("DMDA/DMDACreate2d"))
"""
function DMDACreate2d end

@for_petsc function DMDACreate2d(
    ::$UnionPetscLib,
    comm::MPI.Comm,
    boundary_type_x::DMBoundaryType,
    boundary_type_y::DMBoundaryType,
    stencil_type::DMDAStencilType,
    global_dim_x,
    global_dim_y,
    procs_x,
    procs_y,
    dof_per_node,
    stencil_width,
    points_per_proc_x::Union{Nothing, Vector{$PetscInt}},
    points_per_proc_y::Union{Nothing, Vector{$PetscInt}};
    options...,
)
    opts = Options($petsclib, options...)
    ref_points_per_proc_x = if isnothing(points_per_proc_x)
        C_NULL
    else
        @assert length(points_per_proc_x) == procs_x
        points_per_proc_x
    end
    ref_points_per_proc_y = if isnothing(points_per_proc_y)
        C_NULL
    else
        @assert length(points_per_proc_y) == procs_y
        points_per_proc_y
    end
    da = DM{$PetscLib}(C_NULL, opts)
    with(da.opts) do
        @chk ccall(
            (:DMDACreate2d, $petsc_library),
            PetscErrorCode,
            (
                MPI.MPI_Comm,
                DMBoundaryType,
                DMBoundaryType,
                DMDAStencilType,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                Ptr{$PetscInt},
                Ptr{$PetscInt},
                Ptr{CDM},
            ),
            comm,
            boundary_type_x,
            boundary_type_y,
            stencil_type,
            global_dim_x,
            global_dim_y,
            procs_x,
            procs_y,
            dof_per_node,
            stencil_width,
            ref_points_per_proc_x,
            ref_points_per_proc_y,
            da,
        )
    end
    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if comm == MPI.COMM_SELF
        finalizer(destroy, da)
    end
    return da
end

"""
    DMDACreate2d(
        ::PetscLib
        comm::MPI.Comm,
        boundary_type_x::DMBoundaryType,
        boundary_type_y::DMBoundaryType,
        boundary_type_z::DMBoundaryType,
        stencil_type::DMDAStencilType,
        global_dim_x,
        global_dim_y,
        global_dim_z,
        procs_x,
        procs_y,
        procs_z,
        global_dim_z,
        dof_per_node,
        stencil_width,
        points_per_proc_x::Union{Nothing, Vector{PetscInt}};
        points_per_proc_y::Union{Nothing, Vector{PetscInt}};
        points_per_proc_z::Union{Nothing, Vector{PetscInt}};
        options...
    )

Creates a 3-D distributed array with the options specified using keyword
arguments.

# External Links
$(_doc_external("DMDA/DMDACreate3d"))
"""
function DMDACreate3d end

@for_petsc function DMDACreate3d(
    ::$UnionPetscLib,
    comm::MPI.Comm,
    boundary_type_x::DMBoundaryType,
    boundary_type_y::DMBoundaryType,
    boundary_type_z::DMBoundaryType,
    stencil_type::DMDAStencilType,
    global_dim_x,
    global_dim_y,
    global_dim_z,
    procs_x,
    procs_y,
    procs_z,
    dof_per_node,
    stencil_width,
    points_per_proc_x::Union{Nothing, Vector{$PetscInt}},
    points_per_proc_y::Union{Nothing, Vector{$PetscInt}},
    points_per_proc_z::Union{Nothing, Vector{$PetscInt}};
    options...,
)
    opts = Options($petsclib, options...)
    ref_points_per_proc_x = if isnothing(points_per_proc_x)
        C_NULL
    else
        @assert length(points_per_proc_x) == procs_x
        points_per_proc_x
    end
    ref_points_per_proc_y = if isnothing(points_per_proc_y)
        C_NULL
    else
        @assert length(points_per_proc_y) == procs_y
        points_per_proc_y
    end
    ref_points_per_proc_z = if isnothing(points_per_proc_z)
        C_NULL
    else
        @assert length(points_per_proc_z) == procs_z
        points_per_proc_z
    end
    da = DM{$PetscLib}(C_NULL, opts)
    with(da.opts) do
        @chk ccall(
            (:DMDACreate3d, $petsc_library),
            PetscErrorCode,
            (
                MPI.MPI_Comm,
                DMBoundaryType,
                DMBoundaryType,
                DMBoundaryType,
                DMDAStencilType,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                $PetscInt,
                Ptr{$PetscInt},
                Ptr{$PetscInt},
                Ptr{$PetscInt},
                Ptr{CDM},
            ),
            comm,
            boundary_type_x,
            boundary_type_y,
            boundary_type_z,
            stencil_type,
            global_dim_x,
            global_dim_y,
            global_dim_z,
            procs_x,
            procs_y,
            procs_z,
            dof_per_node,
            stencil_width,
            ref_points_per_proc_x,
            ref_points_per_proc_y,
            ref_points_per_proc_z,
            da,
        )
    end
    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if comm == MPI.COMM_SELF
        finalizer(destroy, da)
    end
    return da
end

"""
    DMDAGetInfo(da::AbstractDM)

Get the info associated with the distributed array `da`.

# External Links
$(_doc_external("DMDA/DMDAGetInfo"))
"""
function DMDAGetInfo end

@for_petsc function DMDAGetInfo(da::AbstractDM{$PetscLib})
    dim = [$PetscInt(0)]
    glo_size = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    procs_per_dim = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    dof_per_node = [$PetscInt(0)]
    stencil_width = [$PetscInt(0)]
    boundary_type = [DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE]
    stencil_type = [DMDA_STENCIL_STAR]
    @chk ccall(
        (:DMDAGetInfo, $petsc_library),
        PetscErrorCode,
        (
            CDM,
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{DMBoundaryType},
            Ref{DMBoundaryType},
            Ref{DMBoundaryType},
            Ref{DMDAStencilType},
        ),
        da,
        dim,
        Ref(glo_size, 1),
        Ref(glo_size, 2),
        Ref(glo_size, 3),
        Ref(procs_per_dim, 1),
        Ref(procs_per_dim, 2),
        Ref(procs_per_dim, 3),
        dof_per_node,
        stencil_width,
        Ref(boundary_type, 1),
        Ref(boundary_type, 2),
        Ref(boundary_type, 3),
        stencil_type,
    )
    return (
        dim = dim[1],
        global_size = glo_size,
        procs_per_dim = procs_per_dim,
        dof_per_node = dof_per_node[1],
        boundary_type = boundary_type,
        stencil_width = stencil_width[1],
        stencil_type = stencil_type[1],
    )
end

"""
    DMDAGetCorners(da::AbstractDM)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(_doc_external("DMDA/DMDAGetCorners"))
"""
function DMDAGetCorners end

@for_petsc function DMDAGetCorners(da::AbstractDM{$PetscLib})
    info = DMDALocalInfo{$PetscInt}()
    corners = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    local_size = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    @chk ccall(
        (:DMDAGetCorners, $petsc_library),
        PetscErrorCode,
        (
            CDM,
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
        ),
        da,
        Ref(corners, 1),
        Ref(corners, 2),
        Ref(corners, 3),
        Ref(local_size, 1),
        Ref(local_size, 2),
        Ref(local_size, 3),
    )
    corners .+= 1
    return (
        lower = corners,
        upper = corners .+ local_size .- $PetscInt(1),
        size = local_size,
    )
end

"""
    DMDAGetGhostCorners(da::AbstractDM)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(_doc_external("DMDA/DMDAGetGhostCorners"))
"""
function DMDAGetGhostCorners end

@for_petsc function DMDAGetGhostCorners(da::AbstractDM{$PetscLib})
    info = DMDALocalInfo{$PetscInt}()
    corners = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    local_size = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    @chk ccall(
        (:DMDAGetGhostCorners, $petsc_library),
        PetscErrorCode,
        (
            CDM,
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
        ),
        da,
        Ref(corners, 1),
        Ref(corners, 2),
        Ref(corners, 3),
        Ref(local_size, 1),
        Ref(local_size, 2),
        Ref(local_size, 3),
    )
    corners .+= 1
    return (
        lower = corners,
        upper = corners .+ local_size .- 1,
        size = local_size,
    )
end
