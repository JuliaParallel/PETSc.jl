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