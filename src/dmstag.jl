abstract type AbstractDMStag{PetscLib} <: AbstractDM{PetscLib} end

mutable struct DMStag{PetscLib} <: AbstractDMStag{PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    age::Int
end

mutable struct DMStagPtr{PetscLib} <: AbstractDMStag{PetscLib}
    ptr::CDM
    age::Snt
    own::Bool
end

import PETSc.LibPETSc:  DMStagStencilLocation, DMType, DMStagStencil


include("dmstag_wrapped.jl")    

"""
    empty(dm::DMStag)

return an uninitialized `DMStag` struct.
"""
Base.empty(dm::DMStag{PetscLib}) where {PetscLib} =
    DMStag{PetscLib}(C_NULL, dm.opts, dm.age)

#abstract type DMStagStencilLocation{PetscLib} end

"""
    DMStag(
        petsclib::PetscLib
        comm::MPI.Comm,
        boundary_type::NTuple{D, DMBoundaryType},
        global_dim::NTuple{D, Integer},
        dof_per_node::NTuple{1 + D, Integer},
        stencil_width::Integer,
        stencil_type = DMSTAG_STENCIL_BOX;
        points_per_proc::Tuple,
        processors::Tuple,
        setfromoptions = true,
        dmsetup = true,
        options...
    )

Creates a D-dimensional distributed staggered array with the options specified
using keyword arguments.

The Tuple `dof_per_node` specifies how many degrees of freedom are at all the
staggerings in the order:
 - 1D: `(vertex, element)`
 - 2D: `(vertex, edge, element)`
 - 3D: `(vertex, edge, face, element)`

If keyword argument `points_per_proc[k] isa Vector{petsclib.PetscInt}` then this
specifies the points per processor in dimension `k`.

If keyword argument `processors[k] isa Integer` then this specifies the number of
processors used in dimension `k`; ignored when `D == 1`.

If keyword argument `setfromoptions == true` then `setfromoptions!` called.

If keyword argument `dmsetup == true` then `setup!` is called.

When `D == 1` the `stencil_type` argument is not required and ignored if
specified.

# External Links
$(_doc_external("DMStag/DMStagCreate1d"))
$(_doc_external("DMStag/DMStagCreate2d"))
$(_doc_external("DMStag/DMStagCreate3d"))
"""
function DMStag(
    petsclib::PetscLib,
    comm::MPI.Comm,
    boundary_type::NTuple{1, DMBoundaryType},
    global_dim::NTuple{1, Integer},
    dof_per_node::NTuple{2, Integer},
    stencil_width::Integer,
    stencil_type = DMSTAG_STENCIL_BOX;
    points_per_proc::Tuple = (nothing,),
    processors = nothing,
    dmsetfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib}
    opts = Options(petsclib; options...)
    dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

    ref_points_per_proc =
        if isnothing(points_per_proc[1]) || points_per_proc[1] == PETSC_DECIDE
            C_NULL
        else
           # @assert points_per_proc[1] isa Array{PetscLib.PetscInt}
            @assert length(points_per_proc[1]) == MPI.Comm_size(comm)
            points_per_proc[1]
        end

    with(dm.opts) do
        LibPETSc.DMStagCreate1d(
            PetscLib,
            comm,
            boundary_type[1],
            global_dim[1],
            dof_per_node[1],
            dof_per_node[2],
            stencil_type,
            stencil_width,
            ref_points_per_proc,
            dm,
        )
    end

    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dm)
    end

    return dm
end


function DMStag(
    petsclib::PetscLib,
    comm::MPI.Comm,
    boundary_type::NTuple{2, DMBoundaryType},
    global_dim::NTuple{2, Integer},
    dof_per_node::NTuple{3, Integer},
    stencil_width::Integer,
    stencil_type = DMSTAG_STENCIL_BOX;
    points_per_proc::Tuple = (nothing, nothing),
    processors::Tuple = (PETSC_DECIDE, PETSC_DECIDE),
    dmsetfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib}
    opts = Options(petsclib; options...)
    dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

    ref_points_per_proc = ntuple(2) do d
        if isnothing(points_per_proc[d]) || points_per_proc[d] == PETSC_DECIDE
            C_NULL
        else
            @assert points_per_proc[d] isa Array{PetscLib.PetscInt}
            @assert length(points_per_proc[d]) == MPI.Comm_size(comm)
            points_per_proc[d]
        end
    end

    with(dm.opts) do
        LibPETSc.DMStagCreate2d(
            PetscLib,
            comm,
            boundary_type[1],
            boundary_type[2],
            global_dim[1],
            global_dim[2],
            processors[1],
            processors[2],
            dof_per_node[1],
            dof_per_node[2],
            dof_per_node[3],
            stencil_type,
            stencil_width,
            ref_points_per_proc[1],
            ref_points_per_proc[2],
            dm,
        )
    end

    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dm)
    end

    return dm
end

function DMStag(
    petsclib::PetscLib,
    comm::MPI.Comm,
    boundary_type::NTuple{3, DMBoundaryType},
    global_dim::NTuple{3, Integer},
    dof_per_node::NTuple{4, Integer},
    stencil_width::Integer,
    stencil_type = DMSTAG_STENCIL_BOX;
    points_per_proc::Tuple = (nothing, nothing, nothing),
    processors::Tuple = (PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE),
    dmsetfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib}
    opts = Options(petsclib; options...)
    dm = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

    ref_points_per_proc = ntuple(3) do d
        if isnothing(points_per_proc[d]) || points_per_proc[d] == PETSC_DECIDE
            C_NULL
        else
            @assert points_per_proc[d] isa Array{PetscLib.PetscInt}
            @assert length(points_per_proc[d]) == MPI.Comm_size(comm)
            points_per_proc[d]
        end
    end

    with(dm.opts) do
        LibPETSc.DMStagCreate3d(
            PetscLib,
            comm,
            boundary_type[1],
            boundary_type[2],
            boundary_type[3],
            global_dim[1],
            global_dim[2],
            global_dim[3],
            processors[1],
            processors[2],
            processors[3],
            dof_per_node[1],
            dof_per_node[2],
            dof_per_node[3],
            dof_per_node[4],
            stencil_type,
            stencil_width,
            ref_points_per_proc[1],
            ref_points_per_proc[2],
            ref_points_per_proc[3],
            dm,
        )
    end

    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dm)
    end

    return dm
end

function DMStag(
    dm::AbstractDMStag{PetscLib},
    dof_per_node::Union{NTuple{2,Int},NTuple{3,Int},NTuple{4,Int}},
    dmsetfromoptions = true,
    dmsetup = true,
    options...,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    opts = Options(petsclib; options...)
    dmnew = DMStag{PetscLib}(C_NULL, opts, petsclib.age)

    s = size(dof_per_node,1)

    dof_per_node_C = [0,0,0,0]

    for (i, value) in enumerate(dof_per_node)
        dof_per_node_C[i] = value
    end


    with(dm.opts) do
        LibPETSc.DMStagCreateCompatibleDMStag(
            PetscLib,
            dm,
            dof_per_node_C[1],
            dof_per_node_C[2],
            dof_per_node_C[3],
            dof_per_node_C[4],
            dmnew,
        )
    end

    dmsetfromoptions && setfromoptions!(dmnew)
    dmsetup && setup!(dmnew)

    comm  = getcomm(dm);

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dmnew)
    end    

    return dmnew
end


Base.size(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm)
globalsize(dm::AbstractDMStag) = DMStagGetGlobalSizes(dm::AbstractDMStag)
boundarytypes(dm::AbstractDMStag)  = DMStagGetBoundaryTypes(dm::AbstractDMStag) 

"""
    getcorners(dm::AbstractDMStag)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. Also included is `nextra` of
the number of extra partial elements in each direction.

# External Links
$(_doc_external("DMDA/DMStagGetCorners"))
"""
function getcorners(dm::AbstractDMStag{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    
    x,y,z,m,n,p,nExtrax,nExtray,nExtraz = DMStagGetCorners(dm)

    corners = [x,y,z]
    local_size = [m,n,p]
    nextra = [nExtrax,nExtray,nExtraz]
 

    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)
    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
        nextra = (nextra...,),
    )
end


"""
    getghostcorners(dm::AbstractDMStag)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(_doc_external("DMDA/DMDAGetGhostCorners"))
"""
function getghostcorners(dm::AbstractDMStag{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt

    x,y,z,m,n,p = DMStagGetGhostCorners(dm)

    corners = [x,y,z]
    local_size = [m,n,p]
    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)
    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


function setuniformcoordinates!(
    dm::AbstractDMStag{PetscLib},
    xyzmin::Union{NTuple{1,Int},NTuple{2,Int},NTuple{3,Int}},
    xyzmax::Union{NTuple{1,Int},NTuple{2,Int},NTuple{3,Int}},
    ) where {PetscLib}
    PetscInt = PetscLib.PetscInt

    xmin = xyzmin[1]
    xmax = xyzmax[1]

    s = size(xyzmin,1)

    ymin = (s > 1) ? xyzmin[2] : PetscInt(0)
    ymax = (s > 1) ? xyzmax[2] : PetscInt(0)

    zmin = (s > 2) ? xyzmin[3] : PetscInt(0)
    zmax = (s > 2) ? xyzmax[3] : PetscInt(0)

    LibPETSc.DMStagSetUniformCoordinatesProduct(
        PetscLib,
        dm,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    )

    return nothing
end

function getcoordinatearray(dm::AbstractDMStag{PetscLib}) where {PetscLib}
    PetscScalar = PetscLib.PetscScalar

    xP = Ref{Ptr{Ptr{PetscScalar}}}(C_NULL)
    yP = Ref{Ptr{Ptr{PetscScalar}}}(C_NULL)
    zP = Ref{Ptr{Ptr{PetscScalar}}}(C_NULL)

    LibPETSc.DMStagGetProductCoordinateArrays(
        PetscLib,
        dm,
        xP,
        yP,
        zP
    )

    corners = getghostcorners(dm)
    mstart  = corners.lower
    s       = corners.size
    
    xP = Base.unsafe_load(xP[],mstart[1])
    xArray = unsafe_wrap(Array, xP, (2,s[1]); own = false)
    xArrayO = OffsetArray(xArray,CartesianIndex(1,mstart[1]):CartesianIndex(2,mstart[1]+s[1]-1))
    if yP[] != (C_NULL)
        yP = Base.unsafe_load(yP[],mstart[2])
        yArray = unsafe_wrap(Array, yP, (2,s[2]); own = false)
        yArrayO = OffsetArray(yArray,CartesianIndex(1,mstart[2]):CartesianIndex(2,mstart[2]+s[2]-1))
    else 
        yArrayO = nothing
    end
    if zP[] != (C_NULL)
        zP = Base.unsafe_load(zP[],mstart[3])
        zArray = unsafe_wrap(Array, zP, (2,s[3]); own = false)
        zArrayO = OffsetArray(zArray,CartesianIndex(1,mstart[3]):CartesianIndex(2,mstart[3]+s[3]-1))
    else
        zArrayO = nothing
    end

    return xArrayO,yArrayO,zArrayO

end


"""
    getentriesperelement(dm::AbstractDMStag)

"""
function getentriesperelement(dm::AbstractDMStag{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    entriesPerElement = [PetscInt(1)]

    LibPETSc.DMStagGetEntriesPerElement(
        PetscLib,
        dm,
        Ref(entriesPerElement,1)
    )

    return entriesPerElement[1]
end


"""
    vecgetarray(dm::AbstractDMStag{PetscLib}, v::AbstractVec{PetscLib})
"""
function vecgetarray(dm::AbstractDMStag{PetscLib}, v::AbstractVec{PetscLib}) where {PetscLib}
    # Note: there is actually no need to call PETSc again, as Julia has the possibility 
    # to wrap an existing array into another one. Our vec already has the array wrapper, 
    # so we reshape that 

    # Extract array from vector. Note: we need to release this by calling 
    # Base.finalize on X1!
    X1       =   unsafe_localarray(v;  write=true, read=true)

    array          =   vecgetarray(dm, X1) 
        
    return array
end


"""
    vecgetarray(dm::AbstractDMStag{PetscLib}, v::Vector)
"""
function vecgetarray(dm::AbstractDMStag{PetscLib}, v::Vector) where {PetscLib}

    entriesPerElement   =   getentriesperelement(dm)
    ghost_corners       =   getghostcorners(dm)
    dim                 =   getdimension(dm) 

    # Dimensions of new array (see the PETSc DMStagVecGetArrayRead routine)
    dim_vec             =   [entriesPerElement; collect(ghost_corners.size[1:dim])]

    # Wrap julia vector to new vector.
    X                   =    Base.view(v,:)
     
    # reshape to correct format
    X                   =   reshape(v, Tuple(dim_vec))
    array               =   PermutedDimsArray(X, Tuple([2:dim+1;1]))   # permute to take care of different array ordering in C/Julia
       
    return array
end  

function getlocationslot(
    dm::AbstractDMStag{PetscLib},
    loc::LibPETSc.DMStagStencilLocation, 
    dof::Integer
    ) where {PetscLib}
    PetscInt     = inttype(PetscLib)
    slot = [PetscInt(1)]

    LibPETSc.DMStagGetLocationSlot(
        PetscLib,
        dm,
        loc,
        dof,
        Ref(slot,1)
    )

    return slot[]
end

DMStagGetLocationSlot(dm, loc, dof)   =  getlocationslot(dm, loc, dof)

"""
    Indices = DMStagGetIndices(dm::DMStag)
    
Return indices of start and end of the central/vertex nodes of a local array built from the input `dm`. 
This takes ghost points into account and helps    

    dm 	        - the DMStag object
    Indices 	- indices of lower and upper range of center and vertex nodes
"""
function DMStagGetIndices end

function DMStagGetIndices(dm::DMStag)
    # In Julia, indices in arrays start @ 1, whereas they can go negative in C
    gc              =   PETSc.getghostcorners(dm);  
    c               =   PETSc.getcorners(dm); 

    # If we have ghosted boundaries, we need to shift the start/end points, as ghosted 
    # boundaries are treated in PETSc with negative numbers, whereas in Julia everything is 1-based

    # NOTE: we have not yet tested this in parallel
    Diff            =   c.lower - gc.lower;
    Start           =   c.lower + Diff;
    End             =   Start + CartesianIndex(c.size) -  CartesianIndex(1,1,1) ;

    # Note that we add the shift for julia/petsc consistency
    shift = 0;
    center = (  x= Start[1]:End[1],
                y= Start[2]:End[2],  
                z= Start[3]:End[3] )

    vertex = (  x= Start[1]:End[1]+1 ,
                y= Start[2]:End[2]+1 ,  
                z= Start[3]:End[3]+1 )

    return (center=center, vertex=vertex)
            
end

