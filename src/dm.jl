import .LibPETSc: AbstractPetscDM, PetscDM, CDM

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscDM{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc DM (null pointer)")
        return
    else
        ty = LibPETSc.DMGetType(PetscLib,v)
        di = LibPETSc.DMGetDimension(PetscLib,v)
        print(io, "PETSc DM $ty object in $di dimensions")
    end
    return nothing
end


function destroy(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) && dm.ptr != C_NULL 
        LibPETSc.DMDestroy(PetscLib, dm)
    end
    dm.ptr = C_NULL
    return nothing
end



"""
   dagetinfo =  getinfo(dm::AbstractPetscDM{PetscLib})
Gets info about a given DM

Output:
===
- `dagetinfo` - named tuple with the following fields:

- `dim` - dimension of the `DMDA` (1, 2, or 3)
- `global_size` - tuple with global dimensions in each direction of the array
- `mpi_proc_size` - tuple with number of MPI processes in each direction
- `dof` - degrees of freedom per node
- `stencil_width` - width of the stencil
- `boundary_type` - tuple with boundary types in each direction
- `stencil_type` - stencil type, either `DMDA_STENCIL_STAR` or `DMDA_STENCIL_BOX`


"""
function getinfo(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

    dim, M, N, P, m, n, p, dof, s, bx, by, bz, st = LibPETSc.DMDAGetInfo(PetscLib, dm)
    global_size   = (M,N,P)
    mpi_proc_size = (m,n,p)
    boundary_type = (bx,by,bz)
    stencil_width = s
    stencil_type  = st
               

	return (;dim,global_size,mpi_proc_size,dof,s,boundary_type,stencil_width,stencil_type)
end

"""
    lower, upper, size = getcorners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.


Calls `LibPETSc.DMDAGetCorners`.
"""
function getcorners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetCorners(PetscLib, dm)
    corners = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]
    
    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end

"""
    lower, upper, size = getghostcorners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size` of the local part of the domain.

Calls `LibPETSc.DMDAGetCorners`.
"""
function getghostcorners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetGhostCorners(PetscLib, dm)
    corners = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]
    
    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


"""
    setup!(dm::DM)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
setup!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetUp(PetscLib, dm)


"""
    setfromoptions!(dm::AbstractPetscDM)

Sets the global options to the `dm`    
# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
setfromoptions!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetFromOptions(PetscLib, dm)



"""
    empty(da::AbstractPetscDM)

return an uninitialized `DMDA` struct.
"""
Base.empty(da::AbstractPetscDM{PetscLib}) where {PetscLib} = PetscDM{PetscLib}(C_NULL, da.age)
