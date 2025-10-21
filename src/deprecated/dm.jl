const CDM = Ptr{Cvoid}
abstract type AbstractDM{PetscLib} end

import PETSc.LibPETSc: DMType, PetscViewer, ISColoringType, DMBlockingType, MPI_Datatype
import PETSc.LibPETSc: PetscObject, DMLabel, PetscMPIInt, PetscSection, DMField, PetscFE
import PETSc.LibPETSc: DMPointLocationType, VecScatter, VecType, MatType, MatFDColoring, DMReorderDefaultFlag,
        MatOrderingType,PetscSF,PetscDS,PetscCopyMode,DMCopyLabelsMode,DMPolytopeType, IS


mutable struct DM{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    age::Int
end

mutable struct DMPtr{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
    age::Int
    own::Bool
end

include("./dm_wrapped.jl")

function destroy(dm::AbstractDM{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) &&
       dm.age == getlib(PetscLib).age &&
       dm.ptr != C_NULL &&
       (!hasfield(typeof(dm), :own) || dm.own)
        LibPETSc.DMDestroy(PetscLib, dm)
    end
    dm.ptr = C_NULL
    return nothing
end



"""
    setfromoptions!(dm::DM, opts=dm.opts)

# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
function setfromoptions!(
    dm::AbstractDM{PetscLib},
    opts::Options = dm.opts,
) where {PetscLib}
    with(opts) do
        LibPETSc.DMSetFromOptions(PetscLib, dm)
    end
end

"""
    setup!(dm::DM, opts=dm.opts)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
function setup!(
    dm::AbstractDM{PetscLib},
    opts::Options = dm.opts,
) where {PetscLib}
    with(opts) do
        LibPETSc.DMSetUp(PetscLib, dm)
    end
end

"""
    gettype(dm::AbstractDM)

Gets type name of the `dm`

# External Links
$(_doc_external("DM/DMGetType"))
"""
gettype(dm::AbstractDM{PetscLib}) where PetscLib = DMGetType(dm)

#function gettype(dm::AbstractDM{PetscLib}) where {PetscLib}
#    t_r = Ref{PETSc.DMType}()
#    LibPETSc.DMGetType(PetscLib, dm, t_r)
#    return unsafe_string(t_r[])
#end

"""
    getdimension(dm::AbstractDM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
getdimension(dm::AbstractDM{PetscLib}) where PetscLib = DMGetDimension(dm)


#function getdimension(dm::AbstractDM{PetscLib}) where {PetscLib}
#    r_dim = Ref{PetscLib.PetscInt}()
#    LibPETSc.DMGetDimension(PetscLib, dm, r_dim)
#    return r_dim[]
#end

"""
    MatAIJ(dm::AbstractDM)

Generates a matrix from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateMatrix"))
"""
function MatAIJ(dm::AbstractDM{PetscLib}) where {PetscLib}
    mat = MatAIJ{PetscLib, PetscLib.PetscScalar}(C_NULL, getlib(PetscLib).age)

    LibPETSc.DMCreateMatrix(PetscLib, dm, mat)

    if MPI.Comm_size(getcomm(mat)) == 1
        finalizer(destroy, mat)
    end

    return mat
end

"""
    DMLocalVec(dm::AbstractDM)

returns a local vector from the `dm` object (has space for ghost).

# External Links
$(_doc_external("DM/DMCreateLocalVector"))
"""
mutable struct DMLocalVec{PetscLib, PetscScalar, DMT <: AbstractDM} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
    age::Int
    dm::DMT
    own::Bool
end
function DMLocalVec(dm::AbstractDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    vec = DMLocalVec{PetscLib, PetscScalar, typeof(dm)}(
        C_NULL,
        petsclib.age,
        dm,
        true,
    )

    LibPETSc.DMCreateLocalVector(PetscLib, dm, vec)

    if MPI.Comm_size(getcomm(vec)) == 1
        finalizer(destroy, vec)
    end

    return vec
end

"""
    DMGlobalVec(dm::AbstractDM)

returns a global vector from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateGlobalVector"))
"""
mutable struct DMGlobalVec{PetscLib, PetscScalar, DMT <: AbstractDM} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
    age::Int
    dm::DMT
    own::Bool
end
function DMGlobalVec(dm::AbstractDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    vec = DMGlobalVec{PetscLib, PetscScalar, typeof(dm)}(
        C_NULL,
        petsclib.age,
        dm,
        true,
    )

    LibPETSc.DMCreateGlobalVector(PetscLib, dm, vec)

    if MPI.Comm_size(getcomm(vec)) == 1
        finalizer(destroy, vec)
    end

    return vec
end

"""
    update!(
        local_vec::DMLocalVec,
        global_vec::AbstractVec,
        mode::InsertMode,
    )

Updates `local_vec` from the `global_vec` with insert `mode`.

Communication and computation can be overlapped with [`updatebegin!`](@ref) and
[`updateend!`](@ref)

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
$(_doc_external("DM/DMGlobalToLocalBegin"))
$(_doc_external("DM/DMGlobalToLocalEnd"))
"""
function update!(
    local_vec::DMLocalVec{PetscLib},
    global_vec::AbstractVec,
    mode::InsertMode,
) where {PetscLib}
    updatebegin!(local_vec, global_vec, mode)
    return updateend!(local_vec, global_vec, mode)
end

"""
    updatebegin!(
        local_vec::DMLocalVec,
        global_vec::AbstractVec,
        mode::InsertMode,
    )

Begin update of `local_vec` from the `global_vec` with insert `mode`.

# External Links
$(_doc_external("DM/DMGlobalToLocalBegin"))
"""
function updatebegin!(
    local_vec::DMLocalVec{PetscLib},
    global_vec::AbstractVec,
    mode::InsertMode,
) where {PetscLib}
    LibPETSc.DMGlobalToLocalBegin(
        PetscLib,
        local_vec.dm,
        global_vec,
        mode,
        local_vec,
    )

    return nothing
end

"""
    updateend!(
        local_vec::DMLocalVec,
        global_vec::AbstractVec,
        mode::InsertMode,
    )

End update of `local_vec` from the `global_vec` with insert `mode`.

# External Links
$(_doc_external("DM/DMGlobalToLocalEnd"))
"""
function updateend!(
    local_vec::DMLocalVec{PetscLib},
    global_vec::AbstractVec,
    mode::InsertMode,
) where {PetscLib}
    LibPETSc.DMGlobalToLocalEnd(
        PetscLib,
        local_vec.dm,
        global_vec,
        mode,
        local_vec,
    )

    return local_vec
end

"""
    update!(
        global_vec::AbstractVec,
        local_vec::DMLocalVec,
        mode::InsertMode,
    )

Updates `global_vec` from the `local_vec` with insert `mode`.

Communication and computation can be overlapped with [`updatebegin!`](@ref) and
[`updateend!`](@ref)

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
$(_doc_external("DM/DMLocalToGlobalBegin"))
$(_doc_external("DM/DMLocalToGlobalEnd"))
"""
function update!(
    global_vec::AbstractVec,
    local_vec::DMLocalVec{PetscLib},
    mode::InsertMode,
) where {PetscLib}
    updatebegin!(global_vec, local_vec, mode)
    return updateend!(global_vec, local_vec, mode)
end

"""
    updatebegin!(
        global_vec::AbstractVec,
        local_vec::DMLocalVec,
        mode::InsertMode,
    )

Begin update of `global_vec` from the `local_vec` with insert `mode`.

# External Links
$(_doc_external("DM/DMLocalToGlobalBegin"))
"""
function updatebegin!(
    global_vec::AbstractVec,
    local_vec::DMLocalVec{PetscLib},
    mode::InsertMode,
) where {PetscLib}
    LibPETSc.DMLocalToGlobalBegin(
        PetscLib,
        local_vec.dm,
        local_vec,
        mode,
        global_vec,
    )

    return nothing
end

"""
    updateend!(
        global_vec::AbstractVec,
        local_vec::DMLocalVec,
        mode::InsertMode,
    )

End update of `global_vec` from the `local_vec` with insert `mode`.

# External Links
$(_doc_external("DM/DMLocalToGlobalEnd"))
"""
function updateend!(
    global_vec::AbstractVec,
    local_vec::DMLocalVec{PetscLib},
    mode::InsertMode,
) where {PetscLib}
    LibPETSc.DMLocalToGlobalEnd(
        PetscLib,
        local_vec.dm,
        local_vec,
        mode,
        global_vec,
    )

    return local_vec
end

"""
    getcoordinateDM(dm::AbstractDM)

Create a `coord_dm` for the coordinates of `dm`.

# External Links
$(_doc_external("DM/DMGetCoordinateDM"))
"""
function getcoordinateDM(dm::AbstractDM{PetscLib}) where {PetscLib}
    coord_dm = empty(dm)
    LibPETSc.DMGetCoordinateDM(PetscLib, dm, coord_dm)

    if gettype(dm) == "stag"
        @assert gettype(coord_dm) == "product"
    else
        @assert gettype(coord_dm) == gettype(dm)
    end

    if MPI.Comm_size(getcomm(coord_dm)) == 1
        finalizer(destroy, coord_dm)
    end

    return coord_dm
end

"""
    coordinatesDMLocalVec(dm::AbstractDM)

Gets a local vector with the coordinates associated with `dm`.

Note that the returned vector is borrowed from the `dm` and is not a new vector.

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function coordinatesDMLocalVec(dm::AbstractDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    coord_vec = DMLocalVec{PetscLib, PetscScalar, typeof(dm)}(
        C_NULL,
        petsclib.age,
        dm,
        false,
    )
    LibPETSc.DMGetCoordinatesLocal(PetscLib, dm, coord_vec)

    return coord_vec
end

"""
    view(dm::AbstractDM, [viewer])

view a `dm` with `viewer`

# External Links
$(_doc_external("DM/DMView"))
"""
function view(
    dm::AbstractDM{PetscLib},
    viewer = LibPETSc.PETSC_VIEWER_STDOUT_(PetscLib, getcomm(dm)),
) where {PetscLib}
    LibPETSc.DMView(PetscLib, dm, viewer)
    return nothing
end
