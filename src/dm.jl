const CDM = Ptr{Cvoid}
abstract type AbstractDM{PetscLib} end

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
function gettype(dm::AbstractDM{PetscLib}) where {PetscLib}
    t_r = Ref{Cstring}()
    LibPETSc.DMGetType(PetscLib, dm, t_r)
    return unsafe_string(t_r[])
end

"""
    getdimension(dm::AbstractDM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
function getdimension(dm::AbstractDM{PetscLib}) where {PetscLib}
    r_dim = Ref{PetscLib.PetscInt}()
    LibPETSc.DMGetDimension(PetscLib, dm, r_dim)
    return r_dim[]
end

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
end
function DMLocalVec(dm::AbstractDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    vec =
        DMLocalVec{PetscLib, PetscScalar, typeof(dm)}(C_NULL, petsclib.age, dm)

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
end
function DMGlobalVec(dm::AbstractDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    vec =
        DMGlobalVec{PetscLib, PetscScalar, typeof(dm)}(C_NULL, petsclib.age, dm)

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

#=
#
# OLD WRAPPERS
#
"""
    view(dm::AbstractDM, viewer::Viewer=ViewerStdout(petsclib, getcomm(dm)))

view a `dm` with `viewer`

# External Links
$(_doc_external("DM/DMView"))
"""
function view(::AbstractDM) end

@for_petsc function view(
    dm::AbstractDM{$PetscLib},
    viewer::AbstractViewer{$PetscLib} = ViewerStdout($petsclib, getcomm(dm)),
)
    @chk ccall(
        (:DMView, $petsc_library),
        PetscErrorCode,
        (CDM, CPetscViewer),
        dm,
        viewer,
    )
    return nothing
end

"""
    getcoordinateDM(dm::AbstractDM)

Create a `coord_dm` for the coordinates of `dm`.

# External Links
$(_doc_external("DM/DMGetCoordinateDM"))
"""
function getcoordinateDM end

@for_petsc function getcoordinateDM(dm::AbstractDM{$PetscLib})
    coord_dm = empty(dm)
    @chk ccall(
        (:DMGetCoordinateDM, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CDM}),
        dm,
        coord_dm,
    )

    # If this fails then the `empty` call above is probably a bad idea!
    if gettype(coord_dm) != "product"
        @assert gettype(dm) == gettype(coord_dm)
    else
        @assert gettype(dm) == "stag"   # product can only be used with stag
    end

    return coord_dm
end

"""
    getcoordinateslocal(dm::AbstractDM)

Gets a local vector with the coordinates associated with `dm`.

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function getcoordinateslocal end

@for_petsc function getcoordinateslocal(dm::AbstractDM{$PetscLib})
    coord_vec = DMLocalVec(C_NULL, dm)
    @chk ccall(
        (:DMGetCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CVec}),
        dm,
        coord_vec,
    )

    return coord_vec
end
=#
