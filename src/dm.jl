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

#=
#
# OLD WRAPPERS
#
"""
    DMLocalVec(v::CVec, dm::AbstractDM)

Container for an PETSc vector we know is "local"

# External Links
$(_doc_external("Vec/Vec"))
"""
mutable struct DMLocalVec{PetscLib, T, T_DM} <: AbstractVec{T}
    ptr::CVec
    dm::T_DM
    function DMLocalVec(ptr, dm::AbstractDM{PetscLib}) where {PetscLib}
        new{PetscLib, scalartype(PetscLib), typeof(dm)}(ptr, dm)
    end
end

"""
    DMGlobalVec(v::CVec, dm::AbstractDM)

Container for an PETSc vector we know is "global"

# External Links
$(_doc_external("Vec/Vec"))
"""
mutable struct DMGlobalVec{PetscLib, T, T_DM} <: AbstractVec{T}
    ptr::CVec
    dm::T_DM
    function DMGlobalVec(ptr, dm::AbstractDM{PetscLib}) where {PetscLib}
        new{PetscLib, scalartype(PetscLib), typeof(dm)}(ptr, dm)
    end
end

# Mainly for DM we do not know the type of, namely ones returned by PETSc
# functions such as `KSPGetDM`
mutable struct PetscDM{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
end

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
    creatematrix(dm::AbstractDM)

Generates a matrix from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateMatrix"))
"""
function creatematrix end

@for_petsc function creatematrix(dm::AbstractDM{$PetscLib})
    mat = Mat{$PetscScalar}(C_NULL)

    @chk ccall(
        (:DMCreateMatrix, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CMat}),
        dm,
        mat,
    )

    return mat
end

"""
    createlocalvector(dm::AbstractDM)

returns a local vector from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateLocalVector"))
"""
function createlocalvector end

@for_petsc function createlocalvector(dm::AbstractDM{$PetscLib})
    vec = DMLocalVec(C_NULL, dm)

    @chk ccall(
        (:DMCreateLocalVector, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CVec}),
        dm,
        vec,
    )

    return vec
end

"""
    createglobalvector(dm::DM; write::Bool = true, read::Bool = true)

returns a global vector from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateGlobalVector"))
"""
function createglobalvector end

@for_petsc function createglobalvector(dm::AbstractDM{$PetscLib})
    vec = DMGlobalVec(C_NULL, dm)

    @chk ccall(
        (:DMCreateGlobalVector, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CVec}),
        dm,
        vec,
    )

    return vec
end

"""
    update!(
        global_vec::DMGlobalVec,
        local_vec::DMLocalVec,
        mode::InsertMode,
    )

Updates `global_vec` from `local_vec` with insert `mode`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function update!( global_vec::DMGlobalVec, local_vec::DMLocalVec, mode::InsertMode) end

@for_petsc function update!(
    global_vec::DMGlobalVec{$PetscLib},
    local_vec::DMLocalVec{$PetscLib},
    mode::InsertMode
)
    @assert local_vec.dm === global_vec.dm

    update_local2global!( global_vec.ptr, local_vec.ptr, mode, local_vec.dm);

    return nothing
end

"""
    update!(
        global_vec::DMGlobalVec,
        local_ptr::CVec,
        mode::InsertMode,
    )

Updates `global_vec` from a pointer to a local vector `local_ptr` with insert `mode`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function update!( global_vec::DMGlobalVec, local_ptr::CVec, mode::InsertMode) end

@for_petsc function update!(
    global_vec::DMGlobalVec{$PetscLib},
    local_ptr::CVec,
    mode::InsertMode
)
    update_local2global!( global_vec.ptr, local_ptr, mode, global_vec.dm);

    return nothing
end

"""
    update!(
        global_ptr::CVec,
        local_vec::DMLocalVec,
        mode::InsertMode,
    )

Updates pointer to global vec `global_ptr` from `local_vec` with insert `mode`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function update!( global_vec::DMGlobalVec, local_ptr::CVec, mode::InsertMode) end

@for_petsc function update!(
    global_ptr::CVec,
    local_vec::DMLocalVec{$PetscLib},
    mode::InsertMode
)
    update_local2global!( global_ptr, local_vec.ptr, mode, local_vec.dm);

    return nothing
end

"""
    update!(
        local_vec::DMLocalVec,
        global_ptr::CVec,
        mode::InsertMode,
    )

Updates `local_vec` from pointer to global vec `global_ptr` with insert `mode`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
function update!(local_vec::DMLocalVec,global_ptr::CVec, mode::InsertMode) end

@for_petsc function update!(
    local_vec::DMLocalVec{$PetscLib},
    global_ptr::CVec,
    mode::InsertMode,
)

    update_global2local!( local_vec.ptr, global_ptr, mode, local_vec.dm);
    return nothing
end

"""
    update!(
        local_ptr::CVec,
        global_vec::DMGlobalVec
        mode::InsertMode,
    )

Updates pointer to local vec `local_ptr` from `global_vec` with insert `mode`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
function update!(local_ptr::CVec, global_vec::DMGlobalVec, mode::InsertMode) end

@for_petsc function update!(
    local_ptr::CVec,
    global_vec::DMGlobalVec{$PetscLib},
    mode::InsertMode,
)

    update_global2local!( local_ptr, global_vec.ptr, mode, global_vec.dm);
    return nothing
end

"""
    update!(
        local_vec::DMLocalVec,
        global_vec::DMGlobalVec,
        mode::InsertMode
    )

Updates `local_vec` from a pointer to a global vec `global_ptr` with insert `mode`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
function update!(local_vec::DMLocalVec, global_vec::DMGlobalVec, mode::InsertMode) end

@for_petsc function update!(
    local_vec::DMLocalVec{$PetscLib},
    global_vec::DMGlobalVec{$PetscLib},
    mode::InsertMode
)
    @assert local_vec.dm === global_vec.dm

    update_global2local!( local_vec.ptr, global_vec.ptr, mode, local_vec.dm);

    return nothing
end

#=
@for_petsc function update!(
    local_ptr::CVec,
    global_vec::DMGlobalVec{$PetscLib},
    mode::InsertMode,
)
    @chk ccall(
        (:DMGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        global_vec.dm,
        global_vec,
        mode,
        local_ptr,
    )

    return nothing
end
=#

"""
    update_local2global!(
        global_ptr::CVec,
        local_ptr::CVec,
        mode::InsertMode,
        dm::AbstractDM
    )

Updates pointer of `global_vec` from pointer of `local_vec` with insert `mode`. 
Both vectors should belong to the same `dm`

This is a low-level routine that is typically called by `update!`
"""
function update_local2global! end

@for_petsc function update_local2global!(
    global_ptr::CVec,
    local_ptr::CVec,
    mode::InsertMode,
    dm::AbstractDM{$PetscLib}
)
    @chk ccall(
        (:DMLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        dm,
        local_ptr,
        mode,
        global_ptr,
    )

    return nothing
end

"""
    update_global2local!(
        local_ptr::CVec,
        global_ptr::CVec,
        mode::InsertMode,
        dm::AbstractDM
    )

Updates pointer to `local` vector from pointer to `global` vector with insert `mode`, assuming that both belong to the same `dm`

This is a low-level routine that is typically called by `update!`

"""
function update_global2local! end

@for_petsc function update_global2local!(
    local_ptr::CVec,
    global_ptr::CVec,
    mode::InsertMode,
    dm::AbstractDM{$PetscLib}
)
    @chk ccall(
        (:DMGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        dm,
        global_ptr,
        mode,
        local_ptr,
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
