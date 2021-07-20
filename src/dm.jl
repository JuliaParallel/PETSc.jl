const CDM = Ptr{Cvoid}
abstract type AbstractDM{PetscLib} end

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
    setup!(da::DM, opts=da.opts)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
function setup! end

@for_petsc function setup!(da::AbstractDM{$PetscLib}, opts::Options = da.opts)
    with(opts) do
        @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (CDM,), da)
    end
end

"""
    setfromoptions!(da::DM, opts=da.opts)

# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
function setfromoptions! end

@for_petsc function setfromoptions!(
    da::AbstractDM{$PetscLib},
    opts::Options = da.opts,
)
    with(opts) do
        @chk ccall(
            (:DMSetFromOptions, $petsc_library),
            PetscErrorCode,
            (CDM,),
            da,
        )
    end
end

@for_petsc begin
    function destroy(da::AbstractDM{$PetscLib})
        finalized($PetscLib) || @chk ccall(
            (:DMDestroy, $petsc_library),
            PetscErrorCode,
            (Ptr{CDM},),
            da,
        )
        da.ptr = C_NULL
        return nothing
    end
end

"""
    getdimension(dm::AbstractDM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
function getdimension(::AbstractDM) end

@for_petsc function getdimension(dm::AbstractDM{$PetscLib})
    dim = Ref{$PetscInt}()

    @chk ccall(
        (:DMGetDimension, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{$PetscInt}),
        dm,
        dim,
    )

    return dim[]
end

"""
    gettype(dm::AbstractDM)

Gets type name of the `dm`

# External Links
$(_doc_external("DM/DMGetType"))
"""
function gettype(::AbstractDM) end

@for_petsc function gettype(dm::AbstractDM{$PetscLib})
    t_r = Ref{Cstring}()
    @chk ccall(
        (:DMGetType, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{Cstring}),
        dm,
        t_r,
    )
    return unsafe_string(t_r[])
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
update!(::DMGlobalVec, ::DMLocalVec, ::InsertMode)

@for_petsc function update!(
    global_vec::DMGlobalVec{$PetscLib},
    local_vec::DMLocalVec{$PetscLib},
    mode::InsertMode,
)
    @assert local_vec.dm === global_vec.dm
    @chk ccall(
        (:DMLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        local_vec.dm,
        local_vec,
        mode,
        global_vec,
    )

    return nothing
end

"""
    update!(
        local_vec::DMLocalVec,
        global_vec::DMGlobalVec,
        mode::InsertMode,
    )

Updates `local_vec` from `global_vec` with insert `mode`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
update!(::DMLocalVec, ::DMGlobalVec, ::InsertMode)

@for_petsc function update!(
    local_vec::DMLocalVec{$PetscLib},
    global_vec::DMGlobalVec{$PetscLib},
    mode::InsertMode,
)
    @assert local_vec.dm === global_vec.dm
    @chk ccall(
        (:DMGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        global_vec.dm,
        global_vec,
        mode,
        local_vec,
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
    @assert gettype(dm) == gettype(coord_dm)

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
