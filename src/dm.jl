const CDM = Ptr{Cvoid}

abstract type AbstractDM{PetscLib} end
mutable struct DM{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    DM{PetscLib}(ptr, opts = Options(PetscLib)) where {PetscLib} =
        new{PetscLib}(ptr, opts)
end

"""
    DMSetUp!(da::DM, opts=da.opts)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
function DMSetUp! end

@for_petsc function DMSetUp!(da::AbstractDM{$PetscLib}, opts::Options = da.opts)
    with(opts) do
        @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (CDM,), da)
    end
end

"""
    DMSetFromOptions!(da::DM, opts=da.opts)

# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
function DMSetFromOptions! end

@for_petsc function DMSetFromOptions!(
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
    function destroy(da::DM{$PetscLib})
        finalized($PetscScalar) || @chk ccall(
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
    DMGetDimension(dm::DM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
function DMGetDimension end

@for_petsc function DMGetDimension(dm::DM{$PetscLib})
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
    gettype(dm::DM)

Gets type name of the `dm`

# External Links
$(_doc_external("DM/DMGetType"))
"""
function gettype(::DM) end

@for_petsc function gettype(dm::DM{$PetscLib})
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
    view(dm::DM, viewer::Viewer=ViewerStdout(petsclib, getcomm(dm)))

view a `dm` with `viewer`

# External Links
$(_doc_external("DM/DMView"))
"""
function view(::DM) end

@for_petsc function view(
    dm::DM{$PetscLib},
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
    DMCreateMatrix(dm::DM)

Generates a matrix from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateMatrix"))
"""
function DMCreateMatrix end

@for_petsc function DMCreateMatrix(dm::DM{$PetscLib})
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
    DMCreateLocalVector(dm::DM)

returns a local vector from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateLocalVector"))
"""
function DMCreateLocalVector end

@for_petsc function DMCreateLocalVector(dm::DM{$PetscLib})
    vec = Vec{$PetscScalar}(C_NULL)

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
    DMCreateGlobalVector(dm::DM; write::Bool = true, read::Bool = true)

returns a global vector from the `dm` object.

# External Links
$(_doc_external("DM/DMCreateGlobalVector"))
"""
function DMCreateGlobalVector end

@for_petsc function DMCreateGlobalVector(dm::DM{$PetscLib})
    vec = Vec{$PetscScalar}(C_NULL)

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
    DMLocalToGlobal!(
        dm::DM
        local_vec::AbstractVec,
        mode::InsertMode,
        global_vec::AbstractVec,
    )

Updates `global_vec` from `local_vec` using the `dm` with insert `mode`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function DMLocalToGlobal! end

@for_petsc function DMLocalToGlobal!(
    dm::DM{$PetscLib},
    local_vec::AbstractVec,
    mode::InsertMode,
    global_vec::AbstractVec,
)
    @chk ccall(
        (:DMLocalToGlobal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        dm,
        local_vec,
        mode,
        global_vec,
    )

    return nothing
end

"""
    DMGlobalToLocal!(
        dm::DM
        global_vec::AbstractVec,
        mode::InsertMode,
        local_vec::AbstractVec,
    )

Updates `local_vec` from `global_vec` using the `dm` with insert `mode`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
function DMGlobalToLocal! end

@for_petsc function DMGlobalToLocal!(
    dm::DM{$PetscLib},
    global_vec::AbstractVec,
    mode::InsertMode,
    local_vec::AbstractVec,
)
    @chk ccall(
        (:DMGlobalToLocal, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, InsertMode, CVec),
        dm,
        global_vec,
        mode,
        local_vec,
    )

    return nothing
end

"""
    DMGetCoordinateDM(
        dm::DM
    )

Create a `DM` for the coordinates of `dm`.

# External Links
$(_doc_external("DM/DMGetCoordinateDM"))
"""
function DMGetCoordinateDM end

@for_petsc function DMGetCoordinateDM(
    dm::DM{$PetscLib},
)
    coord_dm = DM{$PetscLib}(C_NULL)
    @chk ccall(
        (:DMGetCoordinateDM, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CDM}),
        dm,
        coord_dm,
    )

    return coord_dm
end

"""
    DMGetCoordinatesLocal(
        dm::DM
    )

Gets a local vector with the coordinates associated with `dm`.

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function DMGetCoordinatesLocal end

@for_petsc function DMGetCoordinatesLocal(
    dm::DM{$PetscLib},
)
    coord_vec = Vec{$PetscScalar}(C_NULL)
    @chk ccall(
        (:DMGetCoordinatesLocal, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{CVec}),
        dm,
        coord_vec,
    )

    return coord_vec
end
